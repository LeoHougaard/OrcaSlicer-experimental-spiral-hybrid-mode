#include "ContinuousFermat.hpp"

#include "ClipperUtils.hpp"
#include "Exception.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Layer.hpp"
#include "Print.hpp"

#include <algorithm>
#include <array>
#include <boost/log/trivial.hpp>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace Slic3r {
namespace ContinuousFermat {
namespace {

static constexpr double EPS = 1e-9;
static constexpr double PI = 3.141592653589793238462643383279502884;
// Pairwise topology repair is superlinear and becomes counterproductive on
// very dense paths. Above this size, preserve the structurally valid route and
// report any remaining topology defects as geometric advisories.
static constexpr size_t MAX_EXHAUSTIVE_REPAIR_POINTS = 5000;

using PerfClock = std::chrono::steady_clock;

struct PerfCounter
{
    size_t calls { 0 };
    double milliseconds { 0.0 };
};

struct PerfCounters
{
    PerfCounter pair_metrics;
    PerfCounter containment;
};

thread_local PerfCounters perf_counters;

bool perf_timing_enabled()
{
    static const bool enabled = std::getenv("CONTINUOUS_FERMAT_PROFILE") != nullptr;
    return enabled;
}

double elapsed_milliseconds(const PerfClock::time_point start)
{
    return std::chrono::duration<double, std::milli>(PerfClock::now() - start).count();
}

class AggregatePerfTimer
{
public:
    explicit AggregatePerfTimer(PerfCounter &counter) :
        m_counter(counter), m_enabled(perf_timing_enabled()), m_start(m_enabled ? PerfClock::now() : PerfClock::time_point())
    {}
    ~AggregatePerfTimer()
    {
        if (m_enabled) {
            ++m_counter.calls;
            m_counter.milliseconds += elapsed_milliseconds(m_start);
        }
    }

private:
    PerfCounter &m_counter;
    bool m_enabled;
    PerfClock::time_point m_start;
};

class CandidatePerfTimer
{
public:
    CandidatePerfTimer(const size_t area_count, const double spacing_factor) :
        m_area_count(area_count), m_spacing_factor(spacing_factor), m_before(perf_counters), m_start(PerfClock::now())
    {}

    ~CandidatePerfTimer()
    {
        if (!perf_timing_enabled())
            return;
        BOOST_LOG_TRIVIAL(warning) << "Continuous Fermat profile candidate_ms=" << elapsed_milliseconds(m_start)
                                   << " spacing_factor=" << m_spacing_factor << " areas=" << m_area_count
                                   << " contours=" << contour_count << " output_points=" << output_points
                                   << " pair_calls=" << perf_counters.pair_metrics.calls - m_before.pair_metrics.calls
                                   << " pair_ms=" << perf_counters.pair_metrics.milliseconds - m_before.pair_metrics.milliseconds
                                   << " containment_calls=" << perf_counters.containment.calls - m_before.containment.calls
                                   << " containment_ms=" << perf_counters.containment.milliseconds - m_before.containment.milliseconds;
    }

    size_t contour_count { 0 };
    size_t output_points { 0 };

private:
    size_t m_area_count;
    double m_spacing_factor;
    PerfCounters m_before;
    PerfClock::time_point m_start;
};

double maximum_physical_width(const Flow &flow, const double configured_max_line_width)
{
    const double hard_limit = scale_(double(flow.nozzle_diameter()) * 2.0);
    if (!(hard_limit > EPS) || !std::isfinite(hard_limit))
        return 0.0;
    if (configured_max_line_width > 0.0 && std::isfinite(configured_max_line_width))
        return std::min(hard_limit, scale_(configured_max_line_width));

    const double line_height = scale_(double(flow.height()));
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double nominal_cross_section = double(flow.scaled_width()) - rounded_corner_loss;
    return std::min(hard_limit, rounded_corner_loss + nominal_cross_section * 1.60);
}

double maximum_extrusion_multiplier(const Flow &flow, const double configured_max_line_width)
{
    const double line_height = scale_(double(flow.height()));
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double nominal_cross_section = double(flow.scaled_width()) - rounded_corner_loss;
    if (!(nominal_cross_section > EPS))
        return 0.0;
    return (maximum_physical_width(flow, configured_max_line_width) - rounded_corner_loss) /
           nominal_cross_section;
}

double minimum_physical_width(const Flow &flow)
{
    const double rounded_corner_loss = scale_(double(flow.height())) * (1.0 - 0.25 * PI);
    return std::max(scale_(double(flow.nozzle_diameter()) * 0.05), rounded_corner_loss + scale_(0.001));
}

double minimum_extrusion_multiplier(const Flow &flow)
{
    const double rounded_corner_loss = scale_(double(flow.height())) * (1.0 - 0.25 * PI);
    const double nominal_cross_section = double(flow.scaled_width()) - rounded_corner_loss;
    return nominal_cross_section > EPS ?
        (minimum_physical_width(flow) + scale_(0.001) - rounded_corner_loss) / nominal_cross_section : 1.0;
}

struct ContourLoop
{
    size_t level_index { 0 };
    Points points;
    double area { 0.0 };
    double length { 0.0 };
    Point centroid;
    bool closed { true };
    double line_width { 0.0 };
    double flow_multiplier { 1.0 };
};

double point_distance(const Point &a, const Point &b)
{
    return (a - b).cast<double>().norm();
}

double point_distance2(const Point &a, const Point &b)
{
    return (a - b).cast<double>().squaredNorm();
}

Point lerp_point(const Point &a, const Point &b, const double t)
{
    return Point(
        double(a.x()) + (double(b.x()) - double(a.x())) * t,
        double(a.y()) + (double(b.y()) - double(a.y())) * t);
}

double polyline_length(const Points &points, const bool closed)
{
    if (points.size() < 2)
        return 0.0;

    double out = 0.0;
    for (size_t i = 1; i < points.size(); ++i)
        out += point_distance(points[i - 1], points[i]);
    if (closed)
        out += point_distance(points.back(), points.front());
    return out;
}

double point_segment_distance(const Point &p, const Point &a, const Point &b)
{
    const Vec2d ap = (p - a).cast<double>();
    const Vec2d ab = (b - a).cast<double>();
    const double denom = ab.squaredNorm();
    const double t = denom <= EPS ? 0.0 : std::clamp(ap.dot(ab) / denom, 0.0, 1.0);
    const Vec2d q = a.cast<double>() + ab * t;
    return (p.cast<double>() - q).norm();
}

double distance_to_loop(const Point &p, const Points &loop)
{
    if (loop.empty())
        return std::numeric_limits<double>::infinity();

    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < loop.size(); ++i)
        best = std::min(best, point_segment_distance(p, loop[i], loop[(i + 1) % loop.size()]));
    return best;
}

double orient2d(const Point &a, const Point &b, const Point &c)
{
    const Vec2d ab = (b - a).cast<double>();
    const Vec2d ac = (c - a).cast<double>();
    return cross2(ab, ac);
}

bool on_segment(const Point &a, const Point &b, const Point &p)
{
    return std::min(a.x(), b.x()) <= p.x() && p.x() <= std::max(a.x(), b.x()) &&
           std::min(a.y(), b.y()) <= p.y() && p.y() <= std::max(a.y(), b.y()) &&
           std::abs(orient2d(a, b, p)) <= 1e-7;
}

bool segments_intersect(const Point &a, const Point &b, const Point &c, const Point &d)
{
    const double o1 = orient2d(a, b, c);
    const double o2 = orient2d(a, b, d);
    const double o3 = orient2d(c, d, a);
    const double o4 = orient2d(c, d, b);

    if (((o1 > 1e-7 && o2 < -1e-7) || (o1 < -1e-7 && o2 > 1e-7)) &&
        ((o3 > 1e-7 && o4 < -1e-7) || (o3 < -1e-7 && o4 > 1e-7)))
        return true;

    return on_segment(a, b, c) || on_segment(a, b, d) || on_segment(c, d, a) || on_segment(c, d, b);
}

double segment_distance(const Point &a, const Point &b, const Point &c, const Point &d)
{
    if (segments_intersect(a, b, c, d))
        return 0.0;
    return std::min({
        point_segment_distance(a, c, d),
        point_segment_distance(b, c, d),
        point_segment_distance(c, a, b),
        point_segment_distance(d, a, b),
    });
}

std::string layer_diagnostics(const Layer &layer, const ExPolygons &printable_area, const Flow *flow = nullptr)
{
    size_t non_empty_regions = 0;
    for (const LayerRegion *region : layer.regions())
        if (region != nullptr && !region->slices.empty())
            ++non_empty_regions;

    std::ostringstream out;
    out << "layer=" << layer.id() << " print_z=" << layer.print_z << " regions=" << layer.regions().size()
        << " non_empty_regions=" << non_empty_regions << " printable_expolygons=" << printable_area.size();
    if (flow != nullptr)
        out << " flow_width=" << flow->width() << " flow_height=" << flow->height()
            << " flow_spacing=" << unscale<double>(flow->scaled_spacing());
    return out.str();
}

Point polygon_centroid_or_first(const Points &points)
{
    if (points.empty())
        return Point(0, 0);

    double area2 = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        const Point &p = points[i];
        const Point &q = points[(i + 1) % points.size()];
        const double cross = double(p.x()) * double(q.y()) - double(q.x()) * double(p.y());
        area2 += cross;
        cx += (double(p.x()) + double(q.x())) * cross;
        cy += (double(p.y()) + double(q.y())) * cross;
    }
    if (std::abs(area2) <= EPS)
        return points.front();
    return Point(cx / (3.0 * area2), cy / (3.0 * area2));
}

double signed_area(const Points &points)
{
    double area2 = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        const Point &p = points[i];
        const Point &q = points[(i + 1) % points.size()];
        area2 += double(p.x()) * double(q.y()) - double(q.x()) * double(p.y());
    }
    return area2 * 0.5;
}

void append_point(Points &path, const Point &p)
{
    if (path.empty() || path.back() != p)
        path.emplace_back(p);
}

void append_points(Points &path, const Points &points)
{
    for (const Point &point : points)
        append_point(path, point);
}

Point point_at_closed_fraction(const Points &points, double fraction)
{
    if (points.empty())
        return Point(0, 0);
    if (points.size() == 1)
        return points.front();

    fraction = fraction - std::floor(fraction);
    const double total = polyline_length(points, true);
    if (total <= EPS)
        return points.front();

    const double target = total * fraction;
    double walked = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        const Point &p = points[i];
        const Point &q = points[(i + 1) % points.size()];
        const double len = point_distance(p, q);
        if (walked + len >= target)
            return lerp_point(p, q, (target - walked) / std::max(len, EPS));
        walked += len;
    }
    return points.back();
}

double loop_point_index_mod(const Points &points, double index)
{
    const double n = double(points.size());
    index = std::fmod(index, n);
    return index < 0.0 ? index + n : index;
}

Point loop_point_at_index(const Points &points, double index)
{
    if (points.empty())
        return Point(0, 0);

    index = loop_point_index_mod(points, index);
    const size_t base = size_t(std::floor(index)) % points.size();
    const double alpha = index - std::floor(index);
    return lerp_point(points[base], points[(base + 1) % points.size()], alpha);
}

double loop_segment_length(const Points &points, const size_t index)
{
    return point_distance(points[index % points.size()], points[(index + 1) % points.size()]);
}

double loop_nearest_index(const Points &points, const Point &p)
{
    double best_idx = 0.0;
    double best_dist = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < points.size(); ++i) {
        const Point &a = points[i];
        const Point &b = points[(i + 1) % points.size()];
        const Vec2d ab = (b - a).cast<double>();
        const double denom = ab.squaredNorm();
        const double t = denom <= EPS ? 0.0 : std::clamp((p - a).cast<double>().dot(ab) / denom, 0.0, 1.0);
        const Vec2d q = a.cast<double>() + ab * t;
        const double d = (p.cast<double>() - q).squaredNorm();
        if (d < best_dist) {
            best_dist = d;
            best_idx = double(i) + t;
        }
    }

    return loop_point_index_mod(points, best_idx);
}

double loop_furthest_vertex_index(const Points &points, const Point &p)
{
    size_t best_idx = 0;
    double best_dist = -1.0;
    for (size_t i = 0; i < points.size(); ++i) {
        const double d = point_distance2(points[i], p);
        if (d > best_dist) {
            best_dist = d;
            best_idx = i;
        }
    }
    return double(best_idx);
}

double loop_back_by_distance(const Points &points, double index, double distance_along);

double loop_forward_by_distance(const Points &points, double index, double distance_along)
{
    if (distance_along < 0.0)
        return loop_back_by_distance(points, index, -distance_along);

    const size_t n = points.size();
    index = loop_point_index_mod(points, index);
    const Point p = loop_point_at_index(points, index);
    const size_t ceil_idx = size_t(std::ceil(index)) % n;
    double walked = point_distance(p, points[ceil_idx]);
    if (walked >= distance_along && walked > EPS) {
        const double alpha = distance_along / walked;
        return loop_point_index_mod(points, std::ceil(index) * alpha + index * (1.0 - alpha));
    }

    size_t i = ceil_idx;
    for (size_t guard = 0; guard < n + 2; ++guard) {
        const double segment = loop_segment_length(points, i);
        if (walked + segment <= distance_along) {
            walked += segment;
            i = (i + 1) % n;
        } else {
            const double alpha = (distance_along - walked) / std::max(segment, EPS);
            return loop_point_index_mod(points, double(i) + alpha);
        }
    }

    return index;
}

double loop_back_by_distance(const Points &points, double index, double distance_along)
{
    if (distance_along < 0.0)
        return loop_forward_by_distance(points, index, -distance_along);

    const size_t n = points.size();
    index = loop_point_index_mod(points, index);
    const Point p = loop_point_at_index(points, index);
    const size_t floor_idx = size_t(std::floor(index)) % n;
    double walked = point_distance(p, points[floor_idx]);
    if (walked >= distance_along && walked > EPS) {
        const double alpha = distance_along / walked;
        return loop_point_index_mod(points, double(floor_idx) * alpha + index * (1.0 - alpha));
    }

    size_t i = (floor_idx + n - 1) % n;
    for (size_t guard = 0; guard < n + 2; ++guard) {
        const double segment = loop_segment_length(points, i);
        if (walked + segment <= distance_along) {
            walked += segment;
            i = (i + n - 1) % n;
        } else {
            const double alpha = (distance_along - walked) / std::max(segment, EPS);
            return loop_point_index_mod(points, double(i) + 1.0 - alpha);
        }
    }

    return index;
}

bool subset_cycle_param(const double left, const double right, double query, const bool close_left = false, const bool close_right = false)
{
    if (std::abs(query - left) <= 1e-9)
        return close_left;
    if (std::abs(query - right) <= 1e-9)
        return close_right;
    if (std::abs(left - right) <= 1e-9)
        return false;
    if (left < right)
        return left < query && query < right;
    return query > left || query < right;
}

double loop_length_between(const Points &points, double start, double end)
{
    const size_t n = points.size();
    start = loop_point_index_mod(points, start);
    end = loop_point_index_mod(points, end);
    const size_t start_floor = size_t(std::floor(start));
    const size_t end_floor = size_t(std::floor(end));
    if (start_floor == end_floor && end > start)
        return (end - start) * loop_segment_length(points, start_floor);

    double length = (1.0 - start + double(start_floor)) * loop_segment_length(points, start_floor);
    length += (end - double(end_floor)) * loop_segment_length(points, end_floor);
    size_t i = (start_floor + 1) % n;
    while (i != end_floor) {
        length += loop_segment_length(points, i);
        i = (i + 1) % n;
    }
    return length;
}

void append_loop_param(Points &path, const Points &points, const double index)
{
    append_point(path, loop_point_at_index(points, index));
}

Points loop_arc_between_params(const Points &points, double start, double end, const int direction)
{
    if (points.empty())
        return {};

    const size_t n = points.size();
    start = loop_point_index_mod(points, start);
    end = loop_point_index_mod(points, end);

    Points out;
    out.reserve(n + 2);
    append_loop_param(out, points, start);

    if (direction >= 0) {
        size_t i = size_t(std::ceil(start)) % n;
        if (std::abs(double(i) - start) <= 1e-6)
            i = (i + 1) % n;
        for (size_t guard = 0; subset_cycle_param(start, end, double(i), false, false) && guard <= n + 2;
             i = (i + 1) % n, ++guard)
            append_point(out, points[i]);
    } else {
        size_t i = size_t(std::floor(start)) % n;
        if (std::abs(double(i) - start) <= 1e-6)
            i = (i + n - 1) % n;
        for (size_t guard = 0; subset_cycle_param(end, start, double(i), false, false) && guard <= n + 2;
             i = (i + n - 1) % n, ++guard)
            append_point(out, points[i]);
    }

    append_loop_param(out, points, end);
    return out;
}

void append_loop_vertices_forward(Points &path, const Points &points, size_t start_index, size_t stop_index)
{
    const size_t n = points.size();
    size_t i = start_index % n;
    for (size_t guard = 0; i != stop_index % n && guard <= n + 2; ++guard) {
        append_point(path, points[i]);
        i = (i + 1) % n;
    }
}

void append_loop_vertices_backward_until(Points &path, const Points &points, size_t start_index, const double left, const double right)
{
    const size_t n = points.size();
    size_t i = start_index % n;
    for (size_t guard = 0; subset_cycle_param(left, double(i), right, false, false) && guard <= n + 2; ++guard) {
        append_point(path, points[i]);
        i = (i + n - 1) % n;
    }
}

double sampled_loop_gap(const ContourLoop &a, const ContourLoop &b)
{
    double best = std::numeric_limits<double>::infinity();
    const size_t stride_a = std::max<size_t>(1, a.points.size() / 128);
    const size_t stride_b = std::max<size_t>(1, b.points.size() / 128);

    for (size_t i = 0; i < a.points.size(); i += stride_a)
        best = std::min(best, distance_to_loop(a.points[i], b.points));
    for (size_t i = 0; i < b.points.size(); i += stride_b)
        best = std::min(best, distance_to_loop(b.points[i], a.points));

    return best;
}

Points resample_closed_loop(const Points &points, const double target_step)
{
    if (points.size() < 2)
        return points;

    const double length = polyline_length(points, true);
    if (length <= EPS)
        return points;

    const size_t sample_count = std::max<size_t>(16, size_t(std::ceil(length / std::max(target_step, EPS))));
    std::vector<double> segment_lengths;
    segment_lengths.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        segment_lengths.emplace_back(point_distance(points[i], points[(i + 1) % points.size()]));

    Points out;
    out.reserve(sample_count);
    size_t seg_idx = 0;
    double seg_start_distance = 0.0;
    for (size_t sample_idx = 0; sample_idx < sample_count; ++sample_idx) {
        const double target = length * double(sample_idx) / double(sample_count);
        while (seg_idx + 1 < points.size() && seg_start_distance + segment_lengths[seg_idx] < target) {
            seg_start_distance += segment_lengths[seg_idx];
            ++seg_idx;
        }

        const double segment = std::max(segment_lengths[seg_idx], EPS);
        const double t = (target - seg_start_distance) / segment;
        out.emplace_back(lerp_point(points[seg_idx], points[(seg_idx + 1) % points.size()], t));
    }

    return out;
}

Points densify_closed_loop_preserving_vertices(const Points &points, const double max_step)
{
    if (points.size() < 2)
        return points;

    Points simplified;
    simplified.reserve(points.size());
    const double collinear_tolerance = std::max(max_step * 0.002, 1.0);
    for (size_t i = 0; i < points.size(); ++i) {
        const Point &prev = points[(i + points.size() - 1) % points.size()];
        const Point &point = points[i];
        const Point &next = points[(i + 1) % points.size()];
        if (point_distance(prev, next) > EPS && point_segment_distance(point, prev, next) <= collinear_tolerance)
            continue;
        simplified.emplace_back(point);
    }
    const Points &source = simplified.size() >= 3 ? simplified : points;

    Points out;
    out.reserve(source.size());
    for (size_t i = 0; i < source.size(); ++i) {
        const Point &a = source[i];
        const Point &b = source[(i + 1) % source.size()];
        append_point(out, a);

        const double length = point_distance(a, b);
        const size_t parts = size_t(std::ceil(length / std::max(max_step, EPS)));
        for (size_t part = 1; part < parts; ++part)
            append_point(out, lerp_point(a, b, double(part) / double(parts)));
    }

    return out;
}

std::map<size_t, std::vector<ContourLoop>> loops_by_level(const std::vector<ContourLoop> &contours)
{
    std::map<size_t, std::vector<ContourLoop>> grouped;
    for (const ContourLoop &loop : contours)
        grouped[loop.level_index].emplace_back(loop);
    return grouped;
}

std::vector<ContourLoop> filter_ring_medial_overlap(const std::vector<ContourLoop> &contours, const double spacing)
{
    const auto grouped = loops_by_level(contours);
    std::vector<ContourLoop> stable;

    for (const auto &[level, loops] : grouped) {
        if (loops.size() == 2) {
            const ContourLoop *a = &loops[0];
            const ContourLoop *b = &loops[1];
            if (a->area < b->area)
                std::swap(a, b);
            if (sampled_loop_gap(*a, *b) < spacing * 0.98)
                break;
        }
        stable.insert(stable.end(), loops.begin(), loops.end());
    }

    return stable.empty() ? contours : stable;
}

double open_polyline_distance(const Point &p, const Points &points)
{
    if (points.empty())
        return std::numeric_limits<double>::infinity();
    if (points.size() == 1)
        return point_distance(p, points.front());

    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 1; i < points.size(); ++i)
        best = std::min(best, point_segment_distance(p, points[i - 1], points[i]));
    return best;
}

double contour_distance(const Point &p, const ContourLoop &contour)
{
    return contour.closed ? distance_to_loop(p, contour.points) : open_polyline_distance(p, contour.points);
}

struct TerminalAxes
{
    Vec2d center;
    Vec2d major;
    Vec2d minor;
    double min_major { 0.0 };
    double max_major { 0.0 };
    double min_minor { 0.0 };
    double max_minor { 0.0 };
};

bool principal_terminal_axes(const Points &points, TerminalAxes &axes)
{
    if (points.size() < 4)
        return false;

    axes.center = Vec2d::Zero();
    for (const Point &point : points)
        axes.center += point.cast<double>();
    axes.center /= double(points.size());

    double xx = 0.0;
    double xy = 0.0;
    double yy = 0.0;
    for (const Point &point : points) {
        const Vec2d delta = point.cast<double>() - axes.center;
        xx += delta.x() * delta.x();
        xy += delta.x() * delta.y();
        yy += delta.y() * delta.y();
    }
    if (xx + yy <= EPS)
        return false;

    const double angle = 0.5 * std::atan2(2.0 * xy, xx - yy);
    axes.major = Vec2d(std::cos(angle), std::sin(angle));
    axes.minor = Vec2d(-axes.major.y(), axes.major.x());

    auto extrema = [&](const Vec2d &axis) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (const Point &point : points) {
            const double value = (point.cast<double>() - axes.center).dot(axis);
            lo = std::min(lo, value);
            hi = std::max(hi, value);
        }
        return std::pair<double, double>(lo, hi);
    };

    std::tie(axes.min_major, axes.max_major) = extrema(axes.major);
    std::tie(axes.min_minor, axes.max_minor) = extrema(axes.minor);
    if (axes.max_minor - axes.min_minor > axes.max_major - axes.min_major) {
        std::swap(axes.major, axes.minor);
        std::swap(axes.min_major, axes.min_minor);
        std::swap(axes.max_major, axes.max_minor);
    }

    return true;
}

std::vector<ContourLoop> add_terminal_medial_gap_fill(
    const std::vector<ContourLoop> &contours,
    const double line_width,
    const double line_height,
    const double spacing,
    const double max_physical_width)
{
    if (contours.empty())
        return contours;

    size_t last_level = 0;
    for (const ContourLoop &loop : contours)
        last_level = std::max(last_level, loop.level_index);

    size_t terminal_index = 0;
    size_t terminal_count = 0;
    for (size_t i = 0; i < contours.size(); ++i) {
        if (contours[i].level_index == last_level) {
            terminal_index = i;
            ++terminal_count;
        }
    }
    if (terminal_count != 1)
        return contours;

    const ContourLoop &terminal = contours[terminal_index];
    if (!terminal.closed || terminal.points.size() < 3)
        return contours;

    TerminalAxes axes;
    if (!principal_terminal_axes(terminal.points, axes))
        return contours;

    const double major_span = axes.max_major - axes.min_major;
    const double minor_span = axes.max_minor - axes.min_minor;
    if (major_span < spacing * 6.0)
        return contours;
    if (!(line_width * 1.02 < minor_span && minor_span < line_width * 2.05))
        return contours;

    // A closed terminal loop contributes two opposing passes.  Redistribute
    // the remaining band over those two passes instead of adding a third open
    // medial pass: an odd three-pass ladder cannot reconnect to same-side CFS
    // ports without a long retrace.
    const double adaptive_width = (minor_span + line_width) / 2.0;
    if (!(line_width * 1.01 <= adaptive_width && adaptive_width < line_width * 1.55 &&
          adaptive_width <= max_physical_width + EPS))
        return contours;

    const double mid_minor = 0.5 * (axes.min_minor + axes.max_minor);
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double nominal_cross_section = line_width - rounded_corner_loss;
    const double adaptive_cross_section = adaptive_width - rounded_corner_loss;
    if (!(nominal_cross_section > EPS) || !(adaptive_cross_section > nominal_cross_section))
        return contours;

    ContourLoop adjusted = terminal;
    adjusted.points.clear();
    adjusted.points.reserve(terminal.points.size());
    for (const Point &point : terminal.points) {
        const Vec2d rel = point.cast<double>() - axes.center;
        const double major_coord = rel.dot(axes.major);
        const double minor_coord = rel.dot(axes.minor);
        const double side = minor_coord >= mid_minor ? 1.0 : -1.0;
        const Vec2d remapped =
            axes.center + axes.major * major_coord + axes.minor * (mid_minor + side * adaptive_width * 0.5);
        append_point(adjusted.points, Point(remapped.x(), remapped.y()));
    }
    adjusted.points = densify_closed_loop_preserving_vertices(adjusted.points, std::max(spacing * 0.25, 1.0));
    if (signed_area(adjusted.points) < 0.0)
        std::reverse(adjusted.points.begin(), adjusted.points.end());
    adjusted.area = std::abs(signed_area(adjusted.points));
    adjusted.length = polyline_length(adjusted.points, true);
    adjusted.centroid = polygon_centroid_or_first(adjusted.points);
    adjusted.closed = true;
    adjusted.line_width = adaptive_width;
    adjusted.flow_multiplier = adaptive_cross_section / nominal_cross_section;

    std::vector<ContourLoop> out = contours;
    out[terminal_index] = std::move(adjusted);
    std::sort(out.begin(), out.end(), [](const ContourLoop &a, const ContourLoop &b) {
        if (a.level_index != b.level_index)
            return a.level_index < b.level_index;
        return a.area > b.area;
    });
    return out;
}

std::vector<ContourLoop> generate_offset_contours(
    const ExPolygons &printable_area,
    const double line_width,
    const double spacing,
    const size_t precise_wall_levels)
{
    const double first_offset = line_width * 0.5;
    const double min_length = spacing * 2.0;

    std::vector<ContourLoop> contours;
    ExPolygons source = union_ex(printable_area);
    ExPolygons incremental_source;
    for (size_t level = 0; level < 512; ++level) {
        const double offset = first_offset + double(level) * spacing;
        ExPolygons inset;
        if (level <= precise_wall_levels) {
            inset = offset_ex(source, float(-offset), ClipperLib::jtMiter, 3.0);
        } else {
            inset = offset_ex(incremental_source, float(-spacing), ClipperLib::jtMiter, 3.0);
        }
        if (inset.empty())
            break;
        if (level >= precise_wall_levels)
            incremental_source = inset;

        for (const ExPolygon &expoly : inset) {
            auto add_polygon = [&](const Polygon &polygon) {
                if (polygon.points.size() < 3)
                    return;
                ContourLoop loop;
                loop.level_index = level;
                loop.points = level < precise_wall_levels ?
                    densify_closed_loop_preserving_vertices(polygon.points, std::max(spacing * 0.25, 1.0)) :
                    resample_closed_loop(polygon.points, std::max(spacing * 0.65, 1.0));
                loop.area = std::abs(double(polygon.area()));
                loop.length = polyline_length(loop.points, true);
                loop.centroid = polygon_centroid_or_first(loop.points);
                // A small-area inset may still be the only centerline capable
                // of covering a terminal core (for example, a 20 mm square at
                // 1.2 mm line width).  Area alone used to discard that valid
                // loop.  Conversely, perimeter alone retains arbitrarily thin
                // terminal slivers and creates near-retracing strokes.  The
                // hydraulic-diameter proxy distinguishes the two without
                // weakening the unchanged mandatory final gate.
                const double thickness_proxy = loop.length > EPS ? 4.0 * loop.area / loop.length : 0.0;
                if (loop.length >= min_length && thickness_proxy >= spacing * 0.35)
                    contours.emplace_back(std::move(loop));
            };
            add_polygon(expoly.contour);
            for (const Polygon &hole : expoly.holes)
                add_polygon(hole);
        }
    }

    const auto grouped = loops_by_level(contours);
    if (!grouped.empty()) {
        const auto &first = grouped.begin()->second;
        const bool one_hole_ring = first.size() == 2 &&
            point_distance(first[0].centroid, first[1].centroid) < spacing * 2.0;
        if (one_hole_ring)
            contours = filter_ring_medial_overlap(contours, spacing);
    }

    return contours;
}

std::vector<ContourLoop> generate_offset_contours(const ExPolygons &printable_area, const Flow &flow)
{
    return generate_offset_contours(printable_area, double(flow.scaled_width()), double(flow.scaled_spacing()), 3);
}

std::vector<ContourLoop> generate_narrow_feature_contours(
    const ExPolygons &printable_area,
    const Flow &flow,
    const double spacing)
{
    if (printable_area.empty())
        return {};

    const BoundingBox bounds = get_extents(printable_area);
    const double min_span = double(std::min(bounds.max.x() - bounds.min.x(), bounds.max.y() - bounds.min.y()));
    if (!(min_span > EPS))
        return {};

    const double nominal_width = double(flow.scaled_width());
    const double rounded_corner_loss = scale_(double(flow.height())) * (1.0 - 0.25 * PI);
    const double nominal_cross_section = nominal_width - rounded_corner_loss;
    double adaptive_width = std::min(nominal_width, min_span * 0.55);
    adaptive_width = std::max(adaptive_width, minimum_physical_width(flow));

    ExPolygons inset;
    for (size_t attempt = 0; attempt < 8 && inset.empty(); ++attempt) {
        inset = offset_ex(printable_area, float(-adaptive_width * 0.5), ClipperLib::jtMiter, 3.0);
        if (inset.empty())
            adaptive_width *= 0.78;
    }
    if (inset.empty() || adaptive_width + EPS < minimum_physical_width(flow))
        return {};

    const double adaptive_cross_section = adaptive_width - rounded_corner_loss;
    if (!(adaptive_cross_section > EPS) || !(nominal_cross_section > EPS))
        return {};

    std::vector<ContourLoop> contours;
    for (const ExPolygon &expoly : inset) {
        const auto add_loop = [&](const Polygon &polygon) {
            if (polygon.points.size() < 3)
                return;
            ContourLoop loop;
            loop.points = densify_closed_loop_preserving_vertices(
                polygon.points, std::max(std::min(spacing, adaptive_width) * 0.25, 1.0));
            loop.area = std::abs(double(polygon.area()));
            loop.length = polyline_length(loop.points, true);
            loop.centroid = polygon_centroid_or_first(loop.points);
            loop.line_width = adaptive_width;
            loop.flow_multiplier = adaptive_cross_section / nominal_cross_section;
            contours.emplace_back(std::move(loop));
        };
        add_loop(expoly.contour);
        for (const Polygon &hole : expoly.holes)
            add_loop(hole);
    }
    return contours;
}

std::vector<ContourLoop> add_ring_medial_contours(
    const ExPolygons &printable_area,
    std::vector<ContourLoop> contours,
    const Flow &flow)
{
    const auto grouped = loops_by_level(contours);
    if (grouped.empty() || grouped.begin()->second.size() != 2 ||
        point_distance(grouped.begin()->second[0].centroid, grouped.begin()->second[1].centroid) >=
            double(flow.scaled_spacing()) * 2.0)
        return contours;

    Polygons swept;
    for (const ContourLoop &loop : contours) {
        if (!loop.closed || loop.points.size() < 4)
            continue;
        Points closed = loop.points;
        append_point(closed, closed.front());
        Polygons bead = offset(
            Polyline(std::move(closed)),
            float(double(flow.scaled_width()) * 0.5),
            ClipperLib::jtRound,
            SCALED_RESOLUTION,
            ClipperLib::etOpenRound);
        swept.insert(swept.end(), std::make_move_iterator(bead.begin()), std::make_move_iterator(bead.end()));
    }
    if (swept.empty())
        return contours;

    const ExPolygons residual = diff_ex(printable_area, union_(swept));
    if (residual.empty())
        return contours;

    std::vector<ContourLoop> medial = generate_offset_contours(
        residual,
        double(flow.scaled_width()),
        double(flow.scaled_spacing()),
        0);
    if (medial.empty())
        return contours;

    size_t first_medial_level = medial.front().level_index;
    for (const ContourLoop &loop : medial)
        first_medial_level = std::min(first_medial_level, loop.level_index);
    medial.erase(
        std::remove_if(medial.begin(), medial.end(), [first_medial_level](const ContourLoop &loop) {
            return loop.level_index != first_medial_level;
        }),
        medial.end());
    if (medial.size() != 2)
        return contours;

    size_t next_level = 0;
    for (const ContourLoop &loop : contours)
        next_level = std::max(next_level, loop.level_index + 1);
    for (ContourLoop &loop : medial)
        loop.level_index = next_level;
    contours.insert(contours.end(), std::make_move_iterator(medial.begin()), std::make_move_iterator(medial.end()));
    return contours;
}

std::vector<ContourLoop> ordered_loops_for_spiral(const std::vector<ContourLoop> &contours)
{
    const auto grouped = loops_by_level(contours);
    if (grouped.empty())
        return {};

    const size_t first_count = grouped.begin()->second.size();
    if (first_count == 1) {
        std::vector<ContourLoop> ordered;
        for (const auto &[level, loops] : grouped) {
            std::vector<ContourLoop> sorted = loops;
            std::sort(sorted.begin(), sorted.end(), [](const ContourLoop &a, const ContourLoop &b) {
                if (a.area != b.area)
                    return a.area > b.area;
                if (a.centroid.x() != b.centroid.x())
                    return a.centroid.x() < b.centroid.x();
                return a.centroid.y() < b.centroid.y();
            });
            ordered.insert(ordered.end(), sorted.begin(), sorted.end());
        }
        return ordered;
    }

    bool all_two = true;
    for (const auto &[level, loops] : grouped)
        all_two &= loops.size() == 2;
    if (first_count == 2 && all_two) {
        std::vector<ContourLoop> outer_family;
        std::vector<ContourLoop> hole_family;
        for (const auto &[level, loops] : grouped) {
            std::vector<ContourLoop> sorted = loops;
            std::sort(sorted.begin(), sorted.end(), [](const ContourLoop &a, const ContourLoop &b) { return a.area > b.area; });
            outer_family.emplace_back(sorted[0]);
            hole_family.emplace_back(sorted[1]);
        }
        std::reverse(hole_family.begin(), hole_family.end());
        outer_family.insert(outer_family.end(), hole_family.begin(), hole_family.end());
        return outer_family;
    }

    std::vector<ContourLoop> ordered;
    for (const auto &[level, loops] : grouped) {
        std::vector<ContourLoop> sorted = loops;
        std::sort(sorted.begin(), sorted.end(), [](const ContourLoop &a, const ContourLoop &b) { return a.area > b.area; });
        ordered.insert(ordered.end(), sorted.begin(), sorted.end());
    }
    return ordered;
}

Points terminal_medial_band_path(
    const ContourLoop &terminal,
    const ContourLoop &medial,
    const Point &start_hint,
    const Point &end_hint)
{
    if (terminal.line_width <= EPS || medial.points.size() < 2)
        return {};

    const Point &a = medial.points.front();
    const Point &b = medial.points.back();
    const Vec2d axis = (b - a).cast<double>();
    const double axis_length = axis.norm();
    if (axis_length <= EPS)
        return {};

    const Vec2d major = axis / axis_length;
    Vec2d minor(-major.y(), major.x());
    const Vec2d center = (a.cast<double>() + b.cast<double>()) * 0.5;

    if (!terminal.points.empty()) {
        const Point *farthest = &terminal.points.front();
        double best = -1.0;
        for (const Point &point : terminal.points) {
            const double value = std::abs((point.cast<double>() - center).dot(minor));
            if (value > best) {
                best = value;
                farthest = &point;
            }
        }
        if ((farthest->cast<double>() - center).dot(minor) < 0.0)
            minor = -minor;
    }

    auto side_point = [&](const int side, const double level) {
        const Vec2d base = side == 0 ? a.cast<double>() : b.cast<double>();
        const Vec2d p = base + minor * (level * terminal.line_width);
        return Point(p.x(), p.y());
    };

    std::vector<std::pair<double, Points>> candidates;
    for (const int start_side : { 0, 1 }) {
        for (const std::array<double, 3> levels : { std::array<double, 3>{ 1.0, 0.0, -1.0 },
                                                    std::array<double, 3>{ -1.0, 0.0, 1.0 } }) {
            int current_side = start_side;
            Points path;
            for (size_t idx = 0; idx < levels.size(); ++idx) {
                const int next_side = 1 - current_side;
                append_point(path, side_point(current_side, levels[idx]));
                append_point(path, side_point(next_side, levels[idx]));
                if (idx + 1 < levels.size())
                    append_point(path, side_point(next_side, levels[idx + 1]));
                current_side = next_side;
            }
            const double score = point_distance(start_hint, path.front()) + point_distance(path.back(), end_hint);
            candidates.emplace_back(score, std::move(path));
        }
    }

    if (candidates.empty())
        return {};
    return std::min_element(candidates.begin(), candidates.end(), [](const auto &a, const auto &b) { return a.first < b.first; })->second;
}

struct PairMetrics
{
    int crossings { 0 };
    int close_pairs { 0 };
    int bead_overlaps { 0 };
    double min_spacing { std::numeric_limits<double>::infinity() };
    std::string first_bead_overlap;
};

PairMetrics path_pair_metrics(
    const Points &path,
    double spacing,
    const std::vector<double> *segment_widths = nullptr,
    double nominal_line_width = 0.0,
    int crossing_limit = std::numeric_limits<int>::max());

Points build_single_minimum_connected_fermat(
    std::vector<ContourLoop> loops,
    const Point &start_anchor,
    const double spacing,
    const Point &exit_anchor,
    const bool preserve_medial_pockets = true,
    const double seam_port_factor = 1.17)
{
    loops.erase(std::remove_if(loops.begin(), loops.end(), [](const ContourLoop &loop) {
        return loop.points.size() < (loop.closed ? 4 : 2);
    }), loops.end());
    if (loops.empty())
        return {};

    const double seam_port_spacing = spacing * seam_port_factor;
    if (loops.size() == 1) {
        if (!loops.front().closed)
            return loops.front().points;

        const Points &points = loops.front().points;
        const size_t n = points.size();
        const size_t start = size_t(std::round(loop_nearest_index(points, start_anchor))) % n;
        Points path;
        path.reserve(n + 1);
        for (size_t i = 0; i <= n; ++i)
            append_point(path, points[(start + i) % n]);
        return path;
    }

    if (!preserve_medial_pockets) {
        while (loops.size() > 2 && loops.back().length < spacing * 20.0)
            loops.pop_back();
    }

    double in_index = loop_nearest_index(loops.front().points, start_anchor);
    double out_index = loop_nearest_index(loops.front().points, exit_anchor);
    bool out_forward_in = true;
    bool in_run = false;
    bool first_circle = true;
    Points in_branch;
    Points out_branch;
    size_t loop_index = 0;

    while (true) {
        const ContourLoop &loop = loops[loop_index];
        const Points &points = loop.points;
        const size_t n = points.size();

        if (loop.closed && loop_index + 1 < loops.size() && !loops[loop_index + 1].closed &&
            !in_branch.empty() && !out_branch.empty()) {
            Points band_path = terminal_medial_band_path(loop, loops[loop_index + 1], in_branch.back(), out_branch.back());
            if (!band_path.empty()) {
                append_points(in_branch, band_path);
                break;
            }
        }

        if (!loop.closed) {
            if (points.size() < 2)
                break;
            if (in_branch.empty() || out_branch.empty())
                return { points.front(), points.back() };

            const Point &a = points.front();
            const Point &b = points.back();
            const double forward_score = point_distance(in_branch.back(), a) + point_distance(out_branch.back(), b);
            const double reverse_score = point_distance(in_branch.back(), b) + point_distance(out_branch.back(), a);
            if (forward_score <= reverse_score) {
                append_point(in_branch, a);
                append_point(out_branch, b);
                append_point(in_branch, b);
            } else {
                append_point(in_branch, b);
                append_point(out_branch, a);
                append_point(in_branch, a);
            }
            break;
        }

        const double port_spacing = seam_port_spacing;
        const bool circle_small = loop.length < port_spacing * 2.0;

        if (circle_small) {
            const double near_index = loop_point_index_mod(points, 0.5 * (in_index + out_index));
            const Point near_pt = loop_point_at_index(points, near_index);
            const double far_index = loop_furthest_vertex_index(points, near_pt);
            const Point far_pt = loop_point_at_index(points, far_index);

            in_index = double(size_t(std::ceil(near_index)) % n);
            auto small_score = [&](double idx) {
                const Point &p = points[size_t(idx) % n];
                return std::min(point_distance(near_pt, p), point_distance(far_pt, p));
            };

            double best = small_score(in_index);
            for (size_t guard = 0, i = (size_t(in_index) + 1) % n;
                 subset_cycle_param(near_index, far_index, double(i)) && guard <= n + 2;
                 i = (i + 1) % n, ++guard) {
                const double score = small_score(double(i));
                if (score > best) {
                    best = score;
                    in_index = double(i);
                }
            }

            out_index = double(size_t(std::ceil(far_index)) % n);
            best = small_score(out_index);
            for (size_t guard = 0, i = (size_t(out_index) + 1) % n;
                 subset_cycle_param(far_index, near_index, double(i)) && guard <= n + 2;
                 i = (i + 1) % n, ++guard) {
                const double score = small_score(double(i));
                if (score > best) {
                    best = score;
                    out_index = double(i);
                }
            }

            out_forward_in = loop_length_between(points, in_index, out_index) <= loop_length_between(points, out_index, in_index);
        } else {
            if (std::abs(in_index - out_index) < 1e-6) {
                const double out_forward = loop_forward_by_distance(points, in_index, port_spacing);
                const double out_backward = loop_back_by_distance(points, in_index, port_spacing);
                if (first_circle && loop_index + 1 < loops.size()) {
                    const ContourLoop &next_loop = loops[loop_index + 1];
                    if (distance_to_loop(loop_point_at_index(points, out_forward), next_loop.points) <
                        distance_to_loop(loop_point_at_index(points, out_backward), next_loop.points)) {
                        out_index = out_backward;
                        out_forward_in = false;
                    } else {
                        out_index = out_forward;
                        out_forward_in = true;
                    }
                } else if (!in_branch.empty() && !out_branch.empty()) {
                    const Point in0 = loop_point_at_index(points, in_index);
                    const double out_forward = loop_forward_by_distance(points, in_index, port_spacing);
                    const double out_backward = loop_back_by_distance(points, in_index, port_spacing);
                    const Point p_forward = loop_point_at_index(points, out_forward);
                    const Point p_backward = loop_point_at_index(points, out_backward);
                    const Point &prev_in = in_branch.back();
                    const Point &prev_out = out_branch.back();

                    const double forward_score = std::min(point_distance(in0, prev_in) + point_distance(p_forward, prev_out),
                                                           point_distance(p_forward, prev_in) + point_distance(in0, prev_out));
                    const double backward_score = std::min(point_distance(in0, prev_in) + point_distance(p_backward, prev_out),
                                                            point_distance(p_backward, prev_in) + point_distance(in0, prev_out));
                    if (forward_score < backward_score) {
                        if (point_distance(in0, prev_in) + point_distance(p_forward, prev_out) <=
                            point_distance(p_forward, prev_in) + point_distance(in0, prev_out)) {
                            out_index = out_forward;
                            out_forward_in = true;
                        } else {
                            out_index = in_index;
                            in_index = out_forward;
                            out_forward_in = false;
                        }
                    } else {
                        if (point_distance(in0, prev_in) + point_distance(p_backward, prev_out) <=
                            point_distance(p_backward, prev_in) + point_distance(in0, prev_out)) {
                            out_index = out_backward;
                            out_forward_in = false;
                        } else {
                            out_index = in_index;
                            in_index = out_backward;
                            out_forward_in = true;
                        }
                    }
                } else {
                    out_index = out_forward;
                    out_forward_in = true;
                }
            } else {
                double length_io = loop_length_between(points, in_index, out_index);
                double length_oi = loop_length_between(points, out_index, in_index);
                out_forward_in = length_io <= length_oi;
                double arc_length = std::min(length_io, length_oi);
                if (first_circle && arc_length < port_spacing) {
                    out_index = out_forward_in ? loop_forward_by_distance(points, in_index, port_spacing) :
                                                 loop_back_by_distance(points, in_index, port_spacing);
                    length_io = loop_length_between(points, in_index, out_index);
                    length_oi = loop_length_between(points, out_index, in_index);
                    out_forward_in = length_io <= length_oi;
                    arc_length = std::min(length_io, length_oi);
                }

                if (!first_circle) {
                    const double in1_index = out_forward_in ?
                        loop_forward_by_distance(points, in_index, arc_length - port_spacing) :
                        loop_back_by_distance(points, in_index, arc_length - port_spacing);
                    const double out1_index = out_forward_in ?
                        loop_back_by_distance(points, out_index, arc_length - port_spacing) :
                        loop_forward_by_distance(points, out_index, arc_length - port_spacing);

                    const Point in0 = loop_point_at_index(points, in_index);
                    const Point out0 = loop_point_at_index(points, out_index);
                    const Point in1 = loop_point_at_index(points, in1_index);
                    const Point out1 = loop_point_at_index(points, out1_index);
                    const Point &prev_in = in_branch.back();
                    const Point &prev_out = out_branch.back();

                    const double keep_in_score = std::min(point_distance(in0, prev_in) + point_distance(out1, prev_out),
                                                          point_distance(in0, prev_out) + point_distance(out1, prev_in));
                    const double keep_out_score = std::min(point_distance(in1, prev_in) + point_distance(out0, prev_out),
                                                           point_distance(in1, prev_out) + point_distance(out0, prev_in));
                    if (keep_in_score <= keep_out_score) {
                        if (point_distance(in0, prev_in) + point_distance(out1, prev_out) <
                            point_distance(in0, prev_out) + point_distance(out1, prev_in)) {
                            out_index = out1_index;
                        } else {
                            out_index = in_index;
                            in_index = out1_index;
                            out_forward_in = !out_forward_in;
                        }
                    } else {
                        if (point_distance(in1, prev_in) + point_distance(out0, prev_out) <
                            point_distance(in1, prev_out) + point_distance(out0, prev_in)) {
                            in_index = in1_index;
                        } else {
                            in_index = out_index;
                            out_index = in1_index;
                            out_forward_in = !out_forward_in;
                        }
                    }
                }
            }
        }

        if (!first_circle && segments_intersect(
                loop_point_at_index(points, in_index),
                in_branch.back(),
                loop_point_at_index(points, out_index),
                out_branch.back())) {
            std::swap(in_index, out_index);
            out_forward_in = !out_forward_in;
        }

        append_loop_param(in_branch, points, in_index);
        append_loop_param(out_branch, points, out_index);

        if (loop_index + 1 < loops.size()) {
            if (in_run) {
                if (out_forward_in) {
                    const double out_far = loop_forward_by_distance(points, out_index, port_spacing);
                    append_loop_vertices_backward_until(in_branch, points, size_t(std::floor(in_index)) % n, out_index, out_far);
                    append_loop_param(in_branch, points, out_far);
                } else {
                    const double out_back = loop_back_by_distance(points, out_index, port_spacing);
                    for (size_t guard = 0, i = size_t(std::ceil(in_index)) % n;
                         subset_cycle_param(double(i), out_index, out_back, false, false) && guard <= n + 2;
                         i = (i + 1) % n, ++guard)
                        append_point(in_branch, points[i]);
                    append_loop_param(in_branch, points, out_back);
                }
            } else {
                if (out_forward_in) {
                    const double in_back = loop_back_by_distance(points, in_index, port_spacing);
                    for (size_t guard = 0, i = size_t(std::ceil(out_index)) % n;
                         subset_cycle_param(double(i), in_index, in_back, false, false) && guard <= n + 2;
                         i = (i + 1) % n, ++guard)
                        append_point(out_branch, points[i]);
                    append_loop_param(out_branch, points, in_back);
                } else {
                    const double in_far = loop_forward_by_distance(points, in_index, port_spacing);
                    append_loop_vertices_backward_until(out_branch, points, size_t(std::floor(out_index)) % n, in_index, in_far);
                    append_loop_param(out_branch, points, in_far);
                }
            }

            in_run = !in_run;
            ++loop_index;
            const Points &child_points = loops[loop_index].points;
            in_index = loop_nearest_index(child_points, in_branch.back());
            out_index = loop_nearest_index(child_points, out_branch.back());
        } else {
            struct TerminalCandidate
            {
                PairMetrics metrics;
                double arc_length { 0.0 };
                Points path;
            };

            std::vector<TerminalCandidate> candidates;
            for (const int direction : { 1, -1 }) {
                const Points arc = loop_arc_between_params(points, in_index, out_index, direction);
                if (arc.size() < 2)
                    continue;
                const double arc_length = polyline_length(arc, false);
                if (arc_length + EPS < loop.length * 0.5)
                    continue;
                bool bounded_edges = true;
                for (size_t i = 1; i < arc.size(); ++i)
                    bounded_edges &= point_distance(arc[i - 1], arc[i]) <= spacing * 2.0;
                if (!bounded_edges)
                    continue;

                Points candidate = in_branch;
                for (size_t i = 1; i < arc.size(); ++i)
                    append_point(candidate, arc[i]);
                for (auto it = out_branch.rbegin(); it != out_branch.rend(); ++it)
                    append_point(candidate, *it);
                candidates.push_back({ path_pair_metrics(candidate, spacing), arc_length, std::move(candidate) });
            }

            if (candidates.empty())
                return {};
            const auto better = [](const TerminalCandidate &a, const TerminalCandidate &b) {
                if (a.metrics.crossings != b.metrics.crossings)
                    return a.metrics.crossings < b.metrics.crossings;
                if (a.metrics.close_pairs != b.metrics.close_pairs)
                    return a.metrics.close_pairs < b.metrics.close_pairs;
                if (a.metrics.min_spacing != b.metrics.min_spacing)
                    return a.metrics.min_spacing > b.metrics.min_spacing;
                return a.arc_length > b.arc_length;
            };
            return std::min_element(candidates.begin(), candidates.end(), better)->path;
        }

        first_circle = false;
    }

    for (auto it = out_branch.rbegin(); it != out_branch.rend(); ++it)
        append_point(in_branch, *it);
    return in_branch;
}

PairMetrics path_pair_metrics(
    const Points &path,
    const double spacing,
    const std::vector<double> *segment_widths,
    const double nominal_line_width,
    const int crossing_limit)
{
    AggregatePerfTimer perf_timer(perf_counters.pair_metrics);
    PairMetrics metrics;
    if (path.size() < 4)
        return metrics;

    struct Segment { Point a; Point b; };
    std::vector<Segment> segments;
    segments.reserve(path.size() - 1);
    for (size_t i = 1; i < path.size(); ++i)
        segments.push_back({ path[i - 1], path[i] });

    std::vector<double> prefix { 0.0 };
    prefix.reserve(segments.size() + 1);
    for (const Segment &seg : segments)
        prefix.push_back(prefix.back() + point_distance(seg.a, seg.b));

    // Match the documented strict benchmark's 25% spacing tolerance.  Pairs
    // inside the local four-spacing arc window are part of one continuous
    // turn; only non-local centerlines below 75% of nominal are warnings.
    const double min_allowed = spacing * 0.75;
    const bool audit_bead_widths = segment_widths != nullptr && segment_widths->size() == segments.size() && nominal_line_width > EPS;
    const double max_segment_width = audit_bead_widths ?
        *std::max_element(segment_widths->begin(), segment_widths->end()) : 0.0;
    const double search_margin = std::max(min_allowed, max_segment_width);
    const double bin_size = std::max(search_margin * 1.5, 1.0);
    const double local_skip_distance = spacing * 4.0;
    const double bead_local_skip_distance = spacing * 4.0;
    const double allowed_nominal_penetration =
        std::max(0.0, nominal_line_width - spacing) + nominal_line_width * 0.02;
    const double path_length = prefix.back();
    const bool closed_path = point_distance(path.front(), path.back()) <= spacing * 0.20;
    std::unordered_map<unsigned long long, std::vector<size_t>> bins;
    struct BinBounds
    {
        long long x0;
        long long y0;
        long long x1;
        long long y1;
    };
    std::vector<BinBounds> bin_bounds;
    bin_bounds.reserve(segments.size());
    const auto bin_key = [](const long long bx, const long long by) {
        return (static_cast<unsigned long long>(uint32_t(bx)) << 32) |
               static_cast<unsigned long long>(uint32_t(by));
    };

    for (size_t idx = 0; idx < segments.size(); ++idx) {
        const Segment &seg = segments[idx];
        const double min_x = std::min(seg.a.x(), seg.b.x()) - search_margin;
        const double min_y = std::min(seg.a.y(), seg.b.y()) - search_margin;
        const double max_x = std::max(seg.a.x(), seg.b.x()) + search_margin;
        const double max_y = std::max(seg.a.y(), seg.b.y()) + search_margin;
        const long long bx0 = (long long)std::floor(min_x / bin_size);
        const long long by0 = (long long)std::floor(min_y / bin_size);
        const long long bx1 = (long long)std::floor(max_x / bin_size);
        const long long by1 = (long long)std::floor(max_y / bin_size);
        bin_bounds.push_back({ bx0, by0, bx1, by1 });
        for (long long by = by0; by <= by1; ++by)
            for (long long bx = bx0; bx <= bx1; ++bx)
                bins[bin_key(bx, by)].push_back(idx);
    }

    for (size_t i = 0; i < segments.size(); ++i) {
        const Segment &seg = segments[i];
        const BinBounds &bounds = bin_bounds[i];

        for (long long by = bounds.y0; by <= bounds.y1; ++by)
            for (long long bx = bounds.x0; bx <= bounds.x1; ++bx)
                if (auto it = bins.find(bin_key(bx, by)); it != bins.end()) {
                    for (const size_t j : it->second) {
                        if (j <= i)
                            continue;
                        // Expanded segment bounds commonly occupy several
                        // grid cells. Process their pair only in the lowest
                        // shared cell, avoiding a large hash set of already
                        // visited pairs while preserving each exact audit.
                        const BinBounds &other_bounds = bin_bounds[j];
                        if (bx != std::max(bounds.x0, other_bounds.x0) ||
                            by != std::max(bounds.y0, other_bounds.y0))
                            continue;

                        size_t index_distance = j - i;
                        if (closed_path)
                            index_distance = std::min(index_distance, segments.size() - index_distance);
                        if (index_distance <= 1)
                            continue;

                        double path_gap = std::max(0.0, prefix[j] - prefix[i + 1]);
                        if (closed_path) {
                            const double occupied_span = prefix[j + 1] - prefix[i];
                            path_gap = std::min(path_gap, std::max(0.0, path_length - occupied_span));
                        }

                        const double d = segment_distance(seg.a, seg.b, segments[j].a, segments[j].b);
                        if (d <= 1e-7) {
                            // A zero/sub-resolution bridge vertex may split
                            // one physical turn into segments whose indices
                            // are non-adjacent.  Treat only that tiny connected
                            // arc as adjacency; every other intersection,
                            // including a small bow-tie, is a crossing.
                            if (path_gap <= 1e-9)
                                continue;
                            ++metrics.crossings;
                            metrics.min_spacing = std::min(metrics.min_spacing, d);
                            if (metrics.crossings >= crossing_limit)
                                return metrics;
                            continue;
                        }
                        if (path_gap > local_skip_distance) {
                            metrics.min_spacing = std::min(metrics.min_spacing, d);
                            if (d < min_allowed)
                                ++metrics.close_pairs;
                        }
                        if (audit_bead_widths && path_gap > bead_local_skip_distance) {
                            const double radius_sum =
                                0.5 * ((*segment_widths)[i] + (*segment_widths)[j]);
                            if (d < radius_sum - allowed_nominal_penetration) {
                                ++metrics.bead_overlaps;
                                if (metrics.first_bead_overlap.empty()) {
                                    std::ostringstream detail;
                                    detail << "segments " << i << "/" << j << " distance=" << unscale<double>(d)
                                           << " path_gap=" << unscale<double>(path_gap)
                                           << " widths=" << unscale<double>((*segment_widths)[i]) << "/"
                                           << unscale<double>((*segment_widths)[j])
                                           << " a=" << unscale<double>(segments[i].a.x()) << "," << unscale<double>(segments[i].a.y())
                                           << "->" << unscale<double>(segments[i].b.x()) << "," << unscale<double>(segments[i].b.y())
                                           << " b=" << unscale<double>(segments[j].a.x()) << "," << unscale<double>(segments[j].a.y())
                                           << "->" << unscale<double>(segments[j].b.x()) << "," << unscale<double>(segments[j].b.y());
                                    metrics.first_bead_overlap = detail.str();
                                }
                            }
                        }
                    }
                }
    }

    return metrics;
}

Points replace_long_connectors_with_contour_arcs(
    Points path,
    const std::vector<ContourLoop> &contours,
    const double spacing)
{
    for (size_t repair = 0; repair < 8 && path.size() >= 2; ++repair) {
        const PairMetrics baseline = path_pair_metrics(path, spacing);
        Points best_path;
        PairMetrics best_metrics = baseline;
        bool improved = false;

        for (size_t segment_index = 1; segment_index < path.size(); ++segment_index) {
            const Point &a = path[segment_index - 1];
            const Point &b = path[segment_index];
            if (point_distance(a, b) <= spacing * 3.0)
                continue;

            for (const ContourLoop &loop : contours) {
                if (!loop.closed || loop.points.size() < 4)
                    continue;
                const double a_index = loop_nearest_index(loop.points, a);
                const double b_index = loop_nearest_index(loop.points, b);
                if (point_distance(a, loop_point_at_index(loop.points, a_index)) > spacing * 1.5 ||
                    point_distance(b, loop_point_at_index(loop.points, b_index)) > spacing * 1.5)
                    continue;

                for (const int direction : { 1, -1 }) {
                    const Points arc = loop_arc_between_params(loop.points, a_index, b_index, direction);
                    if (arc.size() < 2)
                        continue;
                    bool bounded_edges = true;
                    for (size_t i = 1; i < arc.size(); ++i)
                        bounded_edges &= point_distance(arc[i - 1], arc[i]) <= spacing * 2.0;
                    if (!bounded_edges)
                        continue;

                    Points candidate;
                    candidate.reserve(path.size() + arc.size() + 2);
                    for (size_t i = 0; i < segment_index; ++i)
                        append_point(candidate, path[i]);
                    append_points(candidate, arc);
                    append_point(candidate, b);
                    for (size_t i = segment_index + 1; i < path.size(); ++i)
                        append_point(candidate, path[i]);

                    const PairMetrics metrics = path_pair_metrics(candidate, spacing);
                    const bool better = metrics.crossings < best_metrics.crossings ||
                        (metrics.crossings == best_metrics.crossings && metrics.close_pairs < best_metrics.close_pairs) ||
                        (metrics.crossings == best_metrics.crossings && metrics.close_pairs == best_metrics.close_pairs &&
                         metrics.min_spacing > best_metrics.min_spacing);
                    if (better) {
                        best_path = std::move(candidate);
                        best_metrics = metrics;
                        improved = true;
                    }
                }
            }
        }

        if (!improved)
            break;
        path = std::move(best_path);
    }
    return path;
}

Points complete_outer_boundary_cycle(const Points &path, const ContourLoop &outer_loop, const double spacing)
{
    if (path.size() < 2 || outer_loop.points.size() < 4)
        return path;

    if (path.front() == path.back())
        return path;

    const double start_index = loop_nearest_index(outer_loop.points, path.front());
    const double end_index = loop_nearest_index(outer_loop.points, path.back());

    struct Candidate
    {
        PairMetrics metrics;
        double arc_length { 0.0 };
        Points path;
    };

    std::vector<Candidate> candidates;
    for (const int direction : { 1, -1 }) {
        const double stop_index = start_index;
        Points arc = loop_arc_between_params(outer_loop.points, end_index, stop_index, direction);
        if (arc.empty())
            continue;

        Points candidate;
        candidate.reserve(path.size() + arc.size() + 1);
        append_points(candidate, path);
        // Keep the actual projections onto the outer contour.  Replacing the
        // arc endpoints with the open path endpoints can turn the first arc
        // edge into a long chord across a concavity (and a real extruded
        // crossing) whenever a branch merge leaves an endpoint slightly off
        // the outer contour.
        append_points(candidate, arc);
        append_point(candidate, path.front());

        if (candidate.front() != candidate.back())
            continue;

        candidates.push_back({ path_pair_metrics(candidate, spacing), polyline_length(arc, false), std::move(candidate) });
    }

    if (candidates.empty())
        return path;

    auto better = [](const Candidate &a, const Candidate &b) {
        if (a.metrics.crossings != b.metrics.crossings)
            return a.metrics.crossings < b.metrics.crossings;
        if (a.metrics.close_pairs != b.metrics.close_pairs)
            return a.metrics.close_pairs < b.metrics.close_pairs;
        if (a.metrics.min_spacing != b.metrics.min_spacing)
            return a.metrics.min_spacing > b.metrics.min_spacing;
        return a.arc_length < b.arc_length;
    };
    return std::min_element(candidates.begin(), candidates.end(), better)->path;
}

Points build_best_single_chain_path(
    const std::vector<ContourLoop> &ordered,
    const Point &start_anchor,
    const double spacing,
    const bool preserve_medial_pockets = true,
    const double seam_port_factor = 1.17)
{
    const std::vector<double> exit_fractions { 0.50, 0.67, 0.33, 0.25, 0.75, 0.10, 0.90 };
    Points best_path;
    PairMetrics best_metrics;
    bool have_best = false;

    auto better_than_best = [&](const PairMetrics &metrics) {
        return !have_best || metrics.crossings < best_metrics.crossings ||
            (metrics.crossings == best_metrics.crossings && metrics.close_pairs < best_metrics.close_pairs) ||
            (metrics.crossings == best_metrics.crossings && metrics.close_pairs == best_metrics.close_pairs &&
             metrics.min_spacing > best_metrics.min_spacing);
    };

    // Every retained offset contour represents printable material.  Older
    // code progressively removed inner contours until a clean-looking path
    // appeared, which could return only the outer wall while silently losing
    // almost 90% of a solid layer.  Try routing variants, but never trade away
    // required contours; the mandatory final coverage gate decides whether the
    // complete route is usable.
    for (const double fraction : exit_fractions) {
        const Point exit_anchor = point_at_closed_fraction(ordered.front().points, fraction);
        Points candidate = build_single_minimum_connected_fermat(
            ordered, start_anchor, spacing, exit_anchor, preserve_medial_pockets, seam_port_factor);
        if (candidate.size() < 2)
            continue;
        const PairMetrics metrics = path_pair_metrics(candidate, spacing);
        if (better_than_best(metrics)) {
            best_path = std::move(candidate);
            best_metrics = metrics;
            have_best = true;
        }
        if (metrics.crossings == 0 && metrics.close_pairs == 0)
            break;
    }

    return best_path;
}

struct OpenProjection
{
    double index { 0.0 };
    Point point;
    double distance { std::numeric_limits<double>::infinity() };
};

OpenProjection project_open_polyline(const Points &points, const Point &p)
{
    if (points.size() < 2) {
        const Point q = points.empty() ? Point(0, 0) : points.front();
        return { 0.0, q, point_distance(p, q) };
    }

    OpenProjection best;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        const Point &a = points[i];
        const Point &b = points[i + 1];
        const Vec2d ab = (b - a).cast<double>();
        const double denom = ab.squaredNorm();
        const double t = denom <= EPS ? 0.0 : std::clamp((p - a).cast<double>().dot(ab) / denom, 0.0, 1.0);
        const Point q = lerp_point(a, b, t);
        const double d = point_distance(p, q);
        if (d < best.distance) {
            best.index = double(i) + t;
            best.point = q;
            best.distance = d;
        }
    }
    return best;
}

struct CutOpenPolyline
{
    Points before;
    Points after;
    Point cut;
};

CutOpenPolyline cut_open_polyline(const Points &points, double index)
{
    if (points.empty())
        return {};
    if (points.size() == 1)
        return { points, points, points.front() };

    index = std::clamp(index, 0.0, double(points.size() - 1));
    const size_t base = std::min<size_t>(size_t(std::floor(index)), points.size() - 2);
    const double alpha = index - double(base);
    const Point cut = lerp_point(points[base], points[base + 1], alpha);

    Points before;
    before.reserve(base + 2);
    for (size_t i = 0; i <= base; ++i)
        append_point(before, points[i]);
    append_point(before, cut);

    Points after;
    after.reserve(points.size() - base);
    append_point(after, cut);
    for (size_t i = base + 1; i < points.size(); ++i)
        append_point(after, points[i]);

    return { std::move(before), std::move(after), cut };
}

OpenProjection project_open_polyline_near_arc(
    const Points &points,
    const Point &p,
    const double center_index,
    const double max_arc_length)
{
    if (points.size() < 2 || !std::isfinite(max_arc_length))
        return project_open_polyline(points, p);

    std::vector<double> prefix { 0.0 };
    prefix.reserve(points.size());
    for (size_t i = 1; i < points.size(); ++i)
        prefix.emplace_back(prefix.back() + point_distance(points[i - 1], points[i]));
    const size_t center_segment = std::min<size_t>(size_t(std::floor(center_index)), points.size() - 2);
    const double center_alpha = std::clamp(center_index - double(center_segment), 0.0, 1.0);
    const double center_distance = prefix[center_segment] +
        (prefix[center_segment + 1] - prefix[center_segment]) * center_alpha;

    OpenProjection best;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        const double segment_distance_from_center = center_distance < prefix[i] ?
            prefix[i] - center_distance :
            (center_distance > prefix[i + 1] ? center_distance - prefix[i + 1] : 0.0);
        if (segment_distance_from_center > max_arc_length)
            continue;
        const Point &a = points[i];
        const Point &b = points[i + 1];
        const Vec2d ab = (b - a).cast<double>();
        const double denom = ab.squaredNorm();
        const double t = denom <= EPS ? 0.0 : std::clamp((p - a).cast<double>().dot(ab) / denom, 0.0, 1.0);
        const Point q = lerp_point(a, b, t);
        const double distance = point_distance(p, q);
        if (distance < best.distance) {
            best.index = double(i) + t;
            best.point = q;
            best.distance = distance;
        }
    }
    return std::isfinite(best.distance) ? best : project_open_polyline(points, p);
}

Points merge_child_spiral(
    const Points &parent_path,
    Points child_path,
    const double max_replaced_parent_arc = std::numeric_limits<double>::infinity())
{
    if (parent_path.empty())
        return child_path;
    if (child_path.empty())
        return parent_path;

    OpenProjection start_projection = project_open_polyline(parent_path, child_path.front());
    OpenProjection end_projection = project_open_polyline_near_arc(
        parent_path, child_path.back(), start_projection.index, max_replaced_parent_arc);
    if (start_projection.index > end_projection.index) {
        std::reverse(child_path.begin(), child_path.end());
        start_projection = project_open_polyline(parent_path, child_path.front());
        end_projection = project_open_polyline_near_arc(
            parent_path, child_path.back(), start_projection.index, max_replaced_parent_arc);
    }

    CutOpenPolyline parent_start = cut_open_polyline(parent_path, start_projection.index);
    CutOpenPolyline parent_end = cut_open_polyline(parent_path, end_projection.index);

    Points merged;
    merged.reserve(parent_path.size() + child_path.size() + 4);
    append_points(merged, parent_start.before);
    append_points(merged, child_path);
    append_points(merged, parent_end.after);
    return merged;
}

Points uncross_closed_tour(
    const ExPolygons &printable_area,
    Points path,
    double line_width,
    double spacing);
int count_turnback_violations(const Points &path, double line_width, std::string *first_detail);

Points build_branch_connected_fermat(
    const std::vector<ContourLoop> &contours,
    const ExPolygons &printable_area,
    const Point &start_anchor,
    const double spacing)
{
    const auto grouped = loops_by_level(contours);
    if (grouped.empty())
        return {};

    bool multi_loop = false;
    for (const auto &[level, loops] : grouped)
        multi_loop |= loops.size() > 1;
    if (!multi_loop)
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    // Follow one root contour inward through every offset level. All other
    // contours become explicit branches. This also handles multiply connected
    // input, where the outer boundary and all hole boundaries already coexist
    // at level zero and there is no pre-split trunk.
    std::vector<ContourLoop> trunk;
    std::map<size_t, std::vector<ContourLoop>> branch_levels;
    const ContourLoop *previous_root = nullptr;
    for (const auto &[level, loops] : grouped) {
        if (loops.empty())
            continue;
        size_t root_index = 0;
        if (previous_root == nullptr) {
            for (size_t i = 1; i < loops.size(); ++i)
                if (loops[i].area > loops[root_index].area)
                    root_index = i;
        } else {
            double best_score = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < loops.size(); ++i) {
                const double centroid_distance = point_distance(previous_root->centroid, loops[i].centroid);
                const double growth = std::max(0.0, loops[i].area - previous_root->area) /
                                      std::max(previous_root->area, 1.0);
                const double score = centroid_distance + growth * spacing * 20.0;
                if (score < best_score) {
                    best_score = score;
                    root_index = i;
                }
            }
        }
        trunk.emplace_back(loops[root_index]);
        previous_root = &trunk.back();
        for (size_t i = 0; i < loops.size(); ++i)
            if (i != root_index)
                branch_levels[level].emplace_back(loops[i]);
    }

    if (trunk.empty())
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    std::vector<std::vector<ContourLoop>> chains;
    for (const auto &[level, loops] : branch_levels) {
        std::vector<bool> used(loops.size(), false);
        for (std::vector<ContourLoop> &chain : chains) {
            double best_distance2 = std::numeric_limits<double>::infinity();
            size_t best_index = size_t(-1);
            for (size_t i = 0; i < loops.size(); ++i) {
                if (used[i])
                    continue;
                const double d = point_distance2(chain.back().centroid, loops[i].centroid);
                if (d < best_distance2) {
                    best_distance2 = d;
                    best_index = i;
                }
            }
            const double max_link_distance = std::max(
                spacing * 4.0,
                std::sqrt(std::max(chain.back().area, 1.0) / PI) * 0.5);
            if (best_index != size_t(-1) && best_distance2 <= max_link_distance * max_link_distance) {
                used[best_index] = true;
                chain.emplace_back(loops[best_index]);
            }
        }
        // Never discard a contour merely because the offset topology gained a
        // branch or a previous branch terminated. Every unmatched loop starts
        // a new subtree that will be merged into the root traversal below.
        for (size_t i = 0; i < loops.size(); ++i)
            if (!used[i])
                chains.push_back({ loops[i] });
    }

    const Point trunk_exit = point_at_closed_fraction(trunk.front().points, 0.5);
    Points parent_path = build_single_minimum_connected_fermat(trunk, start_anchor, spacing, trunk_exit);
    if (parent_path.size() < 2)
        return parent_path;

    for (const std::vector<ContourLoop> &chain : chains) {
        if (chain.empty() || chain.front().points.empty())
            continue;

        std::vector<ContourLoop> ordered_chain = chain;
        if (ordered_chain.back().area > ordered_chain.front().area)
            std::reverse(ordered_chain.begin(), ordered_chain.end());

        const size_t stride = std::max<size_t>(1, ordered_chain.front().points.size() / 48);
        Point anchor = ordered_chain.front().points.front();
        size_t anchor_vertex = 0;
        double best_distance = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < ordered_chain.front().points.size(); i += stride) {
            const OpenProjection projection = project_open_polyline(parent_path, ordered_chain.front().points[i]);
            if (projection.distance < best_distance) {
                best_distance = projection.distance;
                anchor = ordered_chain.front().points[i];
                anchor_vertex = i;
            }
        }

        // A branch must enter and leave through a narrow slot at its neck.
        // Using the opposite side of the child contour as the exit creates a
        // long extruded chord back to the trunk and can cross an already
        // printed arm.  Try nearby ports on both sides of the closest neck
        // point and keep the safest complete merge.
        Points best_merged;
        PairMetrics best_metrics;
        int best_turnbacks = std::numeric_limits<int>::max();
        bool have_merged = false;
        for (const int direction : { 1, -1 }) {
            for (const double port_distance : {
                     spacing * 0.5,
                     spacing * 0.75,
                     spacing,
                     spacing * 1.25,
                     spacing * 1.5,
                     spacing * 2.0,
                 }) {
                const double exit_index = direction > 0 ?
                    loop_forward_by_distance(ordered_chain.front().points, double(anchor_vertex), port_distance) :
                    loop_back_by_distance(ordered_chain.front().points, double(anchor_vertex), port_distance);
                const Point child_exit = loop_point_at_index(ordered_chain.front().points, exit_index);
                Points child_path = build_single_minimum_connected_fermat(ordered_chain, anchor, spacing, child_exit);
                if (child_path.size() < 2)
                    continue;

                Points merged = merge_child_spiral(parent_path, std::move(child_path), spacing * 4.0);
                merged = uncross_closed_tour(printable_area, std::move(merged), spacing, spacing);
                const PairMetrics metrics = path_pair_metrics(merged, spacing);
                const int turnbacks = count_turnback_violations(merged, spacing, nullptr);
                const bool better = !have_merged || metrics.crossings < best_metrics.crossings ||
                    (metrics.crossings == best_metrics.crossings && turnbacks < best_turnbacks) ||
                    (metrics.crossings == best_metrics.crossings && turnbacks == best_turnbacks &&
                     metrics.close_pairs < best_metrics.close_pairs) ||
                    (metrics.crossings == best_metrics.crossings && turnbacks == best_turnbacks &&
                     metrics.close_pairs == best_metrics.close_pairs &&
                     metrics.min_spacing > best_metrics.min_spacing);
                if (better) {
                    best_merged = std::move(merged);
                    best_metrics = metrics;
                    best_turnbacks = turnbacks;
                    have_merged = true;
                }
                if (metrics.crossings == 0 && turnbacks == 0 && metrics.close_pairs == 0)
                    break;
            }
            if (have_merged && best_metrics.crossings == 0 && best_turnbacks == 0 && best_metrics.close_pairs == 0)
                break;
        }

        if (have_merged)
            parent_path = std::move(best_merged);
    }

    return parent_path;
}

std::vector<Polygon> collect_holes(const ExPolygons &printable_area)
{
    std::vector<Polygon> holes;
    for (const ExPolygon &expoly : union_ex(printable_area))
        holes.insert(holes.end(), expoly.holes.begin(), expoly.holes.end());
    return holes;
}

ContourLoop build_capsule_loop(const Point &left_center, const Point &right_center, const double radius, const size_t level_index)
{
    static constexpr size_t CAPSULE_SEGMENTS = 72;
    const Vec2d a = left_center.cast<double>();
    const Vec2d b = right_center.cast<double>();
    Vec2d axis = b - a;
    const double axis_length = axis.norm();
    Points points;

    if (axis_length <= EPS) {
        points.reserve(CAPSULE_SEGMENTS);
        for (size_t i = 0; i < CAPSULE_SEGMENTS; ++i) {
            const double angle = 2.0 * PI * double(i) / double(CAPSULE_SEGMENTS);
            points.emplace_back(a.x() + std::cos(angle) * radius, a.y() + std::sin(angle) * radius);
        }
    } else {
        const Vec2d u = axis / axis_length;
        const Vec2d v(-u.y(), u.x());
        const size_t half_segments = CAPSULE_SEGMENTS / 2;
        points.reserve(CAPSULE_SEGMENTS + 4);

        points.emplace_back(a.x() + v.x() * radius, a.y() + v.y() * radius);
        points.emplace_back(b.x() + v.x() * radius, b.y() + v.y() * radius);
        for (size_t i = 1; i <= half_segments; ++i) {
            const double angle = PI * 0.5 - PI * double(i) / double(half_segments);
            const Vec2d p = b + u * (std::cos(angle) * radius) + v * (std::sin(angle) * radius);
            points.emplace_back(p.x(), p.y());
        }

        points.emplace_back(a.x() - v.x() * radius, a.y() - v.y() * radius);
        for (size_t i = 1; i <= half_segments; ++i) {
            const double angle = -PI * 0.5 - PI * double(i) / double(half_segments);
            const Vec2d p = a + u * (std::cos(angle) * radius) + v * (std::sin(angle) * radius);
            points.emplace_back(p.x(), p.y());
        }
    }

    if (signed_area(points) < 0.0)
        std::reverse(points.begin(), points.end());

    ContourLoop loop;
    loop.level_index = level_index;
    loop.points = std::move(points);
    loop.area = std::abs(signed_area(loop.points));
    loop.length = polyline_length(loop.points, true);
    loop.centroid = polygon_centroid_or_first(loop.points);
    return loop;
}

Points build_merged_hole_connected_fermat(
    const std::vector<ContourLoop> &contours,
    const ExPolygons &printable_area,
    const Flow &flow,
    const Point &start_anchor)
{
    const double spacing = double(flow.scaled_spacing());
    const double line_width = double(flow.scaled_width());
    const std::vector<Polygon> holes = collect_holes(printable_area);
    if (holes.size() < 2)
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    const auto grouped = loops_by_level(contours);
    std::vector<ContourLoop> outer;
    const double min_outer_area = spacing * spacing * 20.0;
    for (const auto &[level, loops] : grouped) {
        if (loops.empty())
            continue;
        const ContourLoop *largest = &loops.front();
        for (const ContourLoop &loop : loops)
            if (loop.area > largest->area)
                largest = &loop;
        if (largest->area > min_outer_area)
            outer.emplace_back(*largest);
    }

    std::vector<ContourLoop> stable_outer;
    double previous_area = -1.0;
    for (const ContourLoop &loop : outer) {
        if (previous_area > 0.0 && loop.area < previous_area * 0.2)
            break;
        stable_outer.emplace_back(loop);
        previous_area = loop.area;
    }
    outer = std::move(stable_outer);
    if (outer.empty())
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    std::vector<Point> centers;
    centers.reserve(holes.size());
    for (const Polygon &hole : holes)
        centers.emplace_back(polygon_centroid_or_first(hole.points));

    size_t left_index = 0;
    size_t right_index = 1;
    double farthest = -1.0;
    for (size_t i = 0; i < centers.size(); ++i) {
        for (size_t j = i + 1; j < centers.size(); ++j) {
            const double d = point_distance2(centers[i], centers[j]);
            if (d > farthest) {
                farthest = d;
                left_index = i;
                right_index = j;
            }
        }
    }

    Point left_center = centers[left_index];
    Point right_center = centers[right_index];
    if (right_center.x() < left_center.x())
        std::swap(left_center, right_center);

    double base_radius = 0.0;
    for (const Polygon &hole : holes)
        for (const Point &point : hole.points)
            base_radius = std::max(base_radius, point_segment_distance(point, left_center, right_center));
    base_radius += line_width * 0.5;

    const Vec2d axis = (right_center - left_center).cast<double>();
    const double axis_length = axis.norm();
    if (axis_length <= EPS)
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    const Vec2d u = axis / axis_length;
    const Vec2d v(-u.y(), u.x());
    const Vec2d mid = (left_center.cast<double>() + right_center.cast<double>()) * 0.5;
    double min_perp = std::numeric_limits<double>::infinity();
    double max_perp = -std::numeric_limits<double>::infinity();
    for (const ExPolygon &expoly : union_ex(printable_area)) {
        for (const Point &point : expoly.contour.points) {
            const double projection = (point.cast<double>() - mid).dot(v);
            min_perp = std::min(min_perp, projection);
            max_perp = std::max(max_perp, projection);
        }
    }

    const double max_radius = std::min(max_perp, -min_perp) - line_width;
    if (max_radius <= base_radius)
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    std::vector<ContourLoop> inner;
    for (size_t level = 0; level < outer.size(); ++level) {
        const double radius = base_radius + double(level) * spacing;
        if (radius > max_radius)
            break;
        inner.emplace_back(build_capsule_loop(left_center, right_center, radius, level));
    }

    if (inner.empty())
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    std::vector<ContourLoop> ordered;
    const size_t paired_count = std::min(outer.size(), inner.size());
    ordered.reserve(paired_count * 2);
    ordered.insert(ordered.end(), outer.begin(), outer.begin() + paired_count);
    for (auto it = inner.rbegin(); it != inner.rend(); ++it)
        ordered.emplace_back(*it);

    std::vector<Point> candidate_anchors { start_anchor };
    for (const double fraction : { 0.04, 0.06, 0.10, 0.25, 0.50, 0.75 })
        candidate_anchors.emplace_back(point_at_closed_fraction(ordered.front().points, fraction));

    Points best_path;
    PairMetrics best_metrics;
    bool have_best = false;
    for (const Point &anchor : candidate_anchors) {
        const Point exit_anchor = point_at_closed_fraction(ordered.front().points, 0.5);
        Points candidate = build_single_minimum_connected_fermat(ordered, anchor, spacing, exit_anchor);
        PairMetrics metrics = path_pair_metrics(candidate, spacing);
        if (!have_best || metrics.crossings < best_metrics.crossings ||
            (metrics.crossings == best_metrics.crossings && metrics.close_pairs < best_metrics.close_pairs) ||
            (metrics.crossings == best_metrics.crossings && metrics.close_pairs == best_metrics.close_pairs &&
             metrics.min_spacing > best_metrics.min_spacing)) {
            best_path = std::move(candidate);
            best_metrics = metrics;
            have_best = true;
        }
        if (metrics.crossings == 0 && metrics.close_pairs == 0)
            break;
    }

    return best_path;
}

std::vector<double> open_path_prefix_lengths(const Points &path)
{
    std::vector<double> prefix { 0.0 };
    prefix.reserve(path.size());
    for (size_t i = 1; i < path.size(); ++i)
        prefix.push_back(prefix.back() + point_distance(path[i - 1], path[i]));
    return prefix;
}

Point open_path_point_at_index(const Points &path, double index)
{
    if (path.empty())
        return Point(0, 0);
    if (path.size() == 1)
        return path.front();

    index = std::clamp(index, 0.0, double(path.size() - 1));
    const size_t base = std::min<size_t>(size_t(std::floor(index)), path.size() - 2);
    return lerp_point(path[base], path[base + 1], index - double(base));
}

double open_path_length_at_index(const Points &path, const std::vector<double> &prefix, double index)
{
    if (path.size() < 2)
        return 0.0;

    index = std::clamp(index, 0.0, double(path.size() - 1));
    const size_t base = std::min<size_t>(size_t(std::floor(index)), path.size() - 2);
    return prefix[base] + (index - double(base)) * point_distance(path[base], path[base + 1]);
}

double open_path_index_at_length(const Points &path, const std::vector<double> &prefix, double target)
{
    if (path.size() < 2)
        return 0.0;

    target = std::clamp(target, 0.0, prefix.back());
    const auto it = std::lower_bound(prefix.begin(), prefix.end(), target);
    const size_t upper = size_t(std::distance(prefix.begin(), it));
    const size_t idx = upper == 0 ? 0 : std::min(upper - 1, path.size() - 2);
    const double segment_length = std::max(point_distance(path[idx], path[idx + 1]), EPS);
    return double(idx) + (target - prefix[idx]) / segment_length;
}

bool printable_area_contains(const ExPolygons &printable_area, const Point &point)
{
    for (const ExPolygon &expoly : printable_area)
        if (expoly.contains(point))
            return true;
    return false;
}

int count_containment_violations(const ExPolygons &printable_area, const Points &path, const double /* spacing */)
{
    AggregatePerfTimer perf_timer(perf_counters.containment);
    if (path.size() < 2)
        return 0;

    const Polylines outside = diff_pl(Polylines { Polyline(path) }, printable_area);
    return int(std::count_if(outside.begin(), outside.end(), [](const Polyline &polyline) {
        return polyline.length() > SCALED_EPSILON;
    }));
}

double loop_uncovered_fraction(const ContourLoop &loop, const Points &path, const double spacing)
{
    if (path.empty())
        return 1.0;

    const size_t stride = std::max<size_t>(1, loop.points.size() / 32);
    size_t samples = 0;
    size_t uncovered = 0;
    for (size_t i = 0; i < loop.points.size(); i += stride) {
        ++samples;
        if (project_open_polyline(path, loop.points[i]).distance > spacing * 0.55)
            ++uncovered;
    }

    return samples == 0 ? 0.0 : double(uncovered) / double(samples);
}

std::vector<std::vector<ContourLoop>> uncovered_pocket_groups(
    const std::vector<ContourLoop> &contours,
    const Points &path,
    const double spacing)
{
    std::vector<ContourLoop> uncovered;
    const double max_pocket_area = spacing * spacing * 90.0;
    for (const ContourLoop &loop : contours) {
        if (loop.level_index > 0 && loop.area < max_pocket_area && loop_uncovered_fraction(loop, path, spacing) > 0.65)
            uncovered.emplace_back(loop);
    }

    std::sort(uncovered.begin(), uncovered.end(), [](const ContourLoop &a, const ContourLoop &b) { return a.area > b.area; });

    std::vector<std::vector<ContourLoop>> groups;
    for (const ContourLoop &loop : uncovered) {
        bool inserted = false;
        for (std::vector<ContourLoop> &group : groups) {
            if (point_distance(loop.centroid, group.front().centroid) < spacing * 7.0) {
                group.emplace_back(loop);
                inserted = true;
                break;
            }
        }
        if (!inserted)
            groups.push_back({ loop });
    }
    return groups;
}

std::vector<double> unique_port_candidates(std::vector<std::pair<double, double>> candidates, const size_t limit)
{
    std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    std::vector<double> selected;
    selected.reserve(limit);
    for (const auto &[distance, index] : candidates) {
        bool unique = true;
        for (const double existing : selected) {
            if (std::abs(index - existing) <= 2.0) {
                unique = false;
                break;
            }
        }
        if (!unique)
            continue;
        selected.emplace_back(index);
        if (selected.size() >= limit)
            break;
    }
    return selected;
}

Points try_insert_pocket_group(
    const ExPolygons &printable_area,
    const Points &parent_path,
    const std::vector<ContourLoop> &group,
    const double spacing)
{
    if (parent_path.size() < 4 || group.empty())
        return {};

    const std::vector<double> prefix = open_path_prefix_lengths(parent_path);
    Vec2d centroid = Vec2d::Zero();
    for (const ContourLoop &loop : group)
        centroid += loop.centroid.cast<double>();
    centroid /= double(group.size());

    std::vector<std::pair<double, double>> port_candidates;
    OpenProjection projection = project_open_polyline(parent_path, Point(centroid.x(), centroid.y()));
    port_candidates.emplace_back(projection.distance, projection.index);
    for (const ContourLoop &loop : group) {
        const size_t stride = std::max<size_t>(1, loop.points.size() / 16);
        for (size_t i = 0; i < loop.points.size(); i += stride) {
            projection = project_open_polyline(parent_path, loop.points[i]);
            port_candidates.emplace_back(projection.distance, projection.index);
        }
    }

    const std::vector<ContourLoop> ordered = ordered_loops_for_spiral(group);
    const std::vector<double> half_widths {
        spacing * 0.60,
        spacing * 1.00,
        spacing * 1.35,
        spacing * 2.25,
        spacing * 3.25,
    };
    for (double center_index : unique_port_candidates(std::move(port_candidates), 4)) {
        center_index = std::clamp(center_index, 3.0, double(parent_path.size() - 4));
        const double center_length = open_path_length_at_index(parent_path, prefix, center_index);
        for (const double half_width : half_widths) {
            const double left_index = open_path_index_at_length(parent_path, prefix, center_length - half_width);
            const double right_index = open_path_index_at_length(parent_path, prefix, center_length + half_width);
            const Point left_port = open_path_point_at_index(parent_path, left_index);
            const Point right_port = open_path_point_at_index(parent_path, right_index);
            if (point_distance(left_port, right_port) < spacing * 0.5)
                continue;

            for (const bool reverse_anchors : { false, true }) {
                const Point &start_anchor = reverse_anchors ? right_port : left_port;
                const Point &exit_anchor = reverse_anchors ? left_port : right_port;
                Points child_path = build_single_minimum_connected_fermat(ordered, start_anchor, spacing, exit_anchor, true);

                for (const bool reverse_child : { false, true }) {
                    Points candidate_child = child_path;
                    if (reverse_child)
                        std::reverse(candidate_child.begin(), candidate_child.end());

                    CutOpenPolyline parent_start = cut_open_polyline(parent_path, left_index);
                    CutOpenPolyline parent_end = cut_open_polyline(parent_path, right_index);
                    Points candidate;
                    candidate.reserve(parent_path.size() + candidate_child.size() + 4);
                    append_points(candidate, parent_start.before);
                    append_points(candidate, candidate_child);
                    append_points(candidate, parent_end.after);

                    const PairMetrics metrics = path_pair_metrics(candidate, spacing);
                    if (metrics.crossings == 0 && metrics.close_pairs <= 3 &&
                        count_containment_violations(printable_area, candidate, spacing) == 0)
                        return candidate;
                }
            }
        }
    }

    return {};
}

Points insert_uncovered_pocket_spirals(
    const ExPolygons &printable_area,
    const std::vector<ContourLoop> &contours,
    Points path,
    const double spacing)
{
    for (const std::vector<ContourLoop> &group : uncovered_pocket_groups(contours, path, spacing)) {
        Points candidate = try_insert_pocket_group(printable_area, path, group, spacing);
        if (!candidate.empty())
            path = std::move(candidate);
    }
    return path;
}

double contour_group_area(const std::vector<ContourLoop> &group)
{
    double area = 0.0;
    for (const ContourLoop &loop : group)
        area += loop.area;
    return area;
}

std::vector<std::vector<ContourLoop>> residual_gap_groups(const std::vector<ContourLoop> &contours)
{
    if (contours.empty())
        return {};

    size_t first_level = std::numeric_limits<size_t>::max();
    for (const ContourLoop &loop : contours)
        first_level = std::min(first_level, loop.level_index);

    std::vector<ContourLoop> seeds;
    for (const ContourLoop &loop : contours)
        if (loop.level_index == first_level)
            seeds.emplace_back(loop);

    if (seeds.size() <= 1)
        return { contours };

    std::vector<std::vector<ContourLoop>> groups(seeds.size());
    for (const ContourLoop &loop : contours) {
        size_t best = 0;
        double best_distance = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < seeds.size(); ++i) {
            const double d = point_distance2(loop.centroid, seeds[i].centroid);
            if (d < best_distance) {
                best_distance = d;
                best = i;
            }
        }
        groups[best].emplace_back(loop);
    }

    groups.erase(std::remove_if(groups.begin(), groups.end(), [](const std::vector<ContourLoop> &group) { return group.empty(); }), groups.end());
    return groups;
}

Points try_merge_gap_group(
    const ExPolygons &printable_area,
    const Points &parent_path,
    const std::vector<ContourLoop> &group,
    const double spacing)
{
    if (parent_path.size() < 4 || group.empty())
        return {};

    const std::vector<ContourLoop> ordered = ordered_loops_for_spiral(group);
    if (ordered.empty())
        return {};

    const ContourLoop &first_loop = ordered.front();
    std::vector<std::pair<double, Point>> nearest_samples;
    const size_t stride = std::max<size_t>(1, first_loop.points.size() / 24);
    for (size_t i = 0; i < first_loop.points.size(); i += stride) {
        const OpenProjection projection = project_open_polyline(parent_path, first_loop.points[i]);
        nearest_samples.emplace_back(projection.distance, first_loop.points[i]);
    }
    std::sort(nearest_samples.begin(), nearest_samples.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    std::vector<Point> anchor_candidates;
    for (size_t i = 0; i < std::min<size_t>(4, nearest_samples.size()); ++i)
        anchor_candidates.emplace_back(nearest_samples[i].second);
    for (const double fraction : { 0.0, 0.25, 0.5, 0.75 })
        anchor_candidates.emplace_back(point_at_closed_fraction(first_loop.points, fraction));

    for (const Point &start_anchor : anchor_candidates) {
        for (const double exit_fraction : { 0.5, 0.25, 0.75 }) {
            const Point exit_anchor = point_at_closed_fraction(first_loop.points, exit_fraction);
            Points child_path = build_single_minimum_connected_fermat(ordered, start_anchor, spacing, exit_anchor, true);
            if (child_path.size() < 2)
                continue;

            for (const bool reverse_child : { false, true }) {
                Points candidate_child = child_path;
                if (reverse_child)
                    std::reverse(candidate_child.begin(), candidate_child.end());
                Points candidate = merge_child_spiral(parent_path, std::move(candidate_child), spacing * 4.0);
                const PairMetrics metrics = path_pair_metrics(candidate, spacing);
                if (metrics.crossings == 0 && metrics.close_pairs <= 3 &&
                    count_containment_violations(printable_area, candidate, spacing) == 0)
                    return candidate;
            }
        }
    }

    return {};
}

ExPolygons residual_gap_area(const ExPolygons &printable_area, const Points &path, const Flow &flow, const double spacing)
{
    if (path.size() < 2)
        return {};

    // Residual discovery must use the same nominal bead radius as the final
    // swept-area validator.  Inflating this audit footprint hid the medial
    // band where offset fronts meet (notably on annuli), so no repair contour
    // was generated even though roughly ten percent of the layer remained
    // physically uncovered.
    const double stroke_radius = double(flow.scaled_width()) * 0.5;
    Polygons covered = offset(Polyline(path), float(stroke_radius), ClipperLib::jtRound, SCALED_RESOLUTION, ClipperLib::etOpenRound);
    if (covered.empty())
        return {};

    return diff_ex(printable_area, covered);
}

Points insert_residual_gap_spirals(
    const ExPolygons &printable_area,
    const Flow &flow,
    Points path,
    const double spacing,
    std::vector<ContourLoop> *width_contours)
{
    const ExPolygons residual_area = residual_gap_area(printable_area, path, flow, spacing);
    if (residual_area.empty())
        return path;

    const double residual_line_width = scale_(double(flow.nozzle_diameter()) * 0.85);
    const double line_height = scale_(double(flow.height()));
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double nominal_cross_section = double(flow.scaled_width()) - rounded_corner_loss;
    const double residual_cross_section = residual_line_width - rounded_corner_loss;
    if (!(residual_cross_section > EPS) || !(nominal_cross_section > EPS))
        return path;
    const double residual_spacing = residual_cross_section;
    std::vector<ContourLoop> residual_contours =
        generate_offset_contours(residual_area, residual_line_width, residual_spacing, 0);
    residual_contours.erase(
        std::remove_if(
            residual_contours.begin(),
            residual_contours.end(),
            [&](const ContourLoop &loop) { return loop_uncovered_fraction(loop, path, spacing) <= 0.65; }),
        residual_contours.end());
    if (residual_contours.empty())
        return path;

    if (width_contours != nullptr) {
        for (ContourLoop &loop : residual_contours) {
            loop.line_width = residual_line_width;
            loop.flow_multiplier = residual_cross_section / nominal_cross_section;
        }
        width_contours->insert(width_contours->end(), residual_contours.begin(), residual_contours.end());
    }

    std::vector<std::vector<ContourLoop>> groups = residual_gap_groups(residual_contours);
    std::sort(groups.begin(), groups.end(), [](const auto &a, const auto &b) { return contour_group_area(a) > contour_group_area(b); });

    const size_t max_groups = std::min<size_t>(8, groups.size());
    for (size_t i = 0; i < max_groups; ++i) {
        Points candidate = try_insert_pocket_group(printable_area, path, groups[i], spacing);
        if (candidate.empty())
            candidate = try_merge_gap_group(printable_area, path, groups[i], spacing);
        if (!candidate.empty())
            path = std::move(candidate);
    }

    return path;
}

std::vector<float> segment_extrusion_multipliers_for_path(
    const Points &path,
    const std::vector<ContourLoop> &contours,
    const double line_width,
    const double spacing)
{
    std::vector<float> multipliers(path.size() > 1 ? path.size() - 1 : 0, 1.0f);
    if (multipliers.empty())
        return multipliers;

    std::vector<const ContourLoop*> adaptive;
    for (const ContourLoop &loop : contours)
        if (loop.line_width > EPS && std::abs(loop.flow_multiplier - 1.0) > 1e-3 &&
            std::abs(loop.line_width - line_width) > line_width * 0.01)
            adaptive.emplace_back(&loop);
    if (adaptive.empty())
        return multipliers;

    const double tolerance = std::max(spacing * 0.32, line_width * 0.32);
    for (size_t i = 1; i < path.size(); ++i) {
        const Point &a = path[i - 1];
        const Point &b = path[i];
        const Point mid = lerp_point(a, b, 0.5);
        double best_distance = std::numeric_limits<double>::infinity();
        double flow_multiplier = 1.0;
        for (const ContourLoop *loop : adaptive) {
            const double d = std::max({
                contour_distance(a, *loop),
                contour_distance(mid, *loop),
                contour_distance(b, *loop),
            });
            if (d < best_distance && d <= tolerance) {
                best_distance = d;
                flow_multiplier = loop->flow_multiplier;
            }
        }
        multipliers[i - 1] = float(flow_multiplier);
    }
    return multipliers;
}

std::vector<double> nonlocal_segment_clearances(
    const Points &path,
    const double max_physical_width,
    const double spacing)
{
    const size_t segment_count = path.size() > 1 ? path.size() - 1 : 0;
    std::vector<double> clearances(segment_count, std::numeric_limits<double>::infinity());
    if (segment_count < 3)
        return clearances;

    struct Segment { Point a; Point b; };
    std::vector<Segment> segments;
    segments.reserve(segment_count);
    std::vector<double> prefix { 0.0 };
    prefix.reserve(segment_count + 1);
    for (size_t i = 1; i < path.size(); ++i) {
        segments.push_back({ path[i - 1], path[i] });
        prefix.push_back(prefix.back() + point_distance(path[i - 1], path[i]));
    }

    const double max_width = max_physical_width;
    const double bin_size = std::max(max_width * 1.5, 1.0);
    const double local_skip_distance = spacing * 4.0;
    const double path_length = prefix.back();
    const bool closed_path = path.front() == path.back();
    const auto bin_key = [](const long long bx, const long long by) {
        return (static_cast<unsigned long long>(uint32_t(bx)) << 32) |
               static_cast<unsigned long long>(uint32_t(by));
    };

    std::unordered_map<unsigned long long, std::vector<size_t>> bins;
    for (size_t i = 0; i < segment_count; ++i) {
        const Segment &segment = segments[i];
        const long long bx0 = (long long)std::floor((std::min(segment.a.x(), segment.b.x()) - max_width) / bin_size);
        const long long by0 = (long long)std::floor((std::min(segment.a.y(), segment.b.y()) - max_width) / bin_size);
        const long long bx1 = (long long)std::floor((std::max(segment.a.x(), segment.b.x()) + max_width) / bin_size);
        const long long by1 = (long long)std::floor((std::max(segment.a.y(), segment.b.y()) + max_width) / bin_size);
        for (long long by = by0; by <= by1; ++by)
            for (long long bx = bx0; bx <= bx1; ++bx)
                bins[bin_key(bx, by)].push_back(i);
    }

    std::unordered_set<unsigned long long> checked;
    for (size_t i = 0; i < segment_count; ++i) {
        const Segment &segment = segments[i];
        const long long bx0 = (long long)std::floor((std::min(segment.a.x(), segment.b.x()) - max_width) / bin_size);
        const long long by0 = (long long)std::floor((std::min(segment.a.y(), segment.b.y()) - max_width) / bin_size);
        const long long bx1 = (long long)std::floor((std::max(segment.a.x(), segment.b.x()) + max_width) / bin_size);
        const long long by1 = (long long)std::floor((std::max(segment.a.y(), segment.b.y()) + max_width) / bin_size);
        for (long long by = by0; by <= by1; ++by)
            for (long long bx = bx0; bx <= bx1; ++bx)
                if (auto it = bins.find(bin_key(bx, by)); it != bins.end())
                    for (const size_t j : it->second) {
                        if (j <= i)
                            continue;
                        const unsigned long long pair_key = (unsigned long long(i) << 32) ^ unsigned(j);
                        if (!checked.insert(pair_key).second)
                            continue;
                        size_t index_distance = j - i;
                        if (closed_path)
                            index_distance = std::min(index_distance, segment_count - index_distance);
                        if (index_distance <= 1)
                            continue;
                        double path_gap = std::max(0.0, prefix[j] - prefix[i + 1]);
                        if (closed_path) {
                            const double occupied_span = prefix[j + 1] - prefix[i];
                            path_gap = std::min(path_gap, std::max(0.0, path_length - occupied_span));
                        }
                        if (path_gap <= local_skip_distance)
                            continue;
                        const double distance = segment_distance(segment.a, segment.b, segments[j].a, segments[j].b);
                        if (distance <= max_width) {
                            clearances[i] = std::min(clearances[i], distance);
                            clearances[j] = std::min(clearances[j], distance);
                        }
                    }
    }
    return clearances;
}

void clamp_segment_widths(
    const Points &path,
    const Flow &flow,
    const double max_multiplier,
    std::vector<float> &multipliers)
{
    if (path.size() < 2 || multipliers.size() + 1 != path.size())
        return;

    const double line_width = double(flow.scaled_width());
    const double line_height = scale_(double(flow.height()));
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double cross_section_width = line_width - rounded_corner_loss;
    if (!(cross_section_width > EPS))
        return;
    const double min_multiplier = std::clamp(minimum_extrusion_multiplier(flow), 0.01, 1.0);

    for (size_t i = 0; i < multipliers.size(); ++i) {
        multipliers[i] = float(std::clamp(double(multipliers[i]), min_multiplier, max_multiplier));
    }
}

bool distribute_volume_toward_uncovered_area(
    const ExPolygons &printable_area,
    const Points &path,
    const Flow &flow,
    const double max_multiplier_limit,
    const double max_physical_width,
    const double target_volume,
    const double deposited_volume,
    std::vector<float> &multipliers)
{
    if (path.size() < 2 || multipliers.size() + 1 != path.size() || target_volume <= deposited_volume)
        return false;

    const double line_height = scale_(double(flow.height()));
    const double line_width = double(flow.scaled_width());
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double cross_section_width = line_width - rounded_corner_loss;
    if (!(cross_section_width > EPS))
        return false;

    std::map<coord_t, Polylines> segments_by_width;
    for (size_t i = 1; i < path.size(); ++i) {
        if (path[i - 1] == path[i])
            continue;
        const coord_t width = coord_t(std::llround(
            rounded_corner_loss + cross_section_width * double(multipliers[i - 1])));
        segments_by_width[width].emplace_back(Points { path[i - 1], path[i] });
    }

    Polygons swept;
    for (const auto &[width, segments] : segments_by_width) {
        Polygons group = offset(
            segments,
            float(double(width) * 0.5),
            ClipperLib::jtRound,
            SCALED_RESOLUTION,
            ClipperLib::etOpenRound);
        swept.insert(swept.end(), std::make_move_iterator(group.begin()), std::make_move_iterator(group.end()));
    }
    const ExPolygons uncovered = diff_ex(printable_area, union_(swept));

    std::vector<bool> selected(multipliers.size(), false);
    auto select_near_loop = [&](const Points &points) {
        const size_t stride = std::max<size_t>(1, points.size() / 128);
        for (size_t i = 0; i < points.size(); i += stride) {
            const OpenProjection projection = project_open_polyline(path, points[i]);
            const size_t center = std::min(multipliers.size() - 1, size_t(std::floor(projection.index)));
            selected[center] = true;
        }
    };
    for (const ExPolygon &expoly : uncovered) {
        select_near_loop(expoly.contour.points);
        for (const Polygon &hole : expoly.holes)
            select_near_loop(hole.points);
    }

    const auto boundary_distance = [&](const Point &point) {
        double distance = std::numeric_limits<double>::infinity();
        for (const ExPolygon &expoly : printable_area) {
            distance = std::min(distance, distance_to_loop(point, expoly.contour.points));
            for (const Polygon &hole : expoly.holes)
                distance = std::min(distance, distance_to_loop(point, hole.points));
        }
        return distance;
    };
    std::vector<double> max_multiplier(multipliers.size(), max_multiplier_limit);
    for (size_t i = 0; i < multipliers.size(); ++i) {
        const Point midpoint = lerp_point(path[i], path[i + 1], 0.5);
        const double clearance = std::min({
            boundary_distance(path[i]),
            boundary_distance(midpoint),
            boundary_distance(path[i + 1]),
        });
        const double clearance_multiplier =
            (2.0 * clearance - rounded_corner_loss) / cross_section_width;
        double safe_multiplier = std::min(max_multiplier_limit, clearance_multiplier);
        max_multiplier[i] = std::clamp(safe_multiplier, double(multipliers[i]), max_multiplier_limit);
    }

    const double extra_volume = target_volume - deposited_volume;
    auto capacity = [&](const std::vector<bool> &mask) {
        double out = 0.0;
        for (size_t i = 0; i < multipliers.size(); ++i)
            if (mask[i])
                out += point_distance(path[i], path[i + 1]) * line_height * cross_section_width *
                       (max_multiplier[i] - double(multipliers[i]));
        return out;
    };
    if (capacity(selected) + EPS < extra_volume)
        std::fill(selected.begin(), selected.end(), true);
    const double safe_capacity = capacity(selected);
    if (!(safe_capacity > EPS))
        return false;
    const double distributable_volume = std::min(extra_volume, safe_capacity);

    std::vector<double> adjusted(multipliers.begin(), multipliers.end());
    double remaining = distributable_volume;
    for (size_t pass = 0; pass < multipliers.size() + 1 && remaining > EPS; ++pass) {
        double base_volume = 0.0;
        double max_common_increment = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < adjusted.size(); ++i) {
            if (!selected[i] || adjusted[i] >= max_multiplier[i] - 1e-9)
                continue;
            base_volume += point_distance(path[i], path[i + 1]) * line_height * cross_section_width;
            max_common_increment = std::min(max_common_increment, max_multiplier[i] - adjusted[i]);
        }
        if (!(base_volume > EPS) || !std::isfinite(max_common_increment))
            break;
        const double increment = std::min(remaining / base_volume, max_common_increment);
        double used = 0.0;
        for (size_t i = 0; i < adjusted.size(); ++i) {
            if (!selected[i] || adjusted[i] >= max_multiplier[i] - 1e-9)
                continue;
            const double segment_base = point_distance(path[i], path[i + 1]) * line_height * cross_section_width;
            adjusted[i] += increment;
            used += segment_base * increment;
        }
        if (!(used > EPS))
            break;
        remaining = std::max(0.0, remaining - used);
    }
    if (remaining > std::max(EPS, distributable_volume * 1e-9))
        return false;

    for (size_t i = 0; i < multipliers.size(); ++i)
        multipliers[i] = float(adjusted[i]);
    return true;
}

bool rebalance_widths_toward_uncovered_area(
    const ExPolygons &printable_area,
    const Points &path,
    const Flow &flow,
    const double max_multiplier_limit,
    const double max_physical_width,
    std::vector<float> &multipliers)
{
    if (path.size() < 2 || multipliers.size() + 1 != path.size())
        return false;

    const double line_height = scale_(double(flow.height()));
    const double line_width = double(flow.scaled_width());
    const double spacing = double(flow.scaled_spacing());
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double cross_section_width = line_width - rounded_corner_loss;
    if (!(cross_section_width > EPS))
        return false;

    std::map<coord_t, Polylines> segments_by_width;
    std::vector<double> physical_widths(multipliers.size(), line_width);
    for (size_t i = 0; i < multipliers.size(); ++i) {
        physical_widths[i] = rounded_corner_loss + cross_section_width * double(multipliers[i]);
        if (path[i] != path[i + 1])
            segments_by_width[coord_t(std::llround(physical_widths[i]))].emplace_back(
                Points { path[i], path[i + 1] });
    }
    Polygons swept;
    for (const auto &[width, segments] : segments_by_width) {
        Polygons group = offset(
            segments, float(double(width) * 0.5), ClipperLib::jtRound, SCALED_RESOLUTION, ClipperLib::etOpenRound);
        swept.insert(swept.end(), std::make_move_iterator(group.begin()), std::make_move_iterator(group.end()));
    }
    const ExPolygons uncovered = diff_ex(printable_area, union_(swept));
    const double uncovered_area = std::max(0.0, std::abs(area(uncovered)));
    if (!(uncovered_area > EPS))
        return false;

    std::vector<bool> recipient(multipliers.size(), false);
    auto select_near_loop = [&](const Points &points) {
        const size_t stride = std::max<size_t>(1, points.size() / 128);
        for (size_t i = 0; i < points.size(); i += stride) {
            const OpenProjection projection = project_open_polyline(path, points[i]);
            recipient[std::min(multipliers.size() - 1, size_t(std::floor(projection.index)))] = true;
        }
    };
    for (const ExPolygon &expoly : uncovered) {
        select_near_loop(expoly.contour.points);
        for (const Polygon &hole : expoly.holes)
            select_near_loop(hole.points);
    }

    const auto boundary_distance = [&](const Point &point) {
        double distance = std::numeric_limits<double>::infinity();
        for (const ExPolygon &expoly : printable_area) {
            distance = std::min(distance, distance_to_loop(point, expoly.contour.points));
            for (const Polygon &hole : expoly.holes)
                distance = std::min(distance, distance_to_loop(point, hole.points));
        }
        return distance;
    };
    const double min_physical_width = minimum_physical_width(flow);
    const double min_multiplier = std::clamp(minimum_extrusion_multiplier(flow), 0.01, 1.0);
    const std::vector<double> clearances = nonlocal_segment_clearances(path, max_physical_width, spacing);
    std::vector<double> recipient_capacity(multipliers.size(), 0.0);
    std::vector<double> donor_capacity(multipliers.size(), 0.0);
    double total_recipient_capacity = 0.0;
    double total_donor_capacity = 0.0;
    for (size_t i = 0; i < multipliers.size(); ++i) {
        const double base_volume = point_distance(path[i], path[i + 1]) * line_height * cross_section_width;
        if (!(base_volume > EPS))
            continue;
        if (recipient[i]) {
            const Point midpoint = lerp_point(path[i], path[i + 1], 0.5);
            const double clearance = std::min({
                boundary_distance(path[i]), boundary_distance(midpoint), boundary_distance(path[i + 1]),
            });
            const double boundary_multiplier = (2.0 * clearance - rounded_corner_loss) / cross_section_width;
            const double maximum = std::clamp(
                std::min(max_multiplier_limit, boundary_multiplier), double(multipliers[i]), max_multiplier_limit);
            recipient_capacity[i] = base_volume * (maximum - double(multipliers[i]));
            total_recipient_capacity += recipient_capacity[i];
        } else if (clearances[i] < physical_widths[i] * 0.95) {
            donor_capacity[i] = base_volume * std::max(0.0, double(multipliers[i]) - min_multiplier);
            total_donor_capacity += donor_capacity[i];
        }
    }

    const double requested_volume = uncovered_area * line_height * 0.80;
    const double transferred_volume = std::min({ requested_volume, total_recipient_capacity, total_donor_capacity });
    if (!(transferred_volume > EPS))
        return false;
    const double recipient_fraction = transferred_volume / total_recipient_capacity;
    const double donor_fraction = transferred_volume / total_donor_capacity;
    for (size_t i = 0; i < multipliers.size(); ++i) {
        const double base_volume = point_distance(path[i], path[i + 1]) * line_height * cross_section_width;
        if (!(base_volume > EPS))
            continue;
        double adjusted = double(multipliers[i]);
        if (recipient_capacity[i] > 0.0)
            adjusted += recipient_capacity[i] * recipient_fraction / base_volume;
        if (donor_capacity[i] > 0.0)
            adjusted -= donor_capacity[i] * donor_fraction / base_volume;
        multipliers[i] = float(std::clamp(adjusted, min_multiplier, max_multiplier_limit));
    }
    return true;
}

struct LayerPathResult
{
    Polyline path;
    std::vector<float> extrusion_multipliers;
    bool complexity_limited { false };
};

int count_turnback_violations(const Points &path, const double line_width, std::string *first_detail = nullptr)
{
    if (path.size() < 3)
        return 0;

    const bool closed = path.front() == path.back();
    const size_t unique_points = closed ? path.size() - 1 : path.size();
    if (unique_points < 3)
        return 0;

    const double min_leg = line_width * 0.35;
    static constexpr double MAX_SAFE_TURN_COSINE = -0.7071067811865476; // 135 degrees.
    int violations = 0;
    const size_t first = closed ? 0 : 1;
    const size_t end = closed ? unique_points : unique_points - 1;
    for (size_t i = first; i < end; ++i) {
        const size_t prev = (i + unique_points - 1) % unique_points;
        const size_t next = (i + 1) % unique_points;
        const Vec2d incoming = (path[i] - path[prev]).cast<double>();
        const Vec2d outgoing = (path[next] - path[i]).cast<double>();
        const double incoming_length = incoming.norm();
        const double outgoing_length = outgoing.norm();
        if (incoming_length < min_leg || outgoing_length < min_leg)
            continue;
        const double cosine = std::clamp(incoming.dot(outgoing) / (incoming_length * outgoing_length), -1.0, 1.0);
        if (cosine < MAX_SAFE_TURN_COSINE) {
            if (first_detail != nullptr && first_detail->empty()) {
                std::ostringstream detail;
                const size_t prev2 = (prev + unique_points - 1) % unique_points;
                const size_t next2 = (next + 1) % unique_points;
                detail << "point " << i << " around " << unscale<double>(path[prev2].x()) << ","
                       << unscale<double>(path[prev2].y()) << "->" << unscale<double>(path[prev].x()) << ","
                       << unscale<double>(path[prev].y()) << "->" << unscale<double>(path[i].x()) << ","
                       << unscale<double>(path[i].y()) << "->" << unscale<double>(path[next].x()) << ","
                       << unscale<double>(path[next].y()) << "->" << unscale<double>(path[next2].x()) << ","
                       << unscale<double>(path[next2].y()) << " cosine=" << cosine;
                *first_detail = detail.str();
            }
            ++violations;
        }
    }
    return violations;
}

Points uncross_closed_tour(
    const ExPolygons &printable_area,
    Points path,
    const double line_width,
    const double spacing)
{
    if (path.size() < 5 || path.front() != path.back())
        return path;

    for (size_t repair = 0; repair < 128; ++repair) {
        const PairMetrics baseline = path_pair_metrics(path, spacing);
        if (baseline.crossings == 0)
            break;
        bool improved = false;
        size_t evaluated_touching_candidates = 0;

        const size_t segment_count = path.size() - 1;
        std::vector<double> prefix { 0.0 };
        prefix.reserve(path.size());
        for (size_t i = 0; i < segment_count; ++i)
            prefix.emplace_back(prefix.back() + point_distance(path[i], path[i + 1]));
        for (size_t i = 0; i < segment_count && !improved; ++i) {
            for (size_t j = i + 2; j < segment_count; ++j) {
                if (i == 0 && j + 1 == segment_count)
                    continue;

                // This loop is looking only for intersecting or touching
                // segment pairs. Reject disjoint integer-coordinate bounding
                // boxes before the substantially more expensive orientation
                // and point-to-segment distance calculations. This preserves
                // the exact candidate order and repair result.
                if (std::max(path[i].x(), path[i + 1].x()) < std::min(path[j].x(), path[j + 1].x()) ||
                    std::max(path[j].x(), path[j + 1].x()) < std::min(path[i].x(), path[i + 1].x()) ||
                    std::max(path[i].y(), path[i + 1].y()) < std::min(path[j].y(), path[j + 1].y()) ||
                    std::max(path[j].y(), path[j + 1].y()) < std::min(path[i].y(), path[i + 1].y()))
                    continue;

                const double o1 = orient2d(path[i], path[i + 1], path[j]);
                const double o2 = orient2d(path[i], path[i + 1], path[j + 1]);
                const double o3 = orient2d(path[j], path[j + 1], path[i]);
                const double o4 = orient2d(path[j], path[j + 1], path[i + 1]);
                const bool proper_crossing =
                    ((o1 > 1e-7 && o2 < -1e-7) || (o1 < -1e-7 && o2 > 1e-7)) &&
                    ((o3 > 1e-7 && o4 < -1e-7) || (o3 < -1e-7 && o4 > 1e-7));
                if (!proper_crossing) {
                    if (segment_distance(path[i], path[i + 1], path[j], path[j + 1]) > 1e-7)
                        continue;
                    const double occupied_span = prefix[j + 1] - prefix[i];
                    const double path_gap = std::min(
                        std::max(0.0, prefix[j] - prefix[i + 1]),
                        std::max(0.0, prefix.back() - occupied_span));
                    if (path_gap <= 1e-9)
                        continue;
                    if (++evaluated_touching_candidates > 32)
                        break;
                }

                Points candidate = path;
                std::reverse(candidate.begin() + i + 1, candidate.begin() + j + 1);
                if (count_containment_violations(printable_area, candidate, spacing) != 0)
                    continue;
                const PairMetrics metrics = path_pair_metrics(
                    candidate, spacing, nullptr, 0.0, baseline.crossings);
                if (metrics.crossings < baseline.crossings) {
                    path = std::move(candidate);
                    improved = true;
                    break;
                }

                const std::array<size_t, 4> removable { i, i + 1, j, j + 1 };
                for (const size_t vertex : removable) {
                    if (vertex == 0 || vertex + 1 >= path.size())
                        continue;
                    Points simplified = path;
                    simplified.erase(simplified.begin() + vertex);
                    if (count_containment_violations(printable_area, simplified, spacing) != 0)
                        continue;
                    const PairMetrics simplified_metrics = path_pair_metrics(
                        simplified, spacing, nullptr, 0.0, baseline.crossings);
                    if (simplified_metrics.crossings >= baseline.crossings)
                        continue;
                    path = std::move(simplified);
                    improved = true;
                    break;
                }
                if (improved)
                    break;
                const auto try_detour = [&](const size_t moving_segment, const size_t obstacle_segment) {
                    const Point &oa = path[obstacle_segment];
                    const Point &ob = path[obstacle_segment + 1];
                    const Vec2d obstacle = (ob - oa).cast<double>();
                    const double obstacle_length = obstacle.norm();
                    if (!(obstacle_length > EPS))
                        return false;
                    const Vec2d axis = obstacle / obstacle_length;
                    const Vec2d normal(-axis.y(), axis.x());
                    for (const bool around_end : { false, true }) {
                        const Vec2d endpoint = (around_end ? ob : oa).cast<double>();
                        const Vec2d outward = around_end ? axis : -axis;
                        for (const double side : { -1.0, 1.0 }) {
                            const Vec2d waypoint = endpoint + outward * (spacing * 0.45) +
                                                   normal * (side * spacing * 0.18);
                            Points detoured = path;
                            detoured.insert(
                                detoured.begin() + moving_segment + 1,
                                Point(waypoint.x(), waypoint.y()));
                            if (count_containment_violations(printable_area, detoured, spacing) != 0)
                                continue;
                            const PairMetrics detoured_metrics = path_pair_metrics(
                                detoured, spacing, nullptr, 0.0, baseline.crossings);
                            if (detoured_metrics.crossings >= baseline.crossings)
                                continue;
                            path = std::move(detoured);
                            return true;
                        }
                    }
                    return false;
                };
                if (try_detour(j, i) || try_detour(i, j)) {
                    improved = true;
                    break;
                }
                for (const size_t vertex : removable) {
                    for (size_t radius = 1; radius <= 3; ++radius) {
                        const size_t erase_begin = std::max<size_t>(1, vertex > radius ? vertex - radius : 1);
                        const size_t erase_end = std::min(path.size() - 1, vertex + radius + 1);
                        if (erase_begin >= erase_end ||
                            point_distance(path[erase_begin - 1], path[erase_end]) > spacing * 6.0)
                            continue;
                        Points simplified = path;
                        simplified.erase(simplified.begin() + erase_begin, simplified.begin() + erase_end);
                        if (count_containment_violations(printable_area, simplified, spacing) != 0)
                            continue;
                        const PairMetrics simplified_metrics = path_pair_metrics(
                            simplified, spacing, nullptr, 0.0, baseline.crossings);
                        if (simplified_metrics.crossings >= baseline.crossings)
                            continue;
                        path = std::move(simplified);
                        improved = true;
                        break;
                    }
                    if (improved)
                        break;
                }
                if (improved)
                    break;
            }
        }
        if (!improved)
            break;
    }
    return path;
}

Points remove_short_turnback_spikes(Points path, const double line_width, const double spacing)
{
    // A contour port may occasionally leave a short out-and-back vertex after
    // the in/out branches are joined.  Removing such a vertex is safe only as
    // a candidate operation: it must reduce the conservative turnback count,
    // must not worsen centerline crossings or spacing, and the unchanged final
    // footprint validator still checks containment, coverage, overlap, width,
    // and material volume.
    for (size_t repair = 0; repair < 32 && path.size() >= 5; ++repair) {
        const int baseline_turnbacks = count_turnback_violations(path, line_width);
        if (baseline_turnbacks == 0)
            break;
        const PairMetrics baseline_pairs = path_pair_metrics(path, spacing);
        bool improved = false;

        // Preserve the explicit closure endpoints.  A seam turnback remains a
        // hard validation failure instead of silently rotating the path.
        for (size_t i = 1; i + 1 < path.size(); ++i) {
            const double incoming = point_distance(path[i - 1], path[i]);
            const double outgoing = point_distance(path[i], path[i + 1]);
            if (incoming < line_width * 0.35 || outgoing < line_width * 0.35)
                continue;

            const Vec2d a = (path[i] - path[i - 1]).cast<double>();
            const Vec2d b = (path[i + 1] - path[i]).cast<double>();
            const double cosine = std::clamp(a.dot(b) / (a.norm() * b.norm()), -1.0, 1.0);
            if (cosine >= -0.7071067811865476)
                continue;

            // Removing the apex alone may merely move a shallow reversal to
            // the preceding collinear port vertex.  Also try the two-vertex
            // spike as one semantic repair.
            for (size_t remove_count = 1; remove_count <= 6; ++remove_count) {
                if (i + 1 < remove_count + 1)
                    continue;
                const size_t erase_begin = i + 1 - remove_count;
                if (erase_begin == 0)
                    continue;
                const size_t right = i + 1;
                const double replacement = point_distance(path[erase_begin - 1], path[right]);
                double removed_arc = 0.0;
                for (size_t j = erase_begin; j <= right; ++j)
                    removed_arc += point_distance(path[j - 1], path[j]);
                // Wide pockets may leave a longer radial out-and-back spur at
                // a branch port.  The replacement is still local relative to
                // the deposition spacing, and the unchanged final validator
                // rejects any shortcut that loses containment or coverage.
                if (replacement > spacing * 8.0 || removed_arc <= replacement + EPS)
                    continue;

                Points candidate = path;
                candidate.erase(candidate.begin() + erase_begin, candidate.begin() + i + 1);
                const int candidate_turnbacks = count_turnback_violations(candidate, line_width);
                if (candidate_turnbacks >= baseline_turnbacks)
                    continue;
                const PairMetrics candidate_pairs = path_pair_metrics(candidate, spacing);
                if (candidate_pairs.crossings > baseline_pairs.crossings)
                    continue;

                path = std::move(candidate);
                improved = true;
                break;
            }
            if (!improved) {
                for (size_t forward_count = 2; forward_count <= 6; ++forward_count) {
                    const size_t right = i + forward_count;
                    if (right >= path.size())
                        continue;
                    const double replacement = point_distance(path[i - 1], path[right]);
                    double removed_arc = 0.0;
                    for (size_t j = i; j <= right; ++j)
                        removed_arc += point_distance(path[j - 1], path[j]);
                    if (replacement > spacing * 8.0 || removed_arc <= replacement + EPS)
                        continue;
                    Points candidate = path;
                    candidate.erase(candidate.begin() + i, candidate.begin() + right);
                    if (count_turnback_violations(candidate, line_width) >= baseline_turnbacks)
                        continue;
                    const PairMetrics candidate_pairs = path_pair_metrics(candidate, spacing);
                    if (candidate_pairs.crossings > baseline_pairs.crossings)
                        continue;
                    path = std::move(candidate);
                    improved = true;
                    break;
                }
            }
            if (improved)
                break;
        }

        if (!improved)
            break;
    }
    return path;
}

Points remove_sub_gcode_resolution_segments(Points path)
{
    if (path.size() < 4 || path.front() != path.back())
        return path;

    // XYZ is serialized to 0.001 mm. Two distinct points less than the
    // diagonal of one output cell apart may become the same controller
    // coordinate depending on their phase relative to the rounding grid.
    // Remove those vertices before extrusion metadata is calculated; the
    // unchanged final geometry validator then certifies the simplified path.
    const double min_segment_length = scale_(0.002);
    Points filtered;
    filtered.reserve(path.size());
    filtered.emplace_back(path.front());
    for (size_t i = 1; i + 1 < path.size(); ++i)
        if (point_distance(filtered.back(), path[i]) >= min_segment_length)
            filtered.emplace_back(path[i]);
    while (filtered.size() > 3 && point_distance(filtered.back(), filtered.front()) < min_segment_length)
        filtered.pop_back();
    filtered.emplace_back(filtered.front());
    return filtered;
}

LayerPathResult generate_layer_path_candidate(
    const ExPolygons &printable_area,
    const Flow &flow,
    const double configured_max_line_width,
    const double spacing_factor)
{
    CandidatePerfTimer candidate_perf(printable_area.size(), spacing_factor);
    const double max_physical_width = maximum_physical_width(flow, configured_max_line_width);
    const double max_multiplier_limit = maximum_extrusion_multiplier(flow, configured_max_line_width);
    if (!(max_physical_width > EPS) || !(max_multiplier_limit >= 1.0))
        return {};
    // Plan slightly denser centerlines, then normalize their cross-sections to
    // the target material volume. This distributes residual-width error over
    // the whole fill instead of leaving a terminal hole whose size depends on
    // the offset phase.
    const double spacing = double(flow.scaled_spacing()) * spacing_factor;
    std::vector<ContourLoop> contours = generate_offset_contours(
        printable_area, double(flow.scaled_width()), spacing, 3);
    if (contours.empty())
        contours = generate_narrow_feature_contours(printable_area, flow, spacing);
    contours = add_ring_medial_contours(printable_area, std::move(contours), flow);
    contours = add_terminal_medial_gap_fill(
        contours,
        double(flow.scaled_width()),
        scale_(double(flow.height())),
        spacing,
        max_physical_width);
    candidate_perf.contour_count = contours.size();
    if (contours.empty())
        return {};

    const auto grouped = loops_by_level(contours);
    if (grouped.empty())
        return {};

    const std::vector<ContourLoop> ordered = ordered_loops_for_spiral(contours);
    if (ordered.empty())
        return {};

    bool multi_loop = false;
    for (const auto &[level, loops] : grouped)
        multi_loop |= loops.size() > 1;

    const bool branch_single_island = multi_loop;
    const bool multi_hole = collect_holes(printable_area).size() >= 2;
    const Point start_anchor = point_at_closed_fraction(ordered.front().points, 0.0);

    const ContourLoop *outer_loop = &ordered.front();
    for (const ContourLoop &loop : contours) {
        if (loop.level_index == 0 && loop.area > outer_loop->area)
            outer_loop = &loop;
    }

    const auto finish_path = [&](Points path) {
        const PerfClock::time_point finish_start = PerfClock::now();
        LayerPathResult result;
        if (path.size() < 2)
            return result;
        std::vector<ContourLoop> path_contours = contours;

        const auto profile_phase = [&](const char *name, const PerfClock::time_point start, const PerfCounters &before) {
            if (!perf_timing_enabled())
                return;
            BOOST_LOG_TRIVIAL(warning) << "Continuous Fermat profile phase=" << name
                                       << " ms=" << elapsed_milliseconds(start) << " points=" << path.size()
                                       << " pair_calls=" << perf_counters.pair_metrics.calls - before.pair_metrics.calls
                                       << " pair_ms=" << perf_counters.pair_metrics.milliseconds - before.pair_metrics.milliseconds
                                       << " containment_calls=" << perf_counters.containment.calls - before.containment.calls
                                       << " containment_ms=" << perf_counters.containment.milliseconds - before.containment.milliseconds;
        };

        PerfClock::time_point phase_start = PerfClock::now();
        PerfCounters phase_before = perf_counters;
        path = replace_long_connectors_with_contour_arcs(std::move(path), contours, spacing);
        profile_phase("replace_connectors", phase_start, phase_before);
        phase_start = PerfClock::now();
        phase_before = perf_counters;
        path = complete_outer_boundary_cycle(path, *outer_loop, spacing);
        profile_phase("outer_boundary", phase_start, phase_before);
        if (multi_hole) {
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = insert_uncovered_pocket_spirals(printable_area, contours, std::move(path), spacing);
            profile_phase("uncovered_pockets", phase_start, phase_before);
        }
        if (multi_loop) {
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = insert_residual_gap_spirals(printable_area, flow, std::move(path), spacing, &path_contours);
            profile_phase("residual_gaps", phase_start, phase_before);
        }
        const bool dense_path = path.size() > MAX_EXHAUSTIVE_REPAIR_POINTS;
        result.complexity_limited = dense_path;
        if (!dense_path) {
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = uncross_closed_tour(printable_area, std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("uncross_1", phase_start, phase_before);
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = remove_short_turnback_spikes(std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("turnbacks_1", phase_start, phase_before);
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = uncross_closed_tour(printable_area, std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("uncross_2", phase_start, phase_before);
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = remove_short_turnback_spikes(std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("turnbacks_2", phase_start, phase_before);
        } else if (perf_timing_enabled()) {
            BOOST_LOG_TRIVIAL(warning) << "Continuous Fermat profile phase=topology_repairs_skipped"
                                       << " points=" << path.size()
                                       << " limit=" << MAX_EXHAUSTIVE_REPAIR_POINTS;
        }
        if (dense_path) {
            // The contour planner deliberately oversamples curves for robust
            // connector placement. Once the route is fixed, retain that
            // geometry only to 0.01 mm before calculating per-segment flow.
            // This is far below normal extrusion width and G-code resolution,
            // while avoiding thousands of redundant emitter segments.
            Polyline simplified(path);
            simplified.simplify(scale_(0.01));
            if (simplified.points.size() >= 4 && simplified.points.front() == simplified.points.back())
                path = std::move(simplified.points);

            // Repair the compact route rather than its oversampled planning
            // representation. Pairwise work now sees hundreds of segments
            // instead of many thousands.
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = uncross_closed_tour(printable_area, std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("dense_uncross_1", phase_start, phase_before);
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = remove_short_turnback_spikes(std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("dense_turnbacks_1", phase_start, phase_before);
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = uncross_closed_tour(printable_area, std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("dense_uncross_2", phase_start, phase_before);
            phase_start = PerfClock::now();
            phase_before = perf_counters;
            path = remove_short_turnback_spikes(std::move(path), double(flow.scaled_width()), spacing);
            profile_phase("dense_turnbacks_2", phase_start, phase_before);
        }
        path = remove_sub_gcode_resolution_segments(std::move(path));

        phase_start = PerfClock::now();
        phase_before = perf_counters;
        result.extrusion_multipliers =
            segment_extrusion_multipliers_for_path(path, path_contours, double(flow.scaled_width()), spacing);
        clamp_segment_widths(
            path, flow, max_multiplier_limit, result.extrusion_multipliers);
        if (result.extrusion_multipliers.size() + 1 == path.size()) {
            const double line_height = scale_(double(flow.height()));
            const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
            const double nominal_cross_section = line_height * (double(flow.scaled_width()) - rounded_corner_loss);
            const double solid_volume = std::abs(area(printable_area)) * line_height;
            // Leave a small, explicit amount of the validator's material
            // budget available for footprint-directed seam closure.
            const double target_volume = solid_volume * 1.019;
            double deposited_volume = 0.0;
            for (size_t i = 1; i < path.size(); ++i)
                deposited_volume += point_distance(path[i - 1], path[i]) * nominal_cross_section *
                                    double(result.extrusion_multipliers[i - 1]);

            if (deposited_volume > target_volume) {
                const double min_multiplier = std::clamp(minimum_extrusion_multiplier(flow), 0.01, 1.0);
                const double normalization = target_volume / deposited_volume;
                for (float &multiplier : result.extrusion_multipliers)
                    multiplier = float(std::max(min_multiplier, double(multiplier) * normalization));
                deposited_volume = 0.0;
                for (size_t i = 1; i < path.size(); ++i)
                    deposited_volume += point_distance(path[i - 1], path[i]) * nominal_cross_section *
                                        double(result.extrusion_multipliers[i - 1]);
            }

            if (deposited_volume > EPS && target_volume > deposited_volume) {
                const double correction = target_volume / deposited_volume;
                const double max_multiplier = *std::max_element(
                    result.extrusion_multipliers.begin(), result.extrusion_multipliers.end());
                if (correction * max_multiplier <= max_multiplier_limit)
                    distribute_volume_toward_uncovered_area(
                        printable_area,
                        path,
                        flow,
                        max_multiplier_limit,
                        max_physical_width,
                        target_volume,
                        deposited_volume,
                        result.extrusion_multipliers);
            }
        }
        profile_phase("extrusion_widths", phase_start, phase_before);
        result.path = Polyline(std::move(path));
        candidate_perf.output_points = result.path.points.size();
        if (perf_timing_enabled())
            BOOST_LOG_TRIVIAL(warning) << "Continuous Fermat profile finish_ms=" << elapsed_milliseconds(finish_start)
                                       << " points=" << result.path.points.size() << " contours=" << path_contours.size()
                                       << " multi_loop=" << multi_loop << " multi_hole=" << multi_hole;
        return result;
    };

    if (branch_single_island) {
        const PerfClock::time_point branch_start = PerfClock::now();
        Points branch_path = build_branch_connected_fermat(contours, printable_area, start_anchor, spacing);
        if (perf_timing_enabled())
            BOOST_LOG_TRIVIAL(warning) << "Continuous Fermat profile branch_ms=" << elapsed_milliseconds(branch_start)
                                       << " branch_points=" << branch_path.size() << " contours=" << contours.size();

        // The post-finish alternative is recovery for adaptive planning or a
        // thin multiply-connected feature. At nominal width on broad hole
        // layers, retain the established inexpensive centerline choice;
        // otherwise a harmless footprint defect can make every such layer
        // pay for two complete repair and validation passes.
        double boundary_length = 0.0;
        for (const ExPolygon &expoly : printable_area) {
            boundary_length += expoly.contour.length();
            for (const Polygon &hole : expoly.holes)
                boundary_length += hole.length();
        }
        const double thickness_proxy = boundary_length > EPS ?
            2.0 * std::abs(area(printable_area)) / boundary_length :
            std::numeric_limits<double>::infinity();
        const bool adaptive_planning_width =
            double(flow.scaled_width()) > scale_(double(flow.nozzle_diameter()) * 1.01);
        const bool thin_feature = thickness_proxy <= double(flow.scaled_width()) * 3.0;
        if (!adaptive_planning_width && !thin_feature) {
            if (branch_path.size() > MAX_EXHAUSTIVE_REPAIR_POINTS)
                return finish_path(std::move(branch_path));
            Points chain_path = build_best_single_chain_path(ordered, start_anchor, spacing);
            const PairMetrics branch_metrics = path_pair_metrics(branch_path, spacing);
            const PairMetrics chain_metrics = path_pair_metrics(chain_path, spacing);
            const bool chain_is_safer = !multi_hole && !chain_path.empty() && (branch_path.empty() ||
                chain_metrics.crossings < branch_metrics.crossings ||
                (chain_metrics.crossings == branch_metrics.crossings && chain_metrics.close_pairs < branch_metrics.close_pairs) ||
                (chain_metrics.crossings == branch_metrics.crossings && chain_metrics.close_pairs == branch_metrics.close_pairs &&
                 chain_metrics.min_spacing > branch_metrics.min_spacing));
            return finish_path(chain_is_safer ? std::move(chain_path) : std::move(branch_path));
        }

        LayerPathResult branch_result = finish_path(std::move(branch_path));
        const PathValidation branch_validation = validate_layer_path(
            printable_area,
            flow,
            branch_result.path,
            branch_result.extrusion_multipliers,
            configured_max_line_width);
        if (branch_validation.ok)
            return branch_result;

        const auto has_reason = [](const PathValidation &validation, const char *token) {
            return validation.reason.find(token) != std::string::npos;
        };
        const auto has_safety_defect = [&](const PathValidation &validation) {
            return !validation.closed || validation.containment_violations != 0 || validation.crossings != 0 ||
                   validation.turnback_violations != 0 || has_reason(validation, " outside=") ||
                   has_reason(validation, " redeposition=");
        };
        // Coverage or material underfill is intentionally left to the
        // geometry/derived-width retries. Building a second full route for a
        // safe underfill only doubles the dominant repair and audit work.
        if (!has_safety_defect(branch_validation))
            return branch_result;

        Points chain_path = build_best_single_chain_path(ordered, start_anchor, spacing);
        if (chain_path.empty())
            return branch_result;
        LayerPathResult chain_result = finish_path(std::move(chain_path));
        const PathValidation chain_validation = validate_layer_path(
            printable_area,
            flow,
            chain_result.path,
            chain_result.extrusion_multipliers,
            configured_max_line_width);
        if (chain_validation.ok)
            return chain_result;

        const auto safety_failures = [&](const PathValidation &validation) {
            return int(!validation.closed) + int(validation.containment_violations) + validation.crossings +
                   validation.turnback_violations + int(has_reason(validation, " outside=")) +
                   int(has_reason(validation, " redeposition="));
        };
        const int branch_failures = safety_failures(branch_validation);
        const int chain_failures = safety_failures(chain_validation);
        const bool chain_is_safer = chain_failures < branch_failures ||
            (chain_failures == branch_failures &&
             chain_validation.containment_violations < branch_validation.containment_violations) ||
            (chain_failures == branch_failures &&
             chain_validation.containment_violations == branch_validation.containment_violations &&
             chain_validation.crossings < branch_validation.crossings) ||
            (chain_failures == branch_failures &&
             chain_validation.containment_violations == branch_validation.containment_violations &&
             chain_validation.crossings == branch_validation.crossings &&
             chain_validation.turnback_violations < branch_validation.turnback_violations) ||
            (chain_failures == branch_failures &&
             chain_validation.containment_violations == branch_validation.containment_violations &&
             chain_validation.crossings == branch_validation.crossings &&
             chain_validation.turnback_violations == branch_validation.turnback_violations &&
             chain_validation.outside_ratio < branch_validation.outside_ratio) ||
            (chain_failures == branch_failures &&
             chain_validation.containment_violations == branch_validation.containment_violations &&
             chain_validation.crossings == branch_validation.crossings &&
             chain_validation.turnback_violations == branch_validation.turnback_violations &&
             chain_validation.outside_ratio == branch_validation.outside_ratio &&
             chain_validation.redeposition_ratio < branch_validation.redeposition_ratio);
        return chain_is_safer ? std::move(chain_result) : std::move(branch_result);
    }

    // The default port distance preserves coverage on broad, low-aspect
    // beads.  If its completed variable-width footprint fails physical
    // separation, retry with wider in/out ports and accept the retry only if
    // the complete unchanged contract passes.  This is flow-driven recovery,
    // not a nozzle- or shape-specific exception.
    LayerPathResult primary = finish_path(build_best_single_chain_path(ordered, start_anchor, spacing));
    const PathValidation primary_validation =
        validate_layer_path(
            printable_area, flow, primary.path, primary.extrusion_multipliers, configured_max_line_width);
    if (primary_validation.ok ||
        (primary_validation.spacing_violations <= 3 && primary_validation.bead_overlap_violations == 0))
        return primary;

    LayerPathResult retry = finish_path(build_best_single_chain_path(ordered, start_anchor, spacing, true, 1.25));
    const PathValidation retry_validation = validate_layer_path(
        printable_area, flow, retry.path, retry.extrusion_multipliers, configured_max_line_width);
    if (retry_validation.ok)
        return retry;

    const auto hard_failures = [](const PathValidation &validation) {
        return int(validation.containment_violations) + validation.crossings +
               validation.bead_overlap_violations + validation.turnback_violations +
               std::max(0, validation.spacing_violations - 3);
    };
    const int primary_failures = hard_failures(primary_validation);
    const int retry_failures = hard_failures(retry_validation);
    if (retry_failures < primary_failures ||
        (retry_failures == primary_failures && retry_validation.coverage_ratio > primary_validation.coverage_ratio))
        return retry;
    return primary;
}

LayerPathResult generate_layer_path_result(
    const ExPolygons &printable_area,
    const Flow &flow,
    const double configured_max_line_width)
{
    LayerPathResult nominal = generate_layer_path_candidate(
        printable_area, flow, configured_max_line_width, 1.0);
    const PathValidation nominal_validation = validate_layer_path(
        printable_area, flow, nominal.path, nominal.extrusion_multipliers, configured_max_line_width);
    const bool nominal_is_safe_underfill = nominal_validation.containment_violations == 0 &&
        nominal_validation.crossings == 0 && nominal_validation.turnback_violations == 0 &&
        nominal_validation.outside_ratio <= 0.020 + 1e-9 && nominal_validation.redeposition_ratio <= 0.040 + 1e-9 &&
        nominal_validation.exact_coverage_ratio < 0.980 && nominal_validation.material_ratio < 0.980;
    if (nominal_validation.ok)
        return nominal;
    // Once a structurally valid path exists, only a single geometry-derived
    // width retry is worth delaying preview. Spacing sweeps and iterative
    // rebalance remain recovery tools for candidates that cannot be emitted.
    if (nominal_validation.emittable && (nominal.complexity_limited || !nominal_is_safe_underfill))
        return nominal;

    const auto hard_failures = [](const PathValidation &validation) {
        return int(validation.containment_violations) + validation.crossings + validation.turnback_violations;
    };
    LayerPathResult *best = &nominal;
    const PathValidation *best_validation = &nominal_validation;
    std::vector<LayerPathResult> alternatives;
    std::vector<PathValidation> alternative_validations;
    const bool has_holes = !collect_holes(printable_area).empty();
    const std::vector<double> spacing_factors = has_holes ?
        std::vector<double> { 0.97 } : std::vector<double> { 0.97, 0.94, 1.03 };
    alternatives.reserve(spacing_factors.size() + 11);
    alternative_validations.reserve(spacing_factors.size() + 11);

    // A fixed nominal first inset cannot fill some thin corridors or rings:
    // the retained centerlines are safe, but their footprint and material are
    // both too small. Infer one geometry-consistent bead width from the layer
    // area and the discovered centerline length, then regenerate the complete
    // offset geometry at that width and its corresponding Flow spacing. This
    // is a single derived retry, not a width sweep. The fixed quality metrics
    // remain visible whether or not the retry reaches them.
    if (nominal_is_safe_underfill && nominal.path.points.size() >= 4) {
        const double path_length = polyline_length(nominal.path.points, false);
        const double printable_area_value = std::abs(area(printable_area));
        const double line_height = scale_(double(flow.height()));
        const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
        const double nominal_width = double(flow.scaled_width());
        const double nominal_cross_section = nominal_width - rounded_corner_loss;
        const double max_width = maximum_physical_width(flow, configured_max_line_width);
        if (path_length > EPS && printable_area_value > EPS && nominal_cross_section > EPS && max_width > nominal_width + EPS) {
            const double target_cross_section = printable_area_value * 1.019 / path_length;
            const double derived_width = std::clamp(
                rounded_corner_loss + target_cross_section,
                nominal_width,
                max_width);
            if (derived_width >= nominal_width * 1.01) {
                const Flow planning_flow(
                    float(unscale<double>(derived_width)),
                    flow.height(),
                    flow.nozzle_diameter());
                LayerPathResult width_candidate = generate_layer_path_candidate(
                    printable_area, planning_flow, configured_max_line_width, 1.0);
                const double planning_cross_section = double(planning_flow.scaled_width()) - rounded_corner_loss;
                if (planning_cross_section > EPS) {
                    const double metadata_scale = planning_cross_section / nominal_cross_section;
                    for (float &multiplier : width_candidate.extrusion_multipliers)
                        multiplier = float(double(multiplier) * metadata_scale);
                }
                const PathValidation width_validation = validate_layer_path(
                    printable_area,
                    flow,
                    width_candidate.path,
                    width_candidate.extrusion_multipliers,
                    configured_max_line_width);
                if (perf_timing_enabled()) {
                    BOOST_LOG_TRIVIAL(warning) << "Continuous Fermat profile derived_width_retry"
                                               << " derived_width_mm=" << unscale<double>(derived_width)
                                               << " planning_width_mm=" << unscale<double>(planning_flow.scaled_width())
                                               << " planning_spacing_mm=" << unscale<double>(planning_flow.scaled_spacing())
                                               << " source_path_length_mm=" << unscale<double>(path_length)
                                               << " candidate_path_points=" << width_candidate.path.points.size()
                                               << " ok=" << width_validation.ok
                                               << " closed=" << width_validation.closed
                                               << " exact_coverage=" << width_validation.exact_coverage_ratio
                                               << " sampled_coverage=" << width_validation.coverage_ratio
                                               << " outside=" << width_validation.outside_ratio
                                               << " material=" << width_validation.material_ratio
                                               << " redeposition=" << width_validation.redeposition_ratio
                                               << " containment=" << width_validation.containment_violations
                                               << " crossings=" << width_validation.crossings
                                               << " spacing=" << width_validation.spacing_violations
                                               << " bead_overlaps=" << width_validation.bead_overlap_violations
                                               << " turnbacks=" << width_validation.turnback_violations
                                               << " reason=\"" << width_validation.reason << '"';
                }
                if (width_validation.ok)
                    return width_candidate;

                const auto underfill_deficit = [](const PathValidation &validation) {
                    return std::max(0.0, 0.980 - validation.exact_coverage_ratio) +
                           std::max(0.0, 0.980 - validation.material_ratio);
                };
                const bool safety_no_worse = width_validation.containment_violations == 0 &&
                    width_validation.crossings == 0 && width_validation.turnback_violations == 0 &&
                    width_validation.outside_ratio <= std::max(0.020, nominal_validation.outside_ratio) + 1e-9 &&
                    width_validation.redeposition_ratio <= std::max(0.040, nominal_validation.redeposition_ratio) + 1e-9;
                if (safety_no_worse &&
                    underfill_deficit(width_validation) + 1e-7 < underfill_deficit(nominal_validation)) {
                    alternatives.emplace_back(std::move(width_candidate));
                    alternative_validations.emplace_back(width_validation);
                    best = &alternatives.back();
                    best_validation = &alternative_validations.back();
                }
            }
        }

    }
    if (nominal_validation.emittable)
        return std::move(*best);
    for (const double spacing_factor : spacing_factors) {
        alternatives.emplace_back(generate_layer_path_candidate(
            printable_area, flow, configured_max_line_width, spacing_factor));
        alternative_validations.emplace_back(validate_layer_path(
            printable_area,
            flow,
            alternatives.back().path,
            alternatives.back().extrusion_multipliers,
            configured_max_line_width));
        if (alternative_validations.back().ok)
            return std::move(alternatives.back());

        const PathValidation &candidate_validation = alternative_validations.back();
        if (hard_failures(candidate_validation) < hard_failures(*best_validation) ||
            (hard_failures(candidate_validation) == hard_failures(*best_validation) &&
             candidate_validation.exact_coverage_ratio > best_validation->exact_coverage_ratio) ||
            (hard_failures(candidate_validation) == hard_failures(*best_validation) &&
             candidate_validation.exact_coverage_ratio == best_validation->exact_coverage_ratio &&
             candidate_validation.redeposition_ratio < best_validation->redeposition_ratio)) {
            best = &alternatives.back();
            best_validation = &alternative_validations.back();
        }
    }

    const size_t max_rebalance_iterations = has_holes ? 4 : 10;
    for (size_t iteration = 0; iteration < max_rebalance_iterations; ++iteration) {
        LayerPathResult adjusted = *best;
        if (!rebalance_widths_toward_uncovered_area(
                printable_area,
                adjusted.path.points,
                flow,
                maximum_extrusion_multiplier(flow, configured_max_line_width),
                maximum_physical_width(flow, configured_max_line_width),
                adjusted.extrusion_multipliers))
            break;
        const PathValidation adjusted_validation = validate_layer_path(
            printable_area, flow, adjusted.path, adjusted.extrusion_multipliers, configured_max_line_width);
        if (adjusted_validation.ok)
            return adjusted;
        if (hard_failures(adjusted_validation) > hard_failures(*best_validation) ||
            adjusted_validation.exact_coverage_ratio <= best_validation->exact_coverage_ratio + 1e-7)
            break;
        alternatives.emplace_back(std::move(adjusted));
        alternative_validations.emplace_back(adjusted_validation);
        best = &alternatives.back();
        best_validation = &alternative_validations.back();
    }

    return std::move(*best);
}

} // namespace

Polyline generate_layer_path(
    const ExPolygons &printable_area,
    const Flow &flow,
    const double max_line_width)
{
    return generate_layer_path_result(printable_area, flow, max_line_width).path;
}

GeneratedPath generate_layer_path_with_metadata(
    const ExPolygons &printable_area,
    const Flow &flow,
    const double max_line_width)
{
    LayerPathResult generated = generate_layer_path_result(printable_area, flow, max_line_width);
    return { std::move(generated.path), std::move(generated.extrusion_multipliers) };
}

PathValidation validate_layer_path(
    const ExPolygons &printable_area,
    const Flow &flow,
    const Polyline &path,
    const std::vector<float> &extrusion_multipliers,
    const double max_line_width)
{
    PathValidation validation;
    const Points &points = path.points;
    validation.closed = points.size() >= 4 && points.front() == points.back();
    if (!validation.closed) {
        validation.reason = "path is not exactly closed";
        return validation;
    }

    const double printable_area_value = std::abs(area(printable_area));
    if (!(printable_area_value > EPS) || !std::isfinite(printable_area_value)) {
        validation.reason = "printable area is empty or invalid";
        return validation;
    }

    const double line_width = double(flow.scaled_width());
    const double spacing = double(flow.scaled_spacing());
    if (!(line_width > EPS) || !(spacing > EPS) || !std::isfinite(line_width) || !std::isfinite(spacing)) {
        validation.reason = "extrusion width or spacing is invalid";
        return validation;
    }

    // Offset independent open segments and union the capsules.  Passing a
    // self-touching, exactly closed CFS polyline to ClipperOffset as one open
    // path may collapse alternating loops and drastically under-report the
    // physical footprint.
    if (extrusion_multipliers.size() + 1 != points.size()) {
        validation.reason = "per-segment extrusion metadata cardinality does not match the path";
        return validation;
    }
    const double line_height = scale_(double(flow.height()));
    const double rounded_corner_loss = line_height * (1.0 - 0.25 * PI);
    const double nominal_cross_section = line_width - rounded_corner_loss;
    if (!(nominal_cross_section > EPS)) {
        validation.reason = "nominal extrusion cross-section is invalid";
        return validation;
    }
    const double max_physical_width = maximum_physical_width(flow, max_line_width);
    const double min_physical_width = minimum_physical_width(flow);
    if (!(max_physical_width > EPS) || line_width > max_physical_width + EPS) {
        validation.reason = "nominal extrusion width exceeds the configured or two-nozzle safety limit";
        return validation;
    }
    const double max_multiplier = maximum_extrusion_multiplier(flow, max_line_width);
    if (!(max_multiplier >= 1.0) ||
        !std::all_of(extrusion_multipliers.begin(), extrusion_multipliers.end(), [max_multiplier, &flow](const float multiplier) {
            return std::isfinite(multiplier) && multiplier + 1e-6 >= minimum_extrusion_multiplier(flow) &&
                   multiplier <= max_multiplier + 1e-6;
        })) {
        validation.reason = "per-segment extrusion multiplier is outside the safe range";
        return validation;
    }

    std::map<coord_t, Polylines> segments_by_width;
    std::vector<double> physical_segment_widths;
    physical_segment_widths.reserve(points.size() - 1);
    double deposited_volume = 0.0;
    for (size_t i = 1; i < points.size(); ++i) {
        const double multiplier = double(extrusion_multipliers[i - 1]);
        const coord_t segment_width = coord_t(std::llround(rounded_corner_loss + nominal_cross_section * multiplier));
        physical_segment_widths.emplace_back(double(segment_width));
        if (double(segment_width) + EPS < min_physical_width) {
            validation.reason = "adaptive extrusion width is below the minimum continuous-flow limit";
            return validation;
        }
        if (points[i - 1] == points[i])
            continue;
        deposited_volume += point_distance(points[i - 1], points[i]) * line_height * nominal_cross_section * multiplier;
        if (double(segment_width) > max_physical_width + EPS) {
            validation.reason = "adaptive extrusion width exceeds the configured or two-nozzle safety limit";
            return validation;
        }
        segments_by_width[segment_width].emplace_back(Points { points[i - 1], points[i] });
    }

    Polygons swept;
    for (const auto &[segment_width, segments] : segments_by_width) {
        Polygons group = offset(
            segments,
            float(double(segment_width) * 0.5),
            ClipperLib::jtRound,
            SCALED_RESOLUTION,
            ClipperLib::etOpenRound);
        swept.insert(swept.end(), std::make_move_iterator(group.begin()), std::make_move_iterator(group.end()));
    }
    swept = union_(swept);
    if (swept.empty()) {
        validation.reason = "swept extrusion footprint is empty";
        return validation;
    }
    // From this point the path is closed, its per-segment metadata is
    // complete, every physical width is finite and bounded, and it produces a
    // non-empty swept footprint. The remaining checks describe geometric
    // fidelity. They are important diagnostics, but do not make the motion
    // stream malformed or unbounded.
    validation.emittable = true;

    const ExPolygons uncovered = diff_ex(printable_area, swept);
    const double uncovered_area = std::max(0.0, std::abs(area(uncovered)));
    validation.exact_coverage_ratio = std::clamp(1.0 - uncovered_area / printable_area_value, 0.0, 1.0);

    const ExPolygons swept_union = union_ex(swept);
    const double swept_area_value = std::max(0.0, std::abs(area(swept_union)));
    const double deposited_area_value = deposited_volume / line_height;
    validation.redeposition_ratio = std::max(0.0, deposited_area_value - swept_area_value) / printable_area_value;
    const ExPolygons outside = diff_ex(swept_union, printable_area);
    validation.outside_ratio = std::max(0.0, std::abs(area(outside)) / printable_area_value);
    validation.material_ratio = deposited_volume / (printable_area_value * line_height);

    // Retain a coarse, independent 70-cell sample for diagnostics. It is not
    // an acceptance gate: cell-center sampling has a phase-dependent boundary
    // error and was observed to alternate between pass and failure on
    // coordinate-jittered copies of the same exact cross-section. The exact
    // polygon-area coverage above is the authoritative coverage contract.
    const BoundingBox bounds = get_extents(printable_area);
    const double bounds_width = double(bounds.max.x() - bounds.min.x());
    const double bounds_height = double(bounds.max.y() - bounds.min.y());
    const int coverage_cells = 70;
    const int nx = bounds_width >= bounds_height ? coverage_cells :
        std::max(12, int(std::lround(double(coverage_cells) * bounds_width / std::max(bounds_height, 1.0))));
    const int ny = bounds_height > bounds_width ? coverage_cells :
        std::max(12, int(std::lround(double(coverage_cells) * bounds_height / std::max(bounds_width, 1.0))));
    size_t inside_cells = 0;
    size_t covered_cells = 0;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const Point sample(
                double(bounds.min.x()) + (double(x) + 0.5) * bounds_width / double(nx),
                double(bounds.min.y()) + (double(y) + 0.5) * bounds_height / double(ny));
            if (!printable_area_contains(printable_area, sample))
                continue;
            ++inside_cells;
            if (printable_area_contains(swept_union, sample))
                ++covered_cells;
        }
    }
    validation.coverage_ratio = inside_cells == 0 ? 0.0 : double(covered_cells) / double(inside_cells);

    validation.containment_violations = size_t(count_containment_violations(printable_area, points, spacing));
    const PairMetrics pairs = path_pair_metrics(points, spacing, &physical_segment_widths, line_width);
    validation.crossings = pairs.crossings;
    validation.spacing_violations = pairs.close_pairs;
    validation.bead_overlap_violations = pairs.bead_overlaps;
    std::string first_turnback;
    validation.turnback_violations = count_turnback_violations(points, line_width, &first_turnback);

    // A cell-center sample is not an area integral. Keep it visible in
    // diagnostics and corpus output, while exact polygon clipping alone gates
    // swept-area coverage.
    static constexpr double BASE_MIN_EXACT_COVERAGE = 0.980;
    static constexpr double MAX_OUTSIDE_RATIO = 0.020;
    static constexpr double BASE_MAX_REDEPOSITION_RATIO = 0.040;
    static constexpr double MIN_MATERIAL_RATIO = 0.980;
    static constexpr double MAX_MATERIAL_RATIO = 1.020;
    // Keep the quality targets fixed. Thin geometry is categorized separately
    // below instead of silently relaxing coverage or material thresholds.
    const double min_exact_coverage = BASE_MIN_EXACT_COVERAGE;
    const double max_redeposition_ratio = BASE_MAX_REDEPOSITION_RATIO;
    const double min_material_ratio = MIN_MATERIAL_RATIO;
    std::ostringstream reason;
    if (validation.containment_violations != 0)
        reason << " centerline_containment=" << validation.containment_violations;
    if (validation.crossings != 0)
        reason << " crossings=" << validation.crossings;
    if (validation.redeposition_ratio > max_redeposition_ratio + 1e-9)
        reason << " redeposition=" << validation.redeposition_ratio << ">" << max_redeposition_ratio
               << " bead_overlap_pairs=" << validation.bead_overlap_violations;
    if (validation.turnback_violations != 0)
        reason << " turnbacks=" << validation.turnback_violations << "(" << first_turnback << ")";
    if (validation.exact_coverage_ratio + 1e-9 < min_exact_coverage)
        reason << " exact_coverage=" << validation.exact_coverage_ratio << "<" << min_exact_coverage;
    if (validation.outside_ratio > MAX_OUTSIDE_RATIO + 1e-9)
        reason << " outside=" << validation.outside_ratio << ">" << MAX_OUTSIDE_RATIO;
    if (validation.material_ratio + 1e-9 < min_material_ratio)
        reason << " material=" << validation.material_ratio << "<" << min_material_ratio;
    if (validation.material_ratio > MAX_MATERIAL_RATIO + 1e-9)
        reason << " material=" << validation.material_ratio << ">" << MAX_MATERIAL_RATIO;

    validation.reason = reason.str();
    validation.ok = validation.reason.empty();
    if (!validation.ok) {
        // A half-threshold erosion answers whether a 2.5-line-wide disk can
        // move through the section. This detects globally thin sections and
        // local necks without pre-rejecting thin geometry that generated a
        // satisfactory path.
        const double threshold = 2.5 * line_width;
        const ExPolygons eroded = offset_ex(printable_area, float(-0.5 * threshold));
        const double dust_area = line_width * line_width;
        size_t significant_components = 0;
        size_t eroded_holes = 0;
        for (const ExPolygon &component : eroded) {
            if (std::abs(component.area()) < dust_area)
                continue;
            ++significant_components;
            eroded_holes += component.holes.size();
        }
        size_t source_holes = 0;
        for (const ExPolygon &component : printable_area)
            source_holes += component.holes.size();

        validation.thin_feature_threshold_mm = unscale<double>(threshold);
        if (significant_components == 0)
            validation.thin_feature_class = "section_below_2.5_line_widths";
        else if (significant_components > printable_area.size())
            validation.thin_feature_class = "neck_below_2.5_line_widths";
        else if (eroded_holes < source_holes)
            validation.thin_feature_class = "hole_web_below_2.5_line_widths";
    }
    return validation;
}

bool apply_to_layer(Layer &layer)
{
    const PrintConfig &config = layer.object()->print()->config();
    if (!config.spiral_hybrid_non_crossing.value)
        return false;

    auto fail = [](const std::string &reason) -> bool {
        const std::string message = "Continuous Fermat path generation failed: " + reason;
        BOOST_LOG_TRIVIAL(error) << message;
        throw Slic3r::SlicingError(message);
    };

    if (layer.empty()) {
        BOOST_LOG_TRIVIAL(debug) << "Continuous Fermat active on empty layer " << layer.id() << "; skipping normal fill generation.";
        return true;
    }

    LayerRegion *target_region = nullptr;
    size_t non_empty_regions = 0;
    for (LayerRegion *region : layer.regions()) {
        if (region != nullptr && !region->slices.empty()) {
            ++non_empty_regions;
            if (target_region == nullptr)
                target_region = region;
        }
    }
    if (target_region == nullptr)
        return fail("no layer region has printable slices; " + layer_diagnostics(layer, {}));
    if (non_empty_regions > 1)
        return fail("multiple print regions are not supported in strict continuous mode; " + layer_diagnostics(layer, {}));

    ExPolygons printable_area = layer.lslices.empty() ? to_expolygons(target_region->slices.surfaces) : layer.lslices;
    if (printable_area.empty())
        return fail("printable area is empty; " + layer_diagnostics(layer, printable_area));
    if (printable_area.size() != 1)
        return fail("disconnected islands are not supported in strict continuous mode; " + layer_diagnostics(layer, printable_area));
    const BoundingBox printable_bounds = get_extents(printable_area);
    const coord_t min_printable_span = std::min(
        printable_bounds.max.x() - printable_bounds.min.x(),
        printable_bounds.max.y() - printable_bounds.min.y());
    if (min_printable_span < scale_(0.002)) {
        BOOST_LOG_TRIVIAL(debug) << "Continuous Fermat skipping a degenerate zero-area apex on layer " << layer.id() << '.';
        return true;
    }
    // A topologically connected island may contain holes. Every generated
    // path still passes exact centerline containment and swept-footprint
    // validation, including containment against each hole boundary.

    // Initial-layer line width is a conventional adhesion override.  It
    // would change the centerline phase only on layer zero, so continuous
    // mode deliberately uses the normal outer-wall width on every layer.
    // The GUI mirrors the two settings to make this override explicit, while
    // this core rule also protects loaded projects and non-GUI slicing.
    const Flow flow = target_region->region().flow(
        *layer.object(), frExternalPerimeter, layer.height, false);

    // A usable part has one island per layer and every layer is supported by
    // non-empty geometric overlap with the layer below. Cross-sections may
    // otherwise translate, rotate, taper, grow, shrink, or change outline.
    if (layer.lower_layer != nullptr) {
        ExPolygons lower_area = layer.lower_layer->lslices;
        if (lower_area.empty()) {
            for (const LayerRegion *region : layer.lower_layer->regions())
                if (region != nullptr)
                    append(lower_area, to_expolygons(region->slices.surfaces));
        }
        const ExPolygons support_overlap = intersection_ex(printable_area, lower_area);
        if (support_overlap.empty() || !(std::abs(area(support_overlap)) > EPS))
            return fail("layer has no geometric overlap with the layer below; " +
                        layer_diagnostics(layer, printable_area));

        const LayerRegion *lower_region = nullptr;
        for (const LayerRegion *region : layer.lower_layer->regions())
            if (region != nullptr && !region->slices.empty()) {
                if (lower_region != nullptr)
                    return fail("lower layer has multiple printable regions; " + layer_diagnostics(layer, printable_area));
                lower_region = region;
            }
        if (lower_region == nullptr)
            return fail("lower layer has no printable region; " + layer_diagnostics(layer, printable_area));
        const Flow lower_flow = lower_region->region().flow(
            *layer.lower_layer->object(), frExternalPerimeter, layer.lower_layer->height, false);
        if (flow.scaled_width() != lower_flow.scaled_width() ||
            flow.scaled_spacing() != lower_flow.scaled_spacing() ||
            std::abs(double(flow.height()) - double(lower_flow.height())) > EPSILON ||
            std::abs(double(flow.nozzle_diameter()) - double(lower_flow.nozzle_diameter())) > EPSILON)
            return fail("layer extrusion flow differs from the layer below; " + layer_diagnostics(layer, printable_area, &flow));
    } else if (layer.id() != 0) {
        return fail("non-first layer has no lower model layer; " + layer_diagnostics(layer, printable_area));
    } else if (!std::isfinite(layer.bottom_z()) || std::abs(layer.bottom_z()) > EPSILON) {
        return fail("first layer is not in direct contact with the build plate; " + layer_diagnostics(layer, printable_area));
    }

    const double max_line_width = config.continuous_max_line_width.get_abs_value(flow.nozzle_diameter());
    LayerPathResult generated_path = generate_layer_path_result(printable_area, flow, max_line_width);
    if (generated_path.path.points.size() < 2 && layer.upper_layer == nullptr) {
        BOOST_LOG_TRIVIAL(debug) << "Continuous Fermat omitting a nongeneratable zero-area terminal apex on layer " << layer.id() << '.';
        return true;
    }
    if (generated_path.path.points.size() < 2)
        return fail("generated path is empty or degenerate; " + layer_diagnostics(layer, printable_area, &flow));

    const PathValidation validation =
        validate_layer_path(
            printable_area, flow, generated_path.path, generated_path.extrusion_multipliers, max_line_width);
    if (!validation.emittable) {
        std::ostringstream reason;
        reason << "final path is not structurally emittable:" << validation.reason
               << " coverage=" << validation.coverage_ratio << " exact_coverage=" << validation.exact_coverage_ratio
               << " outside=" << validation.outside_ratio
               << " material=" << validation.material_ratio
               << "; " << layer_diagnostics(layer, printable_area, &flow);
        return fail(reason.str());
    }

    std::string validation_warning;
    if (!validation.ok) {
        std::ostringstream warning;
        warning << "layer=" << layer.id() << " z=" << layer.print_z
                << " geometric_validation=" << validation.reason
                << " coverage=" << validation.coverage_ratio
                << " exact_coverage=" << validation.exact_coverage_ratio
                << " outside=" << validation.outside_ratio
                << " material=" << validation.material_ratio
                << " redeposition=" << validation.redeposition_ratio;
        if (!validation.thin_feature_class.empty())
            warning << " thin_feature=" << validation.thin_feature_class
                    << " thin_threshold_mm=" << validation.thin_feature_threshold_mm;
        validation_warning = warning.str();
        BOOST_LOG_TRIVIAL(debug) << "Continuous Fermat emitted with an advisory validation warning: "
                                 << validation_warning;
        layer.object()->add_continuous_slicing_validation_warning(
            L("Continuous slicing generated one or more layers outside its geometric quality targets. "
              "The paths remain available for preview and export; inspect the preview and the "
              "_CONTINUOUS_FERMAT_VALIDATION_WARNING comments before deciding whether to use them."));
    }

    for (LayerRegion *region : layer.regions()) {
        region->perimeters.clear();
        region->fills.clear();
        region->thin_fills.clear();
    }

    ExtrusionPath extrusion(erExternalPerimeter, flow.mm3_per_mm(), flow.width(), flow.height());
    extrusion.polyline = std::move(generated_path.path);
    extrusion.continuous_fermat_extrusion_multipliers = std::move(generated_path.extrusion_multipliers);
    extrusion.continuous_fermat_validation_warning = std::move(validation_warning);
    extrusion.set_continuous_fermat();
    extrusion.set_reverse();

    auto *collection = new ExtrusionEntityCollection();
    collection->no_sort = true;
    collection->entities.emplace_back(new ExtrusionPath(std::move(extrusion)));
    target_region->perimeters.entities.emplace_back(collection);

    BOOST_LOG_TRIVIAL(info) << "Continuous Fermat replaced normal layer entities: " << layer_diagnostics(layer, printable_area, &flow)
                            << " path_points=" << static_cast<const ExtrusionPath*>(collection->entities.front())->polyline.points.size();
    return true;
}

} // namespace ContinuousFermat
} // namespace Slic3r
