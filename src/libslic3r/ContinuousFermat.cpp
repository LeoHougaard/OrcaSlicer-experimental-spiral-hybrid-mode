#include "ContinuousFermat.hpp"

#include "ClipperUtils.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Layer.hpp"
#include "Print.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace Slic3r {
namespace ContinuousFermat {
namespace {

static constexpr double EPS = 1e-9;
static constexpr double PI = 3.141592653589793238462643383279502884;

struct ContourLoop
{
    size_t level_index { 0 };
    Points points;
    double area { 0.0 };
    double length { 0.0 };
    Point centroid;
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

std::vector<ContourLoop> generate_offset_contours(const ExPolygons &printable_area, const Flow &flow)
{
    const double spacing = double(flow.scaled_spacing());
    const double first_offset = double(flow.scaled_width()) * 0.5;
    const double min_area = spacing * spacing * 8.0;
    const double min_length = spacing * 5.0;

    std::vector<ContourLoop> contours;
    ExPolygons source = union_ex(printable_area);
    for (size_t level = 0; level < 512; ++level) {
        const double offset = first_offset + double(level) * spacing;
        ExPolygons inset = offset_ex(source, float(-offset), ClipperLib::jtMiter, 3.0);
        if (inset.empty())
            break;

        for (const ExPolygon &expoly : inset) {
            auto add_polygon = [&](const Polygon &polygon) {
                if (polygon.points.size() < 4)
                    return;
                ContourLoop loop;
                loop.level_index = level;
                loop.points = resample_closed_loop(polygon.points, std::max(spacing * 0.45, 1.0));
                loop.area = std::abs(double(polygon.area()));
                loop.length = polyline_length(loop.points, true);
                loop.centroid = polygon_centroid_or_first(loop.points);
                if (loop.area >= min_area && loop.length >= min_length)
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

Points build_single_minimum_connected_fermat(
    std::vector<ContourLoop> loops,
    const Point &start_anchor,
    const double spacing,
    const Point &exit_anchor,
    const bool preserve_medial_pockets = true)
{
    loops.erase(std::remove_if(loops.begin(), loops.end(), [](const ContourLoop &loop) { return loop.points.size() < 4; }), loops.end());
    if (loops.empty())
        return {};

    const double port_spacing = spacing * 2.5;
    if (loops.size() == 1) {
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
            const Point &prev_in = in_branch.back();
            const Point &prev_out = out_branch.back();
            size_t medial_index = 0;
            double best = -1.0;
            for (size_t i = 0; i < n; ++i) {
                const double score = std::min(point_distance(points[i], prev_in), point_distance(points[i], prev_out));
                if (score > best) {
                    best = score;
                    medial_index = i;
                }
            }
            if (subset_cycle_param(in_index, out_index, double(medial_index), true, false))
                append_loop_vertices_forward(in_branch, points, size_t(std::ceil(in_index)) % n, size_t(std::floor(out_index)) % n);
            else
                append_loop_vertices_forward(out_branch, points, size_t(std::ceil(out_index)) % n, size_t(std::floor(in_index)) % n);
            break;
        }

        first_circle = false;
    }

    for (auto it = out_branch.rbegin(); it != out_branch.rend(); ++it)
        append_point(in_branch, *it);
    return in_branch;
}

struct PairMetrics
{
    int crossings { 0 };
    int close_pairs { 0 };
    double min_spacing { std::numeric_limits<double>::infinity() };
};

PairMetrics path_pair_metrics(const Points &path, const double spacing)
{
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

    const double min_allowed = spacing * 0.90;
    const double bin_size = std::max(spacing * 1.5, 1.0);
    const double local_skip_distance = spacing * 4.0;
    const double path_length = prefix.back();
    const bool closed_path = point_distance(path.front(), path.back()) <= spacing * 0.20;
    std::unordered_map<long long, std::vector<size_t>> bins;

    auto key = [&](const double x, const double y) {
        const long long bx = (long long)std::floor(x / bin_size);
        const long long by = (long long)std::floor(y / bin_size);
        return (bx << 32) ^ (by & 0xffffffffLL);
    };

    for (size_t idx = 0; idx < segments.size(); ++idx) {
        const Segment &seg = segments[idx];
        const double min_x = std::min(seg.a.x(), seg.b.x()) - min_allowed;
        const double min_y = std::min(seg.a.y(), seg.b.y()) - min_allowed;
        const double max_x = std::max(seg.a.x(), seg.b.x()) + min_allowed;
        const double max_y = std::max(seg.a.y(), seg.b.y()) + min_allowed;
        const long long bx0 = (long long)std::floor(min_x / bin_size);
        const long long by0 = (long long)std::floor(min_y / bin_size);
        const long long bx1 = (long long)std::floor(max_x / bin_size);
        const long long by1 = (long long)std::floor(max_y / bin_size);
        for (long long by = by0; by <= by1; ++by)
            for (long long bx = bx0; bx <= bx1; ++bx)
                bins[(bx << 32) ^ (by & 0xffffffffLL)].push_back(idx);
    }

    std::unordered_set<unsigned long long> checked;
    for (size_t i = 0; i < segments.size(); ++i) {
        const Segment &seg = segments[i];
        const double min_x = std::min(seg.a.x(), seg.b.x()) - min_allowed;
        const double min_y = std::min(seg.a.y(), seg.b.y()) - min_allowed;
        const double max_x = std::max(seg.a.x(), seg.b.x()) + min_allowed;
        const double max_y = std::max(seg.a.y(), seg.b.y()) + min_allowed;
        const long long bx0 = (long long)std::floor(min_x / bin_size);
        const long long by0 = (long long)std::floor(min_y / bin_size);
        const long long bx1 = (long long)std::floor(max_x / bin_size);
        const long long by1 = (long long)std::floor(max_y / bin_size);

        std::vector<size_t> candidates;
        for (long long by = by0; by <= by1; ++by)
            for (long long bx = bx0; bx <= bx1; ++bx)
                if (auto it = bins.find((bx << 32) ^ (by & 0xffffffffLL)); it != bins.end())
                    candidates.insert(candidates.end(), it->second.begin(), it->second.end());
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

        for (const size_t j : candidates) {
            if (j <= i)
                continue;
            size_t index_distance = j - i;
            if (closed_path)
                index_distance = std::min(index_distance, segments.size() - index_distance);
            if (index_distance <= 6)
                continue;

            double path_gap = std::max(0.0, prefix[j] - prefix[i + 1]);
            if (closed_path) {
                const double occupied_span = prefix[j + 1] - prefix[i];
                path_gap = std::min(path_gap, std::max(0.0, path_length - occupied_span));
            }
            if (path_gap <= local_skip_distance)
                continue;
            const unsigned long long pair_key = (unsigned long long(i) << 32) ^ unsigned(j);
            if (!checked.insert(pair_key).second)
                continue;

            const double d = segment_distance(seg.a, seg.b, segments[j].a, segments[j].b);
            metrics.min_spacing = std::min(metrics.min_spacing, d);
            if (d <= 1e-7)
                ++metrics.crossings;
            else if (d < min_allowed)
                ++metrics.close_pairs;
        }
    }

    return metrics;
}

Points complete_outer_boundary_cycle(const Points &path, const ContourLoop &outer_loop, const double spacing)
{
    if (path.size() < 2 || outer_loop.points.size() < 4)
        return path;

    const double target_gap = spacing * 0.12;
    if (point_distance(path.front(), path.back()) <= target_gap * 1.5)
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
        const double stop_index = direction > 0 ?
            loop_back_by_distance(outer_loop.points, start_index, target_gap) :
            loop_forward_by_distance(outer_loop.points, start_index, target_gap);
        Points arc = loop_arc_between_params(outer_loop.points, end_index, stop_index, direction);
        if (arc.size() < 2)
            continue;

        arc.front() = path.back();
        Points candidate;
        candidate.reserve(path.size() + arc.size());
        append_points(candidate, path);
        for (size_t i = 1; i < arc.size(); ++i)
            append_point(candidate, arc[i]);

        if (point_distance(candidate.front(), candidate.back()) > spacing * 0.20)
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
    const bool preserve_medial_pockets = true)
{
    const std::vector<double> exit_fractions { 0.50, 0.67, 0.33, 0.25, 0.75, 0.10, 0.90 };
    Points best_path;
    PairMetrics best_metrics;
    size_t best_loop_count = 0;
    bool have_best = false;

    auto better_than_best = [&](const PairMetrics &metrics, const size_t loop_count) {
        return !have_best || metrics.crossings < best_metrics.crossings ||
            (metrics.crossings == best_metrics.crossings && metrics.close_pairs < best_metrics.close_pairs) ||
            (metrics.crossings == best_metrics.crossings && metrics.close_pairs == best_metrics.close_pairs &&
             loop_count > best_loop_count) ||
            (metrics.crossings == best_metrics.crossings && metrics.close_pairs == best_metrics.close_pairs &&
             loop_count == best_loop_count && metrics.min_spacing > best_metrics.min_spacing);
    };

    for (size_t loop_count = ordered.size(); loop_count >= 1; --loop_count) {
        const std::vector<ContourLoop> candidate_loops(ordered.begin(), ordered.begin() + loop_count);
        for (const double fraction : exit_fractions) {
            const Point exit_anchor = point_at_closed_fraction(candidate_loops.front().points, fraction);
            Points candidate = build_single_minimum_connected_fermat(
                candidate_loops, start_anchor, spacing, exit_anchor, preserve_medial_pockets);
            if (candidate.size() < 2)
                continue;

            PairMetrics metrics = path_pair_metrics(candidate, spacing);
            if (better_than_best(metrics, loop_count)) {
                best_path = std::move(candidate);
                best_metrics = metrics;
                best_loop_count = loop_count;
                have_best = true;
            }
            if (metrics.crossings == 0 && metrics.close_pairs == 0)
                return best_path;
        }

        if (loop_count == ordered.size() && best_metrics.crossings == 0 && best_metrics.close_pairs == 0)
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

Points merge_child_spiral(const Points &parent_path, Points child_path)
{
    if (parent_path.empty())
        return child_path;
    if (child_path.empty())
        return parent_path;

    OpenProjection start_projection = project_open_polyline(parent_path, child_path.front());
    OpenProjection end_projection = project_open_polyline(parent_path, child_path.back());
    if (start_projection.index > end_projection.index) {
        std::reverse(child_path.begin(), child_path.end());
        start_projection = project_open_polyline(parent_path, child_path.front());
        end_projection = project_open_polyline(parent_path, child_path.back());
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

Points build_branch_connected_fermat(const std::vector<ContourLoop> &contours, const Point &start_anchor, const double spacing)
{
    const auto grouped = loops_by_level(contours);
    if (grouped.empty())
        return {};

    bool have_split_level = false;
    size_t split_level = 0;
    for (const auto &[level, loops] : grouped) {
        if (loops.size() > 1) {
            split_level = level;
            have_split_level = true;
            break;
        }
    }

    if (!have_split_level)
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    std::vector<ContourLoop> trunk;
    for (const auto &[level, loops] : grouped) {
        if (level >= split_level)
            break;
        if (loops.size() == 1)
            trunk.emplace_back(loops.front());
    }

    if (trunk.empty())
        return build_best_single_chain_path(ordered_loops_for_spiral(contours), start_anchor, spacing);

    const Point trunk_exit = point_at_closed_fraction(trunk.front().points, 0.5);
    Points parent_path = build_single_minimum_connected_fermat(trunk, start_anchor, spacing, trunk_exit);
    if (parent_path.size() < 2)
        return parent_path;

    std::vector<ContourLoop> split_loops = grouped.at(split_level);
    std::sort(split_loops.begin(), split_loops.end(), [](const ContourLoop &a, const ContourLoop &b) {
        if (a.centroid.x() != b.centroid.x())
            return a.centroid.x() < b.centroid.x();
        return a.centroid.y() < b.centroid.y();
    });

    std::vector<std::vector<ContourLoop>> chains;
    chains.reserve(split_loops.size());
    for (const ContourLoop &loop : split_loops)
        chains.push_back({ loop });

    for (const auto &[level, loops] : grouped) {
        if (level <= split_level)
            continue;

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
            if (best_index != size_t(-1)) {
                used[best_index] = true;
                chain.emplace_back(loops[best_index]);
            }
        }
    }

    for (const std::vector<ContourLoop> &chain : chains) {
        if (chain.empty() || chain.front().points.empty())
            continue;

        const size_t stride = std::max<size_t>(1, chain.front().points.size() / 48);
        Point anchor = chain.front().points.front();
        double best_distance = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < chain.front().points.size(); i += stride) {
            const OpenProjection projection = project_open_polyline(parent_path, chain.front().points[i]);
            if (projection.distance < best_distance) {
                best_distance = projection.distance;
                anchor = chain.front().points[i];
            }
        }

        const Point child_exit = point_at_closed_fraction(chain.front().points, 0.5);
        Points child_path = build_single_minimum_connected_fermat(chain, anchor, spacing, child_exit);
        parent_path = merge_child_spiral(parent_path, std::move(child_path));
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

int count_containment_violations(const ExPolygons &printable_area, const Points &path, const double spacing)
{
    int violations = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        const Point &a = path[i - 1];
        const Point &b = path[i];
        const double length = point_distance(a, b);
        const size_t samples = std::max<size_t>(2, size_t(std::ceil(length / std::max(spacing * 0.25, 1.0))));
        for (size_t sample = 0; sample <= samples; ++sample) {
            const Point p = lerp_point(a, b, double(sample) / double(samples));
            if (!printable_area_contains(printable_area, p))
                ++violations;
        }
    }
    return violations;
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
    const double half_width = spacing * 1.35;
    for (double center_index : unique_port_candidates(std::move(port_candidates), 2)) {
        center_index = std::clamp(center_index, 3.0, double(parent_path.size() - 4));
        const double center_length = open_path_length_at_index(parent_path, prefix, center_index);
        const double left_index = open_path_index_at_length(parent_path, prefix, center_length - half_width);
        const double right_index = open_path_index_at_length(parent_path, prefix, center_length + half_width);
        const Point left_port = open_path_point_at_index(parent_path, left_index);
        const Point right_port = open_path_point_at_index(parent_path, right_index);
        if (point_distance(left_port, right_port) < spacing * 0.5)
            continue;

        for (const bool reverse_anchors : { false, true }) {
            const Point &start_anchor = reverse_anchors ? right_port : left_port;
            const Point &exit_anchor = reverse_anchors ? left_port : right_port;
            Points child_path = build_single_minimum_connected_fermat(
                ordered, start_anchor, spacing, exit_anchor, true);

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
                if (metrics.crossings == 0 && metrics.close_pairs == 0 &&
                    count_containment_violations(printable_area, candidate, spacing) == 0)
                    return candidate;
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

} // namespace

Polyline generate_layer_path(const ExPolygons &printable_area, const Flow &flow)
{
    const std::vector<ContourLoop> contours = generate_offset_contours(printable_area, flow);
    if (contours.empty())
        return {};

    const auto grouped = loops_by_level(contours);
    if (grouped.empty())
        return {};

    const std::vector<ContourLoop> ordered = ordered_loops_for_spiral(contours);
    if (ordered.empty())
        return {};

    const size_t first_level_count = grouped.begin()->second.size();
    bool multi_loop = false;
    for (const auto &[level, loops] : grouped)
        multi_loop |= loops.size() > 1;

    const double spacing = double(flow.scaled_spacing());
    const bool branch_single_island = multi_loop && first_level_count == 1;
    const bool multi_hole = collect_holes(printable_area).size() >= 2;
    const bool one_hole_ring = multi_loop && first_level_count == 2 &&
        point_distance(grouped.begin()->second[0].centroid, grouped.begin()->second[1].centroid) < spacing * 2.0;
    const Point start_anchor = point_at_closed_fraction(ordered.front().points, 0.0);

    Points best_path;
    if (multi_hole) {
        best_path = build_merged_hole_connected_fermat(contours, printable_area, flow, start_anchor);
    } else if (branch_single_island) {
        best_path = build_branch_connected_fermat(contours, start_anchor, spacing);
    } else {
        best_path = build_best_single_chain_path(ordered, start_anchor, spacing, !one_hole_ring);
    }

    if (best_path.size() < 2)
        return {};

    const ContourLoop *outer_loop = &ordered.front();
    for (const ContourLoop &loop : contours) {
        if (loop.level_index == 0 && loop.area > outer_loop->area)
            outer_loop = &loop;
    }
    best_path = complete_outer_boundary_cycle(best_path, *outer_loop, spacing);
    if (multi_hole)
        best_path = insert_uncovered_pocket_spirals(printable_area, contours, std::move(best_path), spacing);

    return Polyline(std::move(best_path));
}

bool apply_to_layer(Layer &layer)
{
    const PrintConfig &config = layer.object()->print()->config();
    if (!config.spiral_mode.value || !config.spiral_hybrid_non_crossing.value)
        return false;

    LayerRegion *target_region = nullptr;
    for (LayerRegion *region : layer.regions()) {
        if (region != nullptr && !region->slices.empty()) {
            target_region = region;
            break;
        }
    }
    if (target_region == nullptr)
        return false;

    ExPolygons printable_area = layer.lslices.empty() ? to_expolygons(target_region->slices.surfaces) : layer.lslices;
    if (printable_area.empty())
        return false;

    const Flow flow = target_region->flow(frExternalPerimeter);
    Polyline path = generate_layer_path(printable_area, flow);
    if (path.points.size() < 2)
        return false;

    for (LayerRegion *region : layer.regions()) {
        region->perimeters.clear();
        region->fills.clear();
        region->thin_fills.clear();
    }

    ExtrusionPath extrusion(erExternalPerimeter, flow.mm3_per_mm(), flow.width(), flow.height());
    extrusion.polyline = std::move(path);

    auto *collection = new ExtrusionEntityCollection();
    collection->entities.emplace_back(new ExtrusionPath(std::move(extrusion)));
    target_region->perimeters.entities.emplace_back(collection);
    return true;
}

} // namespace ContinuousFermat
} // namespace Slic3r
