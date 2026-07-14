#include <catch2/catch_all.hpp>

#include "libslic3r/ContinuousFermat.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Line.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>

using namespace Slic3r;

namespace {

Point p(double x, double y)
{
    return Point::new_scale(x, y);
}

double us(coord_t value)
{
    return unscale<double>(value);
}

Points circle_points(double radius, size_t count, Point center = Point(0, 0), bool ccw = true)
{
    Points points;
    points.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        const double angle = 2.0 * PI * double(i) / double(count);
        const double oriented_angle = ccw ? angle : -angle;
        points.emplace_back(Point::new_scale(
            us(center.x()) + radius * std::cos(oriented_angle),
            us(center.y()) + radius * std::sin(oriented_angle)));
    }
    return points;
}

Points star_points(double outer_radius, double inner_radius, size_t tips)
{
    Points points;
    points.reserve(tips * 2);
    for (size_t i = 0; i < tips * 2; ++i) {
        const double radius = (i % 2 == 0) ? outer_radius : inner_radius;
        const double angle = PI * 0.5 + PI * double(i) / double(tips);
        points.emplace_back(Point::new_scale(radius * std::cos(angle), radius * std::sin(angle)));
    }
    return points;
}

bool segments_cross_nonlocal(const Polyline &path, std::string *message = nullptr)
{
    const Lines lines = path.lines();
    std::vector<double> prefix_length(lines.size() + 1, 0.0);
    for (size_t i = 0; i < lines.size(); ++i)
        prefix_length[i + 1] = prefix_length[i] + lines[i].length();
    for (size_t i = 0; i < lines.size(); ++i) {
        for (size_t j = i + 2; j < lines.size(); ++j) {
            if (i == 0 && j + 1 == lines.size())
                continue;
            const double intervening_length = prefix_length[j] - prefix_length[i + 1];
            // Remapping an offset contour may leave a sub-resolution bridge
            // vertex between two physically adjacent pieces of the same turn.
            if (intervening_length <= scale_(0.01))
                continue;
            Point intersection;
            if (lines[i].intersection(lines[j], &intersection)) {
                if (message != nullptr) {
                    std::ostringstream out;
                    out << "crossing segments " << i << " and " << j << " of " << lines.size() << " at " << us(intersection.x()) << ", " <<
                        us(intersection.y()) << "; segment " << i << ": " << us(lines[i].a.x()) << ", " <<
                        us(lines[i].a.y()) << " -> " << us(lines[i].b.x()) << ", " << us(lines[i].b.y()) <<
                        "; segment " << j << ": " << us(lines[j].a.x()) << ", " << us(lines[j].a.y()) << " -> " <<
                        us(lines[j].b.x()) << ", " << us(lines[j].b.y());
                    out << "; first points:";
                    for (size_t k = 0; k < std::min<size_t>(path.points.size(), 12); ++k)
                        out << " " << k << "=(" << us(path.points[k].x()) << "," << us(path.points[k].y()) << ")";
                    out << "; around second segment:";
                    for (size_t k = j > 4 ? j - 4 : 0; k < std::min(path.points.size(), j + 7); ++k)
                        out << " " << k << "=(" << us(path.points[k].x()) << "," << us(path.points[k].y()) << ")";
                    *message = out.str();
                }
                return true;
            }
        }
    }
    return false;
}

double point_segment_distance(const Point &point, const Point &a, const Point &b)
{
    const Vec2d ap = (point - a).cast<double>();
    const Vec2d ab = (b - a).cast<double>();
    const double denom = ab.squaredNorm();
    const double t = denom <= 1e-9 ? 0.0 : std::clamp(ap.dot(ab) / denom, 0.0, 1.0);
    const Vec2d q = a.cast<double>() + ab * t;
    return (point.cast<double>() - q).norm();
}

double path_distance_to_point(const Polyline &path, const Point &point)
{
    double best = std::numeric_limits<double>::infinity();
    for (const Line &line : path.lines())
        best = std::min(best, point_segment_distance(point, line.a, line.b));
    return best;
}

ExPolygon corpus_shape(const size_t index)
{
    const size_t vertices = 8 + index % 17;
    const double base_radius = 14.0 + double((index * 7) % 13);
    const double aspect_x = 0.70 + 0.10 * double((index * 5) % 7);
    const double aspect_y = 0.70 + 0.10 * double((index * 3) % 7);
    const double modulation = 0.05 + 0.04 * double(index % 10);
    const size_t lobes = 2 + (index * 11) % 7;
    const double phase = 2.0 * PI * double((index * 37) % 101) / 101.0;
    Points points;
    points.reserve(vertices);
    for (size_t i = 0; i < vertices; ++i) {
        const double angle = 2.0 * PI * double(i) / double(vertices);
        const double radius = base_radius * (1.0 + modulation * std::cos(double(lobes) * angle + phase));
        points.emplace_back(p(radius * aspect_x * std::cos(angle), radius * aspect_y * std::sin(angle)));
    }
    return ExPolygon(std::move(points));
}

double segment_distance(const Line &a, const Line &b)
{
    Point intersection;
    if (a.intersection(b, &intersection))
        return 0.0;
    return std::min({
        point_segment_distance(a.a, b.a, b.b),
        point_segment_distance(a.b, b.a, b.b),
        point_segment_distance(b.a, a.a, a.b),
        point_segment_distance(b.b, a.a, a.b),
    });
}

double seam_to_adjacent_spacing(const Polyline &path)
{
    if (path.points.size() < 5)
        return 0.0;
    const Line seam(path.points[path.points.size() - 2], path.points.back());
    double best = std::numeric_limits<double>::infinity();
    const Lines lines = path.lines();
    for (size_t i = 4; i + 5 < lines.size(); ++i)
        best = std::min(best, unscale<double>(segment_distance(seam, lines[i])));
    return best == std::numeric_limits<double>::infinity() ? 0.0 : best;
}

ContinuousFermat::PathValidation validate_with_unit_multipliers(
    const ExPolygons &area,
    const Flow &flow,
    const Polyline &path)
{
    return ContinuousFermat::validate_layer_path(
        area,
        flow,
        path,
        std::vector<float>(path.points.empty() ? 0 : path.points.size() - 1, 1.0f));
}

void require_path_inside(const ExPolygon &area, const Polyline &path)
{
    REQUIRE(path.points.size() > 8);
    for (size_t i = 0; i < path.points.size(); ++i) {
        const Point &point = path.points[i];
        if (std::abs(us(point.x())) >= 100000.0 || std::abs(us(point.y())) >= 100000.0) {
            FAIL("point " << i << " has invalid coordinates: " << us(point.x()) << ", " << us(point.y()));
        }
        if (!area.contains(point)) {
            FAIL("point " << i << " is outside printable area: " << us(point.x()) << ", " << us(point.y()));
        }
    }
}

Polyline require_continuous_fermat_path(const ExPolygon &area)
{
    Flow flow(1.2f, 0.2f, 1.2f);
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ area }, flow);
    const Polyline &path = generated.path;

    require_path_inside(area, path);
    REQUIRE(path.points.front() == path.points.back());
    REQUIRE(seam_to_adjacent_spacing(path) <= double(flow.spacing()) * 1.10);
    std::string crossing;
    const bool has_crossing = segments_cross_nonlocal(path, &crossing);
    INFO(crossing);
    CHECK_FALSE(has_crossing);
    const ContinuousFermat::PathValidation validation =
        ContinuousFermat::validate_layer_path({ area }, flow, path, generated.extrusion_multipliers);
    INFO("coverage=" << validation.coverage_ratio << " exact_coverage=" << validation.exact_coverage_ratio <<
         " outside=" << validation.outside_ratio << " material=" << validation.material_ratio <<
         " containment=" << validation.containment_violations << " crossings=" << validation.crossings <<
         " spacing=" << validation.spacing_violations << " bead_overlaps=" << validation.bead_overlap_violations <<
         " turnbacks=" << validation.turnback_violations <<
         " multiplier_min=" << *std::min_element(generated.extrusion_multipliers.begin(), generated.extrusion_multipliers.end()) <<
         " multiplier_max=" << *std::max_element(generated.extrusion_multipliers.begin(), generated.extrusion_multipliers.end()));
    INFO(validation.reason);
    REQUIRE(validation.ok);
    return path;
}

ContinuousFermat::PathValidation validate_generated_continuous_fermat_path(const ExPolygon &area)
{
    const Flow flow(1.2f, 0.2f, 1.2f);
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ area }, flow);
    return ContinuousFermat::validate_layer_path({ area }, flow, generated.path, generated.extrusion_multipliers);
}

void require_rectangle_outer_wall_corners(
    const Polyline &path,
    const Flow &flow,
    const double min_x,
    const double min_y,
    const double max_x,
    const double max_y)
{
    const double inset = double(flow.width()) * 0.5;
    const double tolerance = scale_(0.02);
    for (const Point &corner : {
             p(min_x + inset, min_y + inset),
             p(max_x - inset, min_y + inset),
             p(max_x - inset, max_y - inset),
             p(min_x + inset, max_y - inset),
         }) {
        REQUIRE(path_distance_to_point(path, corner) <= tolerance);
    }
}

} // namespace

TEST_CASE("Continuous Fermat generates one non-crossing path for a rectangle", "[ContinuousFermat]")
{
    ExPolygon rectangle({
        p(-12.0, -12.0),
        p( 12.0, -12.0),
        p( 12.0,  12.0),
        p(-12.0,  12.0),
    });
    Polyline path = require_continuous_fermat_path(rectangle);
    Flow flow(1.2f, 0.2f, 1.2f);
    require_rectangle_outer_wall_corners(path, flow, -12.0, -12.0, 12.0, 12.0);

    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow);
    REQUIRE(generated.path.points == path.points);
    const ContinuousFermat::PathValidation validation =
        ContinuousFermat::validate_layer_path({ rectangle }, flow, generated.path, generated.extrusion_multipliers);
    INFO("points=" << path.points.size() << " length_mm=" << unscale<double>(path.length()));
    INFO("coverage=" << validation.coverage_ratio << " exact_coverage=" << validation.exact_coverage_ratio <<
         " outside=" << validation.outside_ratio <<
         " material=" << validation.material_ratio << " multiplier_min=" <<
         *std::min_element(generated.extrusion_multipliers.begin(), generated.extrusion_multipliers.end()) <<
         " multiplier_max=" <<
         *std::max_element(generated.extrusion_multipliers.begin(), generated.extrusion_multipliers.end()));
    INFO(validation.reason);
    REQUIRE(validation.ok);
    REQUIRE(validation.closed);
    REQUIRE(validation.coverage_ratio >= 0.990);
    REQUIRE(validation.exact_coverage_ratio >= 0.985);
    REQUIRE(validation.outside_ratio <= 0.020);
    REQUIRE(validation.material_ratio >= 0.980);
    REQUIRE(validation.material_ratio <= 1.020);
    REQUIRE(validation.containment_violations == 0);
    REQUIRE(validation.crossings == 0);
    REQUIRE(validation.spacing_violations <= 3);
    REQUIRE(validation.redeposition_ratio <= 0.035);
    REQUIRE(validation.turnback_violations == 0);
}

TEST_CASE("Continuous Fermat validates a common 0.4 mm nozzle rectangle", "[ContinuousFermat][regression]")
{
    const ExPolygon rectangle({ p(-12.0, -12.0), p(12.0, -12.0), p(12.0, 12.0), p(-12.0, 12.0) });
    const Flow flow(0.42f, 0.2f, 0.4f);
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow);
    const ContinuousFermat::PathValidation validation =
        ContinuousFermat::validate_layer_path({ rectangle }, flow, generated.path, generated.extrusion_multipliers);
    INFO(validation.reason);
    REQUIRE(validation.ok);
}

TEST_CASE("Continuous Fermat exact coverage is authoritative over sampled grid phase", "[ContinuousFermat][regression][sampled-coverage]")
{
    // A narrow C-shaped corridor puts a small exact uncovered band through
    // more 70-cell sample centers than its area fraction represents. This is
    // the deterministic form of the phase alias observed on coordinate-
    // jittered, otherwise identical internet-corpus layers.
    const ExPolygon c_corridor({
        p(9.822906640, 10.455977160), p(9.002727440, 11.015341640),
        p(9.002481640, 11.015448480), p(8.037485800, 11.434940800),
        p(8.037212920, 11.434999520), p(6.965782800, 11.665656280),
        p(6.318945920, 11.698911720), p(6.318781120, 11.698920160),
        p(6.318781120, 12.510739880), p(6.318781120, 12.510946680),
        p(12.443185760, 12.510946680), p(12.444745760, 12.510946680),
        p(12.444745760, 0.021619200), p(12.444745760, 0.018437920),
        p(6.320341160, 0.018437920), p(6.318781120, 0.018437920),
        p(6.318781120, 0.018644720), p(6.318781120, 0.830440200),
        p(6.318924360, 0.830446440), p(6.881202800, 0.854945520),
        p(7.978449360, 1.075709840), p(7.978728840, 1.075766080),
        p(8.987005360, 1.505351280), p(8.987262160, 1.505460720),
        p(9.855236840, 2.100280640), p(9.855457920, 2.100432120),
        p(10.554726320, 2.805135080), p(10.554904440, 2.805314560),
        p(11.114268880, 3.625493760), p(11.114411360, 3.625702680),
        p(11.114518240, 3.625948480), p(11.534010520, 4.590944320),
        p(11.534069280, 4.591217200), p(11.764726040, 5.662647320),
        p(11.773482360, 6.781847960), p(11.773484600, 6.782133080),
        p(11.773428360, 6.782412560), p(11.552664040, 7.879659080),
        p(11.552554640, 7.879915920), p(11.122969400, 8.888192400),
        p(10.528149480, 9.756167080), p(10.527998000, 9.756388160),
        p(10.527818480, 9.756566280), p(9.823115560, 10.455834680),
    });
    const Flow flow(0.4f, 0.2f, 0.4f);
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ c_corridor }, flow, 0.8);
    const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
        { c_corridor }, flow, generated.path, generated.extrusion_multipliers, 0.8);
    INFO("sampled=" << validation.coverage_ratio << " exact=" << validation.exact_coverage_ratio <<
         " reason=" << validation.reason);
    REQUIRE(validation.exact_coverage_ratio >= 0.980);
    REQUIRE(validation.coverage_ratio < 0.980);
    REQUIRE(validation.reason.find(" coverage=") == std::string::npos);
}

TEST_CASE("Continuous Fermat derives a wider geometry phase for a thin ring", "[ContinuousFermat][regression][derived-width]")
{
    ExPolygon thin_ring(
        { p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0) },
        { p(-1.0, -1.0), p(-1.0, 1.0), p(1.0, 1.0), p(1.0, -1.0) });
    const Flow flow(0.4f, 0.2f, 0.4f);
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ thin_ring }, flow, 0.8);
    const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
        { thin_ring }, flow, generated.path, generated.extrusion_multipliers, 0.8);
    INFO("sampled=" << validation.coverage_ratio << " exact=" << validation.exact_coverage_ratio <<
         " material=" << validation.material_ratio << " redeposition=" << validation.redeposition_ratio <<
         " reason=" << validation.reason);
    REQUIRE(validation.ok);
    REQUIRE(validation.exact_coverage_ratio >= 0.980);
    REQUIRE(validation.material_ratio >= 0.980);
    REQUIRE(validation.material_ratio <= 1.020);
    REQUIRE(validation.crossings == 0);
    REQUIRE(validation.turnback_violations == 0);
    REQUIRE(*std::max_element(generated.extrusion_multipliers.begin(), generated.extrusion_multipliers.end()) > 1.0f);
}

TEST_CASE("Continuous Fermat validates a 0.5 mm line on a 0.4 mm nozzle", "[ContinuousFermat][regression]")
{
    const Flow flow(0.5f, 0.2f, 0.4f);
    REQUIRE(flow.spacing() == Catch::Approx(0.4570796f));
    const ExPolygon rectangle({ p(-6.0, -6.0), p(6.0, -6.0), p(6.0, 6.0), p(-6.0, 6.0) });
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow);
    const ContinuousFermat::PathValidation validation =
        ContinuousFermat::validate_layer_path({ rectangle }, flow, generated.path, generated.extrusion_multipliers);
    INFO("coverage=" << validation.coverage_ratio << " exact_coverage=" << validation.exact_coverage_ratio <<
         " outside=" << validation.outside_ratio << " material=" << validation.material_ratio);
    INFO(validation.reason);
    REQUIRE(validation.ok);
    REQUIRE(validation.coverage_ratio >= 0.985);
    REQUIRE(validation.exact_coverage_ratio >= 0.985);
    REQUIRE(validation.outside_ratio <= 0.020);
    REQUIRE(validation.material_ratio >= 0.980);
    REQUIRE(validation.material_ratio <= 1.020);
}

TEST_CASE("Continuous Fermat accepts safely validated convex contours", "[ContinuousFermat][regression]")
{
    const Flow flow(1.2f, 0.2f, 1.2f);
    const std::vector<ExPolygon> shapes {
        ExPolygon(circle_points(24.0, 6)),
        ExPolygon(circle_points(24.0, 96)),
        ExPolygon({
            p(-24.0, -14.0), p(20.0, -14.0), p(24.0, -10.0), p(24.0, 10.0),
            p(20.0, 14.0), p(-24.0, 14.0),
        }),
    };
    for (const ExPolygon &shape : shapes) {
        const ContinuousFermat::GeneratedPath generated =
            ContinuousFermat::generate_layer_path_with_metadata({ shape }, flow);
        const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
            { shape }, flow, generated.path, generated.extrusion_multipliers);
        INFO(validation.reason);
        REQUIRE(validation.ok);
        REQUIRE(validation.crossings == 0);
        REQUIRE(validation.redeposition_ratio <= 0.035);
        REQUIRE(validation.turnback_violations == 0);
    }
}

TEST_CASE("Continuous Fermat rectangle acceptance survey", "[.][ContinuousFermat][survey]")
{
    struct SurveyCase
    {
        const char *name;
        double width;
        double height;
        float line_width;
        float layer_height;
        float nozzle_diameter;
        double max_line_width { 0.0 };
    };

    const std::vector<SurveyCase> cases {
        { "square-6", 6.0, 6.0, 1.2f, 0.2f, 1.2f },
        { "square-8", 8.0, 8.0, 1.2f, 0.2f, 1.2f },
        { "square-10", 10.0, 10.0, 1.2f, 0.2f, 1.2f },
        { "square-12", 12.0, 12.0, 1.2f, 0.2f, 1.2f },
        { "square-16", 16.0, 16.0, 1.2f, 0.2f, 1.2f },
        { "square-20", 20.0, 20.0, 1.2f, 0.2f, 1.2f },
        { "square-24", 24.0, 24.0, 1.2f, 0.2f, 1.2f },
        { "square-30", 30.0, 30.0, 1.2f, 0.2f, 1.2f },
        { "square-40", 40.0, 40.0, 1.2f, 0.2f, 1.2f },
        { "wide-2x", 24.0, 12.0, 1.2f, 0.2f, 1.2f },
        { "wide-4x", 48.0, 12.0, 1.2f, 0.2f, 1.2f },
        { "tall-2x", 12.0, 24.0, 1.2f, 0.2f, 1.2f },
        { "tall-4x", 12.0, 48.0, 1.2f, 0.2f, 1.2f },
        { "line-0.42", 12.0, 12.0, 0.42f, 0.2f, 0.4f },
        { "line-0.50", 12.0, 12.0, 0.50f, 0.2f, 0.4f },
        { "line-0.60", 12.0, 12.0, 0.60f, 0.2f, 0.4f },
        { "line-0.72", 12.0, 12.0, 0.72f, 0.2f, 0.4f },
        { "square-6-wide", 6.0, 6.0, 1.2f, 0.2f, 1.2f, 2.4 },
        { "square-16-wide", 16.0, 16.0, 1.2f, 0.2f, 1.2f, 2.4 },
        { "line-0.72-wide", 12.0, 12.0, 0.72f, 0.2f, 0.4f, 0.8 },
    };

    std::cout << "name,width,height,line_width,nozzle,ok,coverage,exact,outside,material,containment,crossings,spacing,overlaps,turnbacks,min_multiplier,max_multiplier,elapsed_ms,reason\n";
    for (const SurveyCase &survey : cases) {
        const ExPolygon rectangle({
            p(-survey.width * 0.5, -survey.height * 0.5),
            p( survey.width * 0.5, -survey.height * 0.5),
            p( survey.width * 0.5,  survey.height * 0.5),
            p(-survey.width * 0.5,  survey.height * 0.5),
        });
        const Flow flow(survey.line_width, survey.layer_height, survey.nozzle_diameter);
        const auto started = std::chrono::steady_clock::now();
        const ContinuousFermat::GeneratedPath generated =
            ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow, survey.max_line_width);
        const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
            { rectangle }, flow, generated.path, generated.extrusion_multipliers, survey.max_line_width);
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - started).count();
        float min_value = 0.0f;
        float max_value = 0.0f;
        if (!generated.extrusion_multipliers.empty()) {
            const auto [min_multiplier, max_multiplier] =
                std::minmax_element(generated.extrusion_multipliers.begin(), generated.extrusion_multipliers.end());
            min_value = *min_multiplier;
            max_value = *max_multiplier;
        }
        std::cout << survey.name << ',' << survey.width << ',' << survey.height << ',' << survey.line_width << ','
                  << survey.nozzle_diameter << ',' << validation.ok << ',' << validation.coverage_ratio << ','
                  << validation.exact_coverage_ratio << ',' << validation.outside_ratio << ',' << validation.material_ratio << ','
                  << validation.containment_violations << ',' << validation.crossings << ',' << validation.spacing_violations << ','
                  << validation.bead_overlap_violations << ',' << validation.turnback_violations << ',' << min_value << ','
                  << max_value << ',' << elapsed << ',' << validation.reason << '\n';
        CHECK(validation.closed);
        CHECK(generated.extrusion_multipliers.size() + 1 == generated.path.points.size());
    }
}

TEST_CASE("Continuous Fermat adaptive-width effect survey", "[.][ContinuousFermat][width-survey]")
{
    const Flow flow(1.2f, 0.2f, 1.2f);
    std::cout << "size,legacy_ok,wide_ok,legacy_coverage,wide_coverage,legacy_material,wide_material\n";
    for (double size = 4.0; size <= 30.0; size += 0.5) {
        const ExPolygon rectangle({
            p(-size * 0.5, -size * 0.5), p(size * 0.5, -size * 0.5),
            p(size * 0.5, size * 0.5), p(-size * 0.5, size * 0.5),
        });
        const ContinuousFermat::GeneratedPath legacy =
            ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow);
        const ContinuousFermat::PathValidation legacy_validation =
            ContinuousFermat::validate_layer_path({ rectangle }, flow, legacy.path, legacy.extrusion_multipliers);
        const ContinuousFermat::GeneratedPath wide =
            ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow, 2.4);
        const ContinuousFermat::PathValidation wide_validation =
            ContinuousFermat::validate_layer_path({ rectangle }, flow, wide.path, wide.extrusion_multipliers, 2.4);
        if (legacy_validation.ok != wide_validation.ok ||
            std::abs(legacy_validation.exact_coverage_ratio - wide_validation.exact_coverage_ratio) > 1e-6 ||
            std::abs(legacy_validation.material_ratio - wide_validation.material_ratio) > 1e-6) {
            std::cout << size << ',' << legacy_validation.ok << ',' << wide_validation.ok << ','
                      << legacy_validation.exact_coverage_ratio << ',' << wide_validation.exact_coverage_ratio << ','
                      << legacy_validation.material_ratio << ',' << wide_validation.material_ratio << '\n';
        }
        CHECK(legacy_validation.closed);
        CHECK(wide_validation.closed);
    }
}

TEST_CASE("Continuous Fermat topology acceptance survey", "[.][ContinuousFermat][topology-survey]")
{
    struct ShapeCase { const char *name; ExPolygon shape; };
    const std::vector<ShapeCase> cases {
        { "triangle", ExPolygon({ p(-24.0, -18.0), p(24.0, -18.0), p(0.0, 24.0) }) },
        { "hexagon", ExPolygon(circle_points(24.0, 6)) },
        { "circle", ExPolygon(circle_points(24.0, 96)) },
        { "chamfered-rectangle", ExPolygon({
            p(-24.0, -14.0), p(20.0, -14.0), p(24.0, -10.0), p(24.0, 10.0),
            p(20.0, 14.0), p(-24.0, 14.0),
        }) },
        { "l-shape", ExPolygon({
            p(-24.0, -24.0), p(24.0, -24.0), p(24.0, -8.0),
            p(-8.0, -8.0), p(-8.0, 24.0), p(-24.0, 24.0),
        }) },
    };
    const Flow flow(1.2f, 0.2f, 1.2f);
    std::cout << "shape,ok,coverage,exact,outside,material,crossings,spacing,overlaps,turnbacks,reason\n";
    for (const ShapeCase &shape_case : cases) {
        const ContinuousFermat::GeneratedPath generated =
            ContinuousFermat::generate_layer_path_with_metadata({ shape_case.shape }, flow);
        const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
            { shape_case.shape }, flow, generated.path, generated.extrusion_multipliers);
        std::cout << shape_case.name << ',' << validation.ok << ',' << validation.coverage_ratio << ','
                  << validation.exact_coverage_ratio << ',' << validation.outside_ratio << ',' << validation.material_ratio << ','
                  << validation.crossings << ',' << validation.spacing_violations << ',' << validation.bead_overlap_violations << ','
                  << validation.turnback_violations << ',' << validation.reason << '\n';
        CHECK((generated.path.points.empty() || validation.closed));
    }
}

TEST_CASE("Continuous Fermat 100-shape eligibility corpus", "[.][ContinuousFermat][corpus-100]")
{
    const Flow flow(1.2f, 0.2f, 1.2f);
    size_t accepted = 0;
    std::map<std::string, size_t> rejection_categories;
    for (size_t index = 0; index < 100; ++index) {
        const ExPolygon shape = corpus_shape(index);
        const ContinuousFermat::GeneratedPath generated =
            ContinuousFermat::generate_layer_path_with_metadata({ shape }, flow, 2.4);
        const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
            { shape }, flow, generated.path, generated.extrusion_multipliers, 2.4);
        if (validation.emittable) {
            ++accepted;
        } else if (!validation.closed) {
            ++rejection_categories["open"];
        } else if (validation.crossings > 0) {
            ++rejection_categories["crossing"];
        } else if (validation.redeposition_ratio > 0.035) {
            ++rejection_categories["redeposition"];
        } else if (validation.turnback_violations > 0) {
            ++rejection_categories["turnback"];
        } else if (validation.exact_coverage_ratio < 0.985) {
            ++rejection_categories["coverage"];
        } else {
            ++rejection_categories["other"];
        }
        if (!validation.ok)
            std::cout << "corpus_advisory index=" << index << " exact=" << validation.exact_coverage_ratio
                      << " sampled=" << validation.coverage_ratio << " material=" << validation.material_ratio
                      << " redeposition=" << validation.redeposition_ratio << " crossings=" << validation.crossings
                      << " overlap_pairs=" << validation.bead_overlap_violations << " reason=" << validation.reason << '\n';
    }
    std::cout << "corpus_100 accepted=" << accepted;
    for (const auto &[category, count] : rejection_categories)
        std::cout << ' ' << category << '=' << count;
    std::cout << '\n';
    INFO("accepted=" << accepted);
    REQUIRE(accepted == 100);
}

TEST_CASE("Continuous Fermat touching-connector corpus regression", "[.][ContinuousFermat][corpus-touch]")
{
    const Flow flow(1.2f, 0.2f, 1.2f);
    const ExPolygon shape = corpus_shape(79);
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ shape }, flow, 2.4);
    const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
        { shape }, flow, generated.path, generated.extrusion_multipliers, 2.4);
    std::string crossing;
    const bool has_crossing = segments_cross_nonlocal(generated.path, &crossing);
    INFO(crossing);
    INFO(validation.reason);
    REQUIRE_FALSE(has_crossing);
    REQUIRE(validation.ok);
}

TEST_CASE("Continuous Fermat sharp-point and hole corpus", "[.][ContinuousFermat][sharp-hole-corpus]")
{
    ExPolygon annulus(circle_points(15.0, 96));
    annulus.holes.emplace_back(circle_points(5.0, 48, Point(0, 0), false));
    ExPolygon square_hole(
        { p(-15.0, -15.0), p(15.0, -15.0), p(15.0, 15.0), p(-15.0, 15.0) },
        { p(-5.0, -5.0), p(-5.0, 5.0), p(5.0, 5.0), p(5.0, -5.0) });
    ExPolygon two_holes({ p(-10.0, -6.0), p(10.0, -6.0), p(10.0, 6.0), p(-10.0, 6.0) });
    two_holes.holes.emplace_back(circle_points(1.5, 24, p(-3.5, 0.0), false));
    two_holes.holes.emplace_back(circle_points(1.5, 24, p(3.5, 0.0), false));

    const std::vector<std::pair<const char *, ExPolygon>> cases {
        { "triangle-small", ExPolygon({ p(-5.0, -4.0), p(5.0, -4.0), p(0.0, 6.0) }) },
        { "triangle", ExPolygon({ p(-10.0, -8.0), p(10.0, -8.0), p(0.0, 12.0) }) },
        { "acute-wedge", ExPolygon({ p(-15.0, -4.0), p(15.0, 0.0), p(-15.0, 4.0) }) },
        { "diamond", ExPolygon({ p(0.0, -15.0), p(6.0, 0.0), p(0.0, 15.0), p(-6.0, 0.0) }) },
        { "star", ExPolygon(star_points(15.0, 6.0, 7)) },
        { "annulus", std::move(annulus) },
        { "square-hole", std::move(square_hole) },
        { "two-holes", std::move(two_holes) },
    };
    const Flow flow(0.4f, 0.2f, 0.4f);
    size_t accepted = 0;
    size_t evaluated = 0;
    const char *case_filter = std::getenv("CONTINUOUS_FERMAT_SHARP_HOLE_CASE");
    for (const auto &[name, shape] : cases) {
        if (case_filter != nullptr && std::string(case_filter) != name)
            continue;
        ++evaluated;
        const ContinuousFermat::GeneratedPath generated =
            ContinuousFermat::generate_layer_path_with_metadata({ shape }, flow, 0.8);
        const ContinuousFermat::PathValidation validation = ContinuousFermat::validate_layer_path(
            { shape }, flow, generated.path, generated.extrusion_multipliers, 0.8);
        std::string crossing_detail;
        if (validation.crossings > 0)
            segments_cross_nonlocal(generated.path, &crossing_detail);
        std::cout << "sharp_hole_case name=" << name << " ok=" << validation.ok
                  << " exact=" << validation.exact_coverage_ratio << " sampled=" << validation.coverage_ratio
                  << " material=" << validation.material_ratio << " redeposition=" << validation.redeposition_ratio
                  << " crossings=" << validation.crossings << " turnbacks=" << validation.turnback_violations
                  << " reason=" << validation.reason << " crossing_detail=" << crossing_detail << '\n';
        accepted += validation.emittable;
    }
    std::cout << "sharp_hole_corpus accepted=" << accepted << '/' << evaluated << '\n';
    REQUIRE(accepted == evaluated);
}

TEST_CASE("Continuous Fermat safely rejects an undercovered size phase", "[ContinuousFermat][unsupported]")
{
    const Flow flow(1.2f, 0.2f, 1.2f);
    const ExPolygon rectangle({ p(-3.0, -3.0), p(3.0, -3.0), p(3.0, 3.0), p(-3.0, 3.0) });
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow);
    const ContinuousFermat::PathValidation validation =
        ContinuousFermat::validate_layer_path({ rectangle }, flow, generated.path, generated.extrusion_multipliers);
    REQUIRE_FALSE(validation.ok);
    REQUIRE(validation.emittable);
    REQUIRE(validation.crossings == 0);
    REQUIRE(validation.coverage_ratio < 0.980);
}

TEST_CASE("Continuous Fermat distinguishes structural failures from geometric advisories", "[ContinuousFermat][safety]")
{
    const ExPolygon rectangle({
        p(-30.0, -18.0),
        p( 30.0, -18.0),
        p( 30.0,  18.0),
        p(-30.0,  18.0),
    });
    const Flow flow(1.2f, 0.2f, 1.2f);

    SECTION("open path")
    {
        const Polyline open({ p(-20.0, -10.0), p(20.0, -10.0), p(20.0, 10.0), p(-20.0, 10.0) });
        const ContinuousFermat::PathValidation validation = validate_with_unit_multipliers({ rectangle }, flow, open);
        REQUIRE_FALSE(validation.ok);
        REQUIRE_FALSE(validation.emittable);
        REQUIRE_FALSE(validation.closed);
    }

    SECTION("underfilled closed path")
    {
        const Polyline underfilled({ p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0), p(-2.0, -2.0) });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ rectangle }, flow, underfilled);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.emittable);
        REQUIRE(validation.closed);
        REQUIRE(validation.coverage_ratio < 0.985);
    }

    SECTION("self-crossing path")
    {
        const Polyline crossing({ p(-20.0, -10.0), p(20.0, 10.0), p(-20.0, 10.0), p(20.0, -10.0), p(-20.0, -10.0) });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ rectangle }, flow, crossing);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.emittable);
        REQUIRE(validation.crossings > 0);
    }

    SECTION("small local bow-tie crossing")
    {
        const Polyline crossing({ p(0.0, 0.0), p(2.0, 2.0), p(0.0, 2.0), p(2.0, 0.0), p(0.0, 0.0) });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ rectangle }, flow, crossing);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.crossings > 0);
    }

    SECTION("sub-0.01 mm bow-tie crossing")
    {
        const Polyline crossing({ p(0.0, 0.0), p(0.004, 0.004), p(0.0, 0.004), p(0.004, 0.0), p(0.0, 0.0) });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ rectangle }, flow, crossing);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.crossings > 0);
    }

    SECTION("short exact containment excursion")
    {
        ExPolygon area_with_tiny_hole(
            rectangle.contour.points,
            {{ p(0.24, -0.01), p(0.24, 0.01), p(0.26, 0.01), p(0.26, -0.01) }});
        const Polyline through_hole({ p(-1.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(-1.0, 1.0), p(-1.0, 0.0) });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ area_with_tiny_hole }, flow, through_hole);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.containment_violations > 0);
    }

    SECTION("extruded hairpin reversal")
    {
        const Polyline hairpin({
            p(-20.0, -10.0),
            p( 20.0, -10.0),
            p(-19.0, -10.0),
            p(-19.0,  10.0),
            p(-20.0,  10.0),
            p(-20.0, -10.0),
        });
        const ContinuousFermat::PathValidation validation = validate_with_unit_multipliers({ rectangle }, flow, hairpin);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.turnback_violations > 0);
    }

    SECTION("turn above the conservative 135 degree limit")
    {
        const Polyline turnback({
            p(-2.0, 0.0),
            p( 0.0, 0.0),
            p(-1.532089, 1.285575),
            p(-2.0, 2.0),
            p(-2.0, 0.0),
        });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ rectangle }, flow, turnback);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.turnback_violations > 0);
    }

    SECTION("missing extrusion metadata")
    {
        const Polyline path({ p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0), p(-2.0, -2.0) });
        const ContinuousFermat::PathValidation validation =
            ContinuousFermat::validate_layer_path({ rectangle }, flow, path, {});
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.reason.find("metadata cardinality") != std::string::npos);
    }

    SECTION("invalid multiplier on a zero-length segment")
    {
        const Polyline path({
            p(-2.0, -2.0), p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0), p(-2.0, -2.0),
        });
        std::vector<float> multipliers(path.points.size() - 1, 1.0f);
        multipliers.front() = std::numeric_limits<float>::quiet_NaN();
        const ContinuousFermat::PathValidation validation =
            ContinuousFermat::validate_layer_path({ rectangle }, flow, path, multipliers);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.reason.find("outside the safe range") != std::string::npos);
    }

    SECTION("multiplier bounds are exact")
    {
        const Polyline path({ p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0), p(-2.0, -2.0) });
        for (const float invalid : { 0.5999f, 1.6001f, std::numeric_limits<float>::infinity() }) {
            std::vector<float> multipliers(path.points.size() - 1, 1.0f);
            multipliers[1] = invalid;
            const ContinuousFermat::PathValidation validation =
                ContinuousFermat::validate_layer_path({ rectangle }, flow, path, multipliers);
            REQUIRE_FALSE(validation.ok);
        }
    }

    SECTION("configured maximum line width controls adaptive metadata")
    {
        const Flow narrow_flow(0.42f, 0.2f, 0.4f);
        const Polyline path({ p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0), p(-2.0, -2.0) });
        const std::vector<float> multipliers(path.points.size() - 1, 1.8f);
        const ContinuousFermat::PathValidation legacy =
            ContinuousFermat::validate_layer_path({ rectangle }, narrow_flow, path, multipliers);
        REQUIRE_FALSE(legacy.ok);
        REQUIRE(legacy.reason.find("multiplier is outside") != std::string::npos);

        const ContinuousFermat::PathValidation relaxed =
            ContinuousFermat::validate_layer_path({ rectangle }, narrow_flow, path, multipliers, 0.8);
        REQUIRE_FALSE(relaxed.ok); // The deliberately tiny path is still underfilled.
        REQUIRE(relaxed.reason.find("multiplier is outside") == std::string::npos);
    }

    SECTION("adaptive width is bounded by nozzle diameter")
    {
        const Flow small_nozzle_flow(1.2f, 0.2f, 0.4f);
        const Polyline path({ p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0), p(-2.0, -2.0) });
        std::vector<float> multipliers(path.points.size() - 1, 1.0f);
        const ContinuousFermat::PathValidation validation =
            ContinuousFermat::validate_layer_path({ rectangle }, small_nozzle_flow, path, multipliers);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.reason.find("two-nozzle safety limit") != std::string::npos);
    }
}

TEST_CASE("Continuous Fermat categorizes failed geometry at 2.5 line widths", "[ContinuousFermat][thin-feature]")
{
    const Flow flow(0.4f, 0.2f, 0.4f);
    const Polyline diagnostic_path({
        p(-1.0, -1.0), p(1.0, -1.0), p(1.0, 1.0), p(-1.0, 1.0), p(-1.0, -1.0),
    });
    const std::vector<float> multipliers(diagnostic_path.points.size() - 1, 1.0f);

    const auto validate = [&](const ExPolygon &shape) {
        return ContinuousFermat::validate_layer_path({ shape }, flow, diagnostic_path, multipliers);
    };
    const auto annulus = [](const double wall) {
        const double inner = 10.0 - wall;
        return ExPolygon(
            { p(-10.0, -10.0), p(10.0, -10.0), p(10.0, 10.0), p(-10.0, 10.0) },
            { p(-inner, -inner), p(-inner, inner), p(inner, inner), p(inner, -inner) });
    };
    const auto dumbbell = [](const double neck) {
        const double half_neck = 0.5 * neck;
        return ExPolygon({
            p(-12.0, -6.0), p(-3.0, -6.0), p(-3.0, -half_neck),
            p(  3.0, -half_neck), p(3.0, -6.0), p(12.0, -6.0),
            p( 12.0,  6.0), p(3.0,  6.0), p(3.0,  half_neck),
            p( -3.0,  half_neck), p(-3.0, 6.0), p(-12.0, 6.0),
        });
    };

    const ContinuousFermat::PathValidation thin_annulus = validate(annulus(0.96));
    REQUIRE(thin_annulus.emittable);
    REQUIRE_FALSE(thin_annulus.ok);
    REQUIRE(thin_annulus.thin_feature_class == "section_below_2.5_line_widths");
    REQUIRE(thin_annulus.thin_feature_threshold_mm == Catch::Approx(1.0));

    const ContinuousFermat::PathValidation wider_annulus = validate(annulus(1.04));
    REQUIRE(wider_annulus.emittable);
    REQUIRE(wider_annulus.thin_feature_class.empty());

    const ContinuousFermat::PathValidation thin_neck = validate(dumbbell(0.96));
    REQUIRE(thin_neck.emittable);
    REQUIRE(thin_neck.thin_feature_class == "neck_below_2.5_line_widths");

    const ContinuousFermat::PathValidation wider_neck = validate(dumbbell(1.04));
    REQUIRE(wider_neck.emittable);
    REQUIRE(wider_neck.thin_feature_class.empty());
}

TEST_CASE("Continuous Fermat extrusion metadata follows or blocks geometry mutations", "[ContinuousFermat][metadata]")
{
    ExtrusionPath path(erExternalPerimeter, 1.0, 1.2, 0.2);
    path.polyline = Polyline({ p(0.0, 0.0), p(1.0, 0.0), p(3.0, 0.0) });
    path.continuous_fermat_extrusion_multipliers = { 1.1f, 1.2f };
    path.continuous_fermat_validation_warning = "exact_coverage=0.9";
    path.set_continuous_fermat();
    path.set_reverse();

    REQUIRE_FALSE(path.can_reverse());
    path.reverse();
    REQUIRE(path.polyline.points == Points({ p(3.0, 0.0), p(1.0, 0.0), p(0.0, 0.0) }));
    REQUIRE(path.continuous_fermat_extrusion_multipliers == std::vector<float>({ 1.2f, 1.1f }));
    REQUIRE(path.continuous_fermat_validation_warning == "exact_coverage=0.9");

    const Points protected_points = path.polyline.points;
    path.clip_end(scale_(0.5));
    path.simplify(scale_(10.0));
    path.simplify_by_fitting_arc(scale_(10.0));
    REQUIRE(path.polyline.points == protected_points);
    REQUIRE(path.continuous_fermat_extrusion_multipliers.size() + 1 == path.polyline.points.size());
}

TEST_CASE("Continuous Fermat emits a generated square-hole footprint", "[ContinuousFermat][regression]")
{
    ExPolygon square_with_hole(
        {
            p(-32.0, -32.0),
            p( 32.0, -32.0),
            p( 32.0,  32.0),
            p(-32.0,  32.0),
        },
        {
            p(-9.0, -9.0),
            p(-9.0,  9.0),
            p( 9.0,  9.0),
            p( 9.0, -9.0),
        });
    const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(square_with_hole);
    INFO(validation.reason);
    REQUIRE(validation.emittable);
}

TEST_CASE("Continuous Fermat handles or safely rejects concave and branched island shapes", "[ContinuousFermat][regression]")
{
    SECTION("c shape")
    {
        ExPolygon c_shape({
            p(-48.0, -32.0),
            p( 48.0, -32.0),
            p( 48.0, -16.0),
            p(-18.0, -16.0),
            p(-18.0,  16.0),
            p( 48.0,  16.0),
            p( 48.0,  32.0),
            p(-48.0,  32.0),
        });
        const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(c_shape);
        INFO(validation.reason);
        REQUIRE(validation.emittable);
    }

    SECTION("dumbbell")
    {
        ExPolygon dumbbell({
            p(-52.0, -25.0),
            p(-20.0, -25.0),
            p(-20.0,  -8.0),
            p( 20.0,  -8.0),
            p( 20.0, -25.0),
            p( 52.0, -25.0),
            p( 52.0,  25.0),
            p( 20.0,  25.0),
            p( 20.0,   8.0),
            p(-20.0,   8.0),
            p(-20.0,  25.0),
            p(-52.0,  25.0),
        });
        const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(dumbbell);
        INFO(validation.reason);
        REQUIRE(validation.emittable);
    }

    SECTION("star")
    {
        ExPolygon star(star_points(45.0, 22.0, 7));
        const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(star);
        INFO(validation.reason);
        REQUIRE(validation.emittable);
    }
}

TEST_CASE("Continuous Fermat handles or safely rejects hole topologies", "[ContinuousFermat][regression]")
{
    SECTION("annulus")
    {
        ExPolygon annulus(circle_points(42.0, 160));
        annulus.holes.emplace_back(circle_points(16.0, 96, Point(0, 0), false));
        const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(annulus);
        INFO(validation.reason);
        REQUIRE(validation.emittable);
    }

    SECTION("two holes")
    {
        ExPolygon two_holes({
            p(-56.0, -34.0),
            p( 56.0, -34.0),
            p( 56.0,  34.0),
            p(-56.0,  34.0),
        });
        two_holes.holes.emplace_back(circle_points(10.0, 64, p(-21.0, 0.0), false));
        two_holes.holes.emplace_back(circle_points(10.0, 64, p( 21.0, 0.0), false));
        const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(two_holes);
        INFO(validation.reason);
        REQUIRE(validation.emittable);
        REQUIRE(validation.closed);
    }
}
