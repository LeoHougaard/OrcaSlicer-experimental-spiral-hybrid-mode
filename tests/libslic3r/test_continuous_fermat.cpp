#include <catch2/catch_all.hpp>

#include "libslic3r/ContinuousFermat.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Line.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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
    REQUIRE(validation.exact_coverage_ratio >= 0.970);
    REQUIRE(validation.outside_ratio <= 0.020);
    REQUIRE(validation.material_ratio >= 0.980);
    REQUIRE(validation.material_ratio <= 1.020);
    REQUIRE(validation.containment_violations == 0);
    REQUIRE(validation.crossings == 0);
    REQUIRE(validation.spacing_violations <= 3);
    REQUIRE(validation.bead_overlap_violations == 0);
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

TEST_CASE("Continuous Fermat safely rejects an undercovered size phase", "[ContinuousFermat][unsupported]")
{
    const Flow flow(1.2f, 0.2f, 1.2f);
    const ExPolygon rectangle({ p(-10.0, -10.0), p(10.0, -10.0), p(10.0, 10.0), p(-10.0, 10.0) });
    const ContinuousFermat::GeneratedPath generated =
        ContinuousFermat::generate_layer_path_with_metadata({ rectangle }, flow);
    const ContinuousFermat::PathValidation validation =
        ContinuousFermat::validate_layer_path({ rectangle }, flow, generated.path, generated.extrusion_multipliers);
    REQUIRE_FALSE(validation.ok);
    REQUIRE(validation.crossings == 0);
    REQUIRE(validation.bead_overlap_violations == 0);
    REQUIRE(validation.coverage_ratio < 0.990);
}

TEST_CASE("Continuous Fermat mandatory validation rejects unsafe final paths", "[ContinuousFermat][safety]")
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
        REQUIRE_FALSE(validation.closed);
    }

    SECTION("underfilled closed path")
    {
        const Polyline underfilled({ p(-2.0, -2.0), p(2.0, -2.0), p(2.0, 2.0), p(-2.0, 2.0), p(-2.0, -2.0) });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ rectangle }, flow, underfilled);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.closed);
        REQUIRE(validation.coverage_ratio < 0.990);
    }

    SECTION("self-crossing path")
    {
        const Polyline crossing({ p(-20.0, -10.0), p(20.0, 10.0), p(-20.0, 10.0), p(20.0, -10.0), p(-20.0, -10.0) });
        const ContinuousFermat::PathValidation validation =
            validate_with_unit_multipliers({ rectangle }, flow, crossing);
        REQUIRE_FALSE(validation.ok);
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

TEST_CASE("Continuous Fermat extrusion metadata follows or blocks geometry mutations", "[ContinuousFermat][metadata]")
{
    ExtrusionPath path(erExternalPerimeter, 1.0, 1.2, 0.2);
    path.polyline = Polyline({ p(0.0, 0.0), p(1.0, 0.0), p(3.0, 0.0) });
    path.continuous_fermat_extrusion_multipliers = { 1.1f, 1.2f };
    path.set_continuous_fermat();
    path.set_reverse();

    REQUIRE_FALSE(path.can_reverse());
    path.reverse();
    REQUIRE(path.polyline.points == Points({ p(3.0, 0.0), p(1.0, 0.0), p(0.0, 0.0) }));
    REQUIRE(path.continuous_fermat_extrusion_multipliers == std::vector<float>({ 1.2f, 1.1f }));

    const Points protected_points = path.polyline.points;
    path.clip_end(scale_(0.5));
    path.simplify(scale_(10.0));
    path.simplify_by_fitting_arc(scale_(10.0));
    REQUIRE(path.polyline.points == protected_points);
    REQUIRE(path.continuous_fermat_extrusion_multipliers.size() + 1 == path.polyline.points.size());
}

TEST_CASE("Continuous Fermat safely rejects a square hole", "[ContinuousFermat][unsupported]")
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
    REQUIRE_FALSE(validation.ok);
}

TEST_CASE("Continuous Fermat rejects unsupported concave and branched island shapes", "[ContinuousFermat][unsupported]")
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
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.crossings > 0);
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
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.spacing_violations > 3);
    }

    SECTION("star")
    {
        ExPolygon star(star_points(45.0, 22.0, 7));
        const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(star);
        REQUIRE_FALSE(validation.ok);
        REQUIRE(validation.spacing_violations > 3);
        REQUIRE(validation.turnback_violations > 0);
    }
}

TEST_CASE("Continuous Fermat rejects unsupported hole topologies", "[ContinuousFermat][unsupported]")
{
    SECTION("annulus")
    {
        ExPolygon annulus(circle_points(42.0, 160));
        annulus.holes.emplace_back(circle_points(16.0, 96, Point(0, 0), false));
        const ContinuousFermat::PathValidation validation = validate_generated_continuous_fermat_path(annulus);
        REQUIRE_FALSE(validation.ok);
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
        REQUIRE_FALSE(validation.ok);
        REQUIRE_FALSE(validation.closed);
    }
}
