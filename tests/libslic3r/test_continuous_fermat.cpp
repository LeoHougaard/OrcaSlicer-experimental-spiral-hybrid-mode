#include <catch2/catch_all.hpp>

#include "libslic3r/ContinuousFermat.hpp"
#include "libslic3r/ExPolygon.hpp"
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
    for (size_t i = 0; i < lines.size(); ++i) {
        for (size_t j = i + 2; j < lines.size(); ++j) {
            if (i == 0 && j + 1 == lines.size())
                continue;
            Point intersection;
            if (lines[i].intersection(lines[j], &intersection)) {
                if (message != nullptr) {
                    std::ostringstream out;
                    out << "crossing segments " << i << " and " << j << " at " << us(intersection.x()) << ", " <<
                        us(intersection.y()) << "; segment " << i << ": " << us(lines[i].a.x()) << ", " <<
                        us(lines[i].a.y()) << " -> " << us(lines[i].b.x()) << ", " << us(lines[i].b.y()) <<
                        "; segment " << j << ": " << us(lines[j].a.x()) << ", " << us(lines[j].a.y()) << " -> " <<
                        us(lines[j].b.x()) << ", " << us(lines[j].b.y());
                    out << "; first points:";
                    for (size_t k = 0; k < std::min<size_t>(path.points.size(), 12); ++k)
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
    Flow flow(1.2f, 0.2f, 0.4f);
    Polyline path = ContinuousFermat::generate_layer_path({ area }, flow);

    require_path_inside(area, path);
    REQUIRE((path.points.front() - path.points.back()).cast<double>().norm() <= double(flow.scaled_spacing()) * 0.20);
    std::string crossing;
    const bool has_crossing = segments_cross_nonlocal(path, &crossing);
    INFO(crossing);
    REQUIRE_FALSE(has_crossing);
    return path;
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
        p(-30.0, -18.0),
        p( 30.0, -18.0),
        p( 30.0,  18.0),
        p(-30.0,  18.0),
    });
    Polyline path = require_continuous_fermat_path(rectangle);
    Flow flow(1.2f, 0.2f, 0.4f);
    require_rectangle_outer_wall_corners(path, flow, -30.0, -18.0, 30.0, 18.0);
}

TEST_CASE("Continuous Fermat generates one non-crossing path for a square hole", "[ContinuousFermat]")
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
    Polyline path = require_continuous_fermat_path(square_with_hole);

    Flow flow(1.2f, 0.2f, 0.4f);
    const double gap_tolerance = double(flow.scaled_spacing()) * 1.6;
    for (const Point &gap_point : { p(-22.0, 0.0), p(22.0, 0.0) })
        REQUIRE(path_distance_to_point(path, gap_point) <= gap_tolerance);
}

TEST_CASE("Continuous Fermat handles concave and branched island shapes", "[ContinuousFermat]")
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
        require_continuous_fermat_path(c_shape);
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
        require_continuous_fermat_path(dumbbell);
    }

    SECTION("star")
    {
        ExPolygon star(star_points(45.0, 22.0, 7));
        require_continuous_fermat_path(star);
    }
}

TEST_CASE("Continuous Fermat handles island shapes with holes", "[ContinuousFermat]")
{
    SECTION("annulus")
    {
        ExPolygon annulus(circle_points(42.0, 160));
        annulus.holes.emplace_back(circle_points(16.0, 96, Point(0, 0), false));
        require_continuous_fermat_path(annulus);
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
        Polyline path = require_continuous_fermat_path(two_holes);

        Flow flow(1.2f, 0.2f, 0.4f);
        const double pocket_tolerance = double(flow.scaled_spacing()) * 2.2;
        for (const Point &pocket : { p(-39.0, -17.0), p(-39.0, 17.0), p(39.0, -17.0), p(39.0, 17.0) })
            REQUIRE(path_distance_to_point(path, pocket) <= pocket_tolerance);
    }
}
