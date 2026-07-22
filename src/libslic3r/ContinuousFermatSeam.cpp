#include "ContinuousFermatSeam.hpp"

#include "ClipperUtils.hpp"
#include "Geometry.hpp"
#include "Line.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>

namespace Slic3r::ContinuousFermat {
namespace {

struct SeamCandidate
{
    size_t loop_a { 0 };
    size_t loop_b { 0 };
    size_t vertex_a { 0 };
    size_t vertex_b { 0 };
    size_t edge_a { 0 };
    size_t edge_b { 0 };
    Point point_a;
    Point point_b;
    double squared_length { std::numeric_limits<double>::infinity() };
    double tangent_mismatch { std::numeric_limits<double>::infinity() };
};

double squared_distance(const Point &a, const Point &b)
{
    return (a - b).cast<double>().squaredNorm();
}

bool candidate_less(const SeamCandidate &a, const SeamCandidate &b)
{
    return std::tie(
               a.squared_length,
               a.tangent_mismatch,
               a.loop_a,
               a.loop_b,
               a.point_a.x(),
               a.point_a.y(),
               a.point_b.x(),
               a.point_b.y(),
               a.edge_a,
               a.edge_b) <
           std::tie(
               b.squared_length,
               b.tangent_mismatch,
               b.loop_a,
               b.loop_b,
               b.point_a.x(),
               b.point_a.y(),
               b.point_b.x(),
               b.point_b.y(),
               b.edge_a,
               b.edge_b);
}

bool seam_is_contained(const ExPolygon &area, const Point &a, const Point &b)
{
    return a != b && area.contains(Line(a, b));
}

SeamCandidate shortest_visible_seam(
    const ExPolygon &area,
    const std::vector<const Polygon *> &loops,
    const size_t loop_a,
    const size_t loop_b)
{
    const Points &a = loops[loop_a]->points;
    const Points &b = loops[loop_b]->points;
    SeamCandidate nearest;
    std::vector<SeamCandidate> candidates;
    candidates.reserve(2 * a.size() * b.size());
    const auto consider = [&](SeamCandidate candidate) {
        const Vec2d tangent_a = (a[(candidate.edge_a + 1) % a.size()] - a[candidate.edge_a]).cast<double>().normalized();
        const Vec2d tangent_b = (b[(candidate.edge_b + 1) % b.size()] - b[candidate.edge_b]).cast<double>().normalized();
        // Hole traversal is reversed when the annulus is opened into a strip.
        // Among equal-length cuts, align the resulting strip edges so the
        // rectangle-to-ring map does not twist immediately beside the seam.
        candidate.tangent_mismatch = 1.0 + tangent_a.dot(tangent_b);
        candidates.emplace_back(candidate);
        if (candidate_less(candidate, nearest))
            nearest = std::move(candidate);
    };

    // For two disjoint line segments their closest pair contains at least one
    // endpoint. Project every vertex of each loop onto every edge of the
    // other. This finds the actual boundary-to-boundary minimum instead of the
    // old vertex-only approximation, which produced diagonal cuts through
    // otherwise parallel square or chamfered walls.
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j) {
            const size_t next_b = (j + 1) % b.size();
            Point projected;
            const double d2 = Line(b[j], b[next_b]).distance_to_squared(a[i], &projected);
            consider({ loop_a, loop_b, i, j, i, j, a[i], projected, d2 });
        }
    for (size_t j = 0; j < b.size(); ++j)
        for (size_t i = 0; i < a.size(); ++i) {
            const size_t next_a = (i + 1) % a.size();
            Point projected;
            const double d2 = Line(a[i], a[next_a]).distance_to_squared(b[j], &projected);
            consider({ loop_a, loop_b, i, j, i, (j + b.size() - 1) % b.size(), projected, b[j], d2 });
        }
    if (std::isfinite(nearest.squared_length) &&
        seam_is_contained(area, nearest.point_a, nearest.point_b))
        return nearest;

    // Another hole may obstruct the nearest pair. Sort the remaining pairs by
    // length and stop at the first visible one; this is exact, deterministic,
    // and only pays for clipping in the uncommon obstructed case.
    std::sort(candidates.begin(), candidates.end(), candidate_less);
    for (const SeamCandidate &candidate : candidates)
        if (seam_is_contained(area, candidate.point_a, candidate.point_b))
            return candidate;
    return {};
}

class DisjointSet
{
public:
    explicit DisjointSet(const size_t count) : m_parent(count), m_rank(count, 0)
    {
        std::iota(m_parent.begin(), m_parent.end(), size_t(0));
    }

    size_t find(const size_t value)
    {
        if (m_parent[value] != value)
            m_parent[value] = find(m_parent[value]);
        return m_parent[value];
    }

    bool unite(const size_t a, const size_t b)
    {
        size_t root_a = find(a);
        size_t root_b = find(b);
        if (root_a == root_b)
            return false;
        if (m_rank[root_a] < m_rank[root_b])
            std::swap(root_a, root_b);
        m_parent[root_b] = root_a;
        if (m_rank[root_a] == m_rank[root_b])
            ++m_rank[root_a];
        return true;
    }

private:
    std::vector<size_t> m_parent;
    std::vector<size_t> m_rank;
};

bool strictly_cross(const Point &a, const Point &b, const Point &c, const Point &d)
{
    if (a == c || a == d || b == c || b == d)
        return false;
    return int(Geometry::orient(a, b, c)) * int(Geometry::orient(a, b, d)) < 0 &&
           int(Geometry::orient(c, d, a)) * int(Geometry::orient(c, d, b)) < 0;
}

void append_distinct(Points &out, const Point &point);

std::vector<Polygon> insert_seam_points(
    const std::vector<const Polygon *> &source,
    const std::vector<SeamCandidate> &tree)
{
    struct Insertion
    {
        size_t edge;
        Point point;
    };
    std::vector<std::vector<Insertion>> insertions(source.size());
    for (const SeamCandidate &seam : tree) {
        insertions[seam.loop_a].push_back({ seam.edge_a, seam.point_a });
        insertions[seam.loop_b].push_back({ seam.edge_b, seam.point_b });
    }

    std::vector<Polygon> result;
    result.reserve(source.size());
    for (size_t loop_index = 0; loop_index < source.size(); ++loop_index) {
        const Points &ring = source[loop_index]->points;
        Points expanded;
        expanded.reserve(ring.size() + insertions[loop_index].size());
        for (size_t edge = 0; edge < ring.size(); ++edge) {
            append_distinct(expanded, ring[edge]);
            const Point &next = ring[(edge + 1) % ring.size()];
            std::vector<Point> on_edge;
            for (const Insertion &insertion : insertions[loop_index])
                if (insertion.edge == edge && insertion.point != ring[edge] && insertion.point != next)
                    on_edge.emplace_back(insertion.point);
            std::sort(on_edge.begin(), on_edge.end(), [&](const Point &lhs, const Point &rhs) {
                return squared_distance(ring[edge], lhs) < squared_distance(ring[edge], rhs);
            });
            for (const Point &point : on_edge)
                append_distinct(expanded, point);
        }
        result.emplace_back(std::move(expanded));
    }
    return result;
}

size_t point_index(const Polygon &polygon, const Point &point)
{
    const auto found = std::find(polygon.points.begin(), polygon.points.end(), point);
    return found == polygon.points.end() ? 0 : size_t(found - polygon.points.begin());
}

struct TreeEdge
{
    size_t parent { 0 };
    size_t child { 0 };
    size_t parent_vertex { 0 };
    size_t child_vertex { 0 };
};

void append_distinct(Points &out, const Point &point)
{
    if (out.empty() || out.back() != point)
        out.emplace_back(point);
}

Points expand_boundary_tree(
    const size_t node,
    const size_t start_vertex,
    const std::vector<const Polygon *> &loops,
    const std::vector<std::vector<TreeEdge>> &children)
{
    const Points &ring = loops[node]->points;
    if (ring.empty())
        return {};

    std::vector<std::vector<TreeEdge>> at_vertex(ring.size());
    for (const TreeEdge &edge : children[node])
        at_vertex[edge.parent_vertex].emplace_back(edge);
    for (std::vector<TreeEdge> &edges : at_vertex)
        std::sort(edges.begin(), edges.end(), [](const TreeEdge &a, const TreeEdge &b) {
            return std::tie(a.child, a.child_vertex) < std::tie(b.child, b.child_vertex);
        });

    Points result;
    for (size_t offset = 0; offset < ring.size(); ++offset) {
        const size_t vertex = (start_vertex + offset) % ring.size();
        append_distinct(result, ring[vertex]);
        for (const TreeEdge &edge : at_vertex[vertex]) {
            const Points child = expand_boundary_tree(edge.child, edge.child_vertex, loops, children);
            for (const Point &point : child)
                append_distinct(result, point);
            // The return side of the zero-width seam is intentionally retained.
            append_distinct(result, ring[vertex]);
        }
    }
    append_distinct(result, ring[start_vertex]);
    return result;
}

ExPolygons build_planning_domain(const ExPolygon &area, const std::vector<CutSeam> &seams)
{
    if (seams.empty())
        return { area };

    Polylines centerlines;
    centerlines.reserve(seams.size());
    for (const CutSeam &seam : seams)
        centerlines.emplace_back(Points { seam.parent_point, seam.child_point });

    // The slit is a planning boundary, so it must be wider than Clipper's
    // 0.0001 mm cleanup epsilon or a connector may numerically jump across it.
    // Even the last fallback below is only 0.008192 mm: far below extrusion
    // width, while the normal half-line-width inset creates the physical
    // centerline clearance on both sides. The open-square caps cross both
    // incident boundaries so no residual point contact can preserve a hole.
    static constexpr std::array<float, 8> SLIT_RADII {
        256.0f, 512.0f, 1024.0f, 1536.0f, 2048.0f, 3072.0f, 4096.0f, 8192.0f,
    };
    for (const float radius : SLIT_RADII) {
        const Polygons slits = offset(
            centerlines, radius, ClipperLib::jtMiter, 2.0, ClipperLib::etOpenSquare);
        if (slits.empty())
            continue;
        ExPolygons domain = diff_ex(area, slits);
        if (domain.size() == 1 && domain.front().holes.empty())
            return domain;
    }
    return {};
}

} // namespace

HoleRemovalResult remove_holes_topologically(const ExPolygon &area)
{
    HoleRemovalResult result;
    if (area.contour.points.size() < 3) {
        result.reason = "outer contour has fewer than three vertices";
        return result;
    }

    std::vector<const Polygon *> loops;
    loops.reserve(area.holes.size() + 1);
    loops.emplace_back(&area.contour);
    for (const Polygon &hole : area.holes) {
        if (hole.points.size() < 3) {
            result.reason = "hole has fewer than three vertices";
            return result;
        }
        loops.emplace_back(&hole);
    }

    if (loops.size() == 1) {
        Points boundary = area.contour.points;
        boundary.emplace_back(boundary.front());
        result.boundary = Polyline(std::move(boundary));
        result.planning_domain = { area };
        result.ok = true;
        return result;
    }

    std::vector<SeamCandidate> graph_edges;
    graph_edges.reserve(loops.size() * (loops.size() - 1) / 2);
    for (size_t a = 0; a < loops.size(); ++a)
        for (size_t b = a + 1; b < loops.size(); ++b) {
            SeamCandidate candidate = shortest_visible_seam(area, loops, a, b);
            if (std::isfinite(candidate.squared_length))
                graph_edges.emplace_back(candidate);
        }
    std::sort(graph_edges.begin(), graph_edges.end(), candidate_less);

    DisjointSet components(loops.size());
    std::vector<SeamCandidate> tree;
    tree.reserve(loops.size() - 1);
    for (const SeamCandidate &candidate : graph_edges) {
        if (components.find(candidate.loop_a) == components.find(candidate.loop_b))
            continue;
        const Point &a = candidate.point_a;
        const Point &b = candidate.point_b;
        bool crosses_tree = false;
        for (const SeamCandidate &accepted : tree) {
            const Point &c = accepted.point_a;
            const Point &d = accepted.point_b;
            if (strictly_cross(a, b, c, d)) {
                crosses_tree = true;
                break;
            }
        }
        if (crosses_tree || !components.unite(candidate.loop_a, candidate.loop_b))
            continue;
        tree.emplace_back(candidate);
        if (tree.size() + 1 == loops.size())
            break;
    }
    if (tree.size() + 1 != loops.size()) {
        result.reason = "could not connect every hole with non-crossing contained seams";
        return result;
    }

    std::vector<Polygon> expanded_loops = insert_seam_points(loops, tree);
    loops.clear();
    loops.reserve(expanded_loops.size());
    for (const Polygon &loop : expanded_loops)
        loops.emplace_back(&loop);

    struct AdjacentEdge
    {
        size_t other;
        size_t here_vertex;
        size_t other_vertex;
    };
    std::vector<std::vector<AdjacentEdge>> adjacency(loops.size());
    for (const SeamCandidate &edge : tree) {
        const size_t vertex_a = point_index(expanded_loops[edge.loop_a], edge.point_a);
        const size_t vertex_b = point_index(expanded_loops[edge.loop_b], edge.point_b);
        adjacency[edge.loop_a].push_back({ edge.loop_b, vertex_a, vertex_b });
        adjacency[edge.loop_b].push_back({ edge.loop_a, vertex_b, vertex_a });
    }

    std::vector<std::vector<TreeEdge>> children(loops.size());
    std::vector<bool> visited(loops.size(), false);
    const auto orient_tree = [&](const auto &self, const size_t node) -> void {
        visited[node] = true;
        std::sort(adjacency[node].begin(), adjacency[node].end(), [](const AdjacentEdge &a, const AdjacentEdge &b) {
            return std::tie(a.here_vertex, a.other, a.other_vertex) <
                   std::tie(b.here_vertex, b.other, b.other_vertex);
        });
        for (const AdjacentEdge &edge : adjacency[node]) {
            if (visited[edge.other])
                continue;
            children[node].push_back({ node, edge.other, edge.here_vertex, edge.other_vertex });
            result.seams.push_back({
                node,
                edge.other,
                edge.here_vertex,
                edge.other_vertex,
                loops[node]->points[edge.here_vertex],
                loops[edge.other]->points[edge.other_vertex],
                std::sqrt(squared_distance(
                    loops[node]->points[edge.here_vertex], loops[edge.other]->points[edge.other_vertex])),
            });
            self(self, edge.other);
        }
    };
    orient_tree(orient_tree, 0);

    Points boundary = expand_boundary_tree(0, 0, loops, children);
    if (boundary.size() < 4 || boundary.front() != boundary.back()) {
        result.reason = "expanded cut boundary is not exactly closed";
        result.seams.clear();
        return result;
    }
    // The only new geometry is the seam tree, and each seam was containment-
    // checked before selection. Do not re-clip the expanded weak boundary:
    // Clipper classifies some exact input hole edges as outside and may cancel
    // coincident opposite seam sides. All remaining edges are copied verbatim
    // from the already valid input rings.

    result.boundary = Polyline(std::move(boundary));
    result.planning_domain = build_planning_domain(area, result.seams);
    if (result.planning_domain.empty()) {
        result.reason = "could not materialize the zero-width seam tree as a hole-free planning domain";
        result.boundary.clear();
        result.seams.clear();
        return result;
    }
    result.ok = true;
    return result;
}

} // namespace Slic3r::ContinuousFermat
