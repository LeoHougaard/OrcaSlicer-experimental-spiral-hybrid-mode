#!/usr/bin/env python3
"""Standalone contour-tree planner for the continuous path lab.

This module intentionally depends only on Python's standard library plus the
existing dependency-free prototype in tools/continuous_fermat.  The planner here
is not the production algorithm; it is a small, inspectable v2 test bed:

    polygon -> SDF iso-contours -> contour containment tree -> planned DFS route

The important property is that all branches are known before the path is
emitted.  The router visits child subtrees deliberately and only finishes in one
chosen final subtree, which prevents local "spiral until trapped" behavior.
"""

from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence


FERMAT_DIR = Path(__file__).resolve().parents[1] / "continuous_fermat"
if str(FERMAT_DIR) not in sys.path:
    sys.path.insert(0, str(FERMAT_DIR))

import fermat_layer as cf  # type: ignore  # noqa: E402


Point = tuple[float, float]
EPS = 1e-9


@dataclass
class TreeNode:
    id: int
    loop: cf.ContourLoop
    parent: int | None = None
    children: list[int] = field(default_factory=list)


@dataclass
class RouteDebugEdge:
    kind: str
    parent: int
    child: int
    a: Point
    b: Point


@dataclass
class ConnectorPair:
    parent_enter_idx: int
    child_enter_idx: int
    child_return_idx: int
    parent_return_idx: int
    parent_enter_rail: Point
    child_enter_rail: Point


@dataclass(frozen=True)
class BarrierModel:
    name: str
    outer: list[Point]
    holes: list[list[Point]]
    barriers: list[tuple[Point, Point]]
    entry_path: list[Point] = field(default_factory=list)
    barrier_radius: float = 0.0

    def bounds(self, margin: float = 0.0) -> tuple[float, float, float, float]:
        pts = list(self.outer)
        for hole in self.holes:
            pts.extend(hole)
        for a, b in self.barriers:
            pts.extend([a, b])
        return (
            min(p[0] for p in pts) - margin,
            min(p[1] for p in pts) - margin,
            max(p[0] for p in pts) + margin,
            max(p[1] for p in pts) + margin,
        )

    def contains(self, p: Point) -> bool:
        return cf.point_in_polygon(p, self.outer) and not any(cf.point_in_polygon(p, hole) for hole in self.holes)

    def sdf(self, p: Point) -> float:
        boundary_distance = cf.distance_to_loop(p, self.outer)
        for hole in self.holes:
            boundary_distance = min(boundary_distance, cf.distance_to_loop(p, hole))
        for a, b in self.barriers:
            boundary_distance = min(boundary_distance, cf.point_segment_distance(p, a, b) - self.barrier_radius)
        if boundary_distance < 0.0 and self.contains(p):
            return boundary_distance
        return boundary_distance if self.contains(p) else -abs(boundary_distance)


def v_add(a: Point, b: Point) -> Point:
    return (a[0] + b[0], a[1] + b[1])


def v_sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1])


def v_mul(a: Point, scale: float) -> Point:
    return (a[0] * scale, a[1] * scale)


def v_len(a: Point) -> float:
    return math.hypot(a[0], a[1])


def v_unit(a: Point) -> Point:
    length = v_len(a)
    if length <= EPS:
        return (1.0, 0.0)
    return (a[0] / length, a[1] / length)


def v_perp(a: Point) -> Point:
    return (-a[1], a[0])


def append_point(path: list[Point], p: Point) -> None:
    if not path or cf.dist(path[-1], p) > EPS:
        path.append(p)


def append_points(path: list[Point], points: Iterable[Point]) -> None:
    for point in points:
        append_point(path, point)


def closed_without_duplicate(points: Sequence[Point]) -> list[Point]:
    out = list(points)
    if len(out) > 1 and cf.dist(out[0], out[-1]) <= EPS:
        out.pop()
    return out


def loop_vertex_index(points: Sequence[Point], p: Point) -> int:
    if not points:
        return 0
    best = 0
    best_d2 = float("inf")
    for idx, q in enumerate(points):
        d2 = cf.dist2(p, q)
        if d2 < best_d2:
            best = idx
            best_d2 = d2
    return best


def loop_average_step(points: Sequence[Point]) -> float:
    if len(points) < 2:
        return 1.0
    return max(cf.polyline_length(points, closed=True) / len(points), EPS)


def port_gap_points(points: Sequence[Point], spacing: float) -> int:
    return max(2, min(max(2, len(points) // 5), int(round(spacing / loop_average_step(points)))))


def advance_index(points: Sequence[Point], idx: int, direction: int, count: int) -> int:
    if not points:
        return 0
    return (idx + direction * count) % len(points)


def forward_distance_indices(points: Sequence[Point], start: int, idx: int, direction: int) -> int:
    n = len(points)
    if n == 0:
        return 0
    if direction >= 0:
        return (idx - start) % n
    return (start - idx) % n


def arc_points(points: Sequence[Point], start: int, end: int, direction: int) -> list[Point]:
    if not points:
        return []
    n = len(points)
    start %= n
    end %= n
    out = [points[start]]
    idx = start
    guard = 0
    while idx != end and guard <= n + 2:
        idx = (idx + direction) % n
        out.append(points[idx])
        guard += 1
    return out


def point_inside_loop(point: Point, loop: cf.ContourLoop) -> bool:
    try:
        return cf.point_in_polygon(point, loop.points)
    except Exception:
        return False


def ray_segment_t(origin: Point, direction: Point, a: Point, b: Point) -> float | None:
    sx = b[0] - a[0]
    sy = b[1] - a[1]
    denom = direction[0] * sy - direction[1] * sx
    if abs(denom) <= EPS:
        return None
    ax = a[0] - origin[0]
    ay = a[1] - origin[1]
    t = (ax * sy - ay * sx) / denom
    u = (ax * direction[1] - ay * direction[0]) / denom
    if t >= -EPS and -EPS <= u <= 1.0 + EPS:
        return max(0.0, t)
    return None


def ray_polygon_hits(origin: Point, direction: Point, polygon: Sequence[Point]) -> list[float]:
    hits: list[float] = []
    for a, b in zip(polygon, [*polygon[1:], polygon[0]]):
        t = ray_segment_t(origin, direction, a, b)
        if t is None:
            continue
        if not any(abs(t - existing) <= 1e-5 for existing in hits):
            hits.append(t)
    hits.sort()
    return hits


def point_along_ray(origin: Point, direction: Point, t: float) -> Point:
    return (origin[0] + direction[0] * t, origin[1] + direction[1] * t)


def barrier_clear_of_other_holes(
    model: cf.PolygonModel,
    hole_index: int,
    a: Point,
    b: Point,
    samples: int = 40,
) -> bool:
    for i in range(1, samples):
        p = cf.lerp(a, b, i / samples)
        if not cf.point_in_polygon(p, model.outer):
            return False
        for idx, hole in enumerate(model.holes):
            if idx != hole_index and cf.point_in_polygon(p, hole):
                return False
    return True


def barrier_for_hole(model: cf.PolygonModel, hole_index: int) -> tuple[Point, Point] | None:
    hole = model.holes[hole_index]
    center = cf.polygon_centroid(hole)
    directions: list[Point] = [
        (1.0, 0.0),
        (-1.0, 0.0),
        (0.0, 1.0),
        (0.0, -1.0),
        v_unit((1.0, 1.0)),
        v_unit((-1.0, 1.0)),
        v_unit((1.0, -1.0)),
        v_unit((-1.0, -1.0)),
    ]
    candidates: list[tuple[float, Point, Point]] = []

    for direction in directions:
        hole_hits = [t for t in ray_polygon_hits(center, direction, hole) if t > 1e-5]
        outer_hits = [t for t in ray_polygon_hits(center, direction, model.outer) if t > 1e-5]
        if not hole_hits or not outer_hits:
            continue
        start_t = min(hole_hits)
        end_candidates = [t for t in outer_hits if t > start_t + 1e-5]
        if not end_candidates:
            continue
        end_t = min(end_candidates)
        a = point_along_ray(center, direction, start_t)
        b = point_along_ray(center, direction, end_t)
        if barrier_clear_of_other_holes(model, hole_index, a, b):
            candidates.append((end_t - start_t, a, b))

    if not candidates:
        return None
    _, a, b = min(candidates, key=lambda item: item[0])
    return a, b


def hole_boundary_point(center: Point, hole: Sequence[Point], direction: Point) -> Point | None:
    hits = [t for t in ray_polygon_hits(center, direction, hole) if t > 1e-5]
    if not hits:
        return None
    return point_along_ray(center, direction, min(hits))


def outer_boundary_point(center: Point, outer: Sequence[Point], direction: Point, min_t: float = 0.0) -> Point | None:
    hits = [t for t in ray_polygon_hits(center, direction, outer) if t > min_t + 1e-5]
    if not hits:
        return None
    return point_along_ray(center, direction, min(hits))


def boundary_arc_between(points: Sequence[Point], a: Point, b: Point) -> list[Point]:
    if not points:
        return [a, b]

    start_idx = loop_vertex_index(points, a)
    end_idx = loop_vertex_index(points, b)
    candidates: list[list[Point]] = []
    for direction in (1, -1):
        candidate = [a]
        append_points(candidate, arc_points(points, start_idx, end_idx, direction))
        append_point(candidate, b)
        candidates.append(candidate)
    return min(candidates, key=lambda path: cf.polyline_length(path, closed=False))


def shared_hole_arc(model: cf.PolygonModel, a: Point, b: Point) -> list[Point] | None:
    tolerance = 1e-4
    for hole in model.holes:
        if cf.distance_to_loop(a, hole) <= tolerance and cf.distance_to_loop(b, hole) <= tolerance:
            return boundary_arc_between(hole, a, b)
    return None


def barrier_entry_path(model: cf.PolygonModel, barriers: Sequence[tuple[Point, Point]]) -> list[Point]:
    if not barriers:
        return []

    path = [barriers[0][0]]
    append_point(path, barriers[0][1])
    for start, end in barriers[1:]:
        arc = shared_hole_arc(model, path[-1], start)
        if arc:
            append_points(path, arc)
        else:
            append_point(path, start)
        append_point(path, end)
    return path


def barrier_chain_for_holes(model: cf.PolygonModel) -> list[tuple[Point, Point]] | None:
    if len(model.holes) < 2:
        return None

    centers = [cf.polygon_centroid(hole) for hole in model.holes]
    order = sorted(range(len(model.holes)), key=lambda idx: (centers[idx][0], centers[idx][1]))
    barriers: list[tuple[Point, Point]] = []
    last_direction: Point | None = None

    for prev_idx, next_idx in zip(order, order[1:]):
        prev_center = centers[prev_idx]
        next_center = centers[next_idx]
        direction = v_unit(v_sub(next_center, prev_center))
        a = hole_boundary_point(prev_center, model.holes[prev_idx], direction)
        b = hole_boundary_point(next_center, model.holes[next_idx], v_mul(direction, -1.0))
        if a is None or b is None:
            return None
        if not barrier_clear_of_other_holes(model, prev_idx, a, b):
            # The endpoint is allowed to touch the destination hole, but samples
            # through a third hole would make the virtual cut invalid.
            blocked_by_third = False
            for sample in range(1, 40):
                p = cf.lerp(a, b, sample / 40)
                for idx, hole in enumerate(model.holes):
                    if idx not in {prev_idx, next_idx} and cf.point_in_polygon(p, hole):
                        blocked_by_third = True
                        break
                if blocked_by_third:
                    break
            if blocked_by_third:
                return None
        barriers.append((a, b))
        last_direction = direction

    last_idx = order[-1]
    last_center = centers[last_idx]
    exit_direction = last_direction or (1.0, 0.0)
    exit_start = hole_boundary_point(last_center, model.holes[last_idx], exit_direction)
    if exit_start is None:
        return None
    min_exit_t = cf.dist(last_center, exit_start)
    exit_end = outer_boundary_point(last_center, model.outer, exit_direction, min_t=min_exit_t)
    if exit_end is None:
        return None
    if not barrier_clear_of_other_holes(model, last_idx, exit_start, exit_end):
        return None
    barriers.append((exit_start, exit_end))
    return barriers


def barrier_model_for_holes(model: cf.PolygonModel, barrier_radius: float = 0.0) -> tuple[BarrierModel | cf.PolygonModel, list[str]]:
    if not model.holes:
        return model, []

    barriers: list[tuple[Point, Point]] = []
    diagnostics: list[str] = []
    if len(model.holes) >= 2:
        chain = barrier_chain_for_holes(model)
        if chain:
            diagnostics.append(f"v3 inserted a {len(chain)}-segment virtual hole chain for planning.")
            return BarrierModel(
                f"{model.name}_v3_barrier",
                model.outer,
                model.holes,
                chain,
                barrier_entry_path(model, chain),
                barrier_radius,
            ), diagnostics

    for idx, _hole in enumerate(model.holes):
        barrier = barrier_for_hole(model, idx)
        if barrier is None:
            diagnostics.append(f"v3 could not find a clear virtual barrier for hole {idx}.")
            continue
        barriers.append(barrier)

    if not barriers:
        return model, diagnostics

    diagnostics.append(f"v3 inserted {len(barriers)} virtual hole barrier(s) for planning.")
    entry_path = barrier_entry_path(model, barriers) if len(barriers) == 1 else []
    return BarrierModel(f"{model.name}_v3_barrier", model.outer, model.holes, barriers, entry_path, barrier_radius), diagnostics


def build_contour_tree(contours: Sequence[cf.ContourLoop]) -> tuple[list[TreeNode], list[str]]:
    diagnostics: list[str] = []
    ordered = sorted(
        enumerate(contours),
        key=lambda item: (item[1].level_index, -item[1].area, item[1].centroid[0], item[1].centroid[1]),
    )
    nodes = [TreeNode(id=i, loop=loop) for i, loop in enumerate(contours)]

    for original_idx, loop in ordered:
        candidates: list[tuple[int, float, float, int]] = []
        for parent_idx, parent in enumerate(contours):
            if parent_idx == original_idx:
                continue
            if parent.area <= loop.area:
                continue
            if parent.level_index > loop.level_index:
                continue
            if not point_inside_loop(loop.centroid, parent):
                continue
            level_gap = loop.level_index - parent.level_index
            boundary_gap = cf.distance_to_loop(loop.centroid, parent.points)
            candidates.append((level_gap, parent.area, boundary_gap, parent_idx))

        if candidates:
            # Prefer the closest containing loop, then the smallest containing area.
            _, _, selected_boundary_gap, parent_idx = min(candidates, key=lambda item: (item[0], item[1], item[2]))
            nearby: list[tuple[int, float, float, float, int]] = []
            for nearby_idx, nearby_parent in enumerate(contours):
                if nearby_idx == original_idx:
                    continue
                if nearby_parent.area <= loop.area or nearby_parent.level_index > loop.level_index:
                    continue
                level_gap = loop.level_index - nearby_parent.level_index
                if level_gap > 3:
                    continue
                boundary_gap = cf.distance_to_loop(loop.centroid, nearby_parent.points)
                centroid_gap = cf.dist(loop.centroid, nearby_parent.centroid)
                nearby.append((level_gap, boundary_gap, centroid_gap, nearby_parent.area, nearby_idx))
            if nearby:
                level_gap, boundary_gap, centroid_gap, _, nearby_idx = min(nearby, key=lambda item: (item[1], item[0], item[2], item[3]))
                nearby_parent = contours[nearby_idx]
                offset_step = abs(loop.offset - nearby_parent.offset)
                spacing_hint = max(offset_step, loop.offset / max(1, loop.level_index + 1), 1.0)
                if (
                    level_gap >= 1
                    and boundary_gap <= max(spacing_hint * 1.75, 3.0)
                    and boundary_gap < selected_boundary_gap * 0.80
                ):
                    parent_idx = nearby_idx
            nodes[original_idx].parent = parent_idx
            nodes[parent_idx].children.append(original_idx)
            continue

        nearest: list[tuple[int, float, float, float, int]] = []
        for parent_idx, parent in enumerate(contours):
            if parent_idx == original_idx:
                continue
            if parent.area <= loop.area or parent.level_index > loop.level_index:
                continue
            level_gap = loop.level_index - parent.level_index
            if level_gap > 3:
                continue
            boundary_gap = cf.distance_to_loop(loop.centroid, parent.points)
            centroid_gap = cf.dist(loop.centroid, parent.centroid)
            nearest.append((level_gap, centroid_gap, boundary_gap, parent.area, parent_idx))
        if nearest:
            level_gap, centroid_gap, boundary_gap, _, parent_idx = min(nearest, key=lambda item: (item[0], item[1], item[2], item[3]))
            parent = contours[parent_idx]
            offset_step = abs(loop.offset - parent.offset)
            spacing_hint = max(offset_step, loop.offset / max(1, loop.level_index + 1), 1.0)
            adjacent_front = level_gap == 1 and centroid_gap <= max(spacing_hint * 3.0, 8.0)
            nearby_branch = boundary_gap <= max(spacing_hint * 1.75, 3.0)
            if level_gap >= 1 and (adjacent_front or nearby_branch):
                nodes[original_idx].parent = parent_idx
                nodes[parent_idx].children.append(original_idx)

    roots = [node.id for node in nodes if node.parent is None]
    if len(roots) > 1:
        diagnostics.append(
            f"Detected {len(roots)} independent contour root(s); v2 lab routes only the largest root because strict mode forbids islands."
        )

    for node in nodes:
        node.children.sort(
            key=lambda child_id: (
                node.loop.level_index,
                loop_vertex_index(node.loop.points, nodes[child_id].loop.centroid),
                -nodes[child_id].loop.area,
            )
        )

    return nodes, diagnostics


def root_node_id(nodes: Sequence[TreeNode]) -> int | None:
    roots = [node for node in nodes if node.parent is None]
    if not roots:
        return None
    return max(roots, key=lambda node: node.loop.area).id


def children_in_route_order(nodes: Sequence[TreeNode], node_id: int, start_idx: int, direction: int) -> list[int]:
    node = nodes[node_id]
    points = node.loop.points
    return sorted(
        node.children,
        key=lambda child_id: forward_distance_indices(
            points,
            start_idx,
            loop_vertex_index(points, nodes[child_id].loop.centroid),
            direction,
        ),
    )


def connector_pair_for_child(
    parent_points: Sequence[Point],
    child_points: Sequence[Point],
    child_centroid: Point,
    direction: int,
    spacing: float,
) -> ConnectorPair:
    """Create paired rails between a parent loop and child loop.

    The return rail is the straight child_return -> parent_return segment.  The
    entry rail is generated as an exact spacing offset from that segment, then
    snapped back to the nearest contour ports with short lead-in segments.
    """

    center_idx = loop_vertex_index(parent_points, child_centroid)
    gap = port_gap_points(parent_points, spacing)
    parent_enter_idx = advance_index(parent_points, center_idx, -direction, gap)
    parent_return_idx = advance_index(parent_points, center_idx, direction, gap)
    parent_enter = parent_points[parent_enter_idx]
    parent_return = parent_points[parent_return_idx]

    child_return_idx = loop_vertex_index(child_points, parent_return)
    child_return = child_points[child_return_idx]
    rail_axis = v_unit(v_sub(parent_return, child_return))
    normal = v_perp(rail_axis)

    candidates: list[tuple[float, Point, Point, int]] = []
    for sign in (-1.0, 1.0):
        offset = v_mul(normal, sign * spacing)
        parent_enter_rail = v_add(parent_return, offset)
        child_enter_rail = v_add(child_return, offset)
        child_enter_idx = loop_vertex_index(child_points, child_enter_rail)
        score = (
            cf.dist(parent_enter, parent_enter_rail)
            + cf.dist(child_points[child_enter_idx], child_enter_rail)
            + abs(cf.dist(parent_enter, parent_return) - spacing)
            + abs(cf.dist(child_points[child_enter_idx], child_return) - spacing)
        )
        candidates.append((score, parent_enter_rail, child_enter_rail, child_enter_idx))

    _, parent_enter_rail, child_enter_rail, child_enter_idx = min(candidates, key=lambda item: item[0])
    return ConnectorPair(
        parent_enter_idx=parent_enter_idx,
        child_enter_idx=child_enter_idx,
        child_return_idx=child_return_idx,
        parent_return_idx=parent_return_idx,
        parent_enter_rail=parent_enter_rail,
        child_enter_rail=child_enter_rail,
    )


def linear_chain_from(nodes: Sequence[TreeNode], node_id: int) -> list[int]:
    chain = [node_id]
    current = node_id
    while len(nodes[current].children) == 1:
        current = nodes[current].children[0]
        chain.append(current)
    return chain


def is_simple_chain_subtree(nodes: Sequence[TreeNode], node_id: int) -> bool:
    current = node_id
    while True:
        child_count = len(nodes[current].children)
        if child_count == 0:
            return True
        if child_count > 1:
            return False
        current = nodes[current].children[0]


def route_simple_chain(
    nodes: Sequence[TreeNode],
    node_id: int,
    start_anchor: Point,
    exit_anchor: Point,
    spacing: float,
) -> tuple[list[Point], int]:
    """Route a single-child contour chain with one slot corridor.

    This replaces repeated per-level tree connectors.  The chain router opens
    every nested loop along the same two-slot corridor, so the return rail is a
    single line and the entry rail remains one spacing away from it.
    """

    chain = linear_chain_from(nodes, node_id)
    loops = [nodes[idx].loop for idx in chain]
    path, _ = cf.build_single_minimum_connected_fermat(
        loops,
        start_anchor=start_anchor,
        spacing=spacing,
        port_spacing=spacing,
        exit_anchor=exit_anchor,
        preserve_medial_pockets=True,
    )
    if len(path) >= 2 and cf.dist2(path[0], exit_anchor) + cf.dist2(path[-1], start_anchor) < cf.dist2(path[0], start_anchor) + cf.dist2(path[-1], exit_anchor):
        path.reverse()
    end_idx = loop_vertex_index(nodes[node_id].loop.points, path[-1] if path else exit_anchor)
    return path, end_idx


def visible_child_port_index(parent_point: Point, child_points: Sequence[Point], existing_path: Sequence[Point], spacing: float) -> int:
    if not child_points:
        return 0

    segments = list(zip(existing_path, existing_path[1:]))
    best_idx = loop_vertex_index(child_points, parent_point)
    best_score: tuple[int, float, float] | None = None
    stride = max(1, len(child_points) // 160)
    candidate_indices = set(range(0, len(child_points), stride))
    candidate_indices.add(best_idx)

    for idx in candidate_indices:
        candidate = child_points[idx]
        crossings = 0
        min_distance = float("inf")
        for seg_idx, (a, b) in enumerate(segments[:-4]):
            distance = cf.segment_distance(parent_point, candidate, a, b)
            min_distance = min(min_distance, distance)
            if distance <= 1e-7:
                crossings += 1
                if crossings > 4:
                    break
        score = (crossings, cf.dist(parent_point, candidate), -min_distance)
        if best_score is None or score < best_score:
            best_score = score
            best_idx = idx
            if crossings == 0 and cf.dist(parent_point, candidate) <= spacing * 2.5:
                break

    return best_idx


def route_cycle_with_children(
    nodes: Sequence[TreeNode],
    node_id: int,
    start_idx: int,
    direction: int,
    spacing: float,
    debug_edges: list[RouteDebugEdge],
    excluded_child: int | None = None,
    return_idx: int | None = None,
) -> tuple[list[Point], int]:
    """Walk one complete parent loop and splice every non-excluded child.

    The previous router walked only from the parent entry to a preselected exit.
    That made children outside that arc unreachable.  This routine makes the
    invariant explicit: if a node is not the final endpoint chain, all of its
    reachable children are fully routed before the node is allowed to return.
    """

    node = nodes[node_id]
    points = node.loop.points
    if len(points) < 3:
        return list(points), start_idx

    start_idx %= len(points)
    current = start_idx
    gap = port_gap_points(points, spacing)
    path: list[Point] = [points[start_idx]]

    child_order = [
        child_id
        for child_id in children_in_route_order(nodes, node_id, start_idx, direction)
        if child_id != excluded_child
    ]

    for child_id in child_order:
        child = nodes[child_id]
        connector = connector_pair_for_child(points, child.loop.points, child.loop.centroid, direction, spacing)
        before = connector.parent_enter_idx
        after = connector.parent_return_idx

        append_points(path, arc_points(points, current, before, direction))

        if is_simple_chain_subtree(nodes, child_id):
            child_path, child_end = route_simple_chain(
                nodes,
                child_id,
                start_anchor=connector.child_enter_rail,
                exit_anchor=child.loop.points[connector.child_return_idx],
                spacing=spacing,
            )
        else:
            child_path, child_end = route_cycle_with_children(
                nodes,
                child_id,
                connector.child_enter_idx,
                -direction,
                spacing,
                debug_edges,
                return_idx=connector.child_return_idx,
            )
        debug_edges.append(RouteDebugEdge("enter", node_id, child_id, points[before], connector.parent_enter_rail))
        debug_edges.append(RouteDebugEdge("enter", node_id, child_id, connector.parent_enter_rail, connector.child_enter_rail))
        if child_path:
            debug_edges.append(RouteDebugEdge("enter", node_id, child_id, connector.child_enter_rail, child_path[0]))
        append_point(path, connector.parent_enter_rail)
        append_point(path, connector.child_enter_rail)
        if child_path:
            append_point(path, child_path[0])
        append_points(path, child_path)
        if child_path:
            debug_edges.append(RouteDebugEdge("exit", node_id, child_id, child_path[-1], points[after]))
        append_point(path, points[after])
        current = after

    final_return_idx = return_idx if return_idx is not None else advance_index(points, start_idx, -direction, max(gap, 2))
    append_points(path, arc_points(points, current, final_return_idx, direction))
    return path, final_return_idx


def route_returning(
    nodes: Sequence[TreeNode],
    node_id: int,
    start_idx: int,
    exit_idx: int,
    direction: int,
    spacing: float,
    debug_edges: list[RouteDebugEdge],
) -> list[Point]:
    """Route a complete subtree from one boundary port to another.

    This is used for side branches.  Every descendant is entered and exited
    before the parent traversal continues, so these branches cannot become the
    accidental endpoint.
    """
    path, _ = route_cycle_with_children(nodes, node_id, start_idx, direction, spacing, debug_edges)
    return path


def choose_final_child(nodes: Sequence[TreeNode], node_id: int, exit_point: Point) -> int | None:
    children = nodes[node_id].children
    if not children:
        return None
    return min(children, key=lambda child_id: cf.dist2(nodes[child_id].loop.centroid, exit_point))


def route_open(
    nodes: Sequence[TreeNode],
    node_id: int,
    start_idx: int,
    direction: int,
    spacing: float,
    debug_edges: list[RouteDebugEdge],
) -> list[Point]:
    """Route a subtree and intentionally finish in one final leaf."""

    node = nodes[node_id]
    points = node.loop.points
    if len(points) < 3:
        return list(points)

    exit_point = points[start_idx]
    final_child = choose_final_child(nodes, node_id, exit_point)

    if final_child is None:
        path, _ = route_cycle_with_children(nodes, node_id, start_idx, direction, spacing, debug_edges)
        return path

    path, exit_idx = route_cycle_with_children(
        nodes,
        node_id,
        start_idx,
        direction,
        spacing,
        debug_edges,
        excluded_child=final_child,
    )
    exit_point = points[exit_idx]

    child = nodes[final_child]
    child_start = visible_child_port_index(exit_point, child.loop.points, path, spacing)
    debug_edges.append(RouteDebugEdge("final", node_id, final_child, exit_point, child.loop.points[child_start]))
    append_point(path, child.loop.points[child_start])
    append_points(path, route_open(nodes, final_child, child_start, -direction, spacing, debug_edges))
    return path


def sample_contour_miss_fraction(loop: cf.ContourLoop, path: Sequence[Point], spacing: float) -> float:
    if not loop.points:
        return 1.0
    stride = max(1, len(loop.points) // 48)
    missed = 0
    total = 0
    for idx in range(0, len(loop.points), stride):
        total += 1
        _, _, distance = cf.project_open_polyline(path, loop.points[idx])
        if distance > spacing * 0.72:
            missed += 1
    return missed / max(1, total)


def loop_distance(loop: cf.ContourLoop, p: Point) -> float:
    if loop.closed:
        return cf.distance_to_loop(p, loop.points)
    return cf.open_polyline_distance(p, loop.points)


def semantic_segment_sources(path: Sequence[Point], nodes: Sequence[TreeNode], spacing: float) -> list[int | None]:
    loop_bounds: list[tuple[int, cf.ContourLoop, tuple[float, float, float, float]]] = []
    margin = spacing * 0.55
    for node in nodes:
        if not node.loop.points:
            continue
        min_x, min_y, max_x, max_y = cf.loop_bounds(node.loop.points, margin=margin)
        loop_bounds.append((node.id, node.loop, (min_x, min_y, max_x, max_y)))

    sources: list[int | None] = []
    max_distance = spacing * 0.38
    for a, b in zip(path, path[1:]):
        mid = cf.lerp(a, b, 0.5)
        best_id: int | None = None
        best_distance = max_distance
        for node_id, loop, bounds in loop_bounds:
            min_x, min_y, max_x, max_y = bounds
            if not (min_x <= mid[0] <= max_x and min_y <= mid[1] <= max_y):
                continue
            distance = loop_distance(loop, mid)
            if distance < best_distance:
                best_distance = distance
                best_id = node_id
        sources.append(best_id)
    return sources


def semantic_path_pair_metrics(
    path: Sequence[Point],
    nodes: Sequence[TreeNode],
    spacing: float,
    spacing_tolerance: float,
    local_skip: int = 6,
) -> tuple[int, int, float, int]:
    if len(path) < 4:
        return 0, 0, float("inf"), 0

    sources = semantic_segment_sources(path, nodes, spacing)
    min_allowed = spacing * (1.0 - spacing_tolerance)
    bin_size = max(spacing * 1.5, EPS)
    bins: dict[tuple[int, int], list[int]] = {}
    segments = list(zip(path, path[1:]))
    segment_lengths = [cf.dist(a, b) for a, b in segments]
    prefix_lengths = [0.0]
    for length in segment_lengths:
        prefix_lengths.append(prefix_lengths[-1] + length)
    path_length = prefix_lengths[-1]
    local_skip_distance = spacing * 4.0
    closed_path = cf.dist(path[0], path[-1]) <= spacing * 0.20

    def key(x: float, y: float) -> tuple[int, int]:
        return (math.floor(x / bin_size), math.floor(y / bin_size))

    for idx, (a, b) in enumerate(segments):
        min_x = min(a[0], b[0]) - min_allowed
        min_y = min(a[1], b[1]) - min_allowed
        max_x = max(a[0], b[0]) + min_allowed
        max_y = max(a[1], b[1]) + min_allowed
        k0 = key(min_x, min_y)
        k1 = key(max_x, max_y)
        for by in range(k0[1], k1[1] + 1):
            for bx in range(k0[0], k1[0] + 1):
                bins.setdefault((bx, by), []).append(idx)

    checked: set[tuple[int, int]] = set()
    self_intersections = 0
    spacing_violations = 0
    semantic_ignored = 0
    min_nonlocal_spacing = float("inf")

    for i, (a, b) in enumerate(segments):
        min_x = min(a[0], b[0]) - min_allowed
        min_y = min(a[1], b[1]) - min_allowed
        max_x = max(a[0], b[0]) + min_allowed
        max_y = max(a[1], b[1]) + min_allowed
        k0 = key(min_x, min_y)
        k1 = key(max_x, max_y)
        for by in range(k0[1], k1[1] + 1):
            for bx in range(k0[0], k1[0] + 1):
                for j in bins.get((bx, by), []):
                    if j <= i:
                        continue
                    pair = (i, j)
                    if pair in checked:
                        continue
                    checked.add(pair)

                    index_distance = abs(i - j)
                    if closed_path:
                        index_distance = min(index_distance, len(segments) - index_distance)
                    if index_distance <= local_skip:
                        continue
                    path_gap = max(0.0, prefix_lengths[j] - prefix_lengths[i + 1])
                    if closed_path:
                        occupied_span = prefix_lengths[j + 1] - prefix_lengths[i]
                        path_gap = min(path_gap, max(0.0, path_length - occupied_span))
                    if path_gap <= local_skip_distance:
                        continue

                    same_source = sources[i] is not None and sources[i] == sources[j]
                    if same_source:
                        semantic_ignored += 1
                        continue

                    c, d = segments[j]
                    pair_distance = cf.segment_distance(a, b, c, d)
                    min_nonlocal_spacing = min(min_nonlocal_spacing, pair_distance)

                    if pair_distance <= 1e-7:
                        self_intersections += 1
                    elif pair_distance < min_allowed:
                        spacing_violations += 1

    return self_intersections, spacing_violations, min_nonlocal_spacing, semantic_ignored


def containment_violations(model: cf.PolygonModel, path: Sequence[Point], spacing: float) -> int:
    violations = 0
    for a, b in zip(path, path[1:]):
        length = cf.dist(a, b)
        samples = max(2, int(math.ceil(length / max(spacing * 0.25, EPS))))
        for i in range(samples + 1):
            if model.sdf(cf.lerp(a, b, i / samples)) < -spacing * 0.05:
                violations += 1
    return violations


def strict_coverage_audit(
    model: cf.PolygonModel,
    path: Sequence[Point],
    bounds: tuple[float, float, float, float],
    line_width: float,
    cells: int,
    segment_widths: Sequence[float] | None = None,
    sample_limit: int = 2400,
) -> dict[str, Any]:
    """Measure true underfill and non-local overlap separately.

    The older raster metric only asked whether a cell was touched by any bead.
    This audit uses exact cell-center to segment distance, no extra dilation,
    and counts multiple non-local bead visits as overlap instead of letting
    excess material hide an uncovered pocket elsewhere.
    """

    min_x, min_y, max_x, max_y = bounds
    width = max(max_x - min_x, EPS)
    height = max(max_y - min_y, EPS)
    if width >= height:
        nx = max(12, cells)
        ny = max(12, int(round(cells * height / width)))
    else:
        ny = max(12, cells)
        nx = max(12, int(round(cells * width / height)))

    segments = list(zip(path, path[1:]))
    widths = (
        list(segment_widths)
        if segment_widths is not None and len(segment_widths) == len(segments)
        else [line_width] * len(segments)
    )
    max_line_width = max([line_width, *widths]) if widths else line_width
    radius = max_line_width * 0.5
    bin_size = max(max_line_width, max(width / nx, height / ny) * 2.0, EPS)
    bins: dict[tuple[int, int], list[int]] = {}
    prefix_lengths = [0.0]
    for a, b in segments:
        prefix_lengths.append(prefix_lengths[-1] + cf.dist(a, b))

    def key(p: Point) -> tuple[int, int]:
        return (math.floor(p[0] / bin_size), math.floor(p[1] / bin_size))

    for idx, (a, b) in enumerate(segments):
        segment_radius = widths[idx] * 0.5
        k0 = key((min(a[0], b[0]) - segment_radius, min(a[1], b[1]) - segment_radius))
        k1 = key((max(a[0], b[0]) + segment_radius, max(a[1], b[1]) + segment_radius))
        for by in range(k0[1], k1[1] + 1):
            for bx in range(k0[0], k1[0] + 1):
                bins.setdefault((bx, by), []).append(idx)

    def nonlocal_visit_count(indices: Sequence[int]) -> int:
        if not indices:
            return 0
        visits = 1
        last = indices[0]
        for idx in indices[1:]:
            if idx - last <= 6:
                last = idx
                continue
            path_gap = max(0.0, prefix_lengths[idx] - prefix_lengths[last + 1])
            if path_gap > line_width * 3.0:
                visits += 1
            last = idx
        return visits

    inside_count = 0
    covered_count = 0
    underfill_count = 0
    overlap_count = 0
    max_visit_count = 0
    worst_underfill = 0.0
    underfill_flags = bytearray(nx * ny)
    overlap_flags = bytearray(nx * ny)
    underfill_samples: list[list[float]] = []
    overlap_samples: list[list[float]] = []

    def append_sample(samples: list[list[float]], p: Point) -> None:
        if len(samples) < sample_limit:
            samples.append([p[0], p[1]])

    for y in range(ny):
        py = min_y + (y + 0.5) * height / ny
        for x in range(nx):
            p = (min_x + (x + 0.5) * width / nx, py)
            cell_idx = y * nx + x
            if not model.contains(p):
                continue
            inside_count += 1
            candidates = bins.get(key(p), [])
            covering: list[int] = []
            nearest = float("inf")
            nearest_gap = float("inf")
            for idx in candidates:
                a, b = segments[idx]
                distance = cf.point_segment_distance(p, a, b)
                nearest = min(nearest, distance)
                nearest_gap = min(nearest_gap, distance - widths[idx] * 0.5)
                if distance <= widths[idx] * 0.5:
                    covering.append(idx)
            if not covering:
                underfill_count += 1
                underfill_flags[cell_idx] = 1
                if nearest_gap < float("inf"):
                    worst_underfill = max(worst_underfill, nearest_gap)
                append_sample(underfill_samples, p)
                continue
            covered_count += 1
            visits = nonlocal_visit_count(sorted(set(covering)))
            max_visit_count = max(max_visit_count, visits)
            if visits > 1:
                overlap_count += 1
                overlap_flags[cell_idx] = 1
                append_sample(overlap_samples, p)

    def component_stats(flags: bytearray) -> tuple[int, int, list[dict[str, Any]]]:
        seen = bytearray(nx * ny)
        components = 0
        largest = 0
        groups: list[dict[str, Any]] = []
        for idx, flag in enumerate(flags):
            if not flag or seen[idx]:
                continue
            components += 1
            size = 0
            points: list[Point] = []
            stack = [idx]
            seen[idx] = 1
            while stack:
                current = stack.pop()
                size += 1
                cx = current % nx
                cy = current // nx
                points.append((min_x + (cx + 0.5) * width / nx, min_y + (cy + 0.5) * height / ny))
                for nx2, ny2 in ((cx - 1, cy), (cx + 1, cy), (cx, cy - 1), (cx, cy + 1)):
                    if nx2 < 0 or nx2 >= nx or ny2 < 0 or ny2 >= ny:
                        continue
                    neighbor = ny2 * nx + nx2
                    if flags[neighbor] and not seen[neighbor]:
                            seen[neighbor] = 1
                            stack.append(neighbor)
            largest = max(largest, size)
            centroid = (
                sum(p[0] for p in points) / max(1, len(points)),
                sum(p[1] for p in points) / max(1, len(points)),
            )
            groups.append(
                {
                    "size": size,
                    "centroid": [centroid[0], centroid[1]],
                    "bounds": [
                        min(p[0] for p in points),
                        min(p[1] for p in points),
                        max(p[0] for p in points),
                        max(p[1] for p in points),
                    ],
                    "points": [[p[0], p[1]] for p in points[:128]],
                }
            )
        groups.sort(key=lambda item: item["size"], reverse=True)
        return components, largest, groups[:32]

    underfill_components, largest_underfill, underfill_groups = component_stats(underfill_flags)
    overlap_components, largest_overlap, overlap_groups = component_stats(overlap_flags)
    coverage_ratio = covered_count / inside_count if inside_count else 0.0
    underfill_ratio = underfill_count / inside_count if inside_count else 1.0
    overlap_ratio = overlap_count / inside_count if inside_count else 0.0
    return {
        "nx": nx,
        "ny": ny,
        "bounds": [min_x, min_y, max_x, max_y],
        "cellSize": max(width / nx, height / ny),
        "coverageRatio": coverage_ratio,
        "underfillRatio": underfill_ratio,
        "internalOverlapRatio": overlap_ratio,
        "insideCells": inside_count,
        "coveredCells": covered_count,
        "underfillCells": underfill_count,
        "overlapCells": overlap_count,
        "underfillComponents": underfill_components,
        "largestUnderfillComponent": largest_underfill,
        "overlapComponents": overlap_components,
        "largestOverlapComponent": largest_overlap,
        "underfillGroups": underfill_groups,
        "overlapGroups": overlap_groups,
        "maxVisitCount": max_visit_count,
        "worstUnderfillDistance": worst_underfill,
        "underfillSamples": underfill_samples,
        "overlapSamples": overlap_samples,
    }


def retain_small_printable_contours(contours: Sequence[cf.ContourLoop], spacing: float) -> list[cf.ContourLoop]:
    min_area = spacing * spacing * 0.75
    min_length = spacing * 2.5
    return [loop for loop in contours if loop.area >= min_area and loop.length >= min_length]


def segment_projection(a: Point, b: Point, p: Point) -> tuple[float, Point, float]:
    ab = v_sub(b, a)
    denom = ab[0] * ab[0] + ab[1] * ab[1]
    if denom <= EPS:
        return 0.0, a, cf.dist(a, p)
    t = max(0.0, min(1.0, ((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / denom))
    q = cf.lerp(a, b, t)
    return t, q, cf.dist(q, p)


def segment_inside_model(model: cf.PolygonModel, a: Point, b: Point, spacing: float) -> bool:
    length = cf.dist(a, b)
    samples = max(2, int(math.ceil(length / max(spacing * 0.35, EPS))))
    for idx in range(samples + 1):
        if model.sdf(cf.lerp(a, b, idx / samples)) < -spacing * 0.04:
            return False
    return True


def fair_underfill_detours(
    model: cf.PolygonModel,
    path: Sequence[Point],
    line_width: float,
    spacing: float,
    spacing_tolerance: float,
    coverage_cells: int,
    coverage_threshold: float,
    max_detours: int,
) -> tuple[list[Point], int, float]:
    """Bend existing segments through strict-audit underfill cells.

    This is a local repair: a segment A-B becomes A-P-B where P is an uncovered
    cell near that segment. It fills seam-sized gaps without creating separate
    islands or retraced in/out connectors.
    """

    current = list(path)
    inserted = 0
    if len(current) < 2 or max_detours <= 0:
        audit = strict_coverage_audit(model, current, model.bounds(margin=line_width), line_width, coverage_cells)
        return current, inserted, float(audit["coverageRatio"])

    passes = 3
    for _pass in range(passes):
        audit = strict_coverage_audit(model, current, model.bounds(margin=line_width), line_width, coverage_cells)
        if float(audit["coverageRatio"]) >= coverage_threshold:
            return current, inserted, float(audit["coverageRatio"])

        candidates: dict[int, tuple[float, Point]] = {}
        segments = list(zip(current, current[1:]))
        min_allowed = spacing * (1.0 - spacing_tolerance)

        def has_clearance(idx: int, p: Point) -> bool:
            a, b = segments[idx]
            for j, (c, d) in enumerate(segments):
                if abs(j - idx) <= 6:
                    continue
                if cf.segment_distance(a, p, c, d) < min_allowed:
                    return False
                if cf.segment_distance(p, b, c, d) < min_allowed:
                    return False
            return True

        for raw in audit.get("underfillSamples", []):
            p = (float(raw[0]), float(raw[1]))
            best: tuple[float, int, Point] | None = None
            for idx, (a, b) in enumerate(segments):
                length = cf.dist(a, b)
                if length < line_width * 0.55:
                    continue
                t, _projected, distance = segment_projection(a, b, p)
                if t < 0.12 or t > 0.88:
                    continue
                if distance > line_width * 1.10:
                    continue
                extra = cf.dist(a, p) + cf.dist(p, b) - length
                if extra > spacing * 1.8:
                    continue
                if not segment_inside_model(model, a, p, spacing) or not segment_inside_model(model, p, b, spacing):
                    continue
                if not has_clearance(idx, p):
                    continue
                score = distance - extra * 0.15
                if best is None or score > best[0]:
                    best = (score, idx, p)
            if best is None:
                continue
            score, idx, p = best
            existing = candidates.get(idx)
            if existing is None or score > existing[0]:
                candidates[idx] = (score, p)

        if not candidates:
            return current, inserted, float(audit["coverageRatio"])

        remaining = max_detours - inserted
        batch = sorted(candidates.items(), key=lambda item: item[1][0], reverse=True)[:remaining]
        if not batch:
            return current, inserted, float(audit["coverageRatio"])

        for idx, (_score, p) in sorted(batch, key=lambda item: item[0], reverse=True):
            if idx + 1 >= len(current):
                continue
            if cf.dist(current[idx], p) <= EPS or cf.dist(current[idx + 1], p) <= EPS:
                continue
            current.insert(idx + 1, p)
            inserted += 1
            if inserted >= max_detours:
                break

        if inserted >= max_detours:
            break

    audit = strict_coverage_audit(model, current, model.bounds(margin=line_width), line_width, coverage_cells)
    return current, inserted, float(audit["coverageRatio"])


def fair_segment_widths(
    model: cf.PolygonModel,
    path: Sequence[Point],
    line_width: float,
    coverage_cells: int,
    coverage_threshold: float,
    max_width_factor: float = 1.65,
) -> tuple[list[float], int, float]:
    if len(path) < 2:
        return [], 0, 0.0

    segments = list(zip(path, path[1:]))
    widths = [line_width] * len(segments)
    max_width = line_width * max_width_factor
    changed_segments: set[int] = set()

    for _pass in range(4):
        audit = strict_coverage_audit(
            model,
            path,
            model.bounds(margin=max_width),
            line_width,
            coverage_cells,
            segment_widths=widths,
        )
        if float(audit["coverageRatio"]) >= coverage_threshold:
            return widths, len(changed_segments), float(audit["coverageRatio"])

        proposed: dict[int, float] = {}
        cell_pad = float(audit.get("cellSize", 0.0)) * 0.75
        for raw in audit.get("underfillSamples", []):
            p = (float(raw[0]), float(raw[1]))
            best: tuple[float, int, float] | None = None
            for idx, (a, b) in enumerate(segments):
                t, _projection, distance = segment_projection(a, b, p)
                if t < 0.05 or t > 0.95:
                    continue
                required_width = min(max_width, 2.0 * (distance + cell_pad))
                if required_width <= widths[idx] + line_width * 0.02:
                    continue
                if required_width > max_width + EPS:
                    continue
                score = required_width - widths[idx]
                if best is None or score < best[0]:
                    best = (score, idx, required_width)
            if best is None:
                continue
            _score, idx, required_width = best
            proposed[idx] = max(proposed.get(idx, widths[idx]), required_width)

        if not proposed:
            return widths, len(changed_segments), float(audit["coverageRatio"])

        for idx, required_width in proposed.items():
            widths[idx] = max(widths[idx], required_width)
            changed_segments.add(idx)

    audit = strict_coverage_audit(
        model,
        path,
        model.bounds(margin=max_width),
        line_width,
        coverage_cells,
        segment_widths=widths,
    )
    return widths, len(changed_segments), float(audit["coverageRatio"])


def max_segment_length(path: Sequence[Point]) -> float:
    return max((cf.dist(a, b) for a, b in zip(path, path[1:])), default=0.0)


def model_payload(model: cf.PolygonModel) -> dict[str, Any]:
    return {
        "name": model.name,
        "outer": [list(p) for p in model.outer],
        "holes": [[[x, y] for x, y in hole] for hole in model.holes],
        "bounds": list(model.bounds(0.0)),
    }


def contour_payload(nodes: Sequence[TreeNode]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for node in nodes:
        loop = node.loop
        out.append(
            {
                "id": node.id,
                "level": loop.level_index,
                "area": loop.area,
                "centroid": list(loop.centroid),
                "parent": node.parent,
                "children": list(node.children),
                "points": [list(p) for p in loop.points],
            }
        )
    return out


def route_edge_payload(edges: Sequence[RouteDebugEdge]) -> list[dict[str, Any]]:
    return [
        {
            "kind": edge.kind,
            "parent": edge.parent,
            "child": edge.child,
            "a": list(edge.a),
            "b": list(edge.b),
        }
        for edge in edges
    ]


def plan_v2_tree(
    model: cf.PolygonModel,
    options: dict[str, Any],
    planning_model: cf.PolygonModel | BarrierModel | None = None,
    algorithm_name: str = "contour_tree_v2",
    initial_diagnostics: Sequence[str] = (),
    initial_path: Sequence[Point] = (),
) -> dict[str, Any]:
    started = time.perf_counter()
    planning_model = planning_model or model
    grid_cells = int(options.get("grid", 180))
    coverage_cells = int(options.get("coverageGrid", 150))
    line_width = float(options.get("lineWidth", 1.2))
    spacing = float(options.get("spacing", line_width))
    max_levels = int(options.get("maxLevels", 256))
    coverage_threshold = float(options.get("coverageThreshold", 0.90))
    spacing_tolerance = float(options.get("spacingTolerance", 0.25))
    spacing_warning_threshold = int(options.get("spacingWarningThreshold", 3))
    overlap_threshold = float(options.get("overlapThreshold", 0.02))
    residual_repair = bool(options.get("residualRepair", False))
    use_all_contours = bool(options.get("useAllContours", False))
    width_fairing = bool(options.get("widthFairing", False))
    underfill_fairing = bool(options.get("underfillFairing", False))
    max_underfill_detours = int(options.get("maxUnderfillDetours", 96))
    semantic_spacing = bool(options.get("semanticSpacing", False))
    start_fraction = float(options.get("startFraction", 0.0))
    direction_name = str(options.get("direction", "inward"))
    route_winding = -1 if int(float(options.get("routeWinding", 1))) < 0 else 1
    direction = route_winding

    diagnostics: list[str] = list(initial_diagnostics)
    grid = cf.build_sdf_grid(planning_model, grid_cells=grid_cells, margin=line_width * 2.0)
    generated_contours = cf.generate_offset_contours(grid, line_width=line_width, spacing=spacing, max_levels=max_levels)
    standard_printable_contours = cf.filter_printable_contours(generated_contours, spacing=spacing)
    printable_contours = (
        retain_small_printable_contours(generated_contours, spacing=spacing)
        if use_all_contours
        else standard_printable_contours
    )
    if use_all_contours and len(printable_contours) > len(standard_printable_contours):
        diagnostics.append(
            f"v5 retained {len(printable_contours) - len(standard_printable_contours)} extra small iso-contour(s) after relaxed filtering."
        )
    if isinstance(planning_model, BarrierModel) or use_all_contours:
        contours = printable_contours
    else:
        contours = cf.filter_topology_stable_levels(printable_contours, spacing=spacing)
    # Keep the first v2 lab route to closed iso-contours only.  The legacy
    # prototype's terminal medial open-line repair needs a separate graph node
    # type before the contour-tree router can handle it without false crossings.

    nodes, tree_diagnostics = build_contour_tree(contours)
    diagnostics.extend(tree_diagnostics)
    root_id = root_node_id(nodes)

    debug_edges: list[RouteDebugEdge] = []
    path: list[Point] = []
    if root_id is None:
        diagnostics.append("No routable contour tree root was generated.")
    else:
        root = nodes[root_id]
        start_point = initial_path[-1] if initial_path else cf.point_at_closed_fraction(root.loop.points, start_fraction)
        start_idx = loop_vertex_index(root.loop.points, start_point)
        routed_path = route_open(nodes, root_id, start_idx, direction, spacing, debug_edges)
        if initial_path:
            path = list(initial_path)
            append_points(path, routed_path)
        else:
            path = routed_path

    if direction_name == "outward":
        path = list(reversed(path))

    residual_inserted = 0
    residual_candidates = 0
    if residual_repair and path:
        repaired_path, residual_inserted, residual_candidates = cf.insert_residual_gap_spirals(
            model,
            grid,
            path,
            line_width=line_width,
            spacing=spacing,
            max_levels=max_levels,
            spacing_tolerance=spacing_tolerance,
        )
        if residual_candidates:
            diagnostics.append(
                f"v5 generated {residual_candidates} residual gap contour(s) from uncovered material."
            )
        if residual_inserted:
            diagnostics.append(f"v5 inserted {residual_inserted} residual gap spiral(s).")
            path = repaired_path

    segment_widths: list[float] = [line_width] * max(0, len(path) - 1)
    width_fairing_segments = 0
    width_fairing_coverage = 0.0
    if width_fairing and path:
        segment_widths, width_fairing_segments, width_fairing_coverage = fair_segment_widths(
            model,
            path,
            line_width=line_width,
            coverage_cells=coverage_cells,
            coverage_threshold=coverage_threshold,
        )
        if width_fairing_segments:
            diagnostics.append(
                f"v6 widened {width_fairing_segments} segment(s); interim strict coverage {width_fairing_coverage:.3f}."
            )

    fairing_detours = 0
    fairing_coverage = 0.0
    if underfill_fairing and path:
        path, fairing_detours, fairing_coverage = fair_underfill_detours(
            model,
            path,
            line_width=line_width,
            spacing=spacing,
            spacing_tolerance=spacing_tolerance,
            coverage_cells=coverage_cells,
            coverage_threshold=coverage_threshold,
            max_detours=max_underfill_detours,
        )
        if fairing_detours:
            diagnostics.append(
                f"v6 inserted {fairing_detours} local underfill detour(s); interim strict coverage {fairing_coverage:.3f}."
            )
            segment_widths = [line_width] * max(0, len(path) - 1)

    elapsed = time.perf_counter() - started
    semantic_ignored_spacing = 0
    if semantic_spacing:
        crossings, spacing_violations, min_nonlocal_spacing, semantic_ignored_spacing = semantic_path_pair_metrics(
            path,
            nodes,
            spacing=spacing,
            spacing_tolerance=spacing_tolerance,
        )
    else:
        crossings, spacing_violations, min_nonlocal_spacing = cf.path_pair_metrics(
            path,
            spacing=spacing,
            spacing_tolerance=spacing_tolerance,
        )
    contain = containment_violations(model, path, spacing)
    legacy_coverage_ratio, legacy_underfill_ratio, overfill_ratio = cf.raster_coverage(
        model,
        path,
        model.bounds(margin=line_width),
        line_width,
        coverage_cells,
        segment_widths=segment_widths,
    )
    coverage_audit = strict_coverage_audit(
        model,
        path,
        model.bounds(margin=line_width),
        line_width,
        coverage_cells,
        segment_widths=segment_widths,
    )
    coverage_ratio = float(coverage_audit["coverageRatio"])
    underfill_ratio = float(coverage_audit["underfillRatio"])
    internal_overlap_ratio = float(coverage_audit["internalOverlapRatio"])
    missed = [
        {
            "id": node.id,
            "level": node.loop.level_index,
            "missFraction": sample_contour_miss_fraction(node.loop, path, spacing),
        }
        for node in nodes
    ]
    missed_bad = [item for item in missed if item["missFraction"] > 0.25]
    independent_roots = [node.id for node in nodes if node.parent is None]
    routed_root_ok = len(independent_roots) <= 1

    if missed_bad:
        diagnostics.append(f"{len(missed_bad)} contour loop(s) are not sufficiently covered by the routed path.")
    if contain:
        diagnostics.append(f"{contain} sampled connector/path point(s) leave the printable polygon.")
    if crossings:
        diagnostics.append(f"{crossings} non-local path segment pair(s) cross or touch.")
    if spacing_violations:
        diagnostics.append(f"{spacing_violations} non-local path segment pair(s) are below spacing tolerance.")
    if semantic_ignored_spacing:
        diagnostics.append(f"Semantic spacing ignored {semantic_ignored_spacing} same-contour close pair(s).")
    if spacing_violations > spacing_warning_threshold:
        diagnostics.append(
            f"Spacing warnings {spacing_violations} exceed requested threshold {spacing_warning_threshold}."
        )
    if coverage_ratio < coverage_threshold:
        diagnostics.append(f"Strict coverage {coverage_ratio:.3f} is below requested threshold {coverage_threshold:.3f}.")
    if coverage_audit["underfillCells"]:
        diagnostics.append(
            f"Strict audit found {coverage_audit['underfillCells']} underfilled cell(s) in "
            f"{coverage_audit['underfillComponents']} group(s); largest group {coverage_audit['largestUnderfillComponent']} cell(s)."
        )
    if internal_overlap_ratio > overlap_threshold:
        diagnostics.append(
            f"Internal overlap {internal_overlap_ratio:.3f} is above requested threshold {overlap_threshold:.3f}."
        )
    if coverage_audit["overlapCells"]:
        diagnostics.append(
            f"Strict audit found {coverage_audit['overlapCells']} internal overlap cell(s) in "
            f"{coverage_audit['overlapComponents']} group(s)."
        )
    if legacy_coverage_ratio - coverage_ratio > 0.01:
        diagnostics.append(
            f"Legacy raster overestimated coverage by {legacy_coverage_ratio - coverage_ratio:.3f}; strict audit found uncovered cells."
        )

    ok = (
        bool(path)
        and routed_root_ok
        and contain == 0
        and crossings == 0
        and spacing_violations <= spacing_warning_threshold
        and not missed_bad
        and coverage_ratio >= coverage_threshold
        and internal_overlap_ratio <= overlap_threshold
    )

    return {
        "algorithm": algorithm_name,
        "ok": ok,
        "model": model_payload(model),
        "grid": {
            "nx": grid.nx,
            "ny": grid.ny,
            "dx": grid.dx,
            "dy": grid.dy,
            "maxSdf": grid.max_sdf(),
        },
        "contours": contour_payload(nodes),
        "routeEdges": route_edge_payload(debug_edges),
        "path": [list(p) for p in path],
        "metrics": {
            "ok": ok,
            "elapsedSeconds": elapsed,
            "pathPoints": len(path),
            "pathLength": cf.polyline_length(path, closed=False),
            "maxSegment": max_segment_length(path),
            "contourCount": len(contours),
            "treeRoots": len(independent_roots),
            "treeDepth": max((node.loop.level_index for node in nodes), default=0),
            "containmentViolations": contain,
            "selfIntersections": crossings,
            "spacingViolations": spacing_violations,
            "spacingWarningThreshold": spacing_warning_threshold,
            "semanticIgnoredSpacingPairs": semantic_ignored_spacing,
            "minNonlocalSpacing": min_nonlocal_spacing if min_nonlocal_spacing < float("inf") else 0.0,
            "coverageRatio": coverage_ratio,
            "underfillRatio": underfill_ratio,
            "overfillRatio": overfill_ratio,
            "internalOverlapRatio": internal_overlap_ratio,
            "strictInsideCells": int(coverage_audit["insideCells"]),
            "strictUnderfillCells": int(coverage_audit["underfillCells"]),
            "strictOverlapCells": int(coverage_audit["overlapCells"]),
            "underfillComponents": int(coverage_audit["underfillComponents"]),
            "largestUnderfillComponent": int(coverage_audit["largestUnderfillComponent"]),
            "overlapComponents": int(coverage_audit["overlapComponents"]),
            "largestOverlapComponent": int(coverage_audit["largestOverlapComponent"]),
            "legacyCoverageRatio": legacy_coverage_ratio,
            "legacyUnderfillRatio": legacy_underfill_ratio,
            "residualGapContours": residual_candidates,
            "residualGapSpirals": residual_inserted,
            "widthFairingSegments": width_fairing_segments,
            "maxSegmentWidth": max(segment_widths, default=line_width),
            "underfillDetours": fairing_detours,
            "missedContourCount": len(missed_bad),
            "start": list(path[0]) if path else None,
            "end": list(path[-1]) if path else None,
            "direction": direction_name,
            "routeWinding": route_winding,
        },
        "segmentWidths": segment_widths,
        "coverageAudit": coverage_audit,
        "missedContours": missed,
        "diagnostics": diagnostics,
    }


def plan_legacy_cfs(model: cf.PolygonModel, options: dict[str, Any]) -> dict[str, Any]:
    line_width = float(options.get("lineWidth", 1.2))
    spacing = float(options.get("spacing", line_width))
    result = cf.plan_one_layer(
        model,
        grid_cells=int(options.get("grid", 180)),
        line_width=line_width,
        spacing=spacing,
        max_levels=int(options.get("maxLevels", 256)),
        coverage_cells=int(options.get("coverageGrid", 150)),
        coverage_threshold=float(options.get("coverageThreshold", 0.82)),
        spacing_tolerance=float(options.get("spacingTolerance", 0.25)),
        start_fraction=float(options.get("startFraction", 0.0)),
        exit_fraction=float(options.get("exitFraction", 0.5)),
    )
    nodes = [TreeNode(id=i, loop=loop) for i, loop in enumerate(result.contours)]
    payload = {
        "algorithm": "legacy_cfs",
        "ok": result.metrics.ok,
        "model": model_payload(model),
        "grid": {
            "nx": result.grid.nx,
            "ny": result.grid.ny,
            "dx": result.grid.dx,
            "dy": result.grid.dy,
            "maxSdf": result.grid.max_sdf(),
        },
        "contours": contour_payload(nodes),
        "routeEdges": [],
        "path": [list(p) for p in result.path],
        "metrics": {
            "ok": result.metrics.ok,
            "elapsedSeconds": result.metrics.elapsed_seconds,
            "pathPoints": result.metrics.path_points,
            "pathLength": result.metrics.path_length,
            "maxSegment": result.metrics.max_segment,
            "contourCount": result.metrics.contour_count,
            "treeRoots": 0,
            "treeDepth": result.metrics.contour_levels,
            "containmentViolations": result.metrics.containment_violations,
            "selfIntersections": result.metrics.self_intersections,
            "spacingViolations": result.metrics.spacing_violations,
            "minNonlocalSpacing": result.metrics.min_nonlocal_spacing,
            "coverageRatio": result.metrics.coverage_ratio,
            "underfillRatio": result.metrics.underfill_ratio,
            "overfillRatio": result.metrics.overfill_ratio,
            "missedContourCount": 0,
            "start": list(result.path[0]) if result.path else None,
            "end": list(result.path[-1]) if result.path else None,
            "direction": "cycle",
        },
        "missedContours": [],
        "diagnostics": list(result.diagnostics),
    }
    return payload


def plan_v3_tree(model: cf.PolygonModel, options: dict[str, Any]) -> dict[str, Any]:
    planning_model, diagnostics = barrier_model_for_holes(model)
    return plan_v2_tree(
        model,
        options,
        planning_model=planning_model,
        algorithm_name="contour_tree_v3",
        initial_diagnostics=diagnostics,
        initial_path=getattr(planning_model, "entry_path", []),
    )


def retry_score(result: dict[str, Any], coverage_threshold: float) -> tuple[float, ...]:
    metrics = result.get("metrics", {})
    coverage = float(metrics.get("coverageRatio") or 0.0)
    return (
        0.0 if result.get("ok") else 1.0,
        float(metrics.get("containmentViolations") or 0),
        float(metrics.get("selfIntersections") or 0),
        float(metrics.get("missedContourCount") or 0),
        max(0.0, coverage_threshold - coverage),
        float(metrics.get("internalOverlapRatio") or 0.0),
        float(metrics.get("spacingViolations") or 0),
        float(metrics.get("underfillRatio") or 0.0),
        float(metrics.get("overfillRatio") or 0.0),
        float(metrics.get("pathLength") or 0.0),
    )


def retry_option_variants(options: dict[str, Any], max_attempts: int) -> list[dict[str, Any]]:
    base_start = float(options.get("startFraction", 0.0)) % 1.0
    base_winding = -1 if int(float(options.get("routeWinding", 1))) < 0 else 1
    start_offsets = [0.0, 0.25, 0.50, 0.75, 0.125, 0.375, 0.625, 0.875]
    windings = [base_winding, -base_winding]
    variants: list[dict[str, Any]] = []
    seen: set[tuple[float, int]] = set()

    for offset in start_offsets:
        for winding in windings:
            start = (base_start + offset) % 1.0
            key = (round(start, 6), winding)
            if key in seen:
                continue
            seen.add(key)
            variant = dict(options)
            variant["startFraction"] = start
            variant["routeWinding"] = winding
            variants.append(variant)
            if len(variants) >= max_attempts:
                return variants

    return variants


def retry_attempt_summary(index: int, options: dict[str, Any], result: dict[str, Any]) -> dict[str, Any]:
    metrics = result.get("metrics", {})
    return {
        "index": index,
        "ok": bool(result.get("ok")),
        "startFraction": float(options.get("startFraction", 0.0)),
        "routeWinding": int(options.get("routeWinding", 1)),
        "coverageRatio": float(metrics.get("coverageRatio") or 0.0),
        "internalOverlapRatio": float(metrics.get("internalOverlapRatio") or 0.0),
        "containmentViolations": int(metrics.get("containmentViolations") or 0),
        "selfIntersections": int(metrics.get("selfIntersections") or 0),
        "spacingViolations": int(metrics.get("spacingViolations") or 0),
        "missedContourCount": int(metrics.get("missedContourCount") or 0),
        "pathPoints": int(metrics.get("pathPoints") or 0),
    }


def format_retry_summary(summary: dict[str, Any]) -> str:
    return (
        f"v4 attempt {summary['index']}: "
        f"start={summary['startFraction']:.3f}, "
        f"winding={summary['routeWinding']:+d}, "
        f"coverage={summary['coverageRatio']:.3f}, "
        f"overlap={summary['internalOverlapRatio']:.3f}, "
        f"cross={summary['selfIntersections']}, "
        f"missed={summary['missedContourCount']}, "
        f"contain={summary['containmentViolations']}, "
        f"{'pass' if summary['ok'] else 'fail'}."
    )


def plan_v4_retry(model: cf.PolygonModel, options: dict[str, Any]) -> dict[str, Any]:
    started = time.perf_counter()
    coverage_threshold = float(options.get("coverageThreshold", 0.90))
    max_attempts = max(1, min(16, int(options.get("retryAttempts", 4))))
    variants = retry_option_variants(options, max_attempts)

    best_result: dict[str, Any] | None = None
    best_score: tuple[float, ...] | None = None
    best_index = 0
    summaries: list[dict[str, Any]] = []

    for index, variant in enumerate(variants, start=1):
        result = plan_v3_tree(model, variant)
        score = retry_score(result, coverage_threshold)
        summaries.append(retry_attempt_summary(index, variant, result))
        if best_score is None or score < best_score:
            best_score = score
            best_result = result
            best_index = index
        if result.get("ok"):
            break

    if best_result is None:
        return plan_v3_tree(model, options)

    elapsed = time.perf_counter() - started
    selected = best_result
    selected["algorithm"] = "contour_tree_v4"
    selected["attempts"] = summaries
    selected["diagnostics"] = [
        f"v4 tried {len(summaries)} deterministic route variant(s); selected attempt {best_index}.",
        *[format_retry_summary(summary) for summary in summaries],
        *list(selected.get("diagnostics", [])),
    ]
    metrics = selected.get("metrics", {})
    metrics["elapsedSeconds"] = elapsed
    metrics["attemptCount"] = len(summaries)
    metrics["selectedAttempt"] = best_index
    selected["metrics"] = metrics
    return selected


def plan_v5_coverage(model: cf.PolygonModel, options: dict[str, Any]) -> dict[str, Any]:
    planning_model, diagnostics = barrier_model_for_holes(model)
    v5_options = dict(options)
    v5_options["useAllContours"] = True
    v5_options["residualRepair"] = True
    diagnostics = [
        "v5 uses all printable iso-contours, residual gap contours, and strict coverage audit.",
        *diagnostics,
    ]
    return plan_v2_tree(
        model,
        v5_options,
        planning_model=planning_model,
        algorithm_name="contour_tree_v5",
        initial_diagnostics=diagnostics,
        initial_path=getattr(planning_model, "entry_path", []),
    )


def plan_v6_coverage(model: cf.PolygonModel, options: dict[str, Any]) -> dict[str, Any]:
    line_width = float(options.get("lineWidth", 1.2))
    planning_model, diagnostics = barrier_model_for_holes(model, barrier_radius=line_width * 0.5)
    v6_options = dict(options)
    v6_options["useAllContours"] = True
    v6_options["residualRepair"] = True
    v6_options.setdefault("widthFairing", False)
    v6_options["semanticSpacing"] = True
    v6_options["underfillFairing"] = True
    v6_options.setdefault("maxUnderfillDetours", 180)
    diagnostics = [
        "v6 starts from v5, reserves printed barrier corridors, then locally repairs strict-audit underfill.",
        *diagnostics,
    ]
    return plan_v2_tree(
        model,
        v6_options,
        planning_model=planning_model,
        algorithm_name="contour_tree_v6",
        initial_diagnostics=diagnostics,
        initial_path=getattr(planning_model, "entry_path", []),
    )


def plan_model(model: cf.PolygonModel, options: dict[str, Any]) -> dict[str, Any]:
    algorithm = str(options.get("algorithm", "contour_tree_v2"))
    if algorithm == "legacy_cfs":
        return plan_legacy_cfs(model, options)
    if algorithm == "contour_tree_v6":
        return plan_v6_coverage(model, options)
    if algorithm == "contour_tree_v5":
        return plan_v5_coverage(model, options)
    if algorithm in {"contour_tree_v4", "contour_tree_v4_retry"}:
        return plan_v4_retry(model, options)
    if algorithm == "contour_tree_v3":
        return plan_v3_tree(model, options)
    return plan_v2_tree(model, options)


def available_shapes() -> dict[str, cf.PolygonModel]:
    return cf.available_shapes()
