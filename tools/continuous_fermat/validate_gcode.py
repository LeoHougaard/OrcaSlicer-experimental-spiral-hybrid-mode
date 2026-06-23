#!/usr/bin/env python3
"""Validate integrated Continuous Fermat G-code sections.

The checker is intentionally conservative. It scans the full file to maintain
XY/E state, then enforces that marked Continuous Fermat sections contain only a
connected sequence of extruding XY moves.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


BEGIN_MARKER = ";_CONTINUOUS_FERMAT_BEGIN"
END_MARKER = ";_CONTINUOUS_FERMAT_END"
AXIS_RE = re.compile(r"([A-Za-z])([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)")


@dataclass
class MachineState:
    absolute_xy: bool = True
    absolute_e: bool = True
    x: float | None = None
    y: float | None = None
    z: float | None = None
    e: float = 0.0


@dataclass
class SectionState:
    start_line: int
    previous_extrusion_end: tuple[float, float] | None = None
    extrusion_moves: int = 0


@dataclass
class ValidationResult:
    path: Path
    sections: int = 0
    extrusion_moves: int = 0
    violations: list[str] = field(default_factory=list)

    def fail(self, line_no: int, message: str) -> None:
        self.violations.append(f"line {line_no}: {message}")


def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].strip()


def parse_words(code: str) -> tuple[str, dict[str, float]]:
    if not code:
        return "", {}
    parts = code.split(None, 1)
    command = parts[0].upper()
    axes = {axis.upper(): float(value) for axis, value in AXIS_RE.findall(code)}
    return command, axes


def xy_distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def move_end(state: MachineState, axes: dict[str, float]) -> tuple[float | None, float | None, float | None]:
    x = state.x
    y = state.y
    z = state.z
    if "X" in axes:
        x = axes["X"] if state.absolute_xy or state.x is None else state.x + axes["X"]
    if "Y" in axes:
        y = axes["Y"] if state.absolute_xy or state.y is None else state.y + axes["Y"]
    if "Z" in axes:
        z = axes["Z"] if state.absolute_xy or state.z is None else state.z + axes["Z"]
    return x, y, z


def extrusion_delta(state: MachineState, axes: dict[str, float]) -> float:
    if "E" not in axes:
        return 0.0
    return axes["E"] - state.e if state.absolute_e else axes["E"]


def update_state_after_move(state: MachineState, axes: dict[str, float], end: tuple[float | None, float | None, float | None]) -> None:
    state.x, state.y, state.z = end
    if "E" in axes:
        state.e = axes["E"] if state.absolute_e else state.e + axes["E"]


def validate_lines(lines: Iterable[str], path: Path, xy_tolerance: float, e_tolerance: float) -> ValidationResult:
    state = MachineState()
    result = ValidationResult(path=path)
    section: SectionState | None = None

    for line_no, raw in enumerate(lines, start=1):
        if BEGIN_MARKER in raw:
            if section is not None:
                result.fail(line_no, "nested Continuous Fermat begin marker")
            section = SectionState(start_line=line_no)
            result.sections += 1
            continue

        if END_MARKER in raw:
            if section is None:
                result.fail(line_no, "Continuous Fermat end marker without matching begin")
            else:
                if section.extrusion_moves == 0:
                    result.fail(section.start_line, "Continuous Fermat section contains no extrusion moves")
                result.extrusion_moves += section.extrusion_moves
                section = None
            continue

        command, axes = parse_words(strip_comment(raw))
        if not command:
            continue

        if command == "G90":
            state.absolute_xy = True
            continue
        if command == "G91":
            state.absolute_xy = False
            continue
        if command == "M82":
            state.absolute_e = True
            continue
        if command == "M83":
            state.absolute_e = False
            continue

        if command == "G92":
            if "X" in axes:
                state.x = axes["X"]
            if "Y" in axes:
                state.y = axes["Y"]
            if "Z" in axes:
                state.z = axes["Z"]
            if "E" in axes:
                state.e = axes["E"]
            continue

        if section is not None and command in {"G10", "G11", "G22", "G23"}:
            result.fail(line_no, f"firmware retract/unretract command {command} inside Continuous Fermat section")

        if command not in {"G0", "G1", "G2", "G3"}:
            continue

        start_xy = (state.x, state.y)
        end_x, end_y, end_z = move_end(state, axes)
        de = extrusion_delta(state, axes)
        has_xy_command = "X" in axes or "Y" in axes
        has_known_xy = start_xy[0] is not None and start_xy[1] is not None and end_x is not None and end_y is not None
        has_xy_motion = has_xy_command and (
            not has_known_xy or xy_distance((start_xy[0], start_xy[1]), (end_x, end_y)) > xy_tolerance
        )

        if section is not None:
            if "E" in axes and abs(de) > e_tolerance and not has_xy_motion:
                result.fail(line_no, "E-only retract/unretract inside Continuous Fermat section")
            if has_xy_motion and de <= e_tolerance:
                result.fail(line_no, f"non-extruding XY move inside Continuous Fermat section ({command})")
            if has_xy_motion and de < -e_tolerance:
                result.fail(line_no, "retracting XY move inside Continuous Fermat section")
            if has_xy_motion and de > e_tolerance:
                if not has_known_xy:
                    result.fail(line_no, "extrusion move starts from unknown XY position")
                else:
                    start = (start_xy[0], start_xy[1])
                    end = (end_x, end_y)
                    if section.previous_extrusion_end is not None:
                        gap = xy_distance(section.previous_extrusion_end, start)
                        if gap > xy_tolerance:
                            result.fail(line_no, f"extrusion continuity gap {gap:.6f} mm before move")
                    section.previous_extrusion_end = end
                    section.extrusion_moves += 1

        update_state_after_move(state, axes, (end_x, end_y, end_z))

    if section is not None:
        result.fail(section.start_line, "unterminated Continuous Fermat section")
    if result.sections == 0:
        result.violations.append("no Continuous Fermat sections found")
    return result


def validate_file(path: Path, xy_tolerance: float, e_tolerance: float) -> ValidationResult:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        return validate_lines(fh, path, xy_tolerance, e_tolerance)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate ;_CONTINUOUS_FERMAT_* G-code sections.")
    parser.add_argument("gcode", nargs="+", type=Path, help="G-code file(s) to validate")
    parser.add_argument("--xy-tolerance", type=float, default=0.001, help="XY continuity tolerance in mm")
    parser.add_argument("--e-tolerance", type=float, default=1e-7, help="Extrusion delta tolerance")
    args = parser.parse_args()

    failed = False
    for path in args.gcode:
        result = validate_file(path, args.xy_tolerance, args.e_tolerance)
        if result.violations:
            failed = True
            print(f"{path}: FAIL ({result.sections} sections, {result.extrusion_moves} extrusion moves)")
            for violation in result.violations:
                print(f"  {violation}")
        else:
            print(f"{path}: OK ({result.sections} sections, {result.extrusion_moves} extrusion moves)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
