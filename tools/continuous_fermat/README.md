# Continuous Fermat Prototype and G-code Auditor

> **Pre-alpha source-level research only. Do not send prototype or Continuous Fermat output to a printer.** Standard application/CLI export, cached-file reuse, local copy, and upload are hard-blocked. A geometry or G-code audit pass is not a machine-safety certification.

This directory contains a legacy standalone geometry prototype and the independent serialized-G-code auditor. Production path planning is implemented separately in C++; prototype results are not production acceptance evidence.

The geometry prototype is intentionally dependency-free:

1. Built-in test shapes are represented as polygons with optional holes.
2. A signed-distance field is sampled over the shape.
3. Marching squares extracts inward offset contours.
4. The experimental generator builds a cyclic contour-weave spiral from inward offset contours.
5. The prototype path starts and ends on the outside contour.
6. A sampled verifier checks centerline containment, self-intersections, non-adjacent spacing, and rasterized coverage.
7. Scanline fallbacks are not accepted because they are not Fermat/spiral paths and do not preserve the intended outer-wall behavior.

The outer-contour endpoint rule is historical. Current C++ integration scope obtains exact layer handoff from deterministic identical rectangular layers; it does not implement a general inter-layer boundary-port solver.

Run all built-in prototype cases:

```powershell
python tools/continuous_fermat/fermat_layer.py --all --draw-contours
```

Run one prototype case:

```powershell
python tools/continuous_fermat/fermat_layer.py --shape annulus --draw-contours
```

## Independent Serialized-G-code Audit

The application does not expose Continuous Fermat output as printer-ready G-code. A non-serialized, per-`Print` capability used only by C++ regression tests can emit an exact development artifact; no profile, GUI, CLI, 3MF, or application environment setting enables it. Cached/imported files containing Continuous Fermat markers are also blocked even when the mode flag has been disabled.

This Python validator is not invoked automatically. Run a diagnostic audit without exact counts only while investigating a historical or deliberately generated test artifact:

```powershell
python tools/continuous_fermat/validate_gcode.py path\to\model.gcode --allow-unknown-counts
```

Any complete test artifact being evaluated against the strict contract must supply both expected model-layer and section counts:

```powershell
python tools/continuous_fermat/validate_gcode.py path\to\model.gcode `
    --expected-layers 120 --expected-sections 120
```

Before every marked section, the validator requires explicitly established `G21`, `G90`, `M83`, `T0`, `M200 D0`, `M220 S100`, `M221 S100`, `M900 K0`, and `M413 S0` state. Before the first section it also requires an exact all-axis `G28` followed by a positive-target `M109 R...` wait. It enforces a narrow Marlin-2-oriented, millimetre, linear-`G1` contract. It checks:

- standalone, balanced Continuous Fermat markers;
- exact expected layer and section accounting;
- finite parsed motion values, known positive feed, and known positive layer Z;
- positive E on every XY move inside a marked section;
- no Z motion, retract, travel, arc, reset, or disruptive mode change inside a section;
- closed section endpoints;
- no XY/E/retract/reset or disruptive mode change between sections; and
- a finite positive Z increase with the next section beginning at the preceding XY endpoint.

`--expected-layers` recognizes Orca's `;LAYER_CHANGE` and `; CHANGE_LAYER` comments and requires exactly one marked section per layer.

Positive-E XY motion outside a marked section is rejected. The test-only serializer rejects executable file/machine/filament start hooks and pressure advance, resets linear advance with `M900 K0`, homes with `G28`, and waits in both heating and cooling directions with `M109 R...`. `--allow-unmarked-before-first-section` is retained only for diagnosing historical or non-production files with a reviewed purge; it is not compatible with the guarded C++ integration. The exception never applies between or after sections and does not prove that the purge location, temperature, homing state, leveling state, or machine envelope is safe.

After the final marker, XY or E motion, firmware retract, `G92`, and modal changes are rejected. At most one positive Z-only clearance lift of no more than `2 mm` is accepted, followed only by the exact shutdown/progress allowlist. Final retraction and parking are not accepted merely because they occur after the model.

The parser rejects malformed or non-finite values for commands it parses, exact-word violations for normalizers and final shutdown, commandless modal axis lines, arcs, common coordinate transforms, tool/material changes, and unreviewed commands in protected spans. It does not semantically range-check every permitted acceleration, fan, temperature, or progress argument.

## Prototype Outputs

Prototype outputs are written to `build/continuous_fermat/`:

- `*.svg` shows the model boundary, optional raw contours, generated path, green start marker, and purple end marker.
- `*.png` contains a raster preview.
- `*.json` contains prototype metrics for regression comparisons.

## Important Limitations

- The SDF/marching-squares prototype deliberately differs from the production Clipper implementation; results do not transfer between them.
- Branch and multi-hole work is research-only. The current guarded C++ entry point rejects holes, branches, concave shapes, and all non-rectangular geometry.
- A path that looks spiral-like may still revisit a slot, overlap, leave the region, or miss material.
- A prototype pass means only that the sampled prototype checks passed. It does not validate C++ geometry, Orca feature gates, or final G-code.
- A G-code validator pass is a manually collected structural-motion result, not proof of printer compatibility.
- The validator does not recompute C++ containment, crossing, spacing, swept coverage, outside area, segment widths, or material ratio from the serialized file.
- It checks that E and F are finite and positive where required, but does not independently prove that serialized E/F values quantitatively reproduce the C++ per-segment metadata.
- It does not know whether generic `G28` homing, post-home leveling state, the configured heater targets, acceleration/current limits, collision envelope, or firmware extensions are safe for the physical printer. The test-only C++ serializer checks the exact post-offset bed polygon and swept bead, global/tool Z limits, quantized layer material, and final volumetric feeds, but this Python parser does not reproduce those checks or model real deposited-bead response.
