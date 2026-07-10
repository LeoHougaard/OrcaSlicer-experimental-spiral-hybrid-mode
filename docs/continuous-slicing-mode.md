# Continuous Slicing Mode (Experimental)

> **Safety status: pre-alpha source-level research; do not send this mode's output to a printer.** Standard G-code generation, cached-file reuse, local export, and upload are hard-blocked. A C++ validation pass or Python G-code audit does not establish safe kinematics, temperatures, extrusion pressure, build-envelope containment, collision clearance, leveling state, or real bead behavior.

The developer-mode setting **Continuous slicing (pre-alpha; do not print)** (`spiral_hybrid_non_crossing`) selects the `ContinuousFermat` planner when **Spiral vase** (`spiral_mode`) is also enabled. It is no longer a simple-mode printer feature. Despite the historical configuration name, the current implementation is neither a conventional vase spiral nor the older "normal walls plus infill with fewer crossings" experiment.

## What the Current Mode Generates

For each eligible non-empty model layer, the planner:

1. accepts one hole-free axis-aligned rectangular cross-section;
2. builds one exactly closed Fermat-style path from offset contours and connectors;
3. attaches one extrusion/cross-section multiplier to every path segment;
4. validates the complete variable-width extrusion footprint; and
5. only after validation succeeds, removes the normal perimeter, infill, and thin-fill entities and replaces them with one marked `ContinuousFermat` extrusion path.

The marked path is excluded from later simplification, clipping, arc fitting, overhang resampling, small-area flow compensation, resonance-avoidance speed raising, and generic print/role flow multipliers. Those operations would invalidate the path-and-metadata pair that passed validation. The filament flow ratio is not a remaining calibration input in strict mode: it must be exactly `1.0`.

The generated path uses the external-perimeter flow and role for the complete layer. Normal wall-loop, sparse-infill, and top/bottom-fill paths are not retained. This mode attempts a solid fill of every non-empty cross-section; it does not implement a selected sparse-infill percentage.

The metadata multiplier is not a literal width multiplier. Physical width is derived from the rounded-rectangle bead cross-section. A multiplier greater than `1.0` reduces feed rate by its inverse; a multiplier below `1.0` retains the nominal feed rate rather than accelerating.

Each layer path is exactly closed. The C++ emitter requires exact equality of the internal endpoint coordinates on adjacent layers and emits a positive Z-only transition with no E or XY motion. The independent serialized-file validator separately allows its configured XY tolerance, `0.001 mm` by default. The result is one closed extrusion section **per layer**, with extrusion stopped during the Z transition, not a continuously extruding 3D helix through the object.

The classic `SpiralVase` post-processor is bypassed in this mode.

The application does not expose the resulting path as executable printer G-code. `Print::export_gcode`, the lower G-code generator boundary, cached-G-code processing, local-copy finalization, and upload preparation all reject the mode. Cached or imported files containing Continuous Fermat section markers are rejected even if the setting was subsequently disabled. A non-serialized, per-`Print` capability exists only so C++ regression tests can exercise the exact serialized contract; no profile, GUI, CLI, 3MF, or environment setting enables it in the application.

## Enabling It

Both settings must be true:

```ini
spiral_mode = 1
spiral_hybrid_non_crossing = 1
```

The GUI currently enables `spiral_mode` when Continuous slicing is selected, but profiles and CLI jobs should set both explicitly.

These flags enable research slicing only. They do not bypass the printer-ready export block.

### Legacy hybrid settings

The following settings remain in the configuration schema for compatibility with profiles from the earlier hybrid experiment. They are hidden from the current UI and are not inputs to the `ContinuousFermat` geometry or physical-footprint validator:

- `spiral_hybrid_interior_clearance`
- `spiral_hybrid_flow_mode`

Do not use either setting as evidence that an internal clearance or requested flow policy was applied. Segment extrusion is chosen by the planner and constrained by the mandatory validator. Wall-loop, top-shell, and sparse-infill controls likewise do not restore conventional wall/infill semantics after the layer entities are replaced.

## Supported Scope

The C++ integration entry point currently admits only a hole-free, four-vertex, axis-aligned rectangular vertical prism. Every non-empty layer must have exactly the same cross-section and external-perimeter flow as the layer below, and layer 0 must begin at the build plate. Admission is not acceptance: a rectangle is still rejected if its generated path fails any mandatory metric. Holes, annuli, concave or branched outlines, multiple islands, tapers, roofs, bridges, overhangs, floating layers, and all other non-rectangular shapes are unsupported and rejected before entity replacement.

The slicer also rejects the job unless it has all of the following:

- exactly one print object and one instance;
- exactly one print region, material, configured nozzle, and extruder;
- logical extruder 0 mapped to physical tool 0;
- Marlin 2 G-code semantics, with no Bambu-specific output;
- printing by layer and relative extrusion distances;
- the same positive configured nozzle temperature for the first and later layers, within the filament's configured temperature range;
- a filament flow ratio of exactly `1.0`;
- a finite positive fixed maximum volumetric speed, with adaptive volumetric speed disabled;
- a finite plate origin, finite positive global/tool-0 height limits, and a finite Z offset;
- every quantized path segment and its complete variable-width bead footprint inside the configured printable polygon;
- no bed-exclusion polygon or tool-specific printable-area polygon;
- firmware machine-limit emission disabled, so this experimental mode cannot overwrite the controller's stored limits;
- no support material, raft, effective skirt, draft shield, brim, or prime tower;
- no object-cancellation G-code or power-loss recovery;
- no external post-processing, post-export line numbering, pressure equalization, fixed pressure advance, or adaptive pressure advance;
- no smooth timelapse, first-layer scanning, wrapping/clumping detection, or ironing;
- no per-layer custom G-code or calibration mode;
- comments only in file-start, machine-start, filament-start, layer-change, timelapse, extrusion-role, machine-end, and filament-end G-code hooks;
- fan speedup and fan kick-start set to zero; and
- no auxiliary fan, air filtration, or chamber-temperature control output.

Strict mode refuses cached extrusion entities because they cannot be re-certified against the current slices. The layers must be regenerated.

In test-only serialized artifacts, retraction, wipe, and Z-hop profile values are not used for the first strict approach or for strict layer transitions. The first approach is direct XY at the already checked first-layer Z; after the first section, layer changes are emitted directly as strictly increasing Z-only moves. Every emitted layer Z and the optional final lift are bounded by the lower positive value of the global and tool-0 printable-height limits. This does **not** certify the first XY traverse: after generic `G28`, the slicer does not know the machine's physical home XY or prove that the approach segment avoids clips or other hardware. That unresolved machine-specific approach is one reason application export is blocked.

The test-only serializer emits a deterministic pre-model sequence after the normal bed-temperature setup: modal normalization, `T0`, `M900 K0`, `G28` homing, and a cooling-aware blocking `M109 R...` wait at the validated nozzle target. This removes arbitrary profile start motion, persistent Marlin linear advance, and the possibility of beginning extrusion before the hotend has reached the target from either direction. It still does not prove that generic all-axis homing is safe, establish the final physical home position, or establish a valid bed-leveling mesh. Marlin may disable leveling during `G28` depending on firmware configuration.

## Mandatory Per-Layer Validation

The production integration accepts a generated layer only when **every** condition passes:

- the polyline has at least four points and is exactly closed;
- the metadata count exactly equals the segment count;
- every extrusion multiplier is finite and within `0.60` through `1.60`;
- every resulting physical segment width is at least `0.85` times and at most `2.0` times the nozzle diameter, and the nominal width is no greater than `2.0` times the nozzle diameter;
- centerline containment violations are zero;
- non-local self-crossings are zero;
- non-local physical bead-overlap violations are zero;
- non-local spacing warnings are at most 3; a warning is a centerline separation below 75% of nominal spacing after excluding adjacent segments and pairs within a four-spacing local path arc;
- hairpin turnback violations are zero; the check flags turns sharper than 135 degrees when both legs are at least 35% of line width;
- sampled swept-bead coverage is at least `0.990`, using 70 cells along the layer's longest dimension;
- exact swept-area coverage is at least `0.970`;
- swept area outside the printable region is at most `0.020` of printable area; and
- modeled deposited-material ratio is within `0.980` through `1.020`.

These are conjunctive limits: no score or visual result can override a failure. The slicer raises `SlicingError`, and the failing layer's normal entities are not cleared. A failure on a later layer should not be described as rolling back layers already processed in memory. These thresholds are regression criteria for an experimental algorithm, not proof that a bead is printable with a particular nozzle, material, profile, or machine.

## Final G-code Audit

Accepted Continuous Fermat paths are delimited by standalone comments:

```gcode
;_CONTINUOUS_FERMAT_BEGIN
...
;_CONTINUOUS_FERMAT_END
```

The application does not produce a printer-ready file to validate. The Python tool is for historical files and exact artifacts deliberately emitted by the C++ regression test capability; it is not invoked automatically. A manual pass remains useful development evidence, but it is not authorization to print. Both expected model-layer and section counts are mandatory unless the explicitly diagnostic-only `--allow-unknown-counts` flag is used:

```powershell
python tools\continuous_fermat\validate_gcode.py path\to\model.gcode `
    --expected-layers 120 --expected-sections 120
```

A strict section must be preceded by explicit `G21`, `G90`, `M83`, `T0`, `M200 D0`, `M220 S100`, `M221 S100`, `M900 K0`, and `M413 S0` state. The validator also requires an exact all-axis `G28` followed by positive-target `M109 R...` before the first section. It enforces a deliberately narrow Marlin-2-oriented, millimetre, linear-`G1` motion contract. It checks marker and layer accounting, finite parsed values, known positive feed and layer Z, positive extrusion on section XY moves, section closure, modal state, and same-XY Z-only handoffs. It rejects pressure-advance mutations, unsupported protected commands, arcs, commandless modal axis lines, coordinate changes, tool/material changes, and unmarked positive-E XY motion.

The validator retains a diagnostic exception for reviewing historical or non-production files containing a pre-model purge:

```powershell
python tools\continuous_fermat\validate_gcode.py path\to\model.gcode `
    --expected-layers 120 --expected-sections 120 `
    --allow-unmarked-before-first-section
```

The exception is not compatible with the current production gate, which rejects executable start hooks. It applies only before the first section and does not certify the purge location, heater state, homing sequence, or machine envelope. It never permits a skirt, brim, or model extrusion between or after sections.

After the final section, the validator rejects XY or E motion, firmware retract, `G92`, and modal changes. It permits at most one positive Z-only clearance lift of no more than `2 mm`, followed only by its tightly allowlisted shutdown/progress commands. The test-only serializer uses a built-in no-XY shutdown instead of configured executable end scripts.

The validator does **not** re-run the C++ geometry checks or independently prove that serialized E/F values reproduce the validated segment widths and material ratio. The C++ emitter now rechecks non-collapsed quantized XY, the exact transformed bed polygon and swept bead, per-layer quantized-E material ratio, and a downward-rounded volumetric feed cap; export regressions parse the resulting file again. The Python tool does not duplicate those quantitative calculations. It does not prove homing safety, heater and fan suitability, acceleration or current limits, collision clearance, firmware extensions, or real deposited-bead response. It validates exact forms for normalizers and final shutdown commands, but does not semantically certify every permitted acceleration, fan, temperature, or progress argument.

## Current Limitations

- C++ integration shape support is rectangle-only. Branch, hole, annulus, and general-topology code remains research-only and is unreachable through the guarded entry point.
- No rectangle dimension, nozzle, line-width, or layer-height envelope is physically certified. Some otherwise eligible rectangle/flow phases are intentionally rejected by regression tests.
- A layer can pass the numerical thresholds while still being unsuitable for a real material/process combination.
- The material model uses nominal rounded-rectangle bead geometry; it is not a pressure, cooling, corner-flow, or adhesion simulation.
- The complete stroke uses external-perimeter flow/role semantics; bridge, top-surface, bottom-surface, and infill-specific behavior is not represented.
- The Python prototype and continuous-path lab are research tools. Their results are not substituted for the C++ production validator, and their acceptance contracts are not production-equivalent.
- No current test establishes firmware-independent behavior or physical-printer certification.
- Printer-ready export and upload remain intentionally unavailable until an allowlisted, physically qualified machine/firmware contract establishes homing, leveling, a safe heater-wait pose, travel bounds/keep-outs, and the complete approach corridor.

Use classic Spiral vase with **Continuous slicing** disabled when a conventional hollow, single-wall vase path is required. See [Continuous Slicing Development Audit](continuous-slicing-development-audit.md) for the iteration history, evidence, unresolved risks, and research basis.
