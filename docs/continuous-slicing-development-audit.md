# Continuous Slicing Development Audit

## Conclusion

The Continuous slicing mode is **pre-alpha source-level research software, not a production-ready printer mode. Do not print its output.** Its strongest current property is fail-closed rejection at several boundaries: unsupported jobs and shapes are rejected, the complete generated layer path is checked before normal entities are cleared, and the test-only G-code emitter enforces one marked closed stroke per layer with an exact same-XY handoff. Standard application/CLI export, cached-file reuse, local copy, and upload are now hard-blocked.

That hardening is not printer certification. The C++ integration scope has been deliberately narrowed to identical hole-free axis-aligned rectangular layers, and an eligible rectangle may still fail the mandatory metrics. A representative three-layer artifact emitted through a non-serialized C++ test capability passes the serialized-motion, material, and volumetric-cap regression, but the supported dimension/flow matrix is not complete. The independent Python G-code validator remains manual and does not revalidate footprint geometry, extrusion volume, machine limits, homing/leveling safety, or real extrusion behavior. Broader machine-aware evidence remains a release blocker, and physical trials must not begin while the hard export gate remains necessary.

The fixed-grid continuous-path-lab v6 rerun on 2026-07-09 passed 0 of 10 built-in shapes under the metrics implemented by the lab. Its rectangle reached `0.993` sampled coverage but still had 77 bead-overlap violations and 8 turnbacks. That result must not be described as a run of the complete production contract: the lab does not calculate production exact swept-area coverage, outside-area ratio, deposited-material ratio, or nozzle-relative physical-width bounds.

## Development Iterations

| Iteration | Main strength | Main issue found in audit |
| --- | --- | --- |
| Initial spiral-hybrid controls (`6ad097bb0b`, April 2026) | Preserved ordinary Orca wall/infill generation while relaxing classic vase restrictions and enabling crossing avoidance. | It was not a single continuous fill. Walls, infill, and layer transitions could still create multiple extrusion strokes; "non-crossing" was a preference, not a proven postcondition. |
| First documentation (`e9e2825ec3`, `cadc16c929`, June 2026) | Clearly labeled the feature fork-specific and experimental. | The guide described the initial walls-plus-infill behavior and became dangerously stale after the planner was replaced. |
| C++ Connected Fermat integration (`eaac3444e2`) | Added a production-side one-path-per-layer planner using Orca geometry and bypassed classic vase post-processing. | Route selection was heuristic, difficult topology was not decomposed as in the research algorithm, and early acceptance did not prove full physical coverage before replacing normal entities. |
| Residual-gap repairs (`0bcbf93d7c`) | Addressed visible pockets left by nested offset contours. | Local repair heuristics could trade underfill for overlap, turnback, or self-proximity defects; visual improvement was not an adequate safety criterion. |
| Fail-closed integration and markers (`dd3a6d6c6c`) plus startup fix (`7983484532`) | Added slicing errors instead of silently accepting invalid generation, explicit section markers, a standalone G-code checker, and corrected a startup crash. | The first checker trusted markers too much and did not fully police motion outside them; geometry, G-code, and feature-scope criteria were still fragmented. |
| Standalone hybrid prototype (`bf6ecf4354`) | Made path geometry visible and enabled rapid dependency-free experiments and metrics. | It is raster/SDF-based research code, not the Clipper-based production implementation, so a prototype pass cannot validate an Orca export. |
| Continuous-path lab v2-v6 (`b7c290a695`) and experimental v7 work | Added SVG/UI inspection, topology cases, sampled coverage, spacing/crossing/containment and physical-bead metrics, and a repeatable benchmark loop. | The audited v6 baseline passed 0/10 under the lab's implemented conjunction, which is smaller than the production contract. A prior endpoint-dropping v7 experiment regressed simple rectangle coverage to roughly `0.292`; preserving a bad candidate and rejecting it is safer than deleting required geometry. |
| Current hardened worktree | Models segment extrusion and deposited volume, validates the final swept bead, gates unsupported print features and shapes, enforces exact internal layer handoff, checks quantized E/F/XY against material and volumetric limits, hardens the independent G-code parser with adversarial tests, and blocks all standard printer-ready export/upload paths. | The C++ entry point is rectangle-only. Qualifying 1.2 mm and common 0.4 mm-nozzle rectangle fixtures pass while unsafe phases reject, but a machine homing/leveling/approach contract, complete supported-scope matrix, and physical evidence are not established. |

## Strengths of the Current Design

- **Fail-closed layer replacement:** normal perimeter/fill entities for a layer are cleared only after its generated path passes final validation.
- **One explicit path per non-empty layer:** the replacement entity is unsorted, marked as Continuous Fermat, and paired with per-segment metadata.
- **Physical-footprint checks:** acceptance uses variable-width swept beads and modeled material, not centerline appearance alone.
- **Bounded adaptive extrusion:** multipliers are within `[0.60, 1.60]` and physical segment widths within `[0.85, 2.0]` nozzle diameters. Wider segments are slowed by `1 / multiplier`; narrowed segments retain nominal feed.
- **Layer-handoff enforcement:** the emitter allows no ordinary extrusion in a strict job and requires the next Continuous Fermat section to begin at the previous endpoint.
- **Post-validation immutability:** simplification, clipping, arc fitting, generic flow boosts, small-area compensation, and resonance speed raising are bypassed for the validated path.
- **Exact vertical-prism scope:** every admitted layer cross-section and external-perimeter flow must equal the layer below, and the first layer must start at the plate.
- **Restricted feature surface:** common auxiliary strokes, post-processors, and multi-object/material/tool configurations are rejected before generation.
- **Deterministic pre-model state:** executable start hooks and firmware-limit overrides are rejected; the test-only serializer emits its own modal normalization, `T0`, `M900 K0`, all-axis `G28`, and cooling-aware blocking `M109 R...` sequence.
- **Printer-output interlock:** ordinary core export, the lower generator boundary, cached-file processing, GUI finalization, and upload preparation all reject this mode. Cached/imported Continuous Fermat markers are rejected even after the configuration flag is disabled. The regression-only bypass is in-memory, per-`Print`, non-serialized, and unavailable to application settings.
- **Serialized machine-envelope and flow checks:** actual post-offset, quantized segments and their complete variable-width bead footprint must remain inside the exact configured bed polygon; emitted Z is bounded by global/tool-0 height; quantized E is rechecked against target layer material; and integral feeds retain a conservative margin below the quantized volumetric cap.
- **Independent serialized audit:** the manual Python parser checks final marker/accounting/modal/motion structure rather than trusting the C++ entity graph. It does not re-run the C++ footprint contract.

## Unresolved Issues

### Production geometry

- The only admitted shape is a hole-free, four-vertex, axis-aligned rectangular prism. There is no documented rectangle dimension, line-width, nozzle, or layer-height envelope that is guaranteed to pass.
- Exact and sampled coverage expose different geometric failures, but neither predicts real bead formation at sharp width or direction changes.
- An outside-area ratio of `0.020` and up to three spacing warnings are engineering regression thresholds, not universal safety limits. A machine or process may require zero outside deposition and tighter spacing.
- Pairwise crossing, spacing, and bead checks can be costly. Runtime and memory still need documented bounds for every supported input envelope.

### Research topology

- Branch, hole, annulus, concave, and multi-island algorithms remain research-only; the production gate rejects these shapes before planning.
- The research planners use specialized heuristics rather than a complete deterministic decomposition and graph traversal.
- The lab acceptance conjunction is not production-equivalent. In particular, it omits production exact swept coverage, outside area, material ratio, and nozzle-relative physical-width bounds.
- Future topology support must not be enabled merely because a lab path looks good or improves a weighted score. It must pass the complete production geometry and serialized-output contract.

### Extrusion semantics

- The complete solid layer uses external-perimeter flow and role metadata. It does not preserve bridge, surface, infill, cooling, or speed semantics of the paths it replaces.
- Cross-section and external-perimeter flow must match exactly between layers; tapers, roofs, bridges, and overhangs are unsupported rather than merely constrained to supported material.
- The deposited-volume model assumes rounded-rectangle beads. It does not model pressure advance, acceleration-induced corner overfill, die swell, cooling, layer adhesion, or extrusion force.
- The configured filament flow ratio is locked to `1.0`. The C++ emitter compares serialized E totals with target layer material and enforces the volumetric cap, and export tests parse those values again; the manually invoked Python validator does not reproduce those quantitative checks.
- No printer-profile-aware lower speed, heater capacity, pressure, or sustained-flow check certifies widened or narrowed segments.

### Whole-file and machine safety

- Standard printer-ready export and upload are blocked before a destination is created or overwritten. Test-only serialization rejects executable file/machine/filament start hooks, line-number rewriting, pressure advance, bed/tool-area exclusions, and firmware-limit emission. It resets Marlin linear advance with `M900 K0`, homes all axes with `G28`, waits in both heating and cooling directions with `M109 R...`, uses a direct no-Z-hop first approach, and replaces executable configured end scripts with a built-in no-XY shutdown.
- The independent validator is a manually invoked, narrow Marlin-2-oriented motion parser, not a general G-code interpreter or an automatic export gate.
- It permits some pre-first-section setup commands and cannot certify homing behavior, heater/fan suitability, acceleration, current, or collision clearance. The C++ emitter now rejects non-finite or out-of-envelope protected XY, but that is not a complete collision proof.
- In particular, all-axis `G28` leaves the slicer without a certified physical home XY. The direct non-extruding move to the first stroke has bounded Z and a validated destination, but its unknown-start segment is not proven to remain within a collision-free machine travel envelope. A machine-profile-specific homing/approach contract is required before printer use.
- Marlin can also disable bed leveling during `G28` unless firmware-specific restore/enable behavior is configured. The current repository has no structured, immutable machine contract for homing direction/position, persistent offsets, leveling/mesh state, motion bounds, keep-outs, toolhead footprint, certified safe Z, heater-wait pose, or a validated approach corridor. `printable_area` is an extrusion region, not that motion contract.
- A validator pass cannot prove that the selected model, material, nozzle, profile, and machine combination is safe.

## Formal Research Acceptance Contract

No weighted score, visual judgment, or "least bad" candidate may override a failed condition. Every applicable condition is an AND requirement.

### Per-layer C++ acceptance

For every non-empty production model layer:

1. The printable input is exactly one hole-free, four-vertex, axis-aligned rectangular `ExPolygon` in one non-empty region.
2. The first layer starts at the build plate. Every later cross-section and external-perimeter flow exactly equals the layer below.
3. There is exactly one generated polyline with at least four points, and its first and last points are identical.
4. Per-segment metadata cardinality exactly matches the path, and every multiplier is finite and within `[0.60, 1.60]`.
5. Every physical segment width is at least `0.85` and at most `2.0` nozzle diameters; nominal width is also at most `2.0` nozzle diameters.
6. Centerline containment violations: `0`.
7. Non-local self-crossings: `0`.
8. Non-local physical bead-overlap violations: `0`.
9. Non-local spacing warnings: at most `3`; a warning is a centerline distance below `0.75` times nominal spacing after adjacent segments and pairs within a four-spacing local path arc are excluded.
10. Hairpin turnbacks: `0`; the check flags turns sharper than 135 degrees when both legs are at least `0.35` times line width.
11. Sampled swept-bead coverage: at least `0.990`, sampled at 70 cells on the longest bounding-box dimension.
12. Exact swept-area coverage: at least `0.970`.
13. Swept area outside the printable region divided by printable area: at most `0.020`.
14. Modeled deposited material divided by target solid-layer material: from `0.980` through `1.020`.

The supported-scope release matrix must show that intended rectangle/flow combinations pass deterministically and every out-of-scope or unsafe case rejects. Passing arbitrary research-lab shapes is not required to retain rectangle-only scope, but it is mandatory before widening production scope to those shapes.

### Print-job scope

The job has exactly one object, instance, region, material, configured nozzle, and extruder. Logical extruder 0 maps to physical tool 0. It uses Marlin 2, is not Bambu-specific, prints by layer with relative E, has a finite plate origin/Z offset and finite positive global/tool-0 height limits, has equal positive first/later nozzle temperatures within the configured filament range, a filament flow ratio of `1.0`, and a finite positive fixed maximum volumetric speed with adaptive volumetric speed disabled. Bed exclusions and tool-specific printable-area polygons are absent. Firmware machine-limit emission is disabled.

Support, raft, skirt, draft shield, brim, prime tower, object cancellation, power-loss recovery, external post-processing, post-export line numbering, pressure equalization, fixed/adaptive pressure advance, smooth timelapse, first-layer scanning, ironing, wrapping/clumping detection, per-layer custom G-code, and calibration modes are absent. File/machine/filament start hooks, layer/timelapse/role-change hooks, and end hooks contain comments only. Fan speedup and kick-start are zero; auxiliary fan, air filtration, and chamber control are disabled. Strict jobs refuse cached extrusion entities.

### Test-only serialized and independently audited G-code

1. Every protected section begins with explicitly established `G21`, `G90`, `M83`, `T0`, `M200 D0`, `M220 S100`, `M221 S100`, `M900 K0`, and `M413 S0` state; the first also requires exact all-axis `G28` followed by positive-target `M109 R...`.
2. The number of standalone `BEGIN`/`END` sections equals the expected model-layer and expected-section counts, with exactly one section per layer marker.
3. Every section has finite millimetre coordinates, known positive Z and feed, linear `G1` motion only, strictly positive E on every XY extrusion move, no Z motion, and a closed XY endpoint.
4. C++ internal closure and adjacent-layer endpoint equality are exact. The independent serialized validator applies its configured XY tolerance, `0.001 mm` by default.
5. Between sections, XY and E do not move; firmware retract, `G92`, and disruptive modal commands are absent; Z increases finitely; and the next section begins at the preceding endpoint.
6. After the first section begins, no unmarked positive-E XY extrusion occurs.
7. After the final section there is no XY/E/retract/reset/modal motion. At most one positive Z-only lift of no more than `2 mm` is permitted, followed only by tightly allowlisted shutdown/progress commands.
8. The current test-only serialized artifact contains no configured purge/start extrusion; the Python pre-first-section exception is diagnostic-only for historical files and other non-production artifacts.
9. Post-offset quantized XY segments are nonzero, their variable-width swept beads remain inside the exact configured bed polygon, and their Z remains within the global/tool-0 limit. The emitter rechecks quantized-E layer material and bounds final feed against the fixed volumetric limit; the representative export regression independently parses those values.

The Python validator checks these serialized structural conditions manually. It does not recompute geometry, check E/F against segment-width metadata, or establish machine safety.

### Required evidence before changing the stage

- All production geometry tests for the supported rectangle/flow matrix pass deterministically, and every documented unsupported or unsafe geometry rejects.
- Feature-gate tests prove each unsupported source of auxiliary extrusion or state mutation rejects.
- Adversarial G-code tests prove malformed markers, modal-state tricks, non-finite values, arcs/macros, extra sections, unsafe postamble motion, and unmarked extrusion fail closed.
- Representative complete test-only Orca artifacts pass the independently invoked validator with exact expected counts.
- Runtime and memory are bounded for an explicitly documented supported-input envelope.
- Printer-profile-aware checks cover build volume, firmware dialect, motion limits, extrusion rate, temperatures, and reviewed start/end behavior.
- Before production scope is widened, every applicable research shape passes lab and C++ production-equivalent metrics, including exact outside/material/width checks that the lab currently lacks.
- Controlled physical trials use a sacrificial/test machine, conservative profiles, active supervision, emergency-stop readiness, and documented acceptance measurements. Physical trials must not begin while the software-only criteria fail.
- The hard export/upload gate may be removed only for an immutable, allowlisted, physically qualified printer/firmware contract that proves the complete serialized sequence from homing through leveling, heater wait, safe-Z approach, descent, and first stroke.

## Research Assessment and Recommended Direction

The implementation is inspired by Zhao et al.'s Connected Fermat Spirals (CFS), but the paper's result depends on decomposing a connected region into spiral-fillable subregions and connecting those fills by traversing a graph. Merely weaving all offset contours together does not inherit that construction's properties. A future general-topology planner should therefore:

1. build a deterministic decomposition/contour tree for islands, branches, and holes;
2. generate and validate a low-curvature continuous fill inside each subregion;
3. choose ports and graph connections using reserved, printable corridors with exact containment and swept-width constraints;
4. assemble one global Eulerian/tree traversal without deleting difficult contours or adding disconnected patches;
5. optimize only among candidates that already satisfy topology and physical constraints; and
6. run exact swept-bead, material, nozzle-width, and turnback checks on the assembled path before preserving the existing fail-closed C++ and G-code boundaries.

Globally continuous hybrid and Euler-graph work provides useful alternatives for residual regions and complex topology. Those approaches should be evaluated as new planners behind the same contract, not mixed in as an unchecked fallback. A hybrid contour/zigzag stitch may improve difficult pocket coverage, but it is acceptable only if it remains one connected non-crossing stroke and passes every physical metric.

## Primary Sources

- Haisen Zhao et al., [Connected Fermat Spirals for Layered Fabrication](https://haisenzhao.github.io/CFS/), ACM Transactions on Graphics 35(4), SIGGRAPH 2016. [Paper PDF](https://haisenzhao.github.io/CFS/files/fermat_spirals.pdf), DOI [10.1145/2897824.2925958](https://doi.org/10.1145/2897824.2925958).
- Lingwei Xia et al., [Globally continuous hybrid path for extrusion-based additive manufacturing](https://doi.org/10.1016/j.autcon.2022.104175), Automation in Construction 137 (2022), 104175.
- Prashant Gupta, Bala Krishnamoorthy, and Gregory Dreifus, [Continuous toolpath planning in a graphical framework for sparse infill additive manufacturing](https://doi.org/10.1016/j.cad.2020.102880), Computer-Aided Design 127 (2020), 102880. [Author preprint](https://arxiv.org/abs/1908.07452).
- Marlin Firmware, [`G0`/`G1` Linear Move](https://marlinfw.org/docs/gcode/G000-G001.html) and [`M82` E Absolute](https://marlinfw.org/docs/gcode/M082.html), used to define the narrow modal-motion contract checked by the independent validator.
- Marlin Firmware, [`G90`](https://marlinfw.org/docs/gcode/G090.html), [`G91`](https://marlinfw.org/docs/gcode/G091.html), and [`M83`](https://marlinfw.org/docs/gcode/M083.html), which document Marlin's E-axis override reset semantics.
- Marlin Firmware, [`M109` Wait for Hotend Temperature](https://marlinfw.org/docs/gcode/M109.html), where `R` waits while heating or cooling, and [`M900` Linear Advance Factor](https://marlinfw.org/docs/gcode/M900.html), where `K0` disables Linear Advance; these define the test-only serializer's blocking temperature wait and pressure-advance reset.
- Marlin Firmware, [`G28` Auto Home](https://marlinfw.org/docs/gcode/G028.html) and [`M420` Bed Leveling State](https://marlinfw.org/docs/gcode/M420.html), documenting configurable homing/safe-homing behavior and the firmware-dependent leveling state after homing.
- Marlin Firmware, [`M206` Set Home Offsets](https://marlinfw.org/docs/gcode/M206.html) and [`M114` Get Current Position](https://marlinfw.org/docs/gcode/M114.html), documenting persistent coordinate offsets and runtime position reporting that a standalone static G-code file cannot consume as a pre-move safety handshake.
- Klipper, [G-Codes](https://www.klipper3d.org/G-Codes.html), [Command Templates](https://www.klipper3d.org/Command_Templates.html), and the official [`gcode_move.py`](https://github.com/Klipper3d/klipper/blob/master/klippy/extras/gcode_move.py), illustrating its independent XYZ/E modal state and why the current Marlin-2-only mode cannot claim Klipper compatibility.

## Related Repository Evidence

- Operator-facing behavior: [continuous-slicing-mode.md](continuous-slicing-mode.md)
- Production planner and validator: `src/libslic3r/ContinuousFermat.cpp`
- Production emission boundary: `src/libslic3r/GCode.cpp`
- Unsupported-feature gate: `src/libslic3r/Print.cpp`
- Geometry safety tests: `tests/libslic3r/test_continuous_fermat.cpp`
- Print-scope tests: `tests/fff_print/test_print.cpp`
- Strict serialized-G-code audit: `tools/continuous_fermat/validate_gcode.py`
- G-code adversarial tests: `tools/continuous_fermat/test_validate_gcode.py`
- Lab contract and benchmark loop: `tools/continuous_path_lab/V6_LOOP.md`
