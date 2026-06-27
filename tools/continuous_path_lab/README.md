# Continuous Path Lab

Standalone local tester for the experimental one-extrusion slicer path.

Run:

```powershell
python tools\continuous_path_lab\server.py --port 8765
```

Open:

```text
http://127.0.0.1:8765/
```

The lab accepts built-in polygons, simple DXF outlines, polygon JSON, and STL
cross-sections. It generates SDF iso-contours, builds a contour containment
tree, routes the tree with a deliberate branch order, then reports coverage,
crossing, spacing, and containment metrics.

This is intentionally outside OrcaSlicer integration. The current backend is
pure Python and dependency-free so failures are easy to reproduce. The API shape
is meant to survive a later C++ geometry backend:

```text
polygon layer + planner options -> contours + contour tree + path + metrics
```

The `legacy_cfs` algorithm option calls the existing `tools/continuous_fermat`
prototype for comparison.

`contour_tree_v4` keeps v3's hole barriers, then tries a small deterministic
set of start/winding variants and returns the best validated path. It is a lab
recovery layer for testing; v2 and v3 remain available as fixed baselines.

`contour_tree_v5` is the coverage-audit branch. It keeps all printable
iso-contours, attempts residual gap spirals after the first route, and reports a
strict coverage audit that separates true underfill from internal bead overlap.

`contour_tree_v6` starts from v5 and adds verifier-guided repairs: printed hole
barriers reserve a full bead corridor, same-contour arc pairs are treated
semantically by the spacing audit, and existing segments may be locally bent
through strict-audit underfill cells. This keeps the path continuous without
adding separate islands or per-shape patches.
