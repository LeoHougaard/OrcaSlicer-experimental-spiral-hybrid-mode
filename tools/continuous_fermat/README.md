# Continuous Fermat One-Layer Prototype

This directory contains a standalone prototype for debugging strict continuous
slicing before integrating the algorithm into OrcaSlicer.

The prototype is intentionally dependency-free:

1. Built-in test shapes are represented as polygons with optional holes.
2. A signed-distance field is sampled over the shape.
3. Marching squares extracts inward offset contours.
4. The current experimental generator builds a cyclic contour-weave spiral from
   inward offset contours.
5. The accepted path must start and end on the outside contour. This models the
   later Orca integration where the inter-layer solver chooses the actual
   boundary cut point.
6. A verifier checks centerline containment, self-intersections, non-adjacent
   spacing violations, and coverage by rasterizing the extrusion stroke.
7. Scanline fallbacks are not accepted for this mode because they are not
   Fermat/spiral paths and do not preserve outer-wall quality.

Run all built-in cases:

```powershell
python tools/continuous_fermat/fermat_layer.py --all --draw-contours
```

Run one case:

```powershell
python tools/continuous_fermat/fermat_layer.py --shape annulus --draw-contours
```

Outputs are written to `build/continuous_fermat/`:

- `*.svg` shows the model boundary, optional raw contours, generated path,
  green start marker, and purple end marker.
- `*.png` contains a raster preview that can be inspected by image tools.
- `*.json` contains metrics that can be regression-tested.

Important limitations:

- This is a geometry prototype, not production slicer code.
- The SDF/marching-squares contour extraction should be replaced by Orca's
  Clipper offset contours during integration.
- Branch and multi-hole regions require full contour-tree CFS traversal.
- A path that visually looks spiral-like may still fail if it revisits a slot
  point, because that would create an overlap/blob in a real continuous
  extrusion.
- A passing prototype result means the path is one continuous extrusion stroke
  for the sampled geometry. It does not yet prove the final Orca integration has
  disabled all travel-generating G-code features.
