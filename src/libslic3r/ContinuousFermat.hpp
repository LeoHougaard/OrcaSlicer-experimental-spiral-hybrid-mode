#ifndef slic3r_ContinuousFermat_hpp_
#define slic3r_ContinuousFermat_hpp_

#include "ExPolygon.hpp"
#include "Flow.hpp"
#include "Polyline.hpp"

#include <string>
#include <vector>

namespace Slic3r {

class Layer;

namespace ContinuousFermat {

// Postconditions for a generated per-layer path. These metrics are evaluated
// against the physical extrusion footprint before normal Orca layer entities
// are removed. `emittable` covers structural, finite-flow, containment, and
// topology safety requirements; `ok` additionally covers geometric-quality
// targets such as coverage and material balance.
struct PathValidation
{
    bool ok { false };
    bool emittable { false };
    bool closed { false };
    double exact_coverage_ratio { 0.0 };
    double coverage_ratio { 0.0 };
    double outside_ratio { 0.0 };
    double material_ratio { 0.0 };
    double redeposition_ratio { 0.0 };
    size_t containment_violations { 0 };
    int crossings { 0 };
    int spacing_violations { 0 };
    int bead_overlap_violations { 0 };
    int turnback_violations { 0 };
    std::string thin_feature_class;
    double thin_feature_threshold_mm { 0.0 };
    std::string reason;
};

struct GeneratedPath
{
    Polyline path;
    // Per-segment cross-section multiplier.  Wider terminal beads slow down by
    // the inverse factor in G-code so volumetric flow stays bounded.
    std::vector<float> extrusion_multipliers;
};

// Replace a layer's normal region extrusion entities with a single continuous
// Fermat-style extrusion path when continuous slicing mode is active.
bool apply_to_layer(Layer &layer);

// Exposed for focused geometry tests and debug tooling.
// `max_line_width` is a physical bead width in millimetres. Zero preserves the
// legacy 1.60 cross-section multiplier. In all cases the implementation keeps
// the non-negotiable two-nozzle physical-width ceiling.
Polyline generate_layer_path(const ExPolygons &printable_area, const Flow &flow, double max_line_width = 0.0);
GeneratedPath generate_layer_path_with_metadata(
    const ExPolygons &printable_area,
    const Flow &flow,
    double max_line_width = 0.0);

// Validate a completed path using the same fail-closed contract applied by the
// production layer integration.  Exposed for focused safety regression tests.
PathValidation validate_layer_path(
    const ExPolygons &printable_area,
    const Flow &flow,
    const Polyline &path,
    const std::vector<float> &extrusion_multipliers,
    double max_line_width = 0.0);

} // namespace ContinuousFermat
} // namespace Slic3r

#endif
