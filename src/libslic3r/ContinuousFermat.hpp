#ifndef slic3r_ContinuousFermat_hpp_
#define slic3r_ContinuousFermat_hpp_

#include "ExPolygon.hpp"
#include "Flow.hpp"
#include "Polyline.hpp"

namespace Slic3r {

class Layer;

namespace ContinuousFermat {

// Replace a layer's normal region extrusion entities with a single continuous
// Fermat-style extrusion path when continuous slicing mode is active.
bool apply_to_layer(Layer &layer);

// Exposed for focused geometry tests and debug tooling.
Polyline generate_layer_path(const ExPolygons &printable_area, const Flow &flow);

} // namespace ContinuousFermat
} // namespace Slic3r

#endif
