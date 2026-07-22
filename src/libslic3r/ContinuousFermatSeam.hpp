#ifndef slic3r_ContinuousFermatSeam_hpp_
#define slic3r_ContinuousFermatSeam_hpp_

#include "ExPolygon.hpp"
#include "Polyline.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace Slic3r::ContinuousFermat {

// A zero-width topological cut. Loop zero is the outer boundary; loop N + 1
// is input hole N. The two logical sides of this cut share these coordinates,
// but remain distinct edges in HoleRemovalResult::boundary.
struct CutSeam
{
    size_t parent_loop { 0 };
    size_t child_loop { 0 };
    size_t parent_vertex { 0 };
    size_t child_vertex { 0 };
    Point  parent_point;
    Point  child_point;
    double length { 0.0 };
};

// A disk-like topological representation of one connected ExPolygon. The
// boundary is weakly simple: each seam occurs twice in opposite directions.
// It is planning topology, not an extrusion path and not a polygon intended
// for boolean cleanup.
struct HoleRemovalResult
{
    bool                 ok { false };
    Polyline             boundary;
    // Integer polygon engines cannot retain two exactly coincident sides of a
    // zero-width seam. planning_domain is the equivalent computational disk:
    // the seam is widened by the smallest representable amount needed to make
    // every hole part of the outer boundary. Its area loss is negligible and
    // it is used only for centerline planning; validation still uses the
    // untouched printable ExPolygon.
    ExPolygons           planning_domain;
    std::vector<CutSeam> seams;
    std::string          reason;
};

// Connect every hole to the outer contour (or another already connected
// contour) with a deterministic, containment-safe minimum spanning seam tree.
// No material is removed and no finite-width slit is introduced.
HoleRemovalResult remove_holes_topologically(const ExPolygon &area);

} // namespace Slic3r::ContinuousFermat

#endif
