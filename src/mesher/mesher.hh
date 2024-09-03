#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "mesh.hh"

int make_delaunay(TriMesh&);

int triangulate(
    const std::vector<Vec2>&,
    const std::vector<Int2>&,
    TriMesh&,
    std::unordered_map<int, int>       &duplicated_vertices,
    std::vector<std::pair<Int2, Int2>> &intersecting_segments,
    std::vector<std::pair<Int2, int>>  &overlapping_vertices);

int triangulate(
    const std::vector<Vec2>&,
    const std::vector<Int2>&,
    TriMesh&);

int hide_exterior_region(TriMesh&, const std::vector<Vec2> &seeds);

int laplacian_smoothing(TriMesh&, const double step, const int max_num_iter);

int local_CVT_smoothing(TriMesh&, const double step, const int max_num_iter);
