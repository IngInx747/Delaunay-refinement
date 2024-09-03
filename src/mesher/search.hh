#include "mesh.hh"

Fh search_triangle_brute_force(const TriMesh&, const Vec2&);

Fh search_triangle_straight_way(const TriMesh&, const Vec2&, const Fh&);

Fh search_triangle_guided_bfs(const TriMesh&, const Vec2&, const Fh&);

Fh fuzzy_search_triangle_boundary(const TriMesh&, const Vec2&, const double tol);

Fh fuzzy_search_triangle_boundary(const TriMesh&, const Vec2&, const double tol);
