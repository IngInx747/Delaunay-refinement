#ifndef SEARCH2D_HH
#define SEARCH2D_HH

#include "mesh.hh"

////////////////////////////////////////////////////////////////
/// General search
////////////////////////////////////////////////////////////////

Fh search_triangle_brute_force(const TriMesh&, const Vec2&);

Fh search_triangle_zigzag(const TriMesh&, const Vec2&, const Fh &first_face);

Fh search_triangle_linear(const TriMesh&, const Vec2&, const Fh &first_face);

//Fh search_triangle_guided_bfs(const TriMesh&, const Vec2&, const Fh&);

Fh fuzzy_search_triangle_boundary(const TriMesh&, const Vec2&, const double tol);

Fh fuzzy_search_triangle_boundary(const TriMesh&, const Vec2&, const double tol);

////////////////////////////////////////////////////////////////
/// Ray tracer
////////////////////////////////////////////////////////////////

enum class RAY_STATUS: int { MISS, EDGE, VERT };

class RayTracer
{
public:

    RayTracer(const TriMesh &_m): m_(_m), st_(RAY_STATUS::MISS) {}

    // current intersecting edge iff status is EDGE
    const Hh &halfedge_handle() const { return hh_; }

    // current intersecting vertex iff status is VERT
    const Vh &vertex_handle() const { return vh_; }

    // status of last intersection test
    RAY_STATUS status() const { return st_; }

    // search for next intersecting primitive
    void next();

    // initialization utils

    void set_halfedge_handle(const Hh &_hh) { hh_ = _hh; }
    void set_vertex_handle(const Vh &_vh) { vh_ = _vh; }
    void set_status(RAY_STATUS _st) { st_ = _st; }
    void set_u0(const Vec2 &_u) { u0_ = _u; }
    void set_u1(const Vec2 &_u) { u1_ = _u; }
    const TriMesh &mesh() const { return m_; }

protected:

    // background mesh
    const TriMesh &m_;

    // status of intersection test
    RAY_STATUS st_;

    // primitive variant's edge
    Hh hh_;

    // primitive variant's vertex
    Vh vh_;

    // start and end of the path
    Vec2 u0_, u1_;
};

void init(RayTracer&, const Vec2&, const Vec2&);

void init(RayTracer&, const Fh&, const Vec2&);

void init(RayTracer&, const Vh&, const Vec2&);

void init(RayTracer&, const Vh&, const Vh&);

#endif