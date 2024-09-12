#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "pred2D.hh"
#include "distance.hh"
#include "triangle.hh"
#include "segment.hh"
#include "mesh.hh"

using namespace OpenMesh;

////////////////////////////////////////////////////////////////
/// pred2D
////////////////////////////////////////////////////////////////

/// Check if a point(u) is inside the triangle(u0, u1, u2), inclusively
static inline bool is_inside(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u)
{
    const int r0 = orientation(u1, u2, u);
    const int r1 = orientation(u2, u0, u);
    const int r2 = orientation(u0, u1, u);

    return
        r0 == r1 && r0 == r2 || // in triangle (without r0 != 0)
        r0 == 0  && r1 == r2 || // on edge (u1, u2)
        r1 == 0  && r2 == r0 || // on edge (u2, u0)
        r2 == 0  && r0 == r1 || // on edge (u0, u1)
        r0 == 0  && r1 == 0  || // on vert u2
        r1 == 0  && r2 == 0  || // on vert u0
        r2 == 0  && r0 == 0  ;  // on vert u1
}

static inline bool is_intersecting(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
{
    const int ru0 = orientation(v0, v1, u0);
    const int ru1 = orientation(v0, v1, u1);
    const int rv0 = orientation(u0, u1, v0);
    const int rv1 = orientation(u0, u1, v1);

    return
        ru0 * ru1 < 0 && rv0 * rv1 < 0        || // intersecting
        ru0 == 0 && ru1 != 0 && rv0*rv1 <= 0  || // u0 lies on [v0,v1]
        ru1 == 0 && ru0 != 0 && rv0*rv1 <= 0  || // u1 lies on [v0,v1]
        rv0 == 0 && rv1 != 0 && ru0*ru1 <= 0  || // v0 lies on [u0,u1]
        rv1 == 0 && rv0 != 0 && ru0*ru1 <= 0  ;  // v1 lies on [u0,u1]
}

////////////////////////////////////////////////////////////////
/// utils
////////////////////////////////////////////////////////////////

static inline bool is_inside(const TriMesh &mesh, const Fh &fh, const Vec2 &u)
{
    const auto hh = mesh.halfedge_handle(fh);
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle  (mesh.next_halfedge_handle(hh)));
    return is_inside(u0, u1, u2, u);
}

static inline bool is_intersecting(const TriMesh &mesh, const Hh &hh, const Vec2 &u0, const Vec2 &u1)
{
    const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    return is_intersecting(u0, u1, v0, v1);
}

static inline Vec2 centroid(const TriMesh &mesh, const Fh &fh)
{
    const auto hh = mesh.halfedge_handle(fh);
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle  (mesh.next_halfedge_handle(hh)));
    return (u0 + u1 + u2) / 3.;
}

static inline double hausdorff_distance(const TriMesh &mesh, const Hh &hh, const Vec2 &u)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    return hausdorff_distance({ u0,u1 }, u);
}

static inline double hausdorff_distance(const TriMesh &mesh, const Eh &eh, const Vec2 &u)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0)));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (mesh.halfedge_handle(eh, 0)));
    return hausdorff_distance({ u0, u1 }, u);
}

static inline double hausdorff_distance(const TriMesh &mesh, const Fh &fh, const Vec2 &u)
{
    const auto hh = mesh.halfedge_handle(fh);
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle  (mesh.next_halfedge_handle(hh)));
    return hausdorff_distance({ u0,u1,u2 }, u);
}

static inline SEG_LOC locate(const TriMesh &mesh, const Eh &eh, const Vec2 &u)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0)));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (mesh.halfedge_handle(eh, 0)));
    return locate(u0, u1, u);
}

////////////////////////////////////////////////////////////////
/// Triangle search
////////////////////////////////////////////////////////////////

Fh search_triangle_brute_force(const TriMesh &mesh, const Vec2 &u)
{
    for (auto face : mesh.faces()) if (is_inside(mesh, face, u)) return face;
    return Fh {};
}

Fh search_triangle_local_way(const TriMesh &mesh, const Vec2 &u, const Fh &fho)
{
    const int nf = (int)mesh.n_faces();
    Fh fh = fho; Hh hh {};

    for (int iter = 0; iter < nf; ++iter)
    {
        //if (is_inside(mesh, fh, u)) break;
        Fh fhn = fh; // next face handle

        for (auto hdge : mesh.fh_range(fh)) if (hdge != hh)
        if (is_intersecting(mesh, hdge, centroid(mesh, fh), u))
        { fhn = hdge.opp().face(); hh = hdge.opp(); break; }

        if (!fhn.is_valid()) break; // run into boundary
        if (fhn == fh) assert(is_inside(mesh, fh, u));
        if (fhn == fh) break; // inside triangle
        fh = fhn; // go to next triangle
    }

    return fh;
}

Fh search_triangle_guided_bfs(const TriMesh &mesh, const Vec2 &u, const Fh &fho)
{
    std::queue<Fh> frontier {};
    std::unordered_set<Fh> visited {};

    frontier.push(fho);
    visited.insert(fho);

    while (!frontier.empty())
    {
        auto ft = frontier.front(); frontier.pop();
        auto fh = search_triangle_local_way(mesh, u, ft);
        if (is_inside(mesh, fh, u)) return fh;

        // Do not start with a visited face, with previous one instead
        if (visited.count(fh) && mesh.is_boundary(ft)) fh = ft;

        assert(mesh.is_boundary(fh));
        visited.insert(fh);

        auto hdge = make_smart(fh, mesh).halfedge();
        if (!hdge.opp().is_boundary()) hdge = hdge.next();
        if (!hdge.opp().is_boundary()) hdge = hdge.next();
        hdge = hdge.opp(); // the boundary halfedge of this face

        Fh fbs[2] { hdge.next().opp().face(), hdge.prev().opp().face() };
        for (auto fb : fbs) if (!visited.count(fb))
        { frontier.push(fb); visited.insert(fb); }
    }

    return Fh {};
}

////////////////////////////////////////////////////////////////
/// Fuzzy search
////////////////////////////////////////////////////////////////

Fh fuzzy_search_triangle_brute_force(const TriMesh &mesh, const Vec2 &u, const double tol)
{
    Fh fh {};
    double mind { tol };

    for (auto face : mesh.faces())
    {
        const double d = hausdorff_distance(mesh, face, u);
        if (mind > d) { mind = d; fh = face; }
    }

    return fh;
}

Fh fuzzy_search_triangle_boundary(const TriMesh &mesh, const Vec2 &u, const double tol)
{
    Hh hh {};
    double mind { tol };

    for (auto hdge : mesh.halfedges()) if (hdge.is_boundary())
    {
        const double d = hausdorff_distance(mesh, hdge, u);
        if (mind > d) { mind = d; hh = hdge; }
    }

    return mesh.opposite_face_handle(hh);
}

Fh search_triangle_brute_force(const TriMesh&, const Vec2&);

Fh search_triangle_local_way(const TriMesh&, const Vec2&, const Fh&);

Fh search_triangle_guided_bfs(const TriMesh&, const Vec2&, const Fh&);

Fh fuzzy_search_triangle_boundary(const TriMesh&, const Vec2&, const double tol);

Fh fuzzy_search_triangle_boundary(const TriMesh&, const Vec2&, const double tol);
