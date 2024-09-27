#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "pred2D.hh"
#include "distance.hh"
#include "triangle.hh"
#include "segment.hh"
#include "mesh.hh"
#include "search.hh"

using namespace OpenMesh;

////////////////////////////////////////////////////////////////
/// Exact search
////////////////////////////////////////////////////////////////

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
        Fh fi = fh; // next face handle

        for (auto hdge : mesh.fh_range(fh)) if (hdge != hh)
        if (is_intersecting(mesh, hdge, centroid(mesh, fh), u))
        { fi = hdge.opp().face(); hh = hdge.opp(); break; }

        if (!fi.is_valid()) break; // run into boundary
        if (fi == fh) assert(is_inside(mesh, fh, u));
        if (fi == fh) break; // inside triangle
        fh = fi; // go to the next triangle
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

////////////////////////////////////////////////////////////////
/// Path search
////////////////////////////////////////////////////////////////

static inline PLOW_STATUS intersection_info(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, const Hh &hh)
{
    const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hh));

    const auto ii = intersection_info(u0, u1, v0, v1);
    const int ru0 = ii[0];
    const int ru1 = ii[1];
    const int rv0 = ii[2];
    const int rv1 = ii[3];

    return
        (ru0 * ru1 < 0) && (rv0 * rv1 < 0)        ? PLOW_STATUS::EDGE : // intersecting exclusively
        (ru0 * ru1 < 0) && (rv1 == 0 && rv0 != 0) ? PLOW_STATUS::VERT : // v1 on (u0,u1), v0 is not
        (ru0 != 0 && ru1 == 0) && (rv1 == 0)      ? PLOW_STATUS::VERT : // v1 overlaps u1
        PLOW_STATUS::MISS; // no intersection, colinear, v0 lying on (u0,u1), and other cases
}

static inline PLOW_STATUS next_primitive(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Hh hho, Hh &hhc, Vh &vhc)
{
    auto res { PLOW_STATUS::MISS };

    Hh hh0 = mesh.prev_halfedge_handle(hho);
    Hh hh1 = mesh.next_halfedge_handle(hho);

    for (Hh hh : { hh0, hh1 })
    {
        const auto ii = intersection_info(mesh, u0, u1, hh);
        if (ii == PLOW_STATUS::MISS) continue;
        res = ii;

        if (ii == PLOW_STATUS::EDGE)
        {
            hhc = hh;
        }
        else if (ii == PLOW_STATUS::VERT)
        {
            vhc = mesh.to_vertex_handle(hh);
        }
    }

    return res;
}

static inline PLOW_STATUS next_primitive(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Vh vho, Hh &hhc, Vh &vhc)
{
    auto res { PLOW_STATUS::MISS };

    const auto uo = get_xy(mesh, vho);

    for (Hh hh : mesh.voh_range(vho))
    {
        if (!mesh.is_boundary(hh)) // the apex halfedge to v0
        {
            hh = mesh.next_halfedge_handle(hh);
        }
        else // an extra halfedge to the vertex adjacent to v0 but not tested
        {
            hh = mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(mesh.ccw_rotated_halfedge_handle(hh)));
        }

        // use (uo,u1) for intersecting test instead of (u0,u1)
        const auto ii = intersection_info(mesh, uo, u1, hh);
        if (ii == PLOW_STATUS::MISS) continue;
        res = ii;

        if (ii == PLOW_STATUS::EDGE)
        {
            hhc = hh;
        }
        else if (ii == PLOW_STATUS::VERT)
        {
            vhc = mesh.to_vertex_handle(hh);
        }
    }

    return res;
}

void PrimitivePlow::next()
{
    if (st_ == PLOW_STATUS::EDGE)
    {
        hh_ = m_.opposite_halfedge_handle(hh_);
    }
    if (st_ == PLOW_STATUS::EDGE)
    {
        st_ = next_primitive(m_, u0_, u1_, hh_, hh_, vh_);
    }
    else if (st_ == PLOW_STATUS::VERT)
    {
        st_ = next_primitive(m_, u0_, u1_, vh_, hh_, vh_);
    }
}

//////////////// Linear path search initialization ////////////////

static inline TRI_LOC locate(const TriMesh &mesh, const Fh &fh, const Vec2 &u, Hh &hh)
{
    Hh hh0 = mesh.halfedge_handle(fh);
    Hh hh1 = mesh.next_halfedge_handle(hh0);
    Hh hh2 = mesh.next_halfedge_handle(hh1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh2));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const auto loc = exact_locate(u0, u1, u2, u);
    if (loc == TRI_LOC::E0) { hh = hh0; }
    if (loc == TRI_LOC::E1) { hh = hh1; }
    if (loc == TRI_LOC::E2) { hh = hh2; }
    if (loc == TRI_LOC::V0) { hh = hh1; }
    if (loc == TRI_LOC::V1) { hh = hh2; }
    if (loc == TRI_LOC::V2) { hh = hh0; }
    return loc;
}

static inline PLOW_STATUS first_primitive(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, const Fh &fho, Hh &hhc, Vh &vhc)
{
    Hh hho {}; const auto loc = locate(mesh, fho, u0, hho);

    if (loc == TRI_LOC::IN) // u0 is exclusively in the triangle
    {
        for (Hh hh : mesh.fh_range(fho))
        if (intersection_info(mesh, u0, u1, hh) == PLOW_STATUS::MISS)
        { hhc = hh; break; } // then start with any non-intersecting edge
        return PLOW_STATUS::EDGE;
    }

    if ((int)loc & (int)TRI_LOC::VS) // u0 overlaps any vertex of the triangle
    {
        vhc = mesh.to_vertex_handle(hho);
        return PLOW_STATUS::VERT;
    }

    if ((int)loc & (int)TRI_LOC::ES) // u0 lies on any edge of the triangle
    {
        const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hho));
        const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hho));
        const int r = orientation(v0, v1, u1);

        if (r > 0) // u1 is to the left of (v0,v1)
        {
            hhc = hho;
        }
        else if (r < 0) // u1 is to the right of (v0,v1)
        {
            hhc = mesh.opposite_halfedge_handle(hho);
        }
        else if (dot(u1-u0, v1-v0) < 0) // (u0,u1) is anti-linear with (v0,v1)
        {
            hhc = hho;
        }
        else // if (dot(u1-u0, v1-v0) > 0) // (u0,u1) is colinear with (v0,v1)
        {
            hhc = mesh.opposite_halfedge_handle(hho); // even work in boundary case!
        }

        return PLOW_STATUS::EDGE;
    }

    return PLOW_STATUS::MISS;
}

void init(PrimitivePlow &pp, const Fh &fh0, const Vec2 &u0, const Vec2 &u1)
{
    Hh hhc {}; Vh vhc {};
    const auto &mesh = pp.mesh();

    auto st = first_primitive(mesh, u0, u1, fh0, hhc, vhc);
    if (st == PLOW_STATUS::EDGE) { hhc = mesh.opposite_halfedge_handle(hhc); }

    pp.set_halfedge_handle(hhc);
    pp.set_vertex_handle(vhc);
    pp.set_status(st);
    pp.set_u0(u0);
    pp.set_u1(u1);
}

void init(PrimitivePlow &pp, const Fh &fh0, const Vec2 &u1)
{
    const auto &mesh = pp.mesh();
    const auto u0 = centroid(mesh, fh0);

    for (Hh hh : mesh.fh_range(fh0))
    if (intersection_info(mesh, u0, u1, hh) == PLOW_STATUS::MISS)
    { pp.set_halfedge_handle(mesh.opposite_halfedge_handle(hh)); break; }

    pp.set_status(PLOW_STATUS::EDGE);
    pp.set_u0(u0);
    pp.set_u1(u1);
}

void init(PrimitivePlow &pp, const Vh &vh0, const Vh &vh1)
{
    const auto &mesh = pp.mesh();
    const auto u0 = get_xy(mesh, vh0);
    const auto u1 = get_xy(mesh, vh1);
    pp.set_status(PLOW_STATUS::VERT);
    pp.set_vertex_handle(vh0);
    pp.set_u0(u0);
    pp.set_u1(u1);
}
