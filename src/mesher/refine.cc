#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "delaunay.hh"
#include "triangle.hh"
#include "segment.hh"
#include "mesh.hh"
#include "topology.hh"
#include "search.hh"

using namespace OpenMesh;

////////////////////////////////////////////////////////////////
/// Locations
////////////////////////////////////////////////////////////////

static inline TRI_LOC locate(TriMesh &mesh, const Fh &fh, const Vec2 &u, Hh &hh)
{
    constexpr double kEps = 1e-3;
    const auto hh0 = mesh.halfedge_handle(fh);
    const auto hh1 = mesh.next_halfedge_handle(hh0);
    const auto hh2 = mesh.next_halfedge_handle(hh1);
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

////////////////////////////////////////////////////////////////
/// Delaunay
////////////////////////////////////////////////////////////////

//   1  
//  / \ 
// 2---0
//  \ / 
//   3  
static inline bool is_delaunay(const TriMesh &mesh, const Eh &eh)
{
    Hh hh0 = mesh.halfedge_handle(eh, 0);
    Hh hh1 = mesh.halfedge_handle(eh, 1);
    const auto u0 = get_xy(mesh, hh0);
    const auto u1 = get_xy(mesh, mesh.next_halfedge_handle(hh0));
    const auto u2 = get_xy(mesh, hh1);
    const auto u3 = get_xy(mesh, mesh.next_halfedge_handle(hh1));
    return is_delaunay(u0, u1, u2, u3);
}

struct EuclideanDelaunay
{
    inline bool operator()(const TriMesh &mesh, const Eh &eh) const
    { return is_sharp(mesh, eh) || is_delaunay(mesh, eh); } // If true, do not flip
};

static int make_delaunay(TriMesh &mesh, Eh eh)
{
    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    const int max_n_flip = (int)mesh.n_edges();

    const Eh ehs[4] {
        mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh, 0))),
        mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(eh, 0))),
        mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(eh, 1))),
        mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(eh, 1)))
    };

    delaunifier.reset(); delaunifier.enqueue(ehs, 4); int n_flip {};

    for (auto eh = delaunifier.flip(); eh.is_valid() && n_flip < max_n_flip; eh = delaunifier.flip(), ++n_flip) {}

    return n_flip;
}

////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////

static inline void split_edge(TriMesh &mesh, Hh hh, Vh vh)
{
    Eh eh  = mesh.edge_handle(hh);
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);

    mesh.split_edge_copy(eh, vh);

    for (auto hdge : mesh.voh_range(vh))
    if (hdge.to() != vh0)
    if (hdge.to() != vh1)
    { set_sharp(mesh, hdge.edge(), false); }
}

static inline void split_face(TriMesh &mesh, Fh fh, Vh vh)
{
    mesh.split_copy(fh, vh);
}

static inline Fh search_triangle(const TriMesh &mesh, const Vec2 &u, Fh fh = Fh {})
{
    if (!fh.is_valid()) for (Fh fi : mesh.faces()) { fh = fi; break; }

    fh = search_triangle_straight_way(mesh, u, fh);

    return fh;
}

static inline Vec2 circumcenter(const TriMesh &mesh, const Fh &fh)
{
    Hh hh0 = mesh.halfedge_handle(fh);
    Hh hh1 = mesh.next_halfedge_handle(hh0);
    Hh hh2 = mesh.next_halfedge_handle(hh1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh2));
    return circumcenter(u0, u1, u2);
}

////////////////////////////////////////////////////////////////
/// Encroachment
////////////////////////////////////////////////////////////////

static inline double apex_squared_cosine(const double min_angle)
{
    const double cs0 = cos(min_angle); // cos(t)
    const double cs1 = cs0*cs0*2. - 1; // cos(pi - 2t)
    return cs1*cs1;
}

static inline bool is_exterior(const TriMesh &mesh, const Hh &hh)
{
    return mesh.is_boundary(hh) || is_hidden(mesh, mesh.face_handle(hh));
}

static inline bool is_encroached(const TriMesh &mesh, const Hh &hh, const double cs2)
{
    const auto d0 = -get_dxy(mesh, mesh.next_halfedge_handle(hh));
    const auto d1 =  get_dxy(mesh, mesh.prev_halfedge_handle(hh));
    const double d00 = dot(d0, d0); // a^2
    const double d11 = dot(d1, d1); // b^2
    const double d01 = dot(d0, d1); // a.b
    return (d01 < 0.) && (d01*d01 >= d00*d11*cs2); // a.b < 0 and |a.b| > (a*b)|cos(t0)|
}

struct Encroachment
{
    inline bool operator()(const TriMesh &mesh, const Hh &hh) const
    { return !is_exterior(mesh, hh) && is_encroached(mesh, hh, cs2); } // If true, split it

    const double cs2; // (cos(pi - 2t))^2
};

////////////////////////////////////////////////////////////////
/// Refinement
////////////////////////////////////////////////////////////////

struct Primitive { Hh hh; Vh vh0, vh1, vh2; }; // either a triangle or a segment

static inline bool is_triangle(const Primitive &primitive)
{
    return
        primitive.vh0.is_valid() &&
        primitive.vh1.is_valid() &&
        primitive.vh2.is_valid();
}

static inline bool is_segment(const Primitive &primitive)
{
    return
        primitive.vh0.is_valid() &&
        primitive.vh1.is_valid() &&
       !primitive.vh2.is_valid();
}

static inline bool is_triangle_valid(const TriMesh &mesh, const Primitive &primitive)
{
    if (mesh.is_boundary(primitive.hh)) return false;

    Hh hh = primitive.hh;
    if (mesh.to_vertex_handle(hh) != primitive.vh0)
        hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != primitive.vh0)
        hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != primitive.vh0)
        return false;

    hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != primitive.vh1)
        return false;

    hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != primitive.vh2)
        return false;

    return true;
}

static inline bool is_segment_valid(const TriMesh &mesh, const Primitive &primitive)
{
    Hh hh = primitive.hh;
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);

    return
        vh0 == primitive.vh0 && vh1 == primitive.vh1 ||
        vh0 == primitive.vh1 && vh1 == primitive.vh0;
}

static inline Primitive make_segment(const TriMesh &mesh, Hh hh)
{
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);
    return { hh, vh0, vh1, Vh {} };
}

static int split_segments(TriMesh &mesh, const Encroachment &encroached)
{
    const int max_n_iter = 100;

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    std::queue<Primitive> primitives; // to split

    for (Eh eh : mesh.edges()) if (is_sharp(mesh, eh))
    {
        for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
        {
            primitives.push(make_segment(mesh, hh));
        }
    }

    while (!primitives.empty())
    {
        auto primitive = primitives.front(); primitives.pop();

        if (is_segment(primitive) && is_segment_valid(mesh, primitive))
        {
            // find a position for splitting
            // [TODO] handle acute corner
            const auto u = mesh.calc_edge_midpoint(primitive.hh);

            // split the segment
            Vh vh = mesh.new_vertex(u);
            split_edge(mesh, primitive.hh, vh);

            // edges to flip
            Eh ehs[4]; int ne {};
            for (auto hdge : mesh.voh_range(vh))
                if (!hdge.next().edge().is_boundary())
                    ehs[ne++] = hdge.next().edge();

            // maintain Delaunayhood
            delaunifier.reset(); delaunifier.enqueue(ehs, ne); int n_iter {};

            // check encroachment of encountered segment while testing Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid() && n_iter < max_n_iter; eh = delaunifier.next(), ++n_iter)
            if (is_sharp(mesh, eh)) for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { primitives.push(make_segment(mesh, hh)); }

            // check encroachment of two new segments afterwards
            for (Eh eh : mesh.ve_range(vh)) if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { primitives.push(make_segment(mesh, hh)); }
        }
    }

    return 0;
}

static int split_interior(TriMesh &mesh, const Encroachment &encroached)
{
    return 0;
}

////////////////////////////////////////////////////////////////
/// Wrapping up
////////////////////////////////////////////////////////////////

int refine(TriMesh &mesh, const double min_angle)
{
    const double cs2 = apex_squared_cosine(min_angle);

    Encroachment encroached { cs2 };

    int err;

    err = split_segments(mesh, encroached);

    return err;
}