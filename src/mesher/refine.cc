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
/// Circumcenter
////////////////////////////////////////////////////////////////

static inline Vec2 circumcenter(const TriMesh &mesh, Fh fh)
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

    delaunifier.reset(); delaunifier.to_flip(ehs, 4); int n_flip {};

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
    if (!fh.is_valid()) fh = mesh.face_handle(0);

    fh = search_triangle_guided_bfs(mesh, u, fh);

    if (!fh.is_valid()) fh = search_triangle_brute_force(mesh, u);

    return fh;
}

////////////////////////////////////////////////////////////////
/// Refinement
////////////////////////////////////////////////////////////////

struct Primitive { Fh fh; Eh eh; }; // either a triangle or a segment
