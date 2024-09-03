#ifndef MESH_TOPOLOGY_HH
#define MESH_TOPOLOGY_HH

#include <queue>
#include <unordered_map>
#include <unordered_set>

#include "mesh.hh"

////////////////////////////////////////////////////////////////
/// Delaunay
////////////////////////////////////////////////////////////////

namespace OpenMesh
{

inline bool is_rotate_ok(const PolyConnectivity &mesh, const HalfedgeHandle &hh)
{
    // boundary edges cannot be rotated
    if (mesh.is_boundary(mesh.edge_handle(hh))) return false;

    const auto oh = mesh.opposite_halfedge_handle(hh);

    // check if the rotated edge is already present in the mesh
    const auto vh0 = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
    const auto vh1 = mesh.to_vertex_handle(mesh.next_halfedge_handle(oh));

    if (vh0 == vh1) return false; // this is generally a bad sign !!!

    for (PolyConnectivity::ConstVertexVertexIter vvi(mesh, vh0); vvi.is_valid(); ++vvi)
        if (*vvi == vh1) return false;

    return true;
}

inline void rotate(PolyConnectivity &mesh, const HalfedgeHandle &hh)
{
    // CAUTION : Rotating a halfedge may result in
    // a non-manifold mesh, hence check for yourself
    // whether this operation is allowed or not!
    assert(!mesh.is_boundary(mesh.edge_handle(hh)));
    assert(is_rotate_ok(mesh, hh));

    auto ha0 = hh;
    auto hb0 = mesh.opposite_halfedge_handle(hh);

    auto ha1 = mesh.next_halfedge_handle(ha0);
    auto ha2 = mesh.next_halfedge_handle(ha1);
    auto ha3 = mesh.prev_halfedge_handle(ha0);

    auto hb1 = mesh.next_halfedge_handle(hb0);
    auto hb2 = mesh.next_halfedge_handle(hb1);
    auto hb3 = mesh.prev_halfedge_handle(hb0);

    auto va0 = mesh.to_vertex_handle(ha0);
    auto vb0 = mesh.to_vertex_handle(hb0);

    auto va1 = mesh.to_vertex_handle(ha1);
    auto vb1 = mesh.to_vertex_handle(hb1);

    auto fa = mesh.face_handle(ha0);
    auto fb = mesh.face_handle(hb0);

    mesh.set_vertex_handle(ha0, va1);
    mesh.set_vertex_handle(hb0, vb1);

    mesh.set_next_halfedge_handle(ha0, ha2);
    mesh.set_next_halfedge_handle(hb1, ha0);
    mesh.set_next_halfedge_handle(ha3, hb1);

    mesh.set_next_halfedge_handle(hb0, hb2);
    mesh.set_next_halfedge_handle(ha1, hb0);
    mesh.set_next_halfedge_handle(hb3, ha1);

    mesh.set_face_handle(ha1, fb);
    mesh.set_face_handle(hb1, fa);

    mesh.set_halfedge_handle(fa, ha0);
    mesh.set_halfedge_handle(fb, hb0);

    if (mesh.halfedge_handle(va0) == hb0)
        mesh.set_halfedge_handle(va0, ha1);
    if (mesh.halfedge_handle(vb0) == ha0)
        mesh.set_halfedge_handle(vb0, hb1);
}

inline void flip(PolyConnectivity &mesh, const EdgeHandle &eh)
{
    rotate(mesh, mesh.halfedge_handle(eh, 0));
}

inline void flip(TriConnectivity &mesh, const EdgeHandle &eh)
{
    mesh.flip(eh);
}

} // namespace OpenMesh

template <class MeshT, class DelaunayT>
class Delaunifier
{
public:
    Delaunifier(MeshT &mesh, const DelaunayT &is_delaunay)
    : mesh(mesh), is_delaunay(is_delaunay) {}

    void reset();

    void to_flip();
    void to_flip(const Eh[], const int);

    Eh flip();

    size_t n_enqueued() const { return frontier.size(); }

protected:
    MeshT &mesh;
    const DelaunayT &is_delaunay;
    std::deque<Eh> frontier;
    std::unordered_set<Eh> enqueued;
};

template <class MeshT, class DelaunayT>
inline void Delaunifier<MeshT, DelaunayT>::reset()
{
    frontier.clear();
    enqueued.clear();
}

template <class MeshT, class DelaunayT>
inline void Delaunifier<MeshT, DelaunayT>::to_flip()
{
    for (auto edge : mesh.edges())
    if (!edge.is_boundary())
    {
        frontier.push_back(edge);
        enqueued.insert(edge);
    }
}

template <class MeshT, class DelaunayT>
inline void Delaunifier<MeshT, DelaunayT>::to_flip(const Eh ehs[], const int ne)
{
    for (int i = 0; i < ne; ++i)
    if (!mesh.is_boundary(ehs[i]))
    {
        frontier.push_back(ehs[i]);
        enqueued.insert(ehs[i]);
    }
}

template <class MeshT, class DelaunayT>
inline Eh Delaunifier<MeshT, DelaunayT>::flip()
{
    while (!frontier.empty())
    {
        auto ehf = frontier.front(); frontier.pop_front();
        enqueued.erase(ehf);

        if (is_delaunay(mesh, ehf)) continue;

        OpenMesh::flip(mesh, ehf);

        const Eh ehs[4] {
            mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(ehf, 0))),
            mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(ehf, 0))),
            mesh.edge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(ehf, 1))),
            mesh.edge_handle(mesh.prev_halfedge_handle(mesh.halfedge_handle(ehf, 1)))
        };

        for (const auto &eh : ehs)
        if (!mesh.is_boundary(eh))
        if (!enqueued.count(eh))
        {
            frontier.push_back(eh);
            enqueued.insert(eh);
        }

        return ehf;
    }

    return mesh.InvalidEdgeHandle;
}

template <class MeshT, class DelaunayT>
inline Delaunifier<MeshT, DelaunayT> make_delaunifier(MeshT &mesh, const DelaunayT &is_delaunay)
{ return Delaunifier<MeshT, DelaunayT>(mesh, is_delaunay); }

////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////

template <class MeshT>
inline int get_boundaries(const MeshT&, std::vector<Hh>&);

template <class MeshT>
inline void sort_boundaries(const MeshT&, std::vector<Hh>&);

#endif