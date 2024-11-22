#ifndef MESH_DEFINITION_HH
#define MESH_DEFINITION_HH

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Casts.hh>
#include "property.hh"
#include "vector_n.hh"

using Hh = OpenMesh::HalfedgeHandle;
using Vh = OpenMesh::VertexHandle;
using Fh = OpenMesh::FaceHandle;
using Eh = OpenMesh::EdgeHandle;
using Mh = OpenMesh::MeshHandle;

////////////////////////////////////////////////////////////////
/// Definition
////////////////////////////////////////////////////////////////

struct MeshTraits : public OpenMesh::DefaultTraitsDouble
{
    // Default types
    typedef double TexCoord1D;
    typedef Vec2   TexCoord2D;
    typedef Vec3   TexCoord3D;

    // Default attributes
    VertexAttributes   (OpenMesh::Attributes::Status);
    FaceAttributes     (OpenMesh::Attributes::Status);
    EdgeAttributes     (OpenMesh::Attributes::Status);
    HalfedgeAttributes (OpenMesh::Attributes::Status);

    // Customized attributes
    VertexTraits   {};
    FaceTraits     {};
    EdgeTraits     {};
    HalfedgeTraits {};
};

using PolyMesh = OpenMesh::PolyMesh_ArrayKernelT<MeshTraits>;

using TriMesh = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;

////////////////////////////////////////////////////////////////
/// Topology
////////////////////////////////////////////////////////////////

template <class MeshT>
inline bool is_sync(const MeshT &mesh, const Hh &hh)
{
    return mesh.halfedge_handle(mesh.edge_handle(hh), 0) == hh;
}

////////////////////////////////////////////////////////////////
/// Navigate
////////////////////////////////////////////////////////////////

template <class MeshT, class PredicateT>
inline Hh ccw_rotated(const MeshT &mesh, const PredicateT &predicate, const Hh &hh)
{
    Hh hi = hh;
    while ((hi = mesh.ccw_rotated_halfedge_handle(hi)) != hh)
        if (predicate(mesh, hi)) break;
    return hi; // return either next CCW sharp halfedge or itself
}

template <class MeshT, class PredicateT>
inline Hh cw_rotated(const MeshT &mesh, const PredicateT &predicate, const Hh &hh)
{
    Hh hi = hh;
    while ((hi = mesh.cw_rotated_halfedge_handle(hi)) != hh)
        if (predicate(mesh, hi)) break;
    return hi; // return either next CW sharp halfedge or itself
}

template <class MeshT, class PredicateT>
inline Hh next(const MeshT &mesh, const PredicateT &predicate, const Hh &hh)
{
    return cw_rotated(mesh, predicate, mesh.opposite_halfedge_handle(hh));
}

template <class MeshT, class PredicateT>
inline Hh prev(const MeshT &mesh, const PredicateT &predicate, const Hh &hh)
{
    return mesh.opposite_halfedge_handle(ccw_rotated(mesh, predicate, hh));
}

////////////////////////////////////////////////////////////////
/// Flags
////////////////////////////////////////////////////////////////

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Fh &fh) { return mesh.status(fh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Fh &fh, const bool val) { mesh.status(fh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_selected(val); }

template <class MeshT>
inline bool is_marked(const MeshT &mesh, const Hh &hh) { return mesh.status(hh).selected(); }
template <class MeshT>
inline void set_marked(MeshT &mesh, const Hh &hh, const bool val) { mesh.status(hh).set_selected(val); }

template <class MeshT>
inline bool is_sharp(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).feature(); }
template <class MeshT>
inline void set_sharp(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_feature(val); }

template <class MeshT>
inline bool is_sharp(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).feature(); }
template <class MeshT>
inline void set_sharp(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_feature(val); }

template <class MeshT>
inline bool is_fixed(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).locked(); }
template <class MeshT>
inline void set_fixed(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_locked(val); }

template <class MeshT>
inline bool is_fixed(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).locked(); }
template <class MeshT>
inline void set_fixed(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_locked(val); }

template <class MeshT>
inline bool is_hidden(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).hidden(); }
template <class MeshT>
inline void set_hidden(MeshT &mesh, const Vh &vh, const bool val) { mesh.status(vh).set_hidden(val); }

template <class MeshT>
inline bool is_hidden(const MeshT &mesh, const Fh &fh) { return mesh.status(fh).hidden(); }
template <class MeshT>
inline void set_hidden(MeshT &mesh, const Fh &fh, const bool val) { mesh.status(fh).set_hidden(val); }

template <class MeshT>
inline bool is_hidden(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).hidden(); }
template <class MeshT>
inline void set_hidden(MeshT &mesh, const Eh &eh, const bool val) { mesh.status(eh).set_hidden(val); }

template <class MeshT>
inline bool is_deleted(const MeshT &mesh, const Vh &vh) { return mesh.status(vh).deleted(); }

template <class MeshT>
inline bool is_deleted(const MeshT &mesh, const Fh &fh) { return mesh.status(fh).deleted(); }

template <class MeshT>
inline bool is_deleted(const MeshT &mesh, const Eh &eh) { return mesh.status(eh).deleted(); }

template <class MeshT>
inline bool is_deleted(const MeshT &mesh, const Hh &hh) { return mesh.status(hh).deleted(); }

////////////////////////////////////////////////////////////////
/// Geometry
////////////////////////////////////////////////////////////////

template <class MeshT>
inline Vec3 get_xyz(const MeshT &mesh, const Vh &vh)
{ return mesh.point(vh); }

template <class MeshT>
inline void set_xyz(MeshT &mesh, const Vh &vh, const Vec3 &p)
{ mesh.set_point(vh, p); }

template <class MeshT>
inline Vec2 get_xy(const MeshT &mesh, const Vh &vh)
{ const auto p = mesh.point(vh); return { p[0], p[1] }; }

template <class MeshT>
inline void set_xy(MeshT &mesh, const Vh &vh, const Vec2 &u)
{ const auto p = mesh.point(vh); mesh.set_point(vh, { u[0], u[1], p[2] }); }

template <class MeshT>
inline Vec2 get_xy(const MeshT &mesh, const Hh &hh)
{ return get_xy(mesh, mesh.to_vertex_handle(hh)); }

template <class MeshT>
inline void set_xy(MeshT &mesh, const Hh &hh, const Vec2 &u)
{ set_xy(mesh, mesh.to_vertex_handle(hh)); }

template <class MeshT>
inline Vec2 get_dxy(const MeshT &mesh, const Hh &hh)
{ return get_xy(mesh, hh) - get_xy(mesh, mesh.opposite_halfedge_handle(hh)); }

template <class MeshT>
inline Vec3 get_dxyz(const MeshT &mesh, const Hh &hh)
{ return get_xyz(mesh, hh) - get_xyz(mesh, mesh.opposite_halfedge_handle(hh)); }

template <class MeshT>
inline Vec2 get_uv(const MeshT &mesh, const Vh &vh)
{ return mesh.texcoord2D(vh); }

template <class MeshT>
inline void set_uv(MeshT &mesh, const Vh &vh, const Vec2 &u)
{ mesh.set_texcoord2D(vh, u); }

template <class MeshT>
inline Vec2 get_uv(const MeshT &mesh, const Hh &hh)
{ return mesh.texcoord2D(mesh.to_vertex_handle(hh)); }

template <class MeshT>
inline void set_uv(MeshT &mesh, const Hh &hh, const Vec2 &u)
{ mesh.set_texcoord2D(mesh.to_vertex_handle(hh), u); }

template <class MeshT>
inline double get_param1D(const MeshT &mesh, const Hh &hh)
{ return mesh.texcoord1D(hh); }

template <class MeshT>
inline void set_param1D(MeshT &mesh, const Hh &hh, const double t)
{ mesh.set_texcoord1D(hh, t); }

template <class MeshT>
inline Vec2 get_param2D(const MeshT &mesh, const Hh &hh)
{ return mesh.texcoord2D(hh); }

template <class MeshT>
inline void set_param2D(MeshT &mesh, const Hh &hh, const Vec2 &u)
{ mesh.set_texcoord2D(hh, u); }

template <class MeshT>
inline Vec3 get_param3D(const MeshT &mesh, const Hh &hh)
{ return mesh.texcoord3D(hh); }

template <class MeshT>
inline void set_param3D(MeshT &mesh, const Hh &hh, const Vec3 &m)
{ mesh.set_texcoord3D(hh, m); }

////////////////////////////////////////////////////////////////
/// Info
////////////////////////////////////////////////////////////////

template <class MeshT>
inline void print_handle(const MeshT&, const Vh &vh, const int offset = 0, const char *end = "")
{
    printf("(%d)%s", vh.idx() + offset, end);
}

template <class MeshT>
inline void print_handle(const MeshT &mesh, const Hh &hh, const int offset = 0, const char *end = "")
{
    printf("(%d, %d)%s",
        mesh.from_vertex_handle(hh).idx() + offset,
        mesh.to_vertex_handle  (hh).idx() + offset, end);
}

template <class MeshT>
inline void print_handle(const MeshT &mesh, const Eh &eh, const int offset = 0, const char *end = "")
{
    printf("(%d, %d)%s",
        mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0)).idx() + offset,
        mesh.to_vertex_handle  (mesh.halfedge_handle(eh, 0)).idx() + offset, end);
}

inline void print_handle(const TriMesh &mesh, const Fh &fh, const int offset = 0, const char *end = "")
{
    printf("(%d, %d, %d)%s",
        mesh.from_vertex_handle(mesh.halfedge_handle(fh)).idx() + offset,
        mesh.to_vertex_handle  (mesh.halfedge_handle(fh)).idx() + offset,
        mesh.to_vertex_handle  (mesh.next_halfedge_handle(mesh.halfedge_handle(fh))).idx() + offset, end);
}

inline void print_handle(const PolyMesh &mesh, const Fh &fh, const int offset = 0, const char *end = "")
{
    int nh {}; for (Hh hh : mesh.fh_range(fh))
    { printf("%s%d", ((nh++) ? ", " : "("), mesh.to_vertex_handle(hh).idx() + offset); }
    printf(")%s", end);
}

////////////////////////////////////////////////////////////////
/// Constants
////////////////////////////////////////////////////////////////

const char *var_v_index();
const char *var_f_index();
const char *var_e_index();
const char *var_h_index();

const char *var_v_label();
const char *var_e_label();
const char *var_h_label();
const char *var_f_label();

const char *var_m_name();
const char *var_m_path();

#endif