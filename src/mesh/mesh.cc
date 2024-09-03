#include "mesh.hh"
#include "mesh.io.hh"

using namespace OpenMesh;

const char *var_v_index() { return "vert:index"; }
const char *var_f_index() { return "face:index"; }
const char *var_e_index() { return "edge:index"; }
const char *var_h_index() { return "hdge:index"; }

const char *var_v_label() { return "vert:label"; }
const char *var_f_label() { return "face:label"; }
const char *var_e_label() { return "edge:label"; }
const char *var_h_label() { return "hdge:label"; }

const char *var_m_name() { return "mesh:name"; }
const char *var_m_path() { return "mesh:path"; }

////////////////////////////////////////////////////////////////
/// Dump handle
////////////////////////////////////////////////////////////////

int dump_handle(const ArrayKernel &mesh, const Vh &vh, const int offset)
{
    return vh.idx() + offset;
}

std::tuple<int, int> dump_handle(const ArrayKernel &mesh, const Hh &hh, const int offset)
{
    return {
        mesh.from_vertex_handle(hh).idx() + offset,
        mesh.to_vertex_handle  (hh).idx() + offset
    };
}

std::tuple<int, int> dump_handle(const ArrayKernel &mesh, const Eh &eh, const int offset)
{
    return {
        mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0)).idx() + offset,
        mesh.to_vertex_handle  (mesh.halfedge_handle(eh, 0)).idx() + offset
    };
}

std::tuple<int, int, int> dump_handle(const TriMesh &mesh, const Fh &fh, const int offset)
{
    return {
        mesh.from_vertex_handle(mesh.halfedge_handle(fh)).idx() + offset,
        mesh.to_vertex_handle  (mesh.halfedge_handle(fh)).idx() + offset,
        mesh.to_vertex_handle  (mesh.next_halfedge_handle(mesh.halfedge_handle(fh))).idx() + offset
    };
}

std::vector<int> dump_handle(const PolyMesh &mesh, const Fh &fh, const int offset)
{
    std::vector<int> vs {};
    for (Hh hh : mesh.fh_range(fh))
    { vs.push_back(mesh.to_vertex_handle(hh).idx() + offset); }
    return vs;
}

////////////////////////////////////////////////////////////////
/// Dump mesh
////////////////////////////////////////////////////////////////

template <class MeshT>
int dump_mesh_local_path(const MeshT &mesh, const char *local_name)
{
    std::string filename {};

    // Save mesh data under <mesh path>, or current path if path is empty.
    if (hasProperty<Mh, std::string>(mesh, var_m_path()))
    {
        filename.append(getProperty<Mh, std::string>(mesh, var_m_path())());
        filename.append("/");
    }

    filename.append(local_name);

    return save_mesh(mesh, filename.c_str());
}

int dump_mesh(const TriMesh &mesh, const char *local_name)
{
    return dump_mesh_local_path(mesh, local_name);
}

int dump_mesh(const PolyMesh &mesh, const char *local_name)
{
    return dump_mesh_local_path(mesh, local_name);
}

////////////////////////////////////////////////////////////////
/// Dump primitive
////////////////////////////////////////////////////////////////

// [TODO]
