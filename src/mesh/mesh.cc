#include <unordered_map>
#include <unordered_set>
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
static inline std::string get_path(const MeshT &mesh)
{
    std::string filename {};

    if (hasProperty<Mh, std::string>(mesh, var_m_path()))
    {
        filename.append(getProperty<Mh, std::string>(mesh, var_m_path())());
        filename.append("/");
    }

    return filename;
}

int dump_mesh(const TriMesh &mesh, const char *local_name)
{
    return save_mesh(mesh, get_path(mesh).append(local_name).c_str());
}

int dump_mesh(const PolyMesh &mesh, const char *local_name)
{
    return save_mesh(mesh, get_path(mesh).append(local_name).c_str());
}

////////////////////////////////////////////////////////////////
/// Dump primitives
////////////////////////////////////////////////////////////////

int dump_faces(const TriMesh &mesh, const std::vector<Fh> &faces, const char *local_name)
{
    std::vector<Vec2> vs {};
    std::vector<Int3> fs {};

    int nv {};
    std::unordered_set<Fh> fset {};
    std::unordered_map<Vh, int> vids {};

    for (Fh fh : faces) { if (fh.is_valid()) if (!is_deleted(mesh, fh)) { fset.insert(fh);
    for (Vh vh : mesh.fv_range(fh)) if (!vids.count(vh)) { vids[vh] = nv++; } } }

    vs.resize(nv);

    for (const auto &vid : vids)
    {
        Vh vh = vid.first;
        const int id = vid.second;
        vs[id] = get_xy(mesh, vh);
    }

    for (Fh fh : fset)
    {
        Vh vh0 = mesh.from_vertex_handle(mesh.halfedge_handle(fh));
        Vh vh1 = mesh.to_vertex_handle  (mesh.halfedge_handle(fh));
        Vh vh2 = mesh.to_vertex_handle  (mesh.next_halfedge_handle(mesh.halfedge_handle(fh)));
        const int i0 = vids[vh0];
        const int i1 = vids[vh1];
        const int i2 = vids[vh2];
        fs.emplace_back(i0, i1, i2);
    }

    return save_mesh(
        (const double *)vs.data(), (int)vs.size(),
        (const int    *)fs.data(), (int)fs.size(),
        nullptr, 0,
        get_path(mesh).append(local_name).c_str());
}
