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

    delaunifier.reset(); delaunifier.enqueue(ehs, 4);

    int n_flip = delaunifier.flip_all(max_n_flip);

    return n_flip;
}

////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////

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

    return search_triangle_local_way(mesh, u, fh);
}

////////////////////////////////////////////////////////////////
/// Incremental triangulation
////////////////////////////////////////////////////////////////

static inline void add_points(TriMesh &mesh, const std::vector<Vec2> &vs)
{
    assert(mesh.n_vertices() == 0);

    for (const auto &u : vs)
    {
        mesh.new_vertex({ u[0], u[1], 0 });
    }
}

template <class MeshT>
static inline VecN<Vec2, 2> get_range(const MeshT &mesh)
{
    Vec2 bl { +1e20, +1e20 };
    Vec2 ur { -1e20, -1e20 };

    for (Vh vh : mesh.vertices())
    {
        const auto u = get_xy(mesh, vh);
        bl[0] = (u[0] < bl[0]) ? u[0] : bl[0];
        bl[1] = (u[1] < bl[1]) ? u[1] : bl[1];
        ur[0] = (u[0] > ur[0]) ? u[0] : ur[0];
        ur[1] = (u[1] > ur[1]) ? u[1] : ur[1];
    }

    return { bl, ur };
}

static inline void set_domain(TriMesh &mesh)
{
    // Get the bounding box
    const auto blur = get_range(mesh);
    auto bl = blur[0], ur = blur[1];

    // enlarge the box 2 times the original size
    const auto dl = ur - bl;
    const double l = (dl[0]<dl[1]) ? dl[0] : dl[1];
    const Vec2 d { l, l };
    bl -= d;
    ur += d;

    // create a square as the initial mesh
    Vh vh0 = mesh.new_vertex({ bl[0], bl[1], 0 });
    Vh vh1 = mesh.new_vertex({ ur[0], bl[1], 0 });
    Vh vh2 = mesh.new_vertex({ ur[0], ur[1], 0 });
    Vh vh3 = mesh.new_vertex({ bl[0], ur[1], 0 });
    mesh.add_face({ vh0, vh1, vh2 });
    mesh.add_face({ vh2, vh3, vh0 });
}

static int insert_vertices(TriMesh &mesh, std::unordered_map<Vh, Vh> &dups)
{
    const int max_n_flip = 100;

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    int n_new_vertices {};

    Fh fh_last {};

    for (Vh vh : mesh.vertices()) if (mesh.is_isolated(vh))
    {
        const auto u = get_xy(mesh, vh);

        // find the triangle where the point locates
        Fh fh = search_triangle(mesh, u, fh_last);

        // skip any point out of domain
        if (!(fh_last = fh).is_valid()) continue;

        // at which part of the triangle the point locates
        Hh hh {}; const auto loc = locate(mesh, fh, u, hh);

        // skip any vertex with invalid location
        if (loc == TRI_LOC::OUT) continue;

        // skip and record duplicated vertices
        if ((int)loc & (int)TRI_LOC::VS)
        { dups[vh] = mesh.to_vertex_handle(hh); continue; }

        // insert the point into the triangle or onto the edge
        if (loc == TRI_LOC::IN) split_face(mesh, fh, vh);
        else                    split_edge(mesh, hh, vh);

        // edges to flip
        Eh ehs[4]; int ne {};
        for (auto hdge : mesh.voh_range(vh))
            if (!hdge.next().edge().is_boundary())
                ehs[ne++] = hdge.next().edge();

        // maintain Delaunay
        delaunifier.reset(); delaunifier.enqueue(ehs, ne);
        int n_flip = delaunifier.flip_all(max_n_flip);

        ++n_new_vertices;
    }

    return n_new_vertices;
}

////////////////////////////////////////////////////////////////
/// Local search
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
/// Constraints
////////////////////////////////////////////////////////////////

//
//   1   
//  / \  
// 2---0 
//  \ /  
//   3   
//
static inline bool is_flippable(const TriMesh &mesh, const Hh &hh)
{
    Hh hi = mesh.opposite_halfedge_handle(hh);
    const auto u0 = get_xy(mesh, hh);
    const auto u1 = get_xy(mesh, mesh.next_halfedge_handle(hh));
    const auto u2 = get_xy(mesh, hi);
    const auto u3 = get_xy(mesh, mesh.next_halfedge_handle(hi));
    return exact_convex(u0, u1, u2, u3);
}

static inline bool is_intersecting(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, const Hh &hh)
{
    const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto ii = intersection_info(u0, u1, v0, v1);
    const int ru0 = ii[0];
    const int ru1 = ii[1];
    const int rv0 = ii[2];
    const int rv1 = ii[3];
    return (ru0 * ru1 < 0) && (rv0 * rv1 < 0); // exclusively intersecting
}

static inline bool is_intersecting(const TriMesh &mesh, const Vh &vh0, const Vh &vh1, const Hh &hh)
{
    const auto u0 = get_xy(mesh, vh0);
    const auto u1 = get_xy(mesh, vh1);
    return is_intersecting(mesh, u0, u1, hh);
}

static int get_intersections(const TriMesh &mesh, const Vh &vh0, const Vh &vh1, std::vector<Hh> &hhs, std::vector<Vh> &vhs)
{
    PrimitivePlow pp(mesh);

    init(pp, vh0, vh1); // setup the plow

    const int max_n_iter = (int)mesh.n_edges(); int n_iter {};

    for (pp.next(); n_iter < max_n_iter; pp.next(), ++n_iter)
    {
        if (pp.status() == PLOW_STATUS::VERT) // check if v1 is reached
        {
            if (pp.vertex_handle() == vh1) return 0;
        }
        if (pp.status() == PLOW_STATUS::EDGE) // record intersecting edges
        {
            hhs.push_back(pp.halfedge_handle());
        }
        else if (pp.status() == PLOW_STATUS::VERT) // record overlapping vertices
        {
            vhs.push_back(pp.vertex_handle()); // no vertex can lie on (u0,u1) other than v0 and v1
        }
        else if (pp.status() == PLOW_STATUS::MISS) // searching lost in vain, for some reasons
        {
            break;
        }
    }

    return 1;
}

static Hh restore_constraint(TriMesh &mesh, Vh vh0, Vh vh1, std::vector<Hh> &hitscs, std::vector<Vh> &vitscs)
{
    Hh hh_rc = mesh.find_halfedge(vh0, vh1);
    if (hh_rc.is_valid()) return hh_rc;

    std::vector<Hh> hhs {}; // edges intersected with (u0,u1)
    std::vector<Vh> vhs {}; // vertices lying on (u0,u1)

    if (get_intersections(mesh, vh0, vh1, hhs, vhs)) return Hh {};

    // check self-intersection
    bool has_self_intersection {};

    for (Vh vh : vhs) // record overlapping vertices
    {
        vitscs.push_back(vh);
        has_self_intersection = true;
    }

    for (Hh hh : hhs) if (is_sharp(mesh, mesh.edge_handle(hh))) // record intersection pairs
    {
        hitscs.push_back(hh);
        has_self_intersection = true;
    }

    if (has_self_intersection) return Hh {};

    // keep flipping until the edge is recovered
    std::deque<Hh> frontier(hhs.begin(), hhs.end());

    const int max_num_iter = (int)mesh.n_edges();

    for (int iter = 0; !frontier.empty() && iter < max_num_iter; ++iter)
    {
        Hh hh = frontier.front(); frontier.pop_front();

        if (!is_flippable(mesh, hh))
        {
            frontier.push_back(hh);
            continue;
        }

        mesh.flip(mesh.edge_handle(hh));

        if (is_intersecting(mesh, vh0, vh1, hh))
        {
            frontier.push_back(hh);
        }
    }

    hh_rc = mesh.find_halfedge(vh0, vh1);

    return hh_rc;
}

static int restore_constraints(
    TriMesh &mesh,
    const std::vector<Int2> &es,
    const std::unordered_map<Vh, Vh> &dups,
    std::vector<std::pair<Int2, Int2>> &ee_itsc,
    std::vector<std::pair<Int2, int>>  &ev_ovlp)
{
    int missing_edges {};

    for (const auto &vv : es)
    {
        std::vector<Hh> hhs {}; // intersecting edges
        std::vector<Vh> vhs {}; // overlapping vertices

        Vh vh0 = mesh.vertex_handle(vv[0]);
        Vh vh1 = mesh.vertex_handle(vv[1]);

        if (dups.count(vh0)) { vh0 = dups.at(vh0); }
        if (dups.count(vh1)) { vh1 = dups.at(vh1); }

        Hh hhc = restore_constraint(mesh, vh0, vh1, hhs, vhs);

        if (!hhc.is_valid())
        {
            for (Vh vh : vhs) // record overlapping vertices
            {
                ev_ovlp.emplace_back(vv, vh.idx());
            }

            for (Hh hh : hhs) // record intersecting segments
            {
                Vh vh0 = mesh.from_vertex_handle(hh);
                Vh vh1 = mesh.to_vertex_handle  (hh);
                ee_itsc.emplace_back(vv, Int2 { vh0.idx(), vh1.idx() });
            }

            ++missing_edges;
        }
        else
        {
            Eh ehc = mesh.edge_handle(hhc);

            set_sharp(mesh, ehc, true); // mark restored edge as segment

            make_delaunay(mesh, ehc); // maintain Delaunay after restoration
        }
    }

    return missing_edges;
}

////////////////////////////////////////////////////////////////
/// Wrapping up
////////////////////////////////////////////////////////////////

static int triangulate(
    const std::vector<Vec2> &vs,
    const std::vector<Int2> &es,
    TriMesh &mesh,
    std::unordered_map<Vh, Vh> &dups,
    std::vector<std::pair<Int2, Int2>> &ee_itsc,
    std::vector<std::pair<Int2, int>>  &ev_ovlp)
{
    int err {};

    // Copy all vertices from polygon. Better do it at
    // beginning so the vertex is ordered accordingly.
    add_points(mesh, vs);

    // Generate extended domain
    set_domain(mesh);

    // Insert points into the domain
    insert_vertices(mesh, dups);

    // Restore constraint edges in the domain
    err = restore_constraints(mesh, es, dups, ee_itsc, ev_ovlp) != 0;

    return err;
}

int triangulate(
    const std::vector<Vec2> &vs,
    const std::vector<Int2> &es,
    TriMesh &mesh,
    std::unordered_map<int, int>       &vv_ovlp,
    std::vector<std::pair<Int2, Int2>> &ee_itsc,
    std::vector<std::pair<Int2, int>>  &ev_ovlp)
{
    std::unordered_map<Vh, Vh> dups {};
    int err = triangulate(vs, es, mesh, dups, ee_itsc, ev_ovlp);
    for (const auto &vv : dups) vv_ovlp[vv.second.idx()] = vv.first.idx();
    return err;
}

int triangulate(
    const std::vector<Vec2> &vs,
    const std::vector<Int2> &es,
    TriMesh &mesh)
{
    std::unordered_map<Vh, Vh> dups {};
    std::vector<std::pair<Int2, Int2>> ee_itsc {};
    std::vector<std::pair<Int2, int>>  ev_ovlp {};
    return triangulate(vs, es, mesh, dups, ee_itsc, ev_ovlp);
}

////////////////////////////////////////////////////////////////
/// Postprocessing
////////////////////////////////////////////////////////////////

static int hide_exterior_region(TriMesh &mesh, const std::vector<Fh> &fhs)
{
    std::queue<Fh> frontier {};
    std::unordered_set<Fh> visited {};

    for (Fh fh : fhs) if (!visited.count(fh))
    {
        frontier.push(fh);
        visited.insert(fh);
        set_hidden(mesh, fh, true);
    }

    // hide triangles
    while (!frontier.empty())
    {
        Fh fh = frontier.front(); frontier.pop();

        for (Hh hh : mesh.fh_range(fh))
        {
            Eh eh = mesh.edge_handle(hh);
            if (is_sharp(mesh, eh)) continue;
            if (mesh.is_boundary(eh)) continue;

            Fh fi = mesh.opposite_face_handle(hh);
            if (visited.count(fi)) continue;

            frontier.push(fi);
            visited.insert(fi);
            set_hidden(mesh, fi, true);
        }
    }

    // hide edges
    for (Fh fh : mesh.all_faces()) if (is_hidden(mesh, fh))
    {
        for (Eh eh : mesh.fe_range(fh)) if (!is_sharp(mesh, eh))
        {
            set_hidden(mesh, eh, true);
        }
    }

    // hide vertices
    for (Hh hh : mesh.halfedges()) if (mesh.is_boundary(hh))
    {
        set_hidden(mesh, mesh.to_vertex_handle(hh), true);
    }

    return 0;
}

int hide_exterior_region(TriMesh &mesh, const std::vector<Vec2> &seeds)
{
    std::vector<Fh> fhs {};
    Fh fh_last {};

    // start with boundary of the bounding box
    for (Hh hh : mesh.halfedges()) if (mesh.is_boundary(hh))
    {
        fhs.push_back(mesh.opposite_face_handle(hh));
    }

    // user-defined holes
    for (const auto &u : seeds)
    {
        // find the triangle where the point locates
        Fh fh = search_triangle(mesh, u, fh_last);

        // skip any point out of domain
        if (!(fh_last = fh).is_valid()) continue;

        fhs.push_back(fh);
    }

    return hide_exterior_region(mesh, fhs);
}
