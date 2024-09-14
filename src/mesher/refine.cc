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

////////////////////////////////////////////////////////////////
/// Utilities
////////////////////////////////////////////////////////////////

static inline bool is_inside(const TriMesh &mesh, const Fh &fh, const Vec2 &u)
{
    Hh hh0 = mesh.halfedge_handle(fh);
    Hh hh1 = mesh.next_halfedge_handle(hh0);
    Hh hh2 = mesh.next_halfedge_handle(hh1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh2));
    const auto loc = exact_locate(u0, u1, u2, u);
    return loc != TRI_LOC::OUT;
}

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

static inline Vec2 centroid(const TriMesh &mesh, const Fh &fh)
{
    Hh hh0 = mesh.halfedge_handle(fh);
    Hh hh1 = mesh.next_halfedge_handle(hh0);
    Hh hh2 = mesh.next_halfedge_handle(hh1);
    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh0));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hh1));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hh2));
    return (u0 + u1 + u2) / 3.;
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

////////////////////////////////////////////////////////////////
/// Local search
////////////////////////////////////////////////////////////////

static Fh search_primitive(const TriMesh &mesh, const Vec2 &u1, const Fh &fh0, std::vector<Eh> &ehs)
{
    if (is_inside(mesh, fh0, u1)) return fh0;

    PrimitivePlow pp(mesh);

    init(pp, fh0, u1); // setup the plow

    const int max_n_iter = (int)mesh.n_edges(); int n_iter {};

    for (pp.next(); n_iter < max_n_iter; pp.next(), ++n_iter)
    {
        if (pp.status() == PLOW_STATUS::EDGE) // record encroached segments
        {
            Eh eh = mesh.edge_handle(pp.halfedge_handle());
            if (is_sharp(mesh, eh)) { ehs.push_back(eh); }
        }
        else if (pp.status() == PLOW_STATUS::VERT)
        {
            for (Eh eh : mesh.ve_range(pp.vertex_handle()))
            if (is_sharp(mesh, eh)) { ehs.push_back(eh); }
        }

        Fh fh {}; // check if the target is reached

        if (pp.status() == PLOW_STATUS::EDGE)
        {
            fh = mesh.face_handle(pp.halfedge_handle());
        }
        else if (pp.status() == PLOW_STATUS::VERT)
        {
            fh = mesh.face_handle(mesh.halfedge_handle(pp.vertex_handle()));
        }

        if (fh.is_valid() && is_inside(mesh, fh, u1))
        {
            return fh;
        }

        if (pp.status() == PLOW_STATUS::MISS) // searching lost in vain, for some reasons
        {
            break;
        }
    }

    return Fh {};
}

////////////////////////////////////////////////////////////////
/// Refinable events
////////////////////////////////////////////////////////////////

struct Event // either a triangle or a segment
{
    Hh hh;
    Vh vh0, vh1, vh2;
};

template<>
struct std::hash<Event>
{
    size_t operator()(const Event &event) const noexcept
    {
        size_t h0 = hash<Hh>{}(event.hh);
        size_t h1 = hash<Vh>{}(event.vh0);
        size_t h2 = hash<Vh>{}(event.vh1);
        size_t h3 = hash<Vh>{}(event.vh2);
        return ((((h0 ^ h1) << 1) ^ h2) << 1) ^ h3;
    }
};

template<>
struct std::equal_to<Event>
{
    bool operator()(const Event &lhs, const Event &rhs) const noexcept
    {
        return
            lhs.hh  == rhs.hh  &&
            lhs.vh0 == rhs.vh0 &&
            lhs.vh1 == rhs.vh1 &&
            lhs.vh2 == rhs.vh2;
        }
};

static inline bool is_triangle(const TriMesh &mesh, const Event &event)
{
    // check if the event type is triangle
    if (!event.vh2.is_valid()) return false;

    // check if the event was not changed
    if (mesh.is_boundary(event.hh)) return false;

    Hh hh = event.hh;
    if (mesh.to_vertex_handle(hh) != event.vh0)
        hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != event.vh0)
        hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != event.vh0)
        return false;

    hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != event.vh1)
        return false;

    hh = mesh.next_halfedge_handle(hh);
    if (mesh.to_vertex_handle(hh) != event.vh2)
        return false;

    return true;
}

static inline bool is_segment(const TriMesh &mesh, const Event &event)
{
    // check if the event type is segment
    if (event.vh2.is_valid()) return false;

    // check if the event was not changed
    Hh hh = event.hh;
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);

    return vh0 == event.vh0 && vh1 == event.vh1;
}

static inline Event make_triangle(const TriMesh &mesh, Fh fh)
{
    Hh hh = mesh.halfedge_handle(fh);
    Hh hi = mesh.next_halfedge_handle(hh);
    Hh hj = mesh.prev_halfedge_handle(hh);
    Vh vh0 = mesh.to_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle(hi);
    Vh vh2 = mesh.to_vertex_handle(hj);
    return { hh, vh0, vh1, vh2 };
}

static inline Event make_segment(const TriMesh &mesh, Hh hh)
{
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);
    return { hh, vh0, vh1, Vh {} };
}

////////////////////////////////////////////////////////////////
/// Some weird data structures
////////////////////////////////////////////////////////////////

template <class T>
struct unique_vector
{
    inline void push_back(const T &v)
    { if (!vt.count(v)) { vs.push_back(v); vt.insert(v); } }

    inline void pop_back()
    { vt.erase(vs.back()); vs.pop_back(); }

    inline void clear()
    { vs.clear(); vt.clear(); }

    inline const std::vector<T> &vector() const
    { return vs; }

protected:

    std::vector<T> vs;
    std::unordered_set<T> vt;
};

////////////////////////////////////////////////////////////////
/// Segment encroachment
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
    Encroachment(const double min_angle):
    cs2(apex_squared_cosine(min_angle)) {}

    inline bool operator()(const TriMesh &mesh, const Hh &hh) const // If true, split it
    { return !is_exterior(mesh, hh) && is_encroached(mesh, hh, cs2); }

    const double cs2; // (cos(pi - 2t))^2
};

////////////////////////////////////////////////////////////////
/// Triangle quality
////////////////////////////////////////////////////////////////

static inline bool is_bad_triangle(
    const TriMesh &mesh, const Fh &fh,
    const double max_cos2,
    const double max_len2,
    const double max_twoa)
{
    Hh hh = mesh.halfedge_handle(fh);
    Hh hi = mesh.next_halfedge_handle(hh);
    Hh hj = mesh.prev_halfedge_handle(hh);

    //      1      //
    //  0  / \  2  //
    //    /   \    //
    //   2-----0   //
    //      1      //

    const auto u0 = get_xy(mesh, mesh.to_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle(hi));
    const auto u2 = get_xy(mesh, mesh.to_vertex_handle(hj));

    const auto d0 = u2 - u1;
    const auto d1 = u0 - u2;
    const auto d2 = u1 - u0;

    const double dd0 = dot(d0, d0); // len^2 of edge 0
    const double dd1 = dot(d1, d1); // len^2 of edge 1
    const double dd2 = dot(d2, d2); // len^2 of edge 2

    const double dp12 = dot(-d1, d2); // dot product at corner 0
    const double dp20 = dot(-d2, d0); // dot product at corner 1
    const double dp01 = dot(-d0, d1); // dot product at corner 2

    const double cs0 = (dp12*dp12)/(dd1*dd2); // cos^2 of corner 0
    const double cs1 = (dp20*dp20)/(dd2*dd0); // cos^2 of corner 1
    const double cs2 = (dp01*dp01)/(dd0*dd1); // cos^2 of corner 2

    const double twoa = fabs(cross(-d1, d2)); // area*2

    return
        max_cos2 < cs0 ||
        max_cos2 < cs1 ||
        max_cos2 < cs2 ||
        max_len2 < dd0 ||
        max_len2 < dd1 ||
        max_len2 < dd2 ||
        max_twoa < twoa;
}

struct BadTriangle
{
    BadTriangle(const double min_angle, const double max_area, const double max_length):
    max_cos2(cos(min_angle)*cos(min_angle)),
    max_len2(max_length*max_length),
    max_twoa(max_area*2.0) {}

    inline bool operator()(const TriMesh &mesh, const Fh &fh) const // If true, split it
    { return is_bad_triangle(mesh, fh, max_cos2, max_len2, max_twoa); }

    const double max_cos2; // lower bound of corner angle
    const double max_len2; // upper bound of edge length
    const double max_twoa; // upper bound of triangle area
};

////////////////////////////////////////////////////////////////
/// Refinement
////////////////////////////////////////////////////////////////

static inline Vec2 splitting_position(const TriMesh &mesh, Hh hh)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    double t = 0.5;

    /// The concentric circle algorithm is adapted from 'Triangle' lib
    ///   to handle small angle in the input. Ref 'Delaunay Refinement
    ///   Algorithms for Triangular Mesh Generation' [01 J.R.Shewchuk]

    Hh hi = mesh.opposite_halfedge_handle(hh);

    bool is_acute_v0 = 
        !mesh.is_boundary(hh) && is_sharp(mesh, mesh.edge_handle(mesh.prev_halfedge_handle(hh))) ||
        !mesh.is_boundary(hi) && is_sharp(mesh, mesh.edge_handle(mesh.next_halfedge_handle(hi)));

    bool is_acute_v1 = 
        !mesh.is_boundary(hh) && is_sharp(mesh, mesh.edge_handle(mesh.next_halfedge_handle(hh))) ||
        !mesh.is_boundary(hi) && is_sharp(mesh, mesh.edge_handle(mesh.prev_halfedge_handle(hi)));

    if (is_acute_v0 || is_acute_v1)
    {
        const double l = norm(u1 - u0);
        double ep = 1; // nearest power of 2

        // find the ratio that splits the segment most evenly
        while (ep*3.0 < l) { ep *= 2.0; }
        while (ep*1.5 > l) { ep *= 0.5; }

        t = ep / l;
        t = is_acute_v1 ? 1 - t : t;
    }

    return u0*(1-t) + u1*t;
}

static int split_segments(TriMesh &mesh, const Encroachment &encroached)
{
    const int max_n_iter = 100;

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    std::deque<Event> events {}; // to split

    // Before enqueueing encroached segments, free vertices in
    // the diametral circle of the segment should be deleted.
    // [TODO]

    for (Eh eh : mesh.edges()) if (is_sharp(mesh, eh))
    for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
    { events.push_back(make_segment(mesh, hh)); }

    while (!events.empty())
    {
        auto event = events.front(); events.pop_front();

        if (is_segment(mesh, event))
        {
            // find a position for splitting
            const auto u = splitting_position(mesh, event.hh);

            // split the segment
            Vh vo = mesh.new_vertex({ u[0], u[1], 0 });
            split_edge(mesh, event.hh, vo);

            // edges to flip
            Eh eos[4]; int ne {};
            for (auto hdge : mesh.voh_range(vo))
                if (!hdge.next().edge().is_boundary())
                    eos[ne++] = hdge.next().edge();

            // maintain Delaunayhood
            delaunifier.reset(); delaunifier.enqueue(eos, ne); int n_iter {};

            // encountered edges
            unique_vector<Eh> ehs {};

            // record encountered edges while testing Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid() && n_iter < max_n_iter; eh = delaunifier.next(), ++n_iter)
            { ehs.push_back(eh); }

            // check encroachment of encountered segments
            for (Eh eh : ehs.vector()) if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { events.push_back(make_segment(mesh, hh)); }

            // check encroachment of two new segments afterwards
            for (Eh eh : mesh.ve_range(vo)) if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { events.push_back(make_segment(mesh, hh)); }
        }
    }

    return 0;
}

static int split_interior(TriMesh &mesh, const BadTriangle &bad_triangle, const Encroachment &encroached)
{
    const int max_n_iter = 100;

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    std::deque<Event> events {}; // to split

    // check quality of all triangles at initialization
    for (Fh fh : mesh.faces()) if (bad_triangle(mesh, fh))
    { events.push_back(make_triangle(mesh, fh)); }

    while (!events.empty())
    {
        auto event = events.front(); events.pop_front();

        if (is_triangle(mesh, event))
        {
            // find the circumcenter of the problematic triangle
            const auto u = circumcenter(mesh, mesh.face_handle(event.hh));

            // the segments that were run over during searching
            std::vector<Eh> ess {};

            // search for the primitive at which the circumcenter locates
            Fh fho = search_primitive(mesh, u, mesh.face_handle(event.hh), ess);

            if (!ess.empty()) // some segments are in the way hence encroached
            {
                std::reverse(ess.begin(), ess.end());

                for (Eh eh : ess) for (Hh hh : mesh.eh_range(eh)) if (!is_exterior(mesh, hh))
                { events.push_front(make_segment(mesh, hh)); }

                continue; // abort splitting
            }

            // skip any point out of domain (won't happen as some segments must
            // be encroached and splitting must have been aborted.)
            if (!fho.is_valid() || is_hidden(mesh, fho)) continue;

            // at which part of the triangle the point locates
            Hh hho {}; const auto loc = locate(mesh, fho, u, hho);

            // skip any vertex with invalid location (won't happen)
            if (loc == TRI_LOC::OUT) continue;

            // skip any duplicated vertices (hardly happens, or Delaunayhood is
            // not maintained as the circumcircle contains a point in the mesh.)
            if ((int)loc & (int)TRI_LOC::VS) continue;

            // allocate a new vertex
            Vh vho = mesh.new_vertex({ u[0], u[1], 0 });

            // insert the vertex into the triangle or onto the edge
            if (loc == TRI_LOC::IN) split_face(mesh, fho, vho);
            else                    split_edge(mesh, hho, vho);

            // edges to flip
            Eh eos[4]; int ne {};
            for (auto hdge : mesh.voh_range(vho))
                if (!hdge.next().edge().is_boundary())
                    eos[ne++] = hdge.next().edge();

            // maintain Delaunayhood
            delaunifier.reset(); delaunifier.enqueue(eos, ne); int n_iter {};

            // encountered edges
            unique_vector<Eh> ehs {};

            // record encountered edges while testing Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid() && n_iter < max_n_iter; eh = delaunifier.next(), ++n_iter)
            { ehs.push_back(eh); }

            // potentially encroached segments
            std::vector<Hh> hss {};

            // check encroachment over encountered segments
            for (Eh eh : ehs.vector()) if (is_sharp(mesh, eh)) for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { hss.push_back(hh); }

            if (!hss.empty()) // some segments are encroached
            {
                for (Hh hh : hss) { events.push_front(make_segment(mesh, hh)); }

                // undo vertex insertion
                // [TODO]

                // check quality over affected triangles
                // [TODO]

                continue; // abort splitting
            }

            // check quality of the encountered triangles
            for (Eh eh : ehs.vector()) for (Fh fh : mesh.ef_range(eh)) if (fh.is_valid())
            if (!is_hidden(mesh, fh)) if (bad_triangle(mesh, fh))
            { events.push_back(make_triangle(mesh, fh)); }
        }

        else if (is_segment(mesh, event))
        {
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////
/// Wrapping up
////////////////////////////////////////////////////////////////

int refine(TriMesh &mesh, const double min_angle)
{
    Encroachment encroached(min_angle);

    BadTriangle bad_triangle(min_angle, 1e10, 1e5);

    int err {};

    err = split_segments(mesh, encroached);

    err = split_interior(mesh, bad_triangle, encroached);

    return err;
}