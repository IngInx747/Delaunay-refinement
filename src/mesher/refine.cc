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
/// Containers
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

static inline bool is_exterior(const TriMesh &mesh, const Fh &fh)
{
    return is_hidden(mesh, fh);
}

static inline bool is_exterior(const TriMesh &mesh, const Eh &eh)
{
    return is_hidden(mesh, eh);
}

static inline bool is_exterior(const TriMesh &mesh, const Hh &hh)
{
    return mesh.is_boundary(hh) || is_hidden(mesh, mesh.face_handle(hh));
}

static inline void split(TriMesh &mesh, Hh hh, Vh vh)
{
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);

    mesh.split_edge_copy(mesh.edge_handle(hh), vh);

    // OpenMesh copies property to all 4 new edges even 2 of
    // which are not sharp. Unsharp the 2 edges accordingly.
    for (Hh hh : mesh.voh_range(vh))
    if (mesh.to_vertex_handle(hh) != vh0)
    if (mesh.to_vertex_handle(hh) != vh1)
    { set_sharp(mesh, mesh.edge_handle(hh), false); }

    // One of 2 new diagonal edges can be in the exterior
    // region and hence should be hidden.
    for (Eh eh : mesh.ve_range(vh))
    if (is_exterior(mesh, mesh.halfedge_handle(eh, 0)))
    if (is_exterior(mesh, mesh.halfedge_handle(eh, 1)))
    { set_hidden(mesh, eh, true); }
}

static inline void split(TriMesh &mesh, Fh fh, Vh vh)
{
    mesh.split_copy(fh, vh);
}

static inline Vh collapse(TriMesh &mesh, Vh vh)
{
    assert(!mesh.is_boundary(vh));

    Hh hh {}; // to collapse

    for (Hh hi : mesh.voh_range(vh))
    {
        Vh vi = mesh.to_vertex_handle(hi);
        const auto u0 = get_xy(mesh, vi);
        bool convex { true };

        for (Hh hj : mesh.voh_range(vh)) // check n-2 triangles
        {
            Vh vj = mesh.to_vertex_handle(hj);
            Vh vk = mesh.to_vertex_handle(mesh.next_halfedge_handle(hj));
            if (vj == vi || vk == vi) continue;

            const auto u1 = get_xy(mesh, vj);
            const auto u2 = get_xy(mesh, vk);
            const int r = orientation(u0, u1, u2);
            if (r <= 0) { convex = false; break; }
        }

        if (convex) { hh = hi; break; }
    }

    assert(hh.is_valid());

    Vh vi = mesh.to_vertex_handle(hh);

    mesh.collapse(hh);

    return vi;
}

////////////////////////////////////////////////////////////////
/// Refinable events
////////////////////////////////////////////////////////////////

struct Event // either a triangle or a segment
{
    Vh vh0, vh1, vh2;
};

template<>
struct std::hash<Event>
{
    size_t operator()(const Event &event) const noexcept
    {
        size_t h0 = hash<Vh>{}(event.vh0);
        size_t h1 = hash<Vh>{}(event.vh1);
        size_t h2 = hash<Vh>{}(event.vh2);
        return ((h0 ^ h1) << 1) ^ h2;
    }
};

template<>
struct std::equal_to<Event>
{
    bool operator()(const Event &lhs, const Event &rhs) const noexcept
    {
        return
            lhs.vh0 == rhs.vh0 &&
            lhs.vh1 == rhs.vh1 &&
            lhs.vh2 == rhs.vh2;
        }
};

static inline bool is_triangle(const Event &event)
{
    return event.vh2.is_valid();
}

static inline bool is_segment(const Event &event)
{
    return !event.vh2.is_valid();
}

static inline Fh get_triangle(const TriMesh &mesh, const Event &event)
{
    if (is_deleted(mesh, event.vh0) ||
        is_deleted(mesh, event.vh1) ||
        is_deleted(mesh, event.vh2))
        return Fh {};

    Hh hh = mesh.find_halfedge(event.vh0, event.vh1);
    if (!hh.is_valid()) return Fh {};

    Vh vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
    if (vh2 != event.vh2) return Fh {};

    return mesh.face_handle(hh);
}

static inline Hh get_segment(const TriMesh &mesh, const Event &event)
{
    if (is_deleted(mesh, event.vh0) ||
        is_deleted(mesh, event.vh1))
        return Hh {};

    return mesh.find_halfedge(event.vh0, event.vh1);
}

static inline Event make_triangle(const TriMesh &mesh, Fh fh)
{
    Hh hh = mesh.halfedge_handle(fh);
    Hh hi = mesh.next_halfedge_handle(hh);
    Hh hj = mesh.prev_halfedge_handle(hh);
    Vh vh0 = mesh.to_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle(hi);
    Vh vh2 = mesh.to_vertex_handle(hj);
    return { vh0, vh1, vh2 };
}

static inline Event make_segment(const TriMesh &mesh, Hh hh)
{
    Vh vh0 = mesh.from_vertex_handle(hh);
    Vh vh1 = mesh.to_vertex_handle  (hh);
    return { vh0, vh1, Vh {} };
}

////////////////////////////////////////////////////////////////
/// Segment encroachment
////////////////////////////////////////////////////////////////

static inline double apex_squared_cosine(const double min_angle)
{
    const double cs0 = cos(min_angle); // cos(t)
    const double cs1 = cs0*cs0*2. - 1; // cos(pi - 2t)
    return cs1*cs1;
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
    const Hh hh[3] {
        mesh.halfedge_handle(fh),
        mesh.next_halfedge_handle(mesh.halfedge_handle(fh)),
        mesh.prev_halfedge_handle(mesh.halfedge_handle(fh)),
    };

    //      1      //
    //  0  / \  2  //
    //    /   \    //
    //   2-----0   //
    //      1      //

    const Vec2 u[3] {
        get_xy(mesh, mesh.to_vertex_handle(hh[1])),
        get_xy(mesh, mesh.to_vertex_handle(hh[2])),
        get_xy(mesh, mesh.to_vertex_handle(hh[0])),
    };

    const Vec2 d[3] {
        u[2] - u[1],
        u[0] - u[2],
        u[1] - u[0],
    };

    const double twoa = fabs(cross(-d[1], d[2])); // area*2

    // check area upper bound
    if (twoa > max_twoa) return true;

    const double dd[3] {
        dot(d[0], d[0]), // len^2 of edge 0
        dot(d[1], d[1]), // len^2 of edge 1
        dot(d[2], d[2]), // len^2 of edge 2
    };

    // check length upper bound
    if (dd[0] > max_len2 ||
        dd[1] > max_len2 ||
        dd[2] > max_len2 ) return true;

    const double dp[3] {
        dot(-d[1], d[2]), // dot product at corner 0
        dot(-d[2], d[0]), // dot product at corner 1
        dot(-d[0], d[1]), // dot product at corner 2
    };

    // use the shortest edge as base
    double ddb = dd[0]; int i = 0;
    if (ddb > dd[1]) { ddb = dd[1]; i = 1; }
    if (ddb > dd[2]) { ddb = dd[2]; i = 2; }

    int j = (i + 1) % 3; // leg 1
    int k = (i + 2) % 3; // leg 2

    const double ratio = dd[j] / dd[k]; // length ratio of two legs of the corner

    if (is_sharp(mesh, mesh.edge_handle(hh[j])) &&
        is_sharp(mesh, mesh.edge_handle(hh[k])) &&
        fabs(ratio - 1.0) < 1e-3) return false;

    //if (is_sharp(mesh, mesh.edge_handle(hh[j])) ||
    //    is_sharp(mesh, mesh.edge_handle(hh[k])) ) return false;

    const double cs2 = (dp[i]*dp[i])/(dd[j]*dd[k]); // cos^2 of the smallest corner angle

    // check angle lower bound
    return cs2 > max_cos2;
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
        !is_exterior(mesh, hh) && is_sharp(mesh, mesh.edge_handle(mesh.prev_halfedge_handle(hh))) ||
        !is_exterior(mesh, hi) && is_sharp(mesh, mesh.edge_handle(mesh.next_halfedge_handle(hi))) ;

    bool is_acute_v1 = 
        !is_exterior(mesh, hh) && is_sharp(mesh, mesh.edge_handle(mesh.next_halfedge_handle(hh))) ||
        !is_exterior(mesh, hi) && is_sharp(mesh, mesh.edge_handle(mesh.prev_halfedge_handle(hi))) ;

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

    for (Eh eh : mesh.edges()) { if (is_sharp(mesh, eh))
    for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
    { events.push_back(make_segment(mesh, hh)); } }

    while (!events.empty())
    {
        auto event = events.front(); events.pop_front();

        if (is_segment(event))
        {
            // restore the segment from the event
            Hh ho = get_segment(mesh, event);

            // the segment was gone during refining
            if (!ho.is_valid()) continue;

            // find a position for splitting
            const auto u = splitting_position(mesh, ho);

            // split the segment
            Vh vo = mesh.new_vertex({ u[0], u[1], 0 });
            split(mesh, ho, vo);

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
            { if (!is_exterior(mesh, eh)) ehs.push_back(eh); }

            // check encroachment of encountered segments
            for (Eh eh : ehs.vector()) { if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { events.push_back(make_segment(mesh, hh)); } }

            // check encroachment of two new segments afterwards
            for (Eh eh : mesh.ve_range(vo)) { if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { events.push_back(make_segment(mesh, hh)); } }
        }
    }

    return 0;
}

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

        if (is_triangle(event))
        {
            // restore the triangle from the event
            Fh fo = get_triangle(mesh, event);

            // the triangle was gone during refining
            if (!fo.is_valid()) continue;

            // find the circumcenter of the problematic triangle
            const auto u = circumcenter(mesh, fo);

            // the segments that were run over during searching
            std::vector<Eh> ess {};

            // search for the primitive at which the circumcenter locates
            fo = search_primitive(mesh, u, fo, ess);

            if (!ess.empty()) // some segments are in the way hence encroached
            {
                std::reverse(ess.begin(), ess.end());

                for (Eh eh : ess) { for (Hh hh : mesh.eh_range(eh)) if (!is_exterior(mesh, hh))
                { events.push_front(make_segment(mesh, hh)); } }

                //for (Hh hh : mesh.eh_range(ess.front())) if (!is_exterior(mesh, hh))
                //{ events.push_front(make_segment(mesh, hh)); }

                continue; // abort splitting
            }

            // skip any point out of domain (won't happen as some segments must
            // be encroached and splitting must have been aborted.)
            if (!fo.is_valid() || is_exterior(mesh, fo)) continue;

            // at which part of the triangle the point locates
            Hh ho {}; const auto loc = locate(mesh, fo, u, ho);

            // skip any vertex with invalid location (won't happen)
            if (loc == TRI_LOC::OUT) continue;

            // skip any duplicated vertices (hardly happens, or Delaunayhood is
            // not maintained as the circumcircle contains a point in the mesh.)
            if ((int)loc & (int)TRI_LOC::VS) continue;

            // allocate a new vertex
            Vh vo = mesh.new_vertex({ u[0], u[1], 0 });

            // insert the vertex into the triangle or onto the edge
            if (loc == TRI_LOC::IN) split(mesh, fo, vo);
            else                    split(mesh, ho, vo);

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
            { if (!is_exterior(mesh, eh)) ehs.push_back(eh); }

            // potentially encroached segments
            std::vector<Hh> hss {};

            // check encroachment over encountered segments
            for (Eh eh : ehs.vector()) { if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { hss.push_back(hh); } }

            if (!hss.empty()) // some segments are encroached
            {
                for (Hh hh : hss) { events.push_front(make_segment(mesh, hh)); }

                // undo vertex insertion
                Vh vc = collapse(mesh, vo);

                // edges to flip
                std::vector<Eh> ecs {};
                for (Eh eh : mesh.ve_range(vc)) ecs.push_back(eh);

                // maintain Delaunayhood
                delaunifier.reset(); delaunifier.enqueue(ecs.data(), (int)ecs.size()); int n_iter {};

                ehs.clear();

                // record encountered edges while testing Delaunayhood
                for (Eh eh = delaunifier.next(); eh.is_valid() && n_iter < max_n_iter; eh = delaunifier.next(), ++n_iter)
                if (!is_exterior(mesh, eh)) { ehs.push_back(eh); }
            }

            // encountered triangles
            unique_vector<Fh> fhs {};
            for (Eh eh : ehs.vector()) { for (Fh fh : mesh.ef_range(eh))
            if (fh.is_valid()) if (!is_exterior(mesh, fh)) { fhs.push_back(fh); } }

            // check quality of these triangles
            for (Fh fh : fhs.vector()) { if (bad_triangle(mesh, fh))
            { events.push_back(make_triangle(mesh, fh)); } }
        }

        else if (is_segment(event))
        {
            // restore the segment from the event
            Hh ho = get_segment(mesh, event);

            // the segment was gone during refining
            if (!ho.is_valid()) continue;

            // find a position for splitting
            const auto u = splitting_position(mesh, ho);

            // split the segment
            Vh vo = mesh.new_vertex({ u[0], u[1], 0 });
            split(mesh, ho, vo);

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
            { if (!is_exterior(mesh, eh)) ehs.push_back(eh); }

            // check encroachment of encountered segments
            for (Eh eh : ehs.vector()) { if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { events.push_front(make_segment(mesh, hh)); } }

            // check encroachment of two new segments afterwards
            for (Eh eh : mesh.ve_range(vo)) { if (is_sharp(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
            { events.push_front(make_segment(mesh, hh)); } }

            // encountered triangles
            unique_vector<Fh> fhs {};
            for (Eh eh : ehs.vector()) { for (Fh fh : mesh.ef_range(eh))
            if (fh.is_valid()) if (!is_exterior(mesh, fh)) { fhs.push_back(fh); } }

            // check quality of these triangles
            for (Fh fh : fhs.vector()) { if (bad_triangle(mesh, fh))
            { events.push_back(make_triangle(mesh, fh)); } }
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

    mesh.garbage_collection();

    return err;
}