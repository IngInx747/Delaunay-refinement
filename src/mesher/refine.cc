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
/// Predicates
////////////////////////////////////////////////////////////////

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

static inline bool is_segment(const TriMesh &mesh, const Eh &eh)
{
    return is_sharp(mesh, eh);
}

static inline bool is_segment(const TriMesh &mesh, const Hh &hh)
{
    return is_sharp(mesh, mesh.edge_handle(hh));
}

static inline bool is_segment(const TriMesh &mesh, const Vh &vh)
{
    for (Eh eh : mesh.ve_range(vh)) if (is_segment(mesh, eh)) return true;
    return false; // is the vertex on a segment
}

static inline bool is_endian(const TriMesh &mesh, const Vh &vh)
{
    return is_sharp(mesh, vh); // is the vertex from an original segment
}

struct IsSegment
{
    inline bool operator()(const TriMesh &mesh, const Hh &hh) const
    { return is_segment(mesh, hh); }
};

static inline Vh segment_head(const TriMesh &mesh, const Hh &hh)
{
    Hh hi = hh; IsSegment pred {};
    while (!is_endian(mesh, mesh.to_vertex_handle(hi)))
    { if ((hi = next(mesh, pred, hi)) == hh) break; }
    return mesh.to_vertex_handle(hi);
}

static inline Vh segment_tail(const TriMesh &mesh, const Hh &hh)
{
    Hh hi = hh; IsSegment pred {};
    while (!is_endian(mesh, mesh.from_vertex_handle(hi)))
    { if ((hi = prev(mesh, pred, hi)) == hh) break; }
    return mesh.from_vertex_handle(hi);
}

static void mark_endians(TriMesh &mesh)
{
    for (Hh hh : mesh.halfedges()) if (is_sharp(mesh, mesh.edge_handle(hh)))
    { set_sharp(mesh, mesh.to_vertex_handle(hh), true); }
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
    return fuzzy_delaunay(u0, u1, u2, u3);
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

static inline bool is_collapsable(TriMesh &mesh, const Hh &hhc)
{
    Vh vhc = mesh.from_vertex_handle(hhc);
    Vh vh0 = mesh.to_vertex_handle  (hhc);
    const auto u0 = get_xy(mesh, vh0);

    for (Hh hh : mesh.voh_range(vhc)) // check n-2 triangles
    {
        Vh vh1 = mesh.to_vertex_handle(hh);
        Vh vh2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
        if (vh1 == vh0 || vh2 == vh0) continue;
        const auto u1 = get_xy(mesh, vh1);
        const auto u2 = get_xy(mesh, vh2);
        const int r = orientation(u0, u1, u2);
        if (r <= 0) return false;
    }

    return true;
}

static inline Vh collapse(TriMesh &mesh, Vh vh)
{
    assert(!mesh.is_boundary(vh));

    Hh hhc {}; // to collapse

    for (Hh hh : mesh.voh_range(vh)) if (is_collapsable(mesh, hh))
    { hhc = hh; break; }

    assert(hhc.is_valid());

    Vh vi = mesh.to_vertex_handle(hhc);
    mesh.collapse(hhc);
    return vi;
}

static Fh search_primitive(const TriMesh &mesh, const Vec2 &u1, const Fh &fh0, std::vector<Eh> &ehs)
{
    if (is_inside(mesh, fh0, u1)) return fh0;

    PrimitivePlow pp(mesh);

    init(pp, fh0, u1); // setup the plow

    const int max_n_iter = (int)mesh.n_edges(); int n_iter {};

    for (pp.next(); n_iter < max_n_iter; pp.next(), ++n_iter)
    {
        if (pp.status() == PLOW_STATUS::EDGE) // record segments
        {
            Eh eh = mesh.edge_handle(pp.halfedge_handle());
            if (is_segment(mesh, eh)) { ehs.push_back(eh); }
        }
        else if (pp.status() == PLOW_STATUS::VERT)
        {
            for (Eh eh : mesh.ve_range(pp.vertex_handle()))
            if (is_segment(mesh, eh)) { ehs.push_back(eh); }
        }

        Fh fh {}; // check if the target is reached

        if (pp.status() == PLOW_STATUS::EDGE)
        {
            fh = mesh.opposite_face_handle(pp.halfedge_handle());
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

static inline bool is_encroached(const TriMesh &mesh, const Hh &hh)
{
    const auto d0 = -get_dxy(mesh, mesh.next_halfedge_handle(hh));
    const auto d1 =  get_dxy(mesh, mesh.prev_halfedge_handle(hh));
    return (dot(d0, d1) < 0.); // a.b < 0
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

static inline bool is_encroached(const TriMesh &mesh, const Hh &hh, const Vec2 &u, const double cs2)
{
    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto d0 = u0 - u;
    const auto d1 = u1 - u;
    const double d00 = dot(d0, d0); // a^2
    const double d11 = dot(d1, d1); // b^2
    const double d01 = dot(d0, d1); // a.b
    return (d01 < 0.) && (d01*d01 >= d00*d11*cs2); // a.b < 0 and |a.b| > (a*b)|cos(t0)|
}

struct Encroachment
{
    Encroachment(const double min_angle): cs2(apex_squared_cosine(min_angle)) {}

    inline bool operator()(const TriMesh &mesh, const Hh &hh) const // If true, split it
    { return !is_exterior(mesh, hh) && is_encroached(mesh, hh, cs2); }

    inline bool operator()(const TriMesh &mesh, const Hh &hh, const Vec2 &u) const // If true, split it
    { return !is_exterior(mesh, hh) && is_encroached(mesh, hh, u, cs2); }

    const double cs2; // (cos(pi - 2t))^2
};

////////////////////////////////////////////////////////////////
/// Triangle quality
////////////////////////////////////////////////////////////////

static inline bool is_subtending_input_angle(const TriMesh &mesh, const Hh &hh)
{
    IsSegment pred {};

    // try getting two segments between which the base is
    Hh hh0 = next(mesh, pred, hh);
    Hh hh1 = prev(mesh, pred, hh);

    // one or both ends of the base do not lie on segments
    if (!is_segment(mesh, hh0) || !is_segment(mesh, hh1)) return false;

    Vh vh00 = segment_head(mesh, hh0);
    Vh vh01 = segment_tail(mesh, hh0);
    Vh vh10 = segment_head(mesh, hh1);
    Vh vh11 = segment_tail(mesh, hh1);
    Vh vhc {}; // common end

    if (vh01 == vh10) vhc = vh01;
    if (vh00 == vh11) vhc = vh00;

    // two segment not sharing a common end
    if (!vhc.is_valid()) return false;

    const auto u0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto u1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    const auto uc = get_xy(mesh, vhc);
    const auto d0 = u0 - uc;
    const auto d1 = u1 - uc;
    const double rd = dot(d0, d0) / dot(d1, d1);

    // check if two ends lie on a circle centering the common end
    return fabs(rd - 1.0) < 1e-3;
}

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

    const double cs2 = (dp[i]*dp[i])/(dd[j]*dd[k]); // cos^2 of the smallest corner angle

    // check angle lower bound
    if (cs2 > max_cos2) {

    // skip if the smallest edge subtends an input angle
    if (!is_segment(mesh, hh[i]) && is_subtending_input_angle(mesh, hh[i])) return false;

    return true; }

    // the triangle passes all tests, hence no refinement is needed
    return false;
}

struct BadTriangle
{
    BadTriangle(const double min_angle, const double max_length, const double max_area):
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

static inline double calc_priority(const TriMesh &mesh, const Hh &hh)
{
    return mesh.calc_edge_sqr_length(hh);
}

static inline double calc_priority(const TriMesh &mesh, const Fh &fh)
{
    Hh hh0 = mesh.halfedge_handle(fh);
    Hh hh1 = mesh.next_halfedge_handle(hh0);
    Hh hh2 = mesh.prev_halfedge_handle(hh0);
    const double l0 = mesh.calc_edge_sqr_length(hh0);
    const double l1 = mesh.calc_edge_sqr_length(hh1);
    const double l2 = mesh.calc_edge_sqr_length(hh2);
    return (l1<l0) ? (l1<l2 ? l1 : l2) : (l0<l2 ? l0 : l2); // the shortest edge
}

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
        !is_exterior(mesh, hh) && is_segment(mesh, mesh.prev_halfedge_handle(hh)) ||
        !is_exterior(mesh, hi) && is_segment(mesh, mesh.next_halfedge_handle(hi)) ;

    bool is_acute_v1 = 
        !is_exterior(mesh, hh) && is_segment(mesh, mesh.next_halfedge_handle(hh)) ||
        !is_exterior(mesh, hi) && is_segment(mesh, mesh.prev_halfedge_handle(hi)) ;

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
    int n_new_vertices {};

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    // Before enqueueing encroached segments, free vertices in
    // the diametral circle of the segment should be deleted.
    // [TODO]

    std::deque<Event> events {}; // primitives to split

    for (Eh eh : mesh.edges()) { if (is_segment(mesh, eh))
    for (Hh hh : mesh.eh_range(eh)) { if (encroached(mesh, hh))
    { events.push_back(make_segment(mesh, hh)); } } }

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

            // affected edges
            unique_vector<Eh> eas {};

            // edges that are potentially non-Delaunay
            { Eh ehs[4]; int ne {}; for (Hh hh : mesh.voh_range(vo)) if (!mesh.is_boundary(hh))
            { ehs[ne++] = mesh.edge_handle(mesh.next_halfedge_handle(hh)); }
            delaunifier.reset(); delaunifier.enqueue(ehs, ne); }

            // record encountered edges while testing Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next())
            { if (!is_exterior(mesh, eh)) eas.push_back(eh); }

            // check encroachment of two new segments
            for (Eh eh : mesh.ve_range(vo)) { if (is_segment(mesh, eh)) eas.push_back(eh); }

            // check encroachment of encountered segments
            for (Eh eh : eas.vector()) { if (!is_deleted(mesh, eh)) if (is_segment(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) { if (encroached(mesh, hh))
            { events.push_back(make_segment(mesh, hh)); } } }

            ++n_new_vertices;
        }
    }

    return 0;
}

static int split_interior(TriMesh &mesh, const BadTriangle &bad_triangle, const Encroachment &encroached)
{
    struct Entry { Event event; double score; };

    const auto comp = [](const Entry &a, const Entry &b)
    {
        const int ia = is_segment(a.event) ? 1 : 0;
        const int ib = is_segment(b.event) ? 1 : 0;
        return (ia == ib) ? a.score < b.score : ia < ib;
    };

    int n_new_vertices {};

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    std::priority_queue<Entry, std::deque<Entry>, decltype(comp)> events(comp);

    // check quality of all triangles at initialization
    for (Fh fh : mesh.faces()) if (bad_triangle(mesh, fh))
    { events.push({ make_triangle(mesh, fh), calc_priority(mesh, fh) }); }

    while (!events.empty())
    {
        auto event = events.top().event; events.pop();

        if (is_triangle(event))
        {
            // restore the triangle from the event
            Fh fr = get_triangle(mesh, event);

            // the triangle was gone during refining
            if (!fr.is_valid()) continue;

            // find the circumcenter of the problematic triangle
            const auto u = circumcenter(mesh, fr);

            // search for the primitive at which the circumcenter locates
            std::vector<Eh> ess {}; Fh fo = search_primitive(mesh, u, fr, ess);

            if (!ess.empty()) // some segments are in the way hence encroached
            {
                for (Hh hh : mesh.eh_range(ess.front())) { if (!is_exterior(mesh, hh))// if (encroached(mesh, hh, u))
                { events.push({ make_segment(mesh, hh), calc_priority(mesh, hh) }); } }
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

            // affected edges and triangles
            unique_vector<Eh> eas {};
            unique_vector<Fh> fas {};

            // edges that are potentially non-Delaunay
            { Eh ehs[4]; int ne {}; for (Hh hh : mesh.voh_range(vo)) if (!mesh.is_boundary(hh))
            { ehs[ne++] = mesh.edge_handle(mesh.next_halfedge_handle(hh)); }
            delaunifier.reset(); delaunifier.enqueue(ehs, ne); }

            // record encountered edges while maintaining Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next())
            { if (!is_exterior(mesh, eh)) eas.push_back(eh); }

            // check encroachment over encountered segments
            std::vector<Hh> hss {};
            for (Eh eh : eas.vector()) { if (!is_deleted(mesh, eh)) if (is_segment(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) { if (encroached(mesh, hh)) { hss.push_back(hh); } } }

            if (!hss.empty()) // some segments are encroached
            {
                for (Hh hh : hss) { events.push({ make_segment(mesh, hh), calc_priority(mesh, hh) }); }

                // undo vertex insertion
                Vh vc = collapse(mesh, vo);

                // edges that are potentially non-Delaunay
                { std::vector<Eh> ecs {}; for (Eh eh : mesh.ve_range(vc)) { ecs.push_back(eh); }
                delaunifier.reset(); delaunifier.enqueue(ecs.data(), (int)ecs.size()); }

                // record encountered edges while maintaining Delaunayhood
                for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next())
                { if (!is_exterior(mesh, eh)) eas.push_back(eh); }
            }

            // encountered triangles
            for (Eh eh : eas.vector()) if (!is_deleted(mesh, eh)) { for (Fh fh : mesh.ef_range(eh))
            if (fh.is_valid()) if (!is_exterior(mesh, fh)) { fas.push_back(fh); } }

            // check quality of these triangles
            for (Fh fh : fas.vector()) if (!is_deleted(mesh, fh)) { if (bad_triangle(mesh, fh))
            { events.push({ make_triangle(mesh, fh), calc_priority(mesh, fh) }); } }

            ++n_new_vertices;
        }

        else if (is_segment(event))
        {
            // restore the segment from the event
            Hh ho = get_segment(mesh, event);

            // the segment was gone during refining
            if (!ho.is_valid()) continue;

            // affected edges and triangles
            unique_vector<Eh> eas {};
            unique_vector<Fh> fas {};

            // Before enqueueing encroached segments, free vertices in
            // the diametral circle of the segment should be deleted.
            for ( ; ho.is_valid() && is_encroached(mesh, ho); ho = get_segment(mesh, event))
            {
                // the apex vertex
                Vh vh = mesh.to_vertex_handle(mesh.next_halfedge_handle(ho));

                // stop if the apex lies on a segment
                if (is_segment(mesh, vh)) break;

                // remove the apex vertex
                Vh vc = collapse(mesh, vh);

                // edges that are potentially non-Delaunay
                std::vector<Eh> ehs {}; for (Eh eh : mesh.ve_range(vc)) { ehs.push_back(eh); }
                delaunifier.reset(); delaunifier.enqueue(ehs.data(), (int)ehs.size());

                // record edges while maintaining Delaunayhood
                for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next())
                { if (!is_exterior(mesh, eh)) eas.push_back(eh); }

                // record affected triangles
                for (Fh fh : mesh.vf_range(vc)) { if (!is_exterior(mesh, fh)) fas.push_back(fh); }
            }

            // won't happen
            if (!ho.is_valid()) continue;

            // find a position for splitting
            const auto u = splitting_position(mesh, ho);

            // split the segment
            Vh vo = mesh.new_vertex({ u[0], u[1], 0 });
            split(mesh, ho, vo);

            // edges that are potentially non-Delaunay
            { Eh ehs[4]; int ne {}; for (Hh hh : mesh.voh_range(vo)) if (!mesh.is_boundary(hh))
            { ehs[ne++] = mesh.edge_handle(mesh.next_halfedge_handle(hh)); }
            delaunifier.reset(); delaunifier.enqueue(ehs, ne); }

            // record encountered edges while maintaining Delaunayhood
            for (Eh eh = delaunifier.next(); eh.is_valid(); eh = delaunifier.next())
            { if (!is_exterior(mesh, eh)) eas.push_back(eh); }

            // check encroachment of two new segments afterwards
            for (Eh eh : mesh.ve_range(vo)) { if (is_segment(mesh, eh)) eas.push_back(eh); }

            // check encroachment of encountered segments
            for (Eh eh : eas.vector()) { if (!is_deleted(mesh, eh)) if (is_segment(mesh, eh))
            for (Hh hh : mesh.eh_range(eh)) { if (encroached(mesh, hh))
            { events.push({ make_segment(mesh, hh), calc_priority(mesh, hh) }); } } }

            // check quality of affected triangles
            for (Eh eh : eas.vector()) { if (!is_deleted(mesh, eh)) for (Fh fh : mesh.ef_range(eh))
            if (fh.is_valid()) if (!is_exterior(mesh, fh)) { fas.push_back(fh); } }

            for (Fh fh : fas.vector()) { if (!is_deleted(mesh, fh)) if (bad_triangle(mesh, fh))
            { events.push({ make_triangle(mesh, fh), calc_priority(mesh, fh) }); } }

            ++n_new_vertices;
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////
/// Wrapping up
////////////////////////////////////////////////////////////////

int refine(TriMesh &mesh, const double min_angle, const double max_length, const double max_area)
{
    int err {};

    Encroachment encroached(min_angle);

    BadTriangle bad_triangle(min_angle, max_length, max_area);

    mark_endians(mesh);

    split_segments(mesh, encroached);

    split_interior(mesh, bad_triangle, encroached);

    mesh.garbage_collection();

    return err;
}