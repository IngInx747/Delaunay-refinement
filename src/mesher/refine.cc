#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "delaunay.hh"
#include "triangle.hh"
#include "segment.hh"
#include "mesh.hh"
#include "topology.hh"

using namespace OpenMesh;

////////////////////////////////////////////////////////////////
/// Locations
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

////////////////////////////////////////////////////////////////
/// Local search
////////////////////////////////////////////////////////////////

enum { NO_ITSC, INTO_EDGE, ON_VERTEX };

static inline int intersection_info(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Hh hh)
{
    const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hh));

    const auto ii = intersection_info(u0, u1, v0, v1);
    const int ru0 = ii[0];
    const int ru1 = ii[1];
    const int rv0 = ii[2];
    const int rv1 = ii[3];

    return
        (ru0 * ru1 < 0) && (rv0 * rv1 < 0)        ? INTO_EDGE : // intersecting
        (ru0 * ru1 < 0) && (rv1 == 0 && rv0 != 0) ? ON_VERTEX : // v1 lies on (u0,u1)
        (ru0 != 0 && ru1 == 0) && (rv1 == 0)      ? ON_VERTEX : // v1 overlaps u1
        NO_ITSC; // no intersecting, overlapping, v0 lies on (u0,u1), and other cases
}

static inline double intersection_param(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Hh hh)
{
    const auto v0 = get_xy(mesh, mesh.from_vertex_handle(hh));
    const auto v1 = get_xy(mesh, mesh.to_vertex_handle  (hh));
    return intersection_param(u0, u1, v0, v1)[0]; // t in eq: u0 + (u1-u0)*t
}

static inline int next_primitive(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Hh hho, Hh &hhc, Vh &vhc)
{
    int res { NO_ITSC };

    Hh hh0 = mesh.prev_halfedge_handle(hho);
    Hh hh1 = mesh.next_halfedge_handle(hho);

    for (Hh hh : { hh0, hh1 })
    {
        const int ii = intersection_info(mesh, u0, u1, hh);
        if (ii == NO_ITSC) continue;

        res = ii;

        if (ii == INTO_EDGE)
        {
            hhc = mesh.opposite_halfedge_handle(hh);
        }
        else if (ii == ON_VERTEX)
        {
            vhc = mesh.to_vertex_handle(hh);
        }
    }

    return res;
}

static inline int next_primitive(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Vh vho, Hh &hhc, Vh &vhc)
{
    int res { NO_ITSC };
    double tmax {}; // the parameter in line equation: u0 + (u1-u0)*t

    for (Hh hh : mesh.voh_range(vho)) if (!mesh.is_boundary(hh))
    {
        hh = mesh.next_halfedge_handle(hh); // apex edge of v0

        const int ii = intersection_info(mesh, u0, u1, hh);
        if (ii == NO_ITSC) continue;

        const double t = intersection_param(mesh, u0, u1, hh);
        if (tmax >= t) continue; // intersected but earlier

        tmax = t;
        res = ii;

        if (ii == INTO_EDGE)
        {
            hhc = mesh.opposite_halfedge_handle(hh);
        }
        else if (ii == ON_VERTEX)
        {
            vhc = mesh.to_vertex_handle(hh);
        }
    }

    return res;
}

static inline int first_primitive(const TriMesh &mesh, const Vec2 &u0, const Vec2 &u1, Fh fh, Hh &hhc, Vh &vhc)
{
    int res  { NO_ITSC };
    double tmin { 1e20 }; // the parameter in line equation: u0 + (u1-u0)*t

    for (Hh hh : mesh.fh_range(fh))
    {
        const int ii = intersection_info(mesh, u0, u1, hh);
        if (ii == NO_ITSC) continue;

        const double t = intersection_param(mesh, u0, u1, hh);
        if (tmin < t) continue;

        tmin = t; // ideally t_min < 0
        res = ii;

        if (ii == INTO_EDGE)
        {
            hhc = hh;
        }
        else if (ii == ON_VERTEX)
        {
            vhc = mesh.to_vertex_handle(hh);
        }
    }

    return res;
}

static Fh search_primitive(const TriMesh &mesh, const Vec2 &u1, const Fh &fh0, std::vector<Eh> &ehs)
{
    Vh vhc {}; Hh hhc {};

    if (is_inside(mesh, fh0, u1)) return fh0;

    const auto u0 = centroid(mesh, fh0);

    int status = first_primitive(mesh, u0, u1, fh0, hhc, vhc);

    const int max_n_iter = (int)mesh.n_edges(); int n_iter {};

    for (n_iter = 0; n_iter < max_n_iter; ++n_iter)
    {
        if (status == INTO_EDGE) // go to the next primitive
        {
            status = next_primitive(mesh, u0, u1, hhc, hhc, vhc);
        }
        else if (status == ON_VERTEX)
        {
            status = next_primitive(mesh, u0, u1, vhc, hhc, vhc);
        }
        else // searching lost in vain, for some reasons
        {
            break;
        }

        if (status == INTO_EDGE) // record encroached segments
        {
            Eh eh = mesh.edge_handle(hhc); if (is_sharp(mesh, eh))
            {
                ehs.push_back(eh);
            }
        }
        else if (status == ON_VERTEX)
        {
            for (Eh eh : mesh.ve_range(vhc)) if (is_sharp(mesh, eh))
            {
                ehs.push_back(eh);
            }
        }

        Fh fh {}; // check if the target is reached

        if (status == INTO_EDGE)
        {
            fh = mesh.face_handle(hhc);
        }
        else if (status == ON_VERTEX)
        {
            fh = mesh.face_handle(mesh.halfedge_handle(vhc));
        }

        if (fh.is_valid() && is_inside(mesh, fh, u1))
        {
            return fh;
        }
    }

    return Fh {};
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

    inline bool operator()(const TriMesh &mesh, const Hh &hh) const
    { return !is_exterior(mesh, hh) && is_encroached(mesh, hh, cs2); } // If true, split it

    const double cs2; // (cos(pi - 2t))^2
};

////////////////////////////////////////////////////////////////
/// Triangle quality
////////////////////////////////////////////////////////////////

struct BadTriangle
{
    inline bool operator()(const TriMesh &mesh, const Hh &hh) const
    { return false; } // If true, split it
};

////////////////////////////////////////////////////////////////
/// Refinable events
////////////////////////////////////////////////////////////////

struct Primitive // either a triangle or a segment
{
    Hh hh;
    Vh vh0, vh1, vh2;
};

template<>
struct std::hash<Primitive>
{
    size_t operator()(const Primitive &primitive) const noexcept
    {
        size_t h0 = hash<Hh>{}(primitive.hh);
        size_t h1 = hash<Vh>{}(primitive.vh0);
        size_t h2 = hash<Vh>{}(primitive.vh1);
        size_t h3 = hash<Vh>{}(primitive.vh2);
        return ((((h0 ^ h1) << 1) ^ h2) << 1) ^ h3;
    }
};

template<>
struct std::equal_to<Primitive>
{
    bool operator()(const Primitive &lhs, const Primitive &rhs) const noexcept
    {
        return
            lhs.hh  == rhs.hh  &&
            lhs.vh0 == rhs.vh0 &&
            lhs.vh1 == rhs.vh1 &&
            lhs.vh2 == rhs.vh2;
        }
};

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

    std::queue<Primitive> primitives; // to split

    // Before enqueueing encroached segments, free vertices in
    // the diametral circle of the segment should be deleted.
    // [TODO]

    for (Eh eh : mesh.edges()) if (is_sharp(mesh, eh))
    for (Hh hh : mesh.eh_range(eh)) if (encroached(mesh, hh))
    { primitives.push(make_segment(mesh, hh)); }

    while (!primitives.empty())
    {
        auto primitive = primitives.front(); primitives.pop();

        if (is_segment(primitive) && is_segment_valid(mesh, primitive))
        {
            // find a position for splitting
            const auto u = splitting_position(mesh, primitive.hh);

            // split the segment
            Vh vh = mesh.new_vertex({ u[0], u[1], 0 });
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

static int split_interior(TriMesh &mesh, const BadTriangle &bad_triangle, const Encroachment &encroached)
{
    const int max_n_iter = 100;

    auto delaunifier = make_delaunifier(mesh, EuclideanDelaunay {});

    return 0;
}

////////////////////////////////////////////////////////////////
/// Wrapping up
////////////////////////////////////////////////////////////////

int refine(TriMesh &mesh, const double min_angle)
{
    Encroachment encroached(min_angle);

    BadTriangle bad_triangle;

    int err;

    err = split_segments(mesh, encroached);

    err = split_interior(mesh, bad_triangle, encroached);

    return err;
}