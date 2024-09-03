#ifndef EUCLIDEAN_SEGMENT_HH
#define EUCLIDEAN_SEGMENT_HH

#include "math_def.hh"
#include "vector_n.hh"
#include "pred2D.hh"

////////////////////////////////////////////////////////////////
/// Segments intersection
////////////////////////////////////////////////////////////////

inline Int4 intersection_info(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
{
    const int ru0 = orientation(v0, v1, u0);
    const int ru1 = orientation(v0, v1, u1);
    const int rv0 = orientation(u0, u1, v0);
    const int rv1 = orientation(u0, u1, v1);
    return { ru0, ru1, rv0, rv1 };
}

inline Vec2 intersection_param(const Vec2 &u0, const Vec2 &u1, const Vec2 &v0, const Vec2 &v1)
{
    const double dt = determinant(u1 - u0, v1 - v0);
    const double t0 = determinant(v0 - u0, v1 - v0);
    const double t1 = determinant(v0 - u0, u1 - u0);
    return Vec2 { t0, t1 } / dt;
}

////////////////////////////////////////////////////////////////
/// Locations
////////////////////////////////////////////////////////////////

//
//  O0      IN      O1 
// ---- V0 ---- V1 ----
//
enum class SEG_LOC : int
{
    IN,
    V0,
    V1,
    O0,
    O1
};

inline SEG_LOC locate(const Vec2 &u0, const Vec2 &u1, const Vec2 &u)
{
    const auto d0 = u - u0;
    const auto d1 = u - u1;
    const auto e = u1 - u0;
    const double w0 = dot(d0, e);
    const double w1 = dot(d1,-e);

    return
        w0 < 0  ? SEG_LOC::O0 :
        w1 < 0  ? SEG_LOC::O1 :
        w0 == 0 ? SEG_LOC::V0 :
        w1 == 0 ? SEG_LOC::V1 :
                  SEG_LOC::IN ;
}

#endif