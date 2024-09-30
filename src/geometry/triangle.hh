#ifndef EUCLIDEAN_TRIANGULAR_HH
#define EUCLIDEAN_TRIANGULAR_HH

#include "math_def.hh"
#include "vector_n.hh"
#include "pred2D.hh"

////////////////////////////////////////////////////////////////
/// Barycentric coordinates
////////////////////////////////////////////////////////////////

inline Vec3 barycentric_coordinate(const Vec2 &a, const Vec2 &b, const Vec2 &c, const Vec2 &p)
{
    const auto ab = b - a;
    const auto bc = c - b;
    const auto ca = a - c;
    const auto ap = p - a;
    const auto bp = p - b;
    const auto cp = p - c;
    const double ma = cross(bc, cp);
    const double mb = cross(ca, ap);
    const double mc = cross(ab, bp);
    const Vec3 phi { ma, mb, mc };
    return phi / (ma + mb + mc);
}

////////////////////////////////////////////////////////////////
/// Locations
////////////////////////////////////////////////////////////////

//
//     V2
//    /  \
//  E1    E0
//  /      \
// V0--E2--V1
//
enum class TRI_LOC : int
{
    OUT = 0x00,
    IN  = 0x01,
    E0  = 0x02,
    E1  = 0x04,
    E2  = 0x08,
    V0  = 0x10,
    V1  = 0x20,
    V2  = 0x40,
    ES  = E0 | E1 | E2,
    VS  = V0 | V1 | V2,
};

inline TRI_LOC fuzzy_locate(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u)
{
    constexpr double kEps = 1e-3;
    const auto bc = barycentric_coordinate(u0, u1, u2, u);

    return
        fabs(bc[1]) < kEps && fabs(bc[2]) < kEps ? TRI_LOC::V0 :
        fabs(bc[2]) < kEps && fabs(bc[0]) < kEps ? TRI_LOC::V1 :
        fabs(bc[0]) < kEps && fabs(bc[1]) < kEps ? TRI_LOC::V2 :
        fabs(bc[0]) < kEps ? TRI_LOC::E0 :
        fabs(bc[1]) < kEps ? TRI_LOC::E1 :
        fabs(bc[2]) < kEps ? TRI_LOC::E2 :
        bc[0]>0 && bc[1]>0 && bc[2]>0 ? TRI_LOC::IN :
        TRI_LOC::OUT ;
}

inline TRI_LOC exact_locate(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u)
{
    const int r0 = orientation(u1, u2, u);
    const int r1 = orientation(u2, u0, u);
    const int r2 = orientation(u0, u1, u);

    return
        r0 == r1 && r0 == r2 ? TRI_LOC::IN : // in triangle
        r0 == 0  && r1 == r2 ? TRI_LOC::E0 : // on edge (u1, u2)
        r1 == 0  && r2 == r0 ? TRI_LOC::E1 : // on edge (u2, u0)
        r2 == 0  && r0 == r1 ? TRI_LOC::E2 : // on edge (u0, u1)
        r0 == 0  && r1 == 0  ? TRI_LOC::V2 : // on vert u2
        r1 == 0  && r2 == 0  ? TRI_LOC::V0 : // on vert u0
        r2 == 0  && r0 == 0  ? TRI_LOC::V1 : // on vert u1
        TRI_LOC::OUT ;
}

#endif