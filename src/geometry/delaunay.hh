#ifndef EUCLIDEAN_DELAUNAY_HH
#define EUCLIDEAN_DELAUNAY_HH

#include "math_def.hh"
#include "vector_n.hh"
#include "pred2D.hh"

//   1  
//  / \ 
// 2---0
//  \ / 
//   3  
inline bool fuzzy_delaunay(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u3)
{
    return incircle(u0, u1, u2, u3, true) <= 0;
}

//   1  
//  / \ 
// 2---0
//  \ / 
//   3  
inline bool exact_delaunay(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u3)
{
    return incircle(u0, u1, u2, u3) <= 0; // equivalent to incircle(u2, u3, u0, u1)
}

//   1  
//  / \ 
// 2---0
//  \ / 
//   3  
inline bool fuzzy_convex(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u3)
{
    constexpr double kEps = 1e-10;
    const auto d1 = normalize(u0 - u1);
    const auto d2 = normalize(u2 - u1);
    const auto d3 = normalize(u1 - u3);
    const double c0 = cross(d3, d1);
    const double c1 = cross(d3, d2);
    return abs(c0)<kEps ? false : abs(c1)<kEps ? false : c0*c1<0;
}

//
//   1   
//  / \  
// 2---0 
//  \ /  
//   3   
//
inline bool exact_convex(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2, const Vec2 &u3)
{
    const int r0 = orientation(u3, u0, u1);
    const int r1 = orientation(u0, u1, u2);
    const int r2 = orientation(u1, u2, u3);
    const int r3 = orientation(u2, u3, u0);
    return r0 > 0 && r1 > 0 && r2 > 0 && r3 > 0 ||
           r0 < 0 && r1 < 0 && r2 < 0 && r3 < 0 ;
}

inline Vec2 circumcenter(const Vec2 &u0, const Vec2 &u1, const Vec2 &u2)
{
    const auto d0 = u1 - u0;
    const auto d1 = u2 - u0;
    const double a = cross(d0, d1);
    const double d00 = dot(d0, d0);
    const double d11 = dot(d1, d1);
    const Vec2 uc { d1[1]*d00 - d0[1]*d11, d0[0]*d11 - d1[0]*d00 };
    return uc / (a*2.) + u0;
}

#endif