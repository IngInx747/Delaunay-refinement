#ifndef GEOMETRY_PREDICATE2D_HH
#define GEOMETRY_PREDICATE2D_HH

#include <predicates.h>
#include "math_def.hh"
#include "vector_n.hh"

enum class ORT_2D : int { CCW = 1, CW = -1, LINEAR = 0 };

inline double determinant(const Vec2 &a, const Vec2 &b, const Vec2 &c)
{
    return orient2d(a.data(), b.data(), c.data());
}

inline int orientation(const Vec2 &a, const Vec2 &b, const Vec2 &c)
{
    return sign(orient2d(a.data(), b.data(), c.data())); // -1:CW, +1:CCW, 0:LINEAR
}

inline double determinant(const Vec2 &d0, const Vec2 &d1)
{
    return determinant({ 0,0 }, d0, d1);
}

#endif