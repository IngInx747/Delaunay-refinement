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

inline double determinant(const Vec2 &d0, const Vec2 &d1)
{
    return determinant({ 0,0 }, d0, d1);
}

inline int orientation(const Vec2 &a, const Vec2 &b, const Vec2 &c)
{
    return sign(orient2d(a.data(), b.data(), c.data())); // +1:CCW, -1:CW, 0:LINEAR
}

inline int orientation(const Vec2 &a, const Vec2 &b, const Vec2 &c, bool fuzzy)
{
    return sign(orient2dfast(a.data(), b.data(), c.data())); // +1:CCW, -1:CW, 0:LINEAR
}

inline int incircle(const Vec2 &a, const Vec2 &b, const Vec2 &c, const Vec2 &d)
{
    return sign(incircle(a.data(), b.data(), c.data(), d.data())); // +1:IN, -1:OUT, 0:ON
}

inline int incircle(const Vec2 &a, const Vec2 &b, const Vec2 &c, const Vec2 &d, bool fuzzy)
{
    return sign(incirclefast(a.data(), b.data(), c.data(), d.data())); // +1:IN, -1:OUT, 0:ON
}

#endif