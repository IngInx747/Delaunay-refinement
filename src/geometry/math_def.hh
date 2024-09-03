#ifndef EUCLIDEAN_MATH_DEF_HH
#define EUCLIDEAN_MATH_DEF_HH

#include <cmath>

template <typename T> constexpr T pi()
{
    return (T)(3.14159265358979323846264338327950288);
}

template <typename T> inline int sign(T val)
{
    return (T(0) < val) - (val < T(0)); // <0: -1, >0: +1, =0: 0
}

template <typename T> inline T lerp(const T &u0, const T &u1, const double t)
{
    return u0*(1-t) + u1*t;
}

template <typename T> inline T round(T val, T mtp)
{
    return mtp * std::round(val / mtp);
}

template <typename T> inline T radian(T deg)
{
    return (T)(deg * pi<T>() / (T)(180.));
}

template <typename T> inline T degree(T rad)
{
    return (T)(rad * (T)(180.) / pi<T>());
}

inline double cosine(double a, double b, double c)
{
    const double cs = (a*a + b*b - c*c)/(a*b*2);
    return (cs < -1) ? -1 : (cs > 1) ? 1 : cs;
}

inline double cos2cot(double cs)
{
    return cs / sqrt(1.0 - cs * cs);
}

#endif