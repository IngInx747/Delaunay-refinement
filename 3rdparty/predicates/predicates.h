#ifndef PREDICATES_H
#define PREDICATES_H

#ifdef __cplusplus
extern "C" {
#endif

/// Initialize the variables used for exact arithmetic.
void exactinit();

/// Return a positive value if the points pa, pb, and pc occur in
/// counterclockwise order; a negative value if they occur in clockwise order;
/// and zero if they are collinear.  The result is also a rough approximation
/// of twice the signed area of the triangle defined by the three points.
double orient2d(const double*, const double*, const double*);

double orient2dfast(const double*, const double*, const double*);

/// Return a positive value if the point pd lies below the plane passing
/// through pa, pb, and pc; "below" is defined so that pa, pb, and pc appear
/// in counterclockwise order when viewed from above the plane; a negative
/// value if pd lies above the plane; and zero if the points are coplanar.
/// The result is also a rough approximation of six times the signed volume
/// of the tetrahedron defined by the four points.
double orient3d(const double*, const double*, const double*, const double*);

double orient3dfast(const double*, const double*, const double*, const double*);

/// Return a positive value if the point pd lies inside the circle passing
/// through pa, pb, and pc; a negative value if it lies outside; and zero if
/// the four points are cocircular.  The points pa, pb, and pc must be in
/// counterclockwise order, or the sign of the result will be reversed.
double incircle(const double*, const double*, const double*, const double*);

double incirclefast(const double*, const double*, const double*, const double*);

/// Return a positive value if the point pe lies inside the sphere passing
/// through pa, pb, pc, and pd; a negative value if it lies outside; and zero
/// if the five points are cospherical.  The points pa, pb, pc, and pd must be
/// ordered so that they have a positive orientation (as defined by orient3d())
/// or the sign of the result will be reversed.
double insphere(const double*, const double*, const double*, const double*, const double*);

double inspherefast(const double*, const double*, const double*, const double*, const double*);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // PREDICATES_H