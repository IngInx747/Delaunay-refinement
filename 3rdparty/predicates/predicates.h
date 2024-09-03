#ifndef PREDICATES_H
#define PREDICATES_H

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Initializes the exact predicates.
/// This function has to be called before calling ::orient2d.
void exactinit();

/// @brief Computes the orientation of the supplied points.
/// This function only returns correct values if exactinit() was called before.
/// @param pa, pb, pc Arrays containing the two coordinates of a point, each.
/// @return A value the sign of which reflects the orientation of the three points.
double orient2d(const double*, const double*, const double*);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // PREDICATES_H