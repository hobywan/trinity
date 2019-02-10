#ifndef NUMERIC_H
#define NUMERIC_H
/* ------------------------------------ */
#include "header.h"
/* ------------------------------------ */
namespace trinity { namespace numeric {

void eigenDecomposeTensor(const double* M, double* val, double* vec1, double* vec2);
void interpolateTensor(const double* M, double* R, int n);
void kroneckerProduct(const double* u1, const double* u2, double* M);

double computeQuality(const double* pa, const double* pb, const double* pc,
                      const double* Ma, const double* Mb, const double* Mc);
void computeSteinerPoint(const double* pa, const double* pb,
                         const double* Ma, const double* Mb, double* p, double* M);

double approxRiemannDist(const double* pa, const double* pb, const double* M);
double approxRiemannDist(const double* pa, const double* pb,
                         const double* Ma, const double* Mb);
void approxRiemannCircum(const double* pa, const double* pb, const double* pc,
                         const double* m, double* p);

}}

#endif
