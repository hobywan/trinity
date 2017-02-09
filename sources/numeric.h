#ifndef NUMERIC_H
#define NUMERIC_H
/* ------------------------------------ */
#include "header.h"
/* ------------------------------------ */
namespace trigen {
  namespace numeric {
    
    void tensor_eigen_decomp(const double* M, double* val, double* vec1, double* vec2);
    void tensor_interpolate(const double* M, double* R, int n);
    void kronecker_product(const double* u1, const double* u2, double* M);

    double quality(const double* pa, const double* pb, const double* pc,
                   const double* Ma, const double* Mb, const double* Mc);
    void steiner_point(const double* pa, const double* pb,
                       const double* Ma, const double* Mb, double* p, double* M);      

    double riemannian_distance(const double* pa, const double* pb, const double* M);
    double riemannian_distance(const double* pa, const double* pb,
                               const double* Ma, const double* Mb); 
    void riemannian_circum(const double* pa, const double* pb, const double* pc,
                           const double* M, double* O);
  }
}

#endif
