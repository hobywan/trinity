/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "mesh.h"
#include "numeric.h"
#include "hessian.h"
/* ------------------------------------ */
namespace trinity {

class Metrics {

public:

   Metrics(Mesh* input, float targ_f, int p_norm, double h_min, double h_max);
  ~Metrics();

  void computeTensorField(Stats* tot);
  void clear();

private:
  // steps
  void computeHessianField();
  void localFieldNormalize();
  void calculComplexity();
  void globalFieldNormalize();
  // kernels
  void calculGradient(int index);
  void calculHessian(int index);

  //
  Mesh* mesh;

  int target;
  int p_norm;
  int chunk;
  int& nb_nodes;
  int& nb_elems;
  int& nb_cores;
  int& verbose;
  int& iter;
  int& rounds;

  double* nabla;
  double* solut;
  double* tens;
  Patch* stenc;

  double h_min;
  double h_max;
  double scale_fact;
  double scale_exp;
  double lambda_min;
  double lambda_max;
  double complexity;

  // timers and stats
  Time start;
  void init();
  void recap(Stats* tot);

};
}
