/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "mesh.h"
#include "numeric.h"
#include "hessian.h"
/* ------------------------------------ */
namespace trinity {

class metric_t {

public:

  metric_t(mesh_t* input, float targ_f, int p_norm, double h_min, double h_max);
  ~metric_t();

  void multi_scale_field(stats_t* tot);
  void clear();

private:
  // steps
  void compute_hessian_field();
  void local_normalization();
  void compute_complexity();
  void global_normalization();
  // kernels
  void calcul_gradient(int index);
  void calcul_hessian(int index);

  //
  mesh_t* mesh;

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
  patch_t* stenc;

  double h_min;
  double h_max;
  double scale_fact;
  double scale_exp;
  double lambda_min;
  double lambda_max;
  double complexity;

  // timers and stats
  time_t start;
  void init();
  void recap(stats_t* tot);

};
}
