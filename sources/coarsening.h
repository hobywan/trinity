/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "mesh.h"
#include "sync.h"
#include "numeric.h"
#include "partition.h"
#include "indep.h"
/* ------------------------------------ */
namespace trinity {

class coarse_t {

public:

  coarse_t(mesh_t* input, partit_t* algo);
  ~coarse_t();

  void run(stats_t* tot);

private:
  // steps
  void preprocess();
  void filtering(std::vector<int>* heap);
  void extract_primal();
  void kernel();

  //
  mesh_t* mesh;
  graph_t primal;
  partit_t* heuris;

  // tasklists
  int* target;
  int* filter;
  int* indep;
  int* off;
  char* fixes;
  char* activ;

  // counters
  int& cores;
  int& nb_nodes;
  int& nb_elems;
  int& nb_indep;
  int& verbose;
  int& iter;
  int& rounds;
  int depth;
  int nb_activ;
  int nb_tasks;

  // kernels
  void identify(int id);
  void collapse(int i, int j);

  // for profiling only
  indep_t* alg;

  // timers and stats
  time_t start;
  time_t round;
  time_t tic;
  //
  void init();
  void save_stat(int level, int* stat, int* form);
  void show_stat(int level, int* form);
  void recap(int* time, int* stat, int* form, stats_t* tot);
};
}
