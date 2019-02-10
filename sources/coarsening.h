/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "mesh.h"
#include "sync.h"
#include "numeric.h"
#include "partition.h"
/* ------------------------------------ */
namespace trinity {

class Coarse {

public:

   Coarse() = default;
   Coarse(Mesh* input, Partit* algo);
  ~Coarse();

  void run(Stats* tot);

private:
  // steps
  void preProcess();
  void filterPoints(std::vector<int>* heap);
  void extractPrimalGraph();
  void processPoints();

  //
  Mesh* mesh;
  Graph primal;
  Partit* heuris;

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

  // timers and stats
  Time start;
  Time round;
  Time tic;
  //
  void init();
  void saveStat(int level, int* stat, int* form);
  void showStat(int level, int* form);
  void recap(int* time, int* stat, int* form, Stats* tot);
};
}
