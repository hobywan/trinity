/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "mesh.h"
#include "Table.h"
#include "numeric.h"
#include "matching.h"
/* ------------------------------------ */
namespace trinity {

class Swap {

public:

   Swap(Mesh* input);
  ~Swap();

  void run(Stats* tot);

private:
  // steps
  void cacheQuality();
  void filterElems(std::vector<int>* heap);
  void extractDualGraph();
  void processFlips();

  //
  Mesh* mesh;
  Graph dual;
  Match heuris;

  // tasklist
  int* map;
  int* tasks;
  int* off;
  char* fixes;
  char* activ;
  double* qualit;

  // counters
  int& cores;
  int& nb_nodes;
  int& nb_elems;
  int& verbose;
  int& iter;
  int& rounds;
  int depth;
  int nb_activ;
  int nb_tasks;
  int nb_comms;
  int count;  // for profiling only

  // processFlips
  int swap(int k1, int k2, int idx);

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
