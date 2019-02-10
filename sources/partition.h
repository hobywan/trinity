#ifndef PARTITION_H
#define PARTITION_H
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
#include "tools.h"
#include "mesh.h"
#include "rmat.h"
/* ------------------------------------ */
namespace trinity {

class Partit {

  friend class Coarse;
  friend class Smooth;

public:

   Partit() = default;
   Partit(int max_size, int max_parts);
  ~Partit();

  void extractIndepSet(const Graph& graph, int nb);
  void extractColoring(const Mesh* mesh);
  void extractPartition(const Mesh* mesh);

  // -------------------------
  // implementation in 'coloring.cpp'
  void processBenchmark(int nb_rounds);
  void colorGraph_Catalyurek(RMAT* graph);
  void colorGraph_Gebremedhin(RMAT* graph);
  void colorGraph_Rokos(RMAT* graph);
  void colorGraph_MonteCarlo(RMAT* graph);
  // -------------------------

private:

  int remain[2];    // nb of tasks per worklist
  int size;         // max storage capacity (>= nb_nodes)
  int max_p;        // max allowed number of parts
  int cores;
  int parts;        // nb of colors/parts
  int defect;       // nb of defective vertices
  int rounds;       // nb of iterations

  //
  int* off;
  int* card;
  int* mapping;
  int** subset;
  int** tasks;     // 2 tasklists

  // kernels
  void reset();

  // TODO use thread-local storage for 'forbidden' and 'conflicts'
  void pseudoColor(const Graph& graph, std::vector<int>& forbidden, int id);
  bool detectErrors(const Graph& graph, std::vector<int>& conflicts, int id);
  void reduceMaxColor(int nb_nodes);
};
}
#endif
