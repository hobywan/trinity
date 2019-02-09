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

class partit_t {

  friend class coarse_t;
  friend class smooth_t;
  friend class indep_t;

public:

  partit_t();
  partit_t(int max_size, int max_parts);
  ~partit_t();

  void indep_subset(const graph_t& graph, int nb);
  void catalyurek(const mesh_t* mesh);
  void partitioning(const mesh_t* mesh);

  // -------------------------
  // implementation in 'coloring.cpp'
  void process_benchmark(int nb_rounds); // done
  void catalyurek(RMAT* graph);    // done
  void gebremedhin(RMAT* graph);   // done
  void rokos_gorman(RMAT* graph);  // done
  void monte_carlo(RMAT* graph);
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
  void flush();

  // TODO use thread-local storage for 'forbidden' and 'conflicts'
  void pseudo_color(const graph_t& graph, std::vector<int>& forbidden, int id);
  bool detect_error(const graph_t& graph, std::vector<int>& conflicts, int id);
  void reduce_maxcol(int nb_nodes);
};
}
#endif
