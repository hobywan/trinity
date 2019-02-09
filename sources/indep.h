/* ------------------------------------ */
#pragma once
/* ------------------------------------ */
#include "mesh.h"
#include "tools.h"
#include "random_engine.h"
#include "partition.h"
/* ------------------------------------ */
// TODO: one engine per thread
namespace trinity {

class indep_t {

  friend class coarse_t;

public:

  indep_t(mesh_t* input, partit_t* algo);
  ~indep_t();

  void luby();
  void metivier(const graph_t& graph, int nb);
  void flush();
  void verif(const graph_t& graph);

private:

  mesh_t* mesh;
  int cores;
  int remain;
  int rounds;
  int nb_nodes;
  int& nb_indep;
  int* subset;
  char* activ;
  char* added;
  double* weight;
};
} // namespace trinity
