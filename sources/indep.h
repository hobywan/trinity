/* ------------------------------------ */
#pragma once
/* ------------------------------------ */
#include "mesh.h"
#include "tools.h"
#include "random_engine.h"
/* ------------------------------------ */
namespace trigen {
  
  class indep_t {

    friend class coarse_t;

public:

     indep_t(mesh_t* input);
    ~indep_t();

    int* luby();
    int* metivier(const graph_t& graph, int nb);

    void flush();
    void verif(const graph_t& graph);

private:

    mesh_t* mesh;

    // todo: one engine per thread

    int cores;
    int remain; 
    int rounds;
    int indep;
    int nb_nodes;

    int* subset;
    double* weight;
    char* activ;
    char* added;
  };
}
