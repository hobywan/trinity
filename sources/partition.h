#ifndef PARTITION_H
#define PARTITION_H
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
#include "tools.h"
#include "mesh.h"
/* ------------------------------------ */
namespace trigen {

  class partit_t {

    friend class coarse_t;
    friend class smooth_t;

public:    


     partit_t();
     partit_t(int max_size, int max_parts);
    ~partit_t();

    void indep_subset(const graph_t& graph, int nb);
    void catalyurek  (const mesh_t* mesh);
    void partitioning(const mesh_t* mesh);

private:
    
    int remain[2];    // nb of tasks per worklist
    int size;         // max storage capacity (>= nb_nodes)
    int max_p;        // max allowed number of parts
    int cores;
    int parts;        // nb of colors/parts
    int defect;       // nb of defective vertices
    int rounds;       // nb of iterations
        
    //
    int*  off;
    int*  card;
    int*  mapping;
    int** subset;
    int** tasks;     // 2 tasklists

    // kernels
    void flush();
  };      
}
#endif
