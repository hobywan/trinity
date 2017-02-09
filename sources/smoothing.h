/* ------------------------------------*/ 
#pragma once
/* ------------------------------------*/ 
#include "mesh.h"
#include "hashtable.h"
#include "numeric.h"
#include "partition.h"
/* ------------------------------------ */
namespace trigen {

  class smooth_t {

public:

     smooth_t(mesh_t* input, partit_t* algo, int level);
    ~smooth_t();

    void run(stats_t* tot);

private:

    // steps
    void preprocess();
    void cache_qualit();
    void kernel();

    mesh_t*   mesh;
    partit_t* heuris;
    //
    char* activ;
    double* qualit;

    // counters
    int& cores;
    int& nb_nodes;
    int& nb_elems;
    int& verbose;
    int& iter;
    int& rounds;
    int  depth;
    int  nb_tasks;
    int  nb_comms;

    // kernel
    int laplacian(int id);

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
