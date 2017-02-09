/* ------------------------------------*/ 
#pragma once
/* ------------------------------------*/ 
#include "mesh.h"
#include "hashtable.h"
#include "numeric.h"
#include "matching.h"
/* ------------------------------------ */
namespace trigen {

  class swap_t {

public:    

     swap_t(mesh_t* input);
    ~swap_t();

    void run(stats_t* tot);

private:
    // steps
    void cache_qualit();
    void filter(std::vector<int>* heap);
    void extract_dual();
    void kernel();

    // 
    mesh_t* mesh;
    graph_t dual;
    match_t heuris;

    // tasklist
    int*  map;
    int*  tasks;
    int*  off;
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
    int  depth;
    int  nb_activ;
    int  nb_tasks;
    int  nb_comms;
    
    // kernel
    int swap(int k1, int k2, int idx);
        
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
