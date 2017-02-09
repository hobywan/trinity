/* ------------------------------------*/ 
#pragma once
/* ------------------------------------*/ 
#include "tools.h"
#include "mesh.h"
#include "hashtable.h"
#include "numeric.h"
/* ------------------------------------ */
namespace trigen {

  class refine_t {

public:    

     refine_t(mesh_t* input, int level);
    ~refine_t();

    void run(stats_t* tot);

private:    

    // steps
    void preprocess(std::vector<int>* heap);
    void filtering (std::vector<int>* heap);
    void compute_steiner();
    void kernel(int tid);
    void dissect(int id, int* offset);

    mesh_t* mesh;             // input mesh
    hashtable<int> steiner;   // mapping: vi -> (vj,s)


    int*  index;              // offset for elems insertion
    int*  edges;              // tasklist for steiner point calc. step
    int*  elems;              // tasklist for kernel step
    int*  off;                // offset for tasklist reduction
    char* activ;              // active elems marking
    char* pattern;            // pattern for each elem

    // counters
    int shift;  
    int nb_adds;
    int nb_split;
    int nb_eval;
    int nb_tasks;  
    int nb_stein;
    //
    int& cores;                
    int& nb_nodes;
    int& nb_elems;
    int& verbose;
    int& iter;
    int& rounds;
    int  depth;                // max refinement level
    int  old_node;
    int  old_elem;
    int  cur_elem;
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

