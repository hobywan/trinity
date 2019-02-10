/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "tools.h"
#include "mesh.h"
#include "Table.h"
#include "numeric.h"
/* ------------------------------------ */
namespace trinity {

  class Refine {

public:

     Refine(Mesh* input, int level);
    ~Refine();

    void run(Stats* tot);

private:

    // steps
    void preProcess(std::vector<int>* heap);
    void filterElems(std::vector<int>* heap);
    void computeSteinerPoints();
    void processElems(int tid);
    void cutElem(int id, int* offset);

    Mesh* mesh;             // input mesh
    Table<int> steiner;   // mapping: vi -> (vj,s)


    int*  index;              // offset for elems insertion
    int*  edges;              // tasklist for steiner point calc. step
    int*  elems;              // tasklist for processFlips step
    int*  off;                // offset for tasklist reduction
    char* activ;              // active elems marking
    char* pattern;            // pattern for each elem

    // counters
    int shift;
    int nb_adds;
    int nb_split;
    int nb_eval;
    int nb_tasks;
    int nb_steiner;
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

