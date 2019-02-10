/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "mesh.h"
#include "Table.h"
#include "numeric.h"
#include "partition.h"
/* ------------------------------------ */
namespace trinity {

  class Smooth {

public:

     Smooth(Mesh* input, Partit* algo, int level);
    ~Smooth();

    void run(Stats* tot);

private:

    // steps
    void preProcess();
    void cacheQuality();
    void movePoints();

    Mesh*   mesh;
    Partit* heuris;
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

    // processFlips
    int moveLaplacian(int id);

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
