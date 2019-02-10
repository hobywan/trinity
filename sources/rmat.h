#pragma once
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
#include "tools.h"
/* ------------------------------------ */
namespace trinity {

struct RMAT {

  Time  start;
  int   nb_nodes;
  int   nb_edges;
  int   nb_rounds;
  int   nb_error;
  int   nb_color;
  int   deg_max;
  int   deg_avg;
  float ratio;

  // adjacency list
  Graph graph;

   RMAT();
  ~RMAT();

  // utils
  void reset();
  void load(std::string path);
  void info(std::string name);
  void saveChrono();
  int elapsed();
};
}
