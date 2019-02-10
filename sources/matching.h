#pragma once
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
#include "tools.h"
#include "mesh.h"
/* ------------------------------------ */
#define MATCHED -2
/* ------------------------------------ */
namespace trinity {

class Match {

  friend class Swap;

public:

   Match();
  ~Match();

  void init(int capa, int* map, int* idx);
  int* computeKarpSipser(const Graph& graph, int nb);
  int* computePothenFan(const Graph& graph, int nb);
  int  getRatio(const Graph& graph, int nb, int* count);

private:
  int size;         // max number of nodes (capacity)
  int depth;        // max depth (unused)
  int cores;        // number of cores used
  int card[2];      // tasklists cardinalities

  //
  int* matched;    // matched vertex pairs
  char* visited;    // flag for DFS
  char* degree;     // active vertex degree
  int* mapping;    // mapping: node id -> index in G
  int* off;        // offset for prefix sum
  int** tasks;      // tasklists (used only for general sparse trinity)

  bool path;

  void flush();
  void matchAndUpdate(int id, const Graph& graph, std::stack<int>* stack);
  bool DFS_lookAhead(int id, const Graph& graph, std::stack<int>* stack);

};
} // namespace trinity
