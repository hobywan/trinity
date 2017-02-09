#ifndef MATCHING_H
#define MATCHING_H
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
#include "tools.h"
#include "mesh.h"
/* ------------------------------------ */
namespace trigen {
  
  class match_t {

    friend class swap_t;

public:

     match_t();
    ~match_t();
    
    void init(int capa, int* map, int* idx);

    int* karp_sipser  (const graph_t& graph, int nb);  
    int* pothen_fan   (const graph_t& graph, int nb);
    int* hopcroft_karp(const graph_t& graph, int nb);

private:
    int size;         // max number of nodes (capacity)
    int depth;        // max depth (unused)
    int cores;        // number of cores used
    int card[2];      // tasklists cardinalities

    //
    int*  matched;    // matched vertex pairs
    char* visited;    // flag for DFS
    char* degree;     // active vertex degree 
    int*  mapping;    // mapping: node id -> index in G
    int*  off;        // offset for prefix sum
    int** tasks;      // tasklists (used only for general sparse graphs)

    bool path;

    void flush();
    void match_and_update(int id, const graph_t& graph, std::stack<int>* stack);
    bool DFS_look_ahead  (int id, const graph_t& graph, std::stack<int>* stack);

  };
}  
#endif  
