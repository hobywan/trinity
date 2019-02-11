/*
 *                          'matching.h'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *               Copyright (c) 2016 Hoby Rakotoarivelo.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
/* --------------------------------------------------------------------------- */
#include "header.h"
#include "timer.h"
#include "tools.h"
#include "mesh.h"
/* --------------------------------------------------------------------------- */
#define MATCHED -2
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
class Match {

  friend class Swap;

public:

  // rule of five
  Match();
  Match(const Match& other) = delete;
  Match& operator=(Match other) = delete;
  Match(Match&& other) noexcept = delete;
  Match& operator=(Match&& other) noexcept = delete;
  ~Match();


  void initialize(size_t capa, int* mapping, int* index);
  int* computeGreedyMatching(const Graph& graph, int nb);
  int* localSearchBipartite(const Graph& graph, int nb);
  int  getRatio(const Graph& graph, int nb, int* count);

private:

  // kernels
  void reset();
  void matchAndUpdate(int id, const Graph& graph, std::stack<int>* stack);
  bool lookAheadDFS(int id, const Graph& graph, std::stack<int>* stack);

  struct {
    int  cores;        // number of cores used
    bool found;
  } param;

  struct {
    int  capa;         // max number of nodes (capacity)
    int  depth;        // max depth
  } max;

  struct {
    int*  matched;    // matched vertex pairs
    int** lists;      // tasklists (used only for general sparse trinity)
    int*  mapping;    // mapping: node id -> index in G
    int   cardin[2];      // tasklists cardinalities
  } task;

  struct {
    char* visited;    // flag for DFS
    char* degree;     // active vertex degree
    int*  off;        // offset for prefix sum
  } sync;
};
/* --------------------------------------------------------------------------- */
} // namespace trinity
