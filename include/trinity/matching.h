/*
 *                          'matching.h'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *                Copyright 2016, Hoby Rakotoarivelo
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once
/* -------------------------------------------------------------------------- */
#include "header.h"
#include "timer.h"
#include "tools.h"
#include "mesh.h"
/* -------------------------------------------------------------------------- */
#define MATCHED -2
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
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
  void matchAndUpdate(int vertex, const Graph& graph, std::stack<int>* stack);
  bool lookAheadDFS(int vertex, const Graph& graph, std::stack<int>* stack);

  struct {
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
/* -------------------------------------------------------------------------- */
} // namespace trinity
