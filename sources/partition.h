/*
 *                          'partition.h'
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
#include "rmat.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
class Partit {

  friend class Coarse;
  friend class Smooth;

public:

  // rule of five
  Partit() = delete;
  Partit(const Partit& other) = delete;
  Partit& operator=(Partit other) = delete;
  Partit(Partit&& other) noexcept = delete;
  Partit& operator=(Partit&& other) noexcept = delete;
  Partit(int max_size, int max_parts);
  ~Partit();

  void extractIndepSet(const Graph& graph, int nb);
  void extractColoring(const Mesh* mesh);
  void extractPartition(const Mesh* mesh);

  // -------------------------
  // implementation in 'coloring.cpp'
  void processBenchmark(int nb_rounds);
  void colorGraph_Catalyurek(RMAT* graph);
  void colorGraph_Gebremedhin(RMAT* graph);
  void colorGraph_Rokos(RMAT* graph);
  // -------------------------

private:

  int remain[2];    // nb of tasks per worklist
  int size;         // max storage capacity (>= nb_nodes)
  int max_p;        // max allowed number of parts
  int cores;
  int parts;        // nb of colors/parts
  int defect;       // nb of defective vertices
  int rounds;       // nb of iterations

  //
  int* off;
  int* card;
  int* mapping;
  int** subset;
  int** tasks;     // 2 tasklists

  // kernels
  void reset();

  // TODO use thread-local storage for 'forbidden' and 'conflicts'
  void pseudoColor(const Graph& graph, std::vector<int>& forbidden, int id);
  bool detectErrors(const Graph& graph, std::vector<int>& conflicts, int id);
  void reduceMaxColor(int nb_nodes);
};
} // namespace trinity
