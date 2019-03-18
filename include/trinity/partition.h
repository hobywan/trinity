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
  ~Partit();

  Partit(int max_graph_size, int max_part_size);
  Partit(Mesh const* mesh, int max_part_size);

  void extractIndepSet(const Graph& graph, int nb);
  void extractColoring(const Mesh* mesh);
  void extractPartition(const Mesh* mesh);

private:

  void reset();

  struct {
    int capa;         // max storage capacity (>= nb_nodes)
    int part;         // max allowed number of parts
  } max;

  struct {
    int*  cardin;
    int*  mapping;
    int** subset;
    int** lists;     // 2 tasklists
  } task;

  struct {
    int remain[2];    // nb of tasks per worklist
    int cores;
    int parts;        // nb of colors/parts
    int defect;       // nb of defective vertices
    int rounds;       // nb of iterations
  } nb;

  struct { int* offset; } sync;
};
/* --------------------------------------------------------------------------- */
} // namespace trinity
