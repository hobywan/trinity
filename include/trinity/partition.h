/*
 *                          'partition.h'
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
#include "rmat.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
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
/* -------------------------------------------------------------------------- */
} // namespace trinity
