/*
 *                          'mesh.h'
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
#include "sync.h"
#include "numeric.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
class Mesh {

  friend class Metrics;
  friend class Refine;
  friend class Coarse;
  friend class Swap;
  friend class Smooth;
  friend class Partit;

public:

  // rule of five
  Mesh() = delete;
  Mesh(const Mesh& other) = delete;
  Mesh& operator=(Mesh other) = delete;
  Mesh(Mesh&& other) noexcept = delete;
  Mesh& operator=(Mesh&& other) noexcept = delete;
  Mesh(int size[2], int bucket, int depth, int verbosity, int rounds);
  ~Mesh();

  // global routines
  void reallocMemory();
  void doFirstTouch();
  void compressStorage();  // TODO
  void reorderCells();     // TODO
  void extractPrimalGraph();
  void extractDualGraph(Graph* dual) const;

  // repair topology
  void rebuildTopology();
  void fixTagged();
  void fixAll();
  void initActivElems();
  bool verifyTopology() const;

  // I/O
  void load(const std::string& path, const std::string& solu);
  void store(const std::string& path) const;
  void storePrimalGraph(const std::string& path) const;

  // topological queries
  Patch getVicinity(int id, int deg) const;
  const int* getElem(int i) const;
  const int* getElemCoord(int id, double* p) const;
  int getElemNeigh(int id, int i, int j) const;
  bool isActiveNode(int i) const;
  bool isActiveElem(int i) const;
  bool isBoundary(int i) const;
  bool isCorner(int i) const;
  int getCapaNode() const;
  int getCapaElem() const;

  // topological updates
  void replaceElem(int id, const int* v);
  void eraseElem(int id);
  void updateStencil(int i, int t);
  void updateStencil(int i, const std::initializer_list<int>& t);
  void copyStencil(int i, int j, int nb_rm);

  // deferred updates (for perfs comparison)
#if DEFER_UPDATES
  void initUpdates();
  void commitUpdates();
  void resetUpdates();
  void deferredAppend(int tid, int i, int t);
  void deferredAppend(int tid, int i, const std::initializer_list<int>& t);
  void deferredRemove(int tid, int i, int t);
#endif

  // geometrical queries
  double computeLength(int i, int j) const;
  double computeQuality(int i) const;
  double computeQuality(const int* t) const;
  bool isCounterclockwise(const int* elem) const;
  void computeSteinerPoint(int edge_i, int edge_j, double* point, double* metric) const;
  void computeQuality(double quality[3]);

private:

  struct {
    int nodes;                        // nb of active + deleted nodes
    int elems;                        // nb of active + deleted elements
    int cores;                        // nb of logical cores used
    int activ_elem;                   // actual number of elements
  } nb;

  struct {
    std::vector<int> elems;           // elem vertices, stride=3
    Graph stenc;                      // nodal incident elements
    Graph vicin;                      // nodal adjacent nodes
  } topo;                             // could use a set but performance suffers

  struct {
    std::vector<double> points;       // nodal coordinates, stride=2
    std::vector<double> tensor;       // nodal metric tensor, stride=3
    std::vector<double> solut;        // solution field for metric computation
    std::vector<double> qualit;       // cache elem quality to speedup kernels
  } geom;

  struct {
    int depth;                        // refine|smooth max depth
    int verb;                         // verbosity level
    int iter;                         // current iteration
    int rounds;                       // remesh rounds
  } param;

  struct {
    int    scale;                     // capacity scale factor
    size_t bucket;                    // bucket capacity
    size_t node;                      // node capacity
    size_t elem;                      // elem capacity
  } capa;

  struct {
    int* deg      = nullptr;          // nodal degree
    int* off      = nullptr;          // offset per thread for prefix sum
    char* fixes   = nullptr;          // nodal marker for topology fixes
    char* activ   = nullptr;          // nodal marker for kernel propagation
    uint8_t* tags = nullptr;          // nodal marker for boundary|corner
    Time tic;                         // time point for profiling purposes
  } sync;

#if DEFER_UPDATES
  struct Updates {
    std::vector<int> add;
    std::vector<int> rem;
  };
  // Table for pending updates
  std::vector<std::vector<Updates>> deferred;
  // scaling factor
  const int def_scale_fact = 32;
#endif
};
} // namespace
