/*
 *                          'mesh.h'
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
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
#include "sync.h"
#include "numeric.h"
/* ------------------------------------ */
namespace trinity {
/* ------------------------------------ */
class Mesh {

  friend class Metrics;
  friend class Refine;
  friend class Coarse;
  friend class Swap;
  friend class Smooth;
  friend class Indep;
  friend class Partit;

public:

   Mesh() = default;
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
  void loadFrom(const std::string& path, const std::string& solu);
  void storeTo(const std::string& path) const;
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
#ifdef DEFERRED_UPDATES
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
  bool isCounterclockwise(const int* t) const;
  void computeSteinerPoint(int i, int j, double* p, double* m) const;
  void computeQuality(double* q);

private:

  struct {
    int nodes;                        // nb of active+deleted nodes
    int elems;                        // nb of active+deleted elems
    int cores;                        // nb of logical cores used
    int activ_elem;                   // actual nb of elems
  } nb;

  struct {
    std::vector<int> elems;
    Graph stenc;                      // nodal incident elems
    Graph vicin;                      // nodal adjacent nodes
  } topo;                             // could use a std::set but performance suffers

  struct {
    std::vector<double> points;       // nodal coordinates, stride=2
    std::vector<double> tensor;       // nodal metric tensor, stride=3
    std::vector<double> solut;        // solut value (for metric field calculation)
    std::vector<double> qualit;       // cache elem quality to speedup kernels
  } geom;

  struct {
    int depth;                        // refine|smooth max depth
    int verb;                         // verbosity level
    int iter;                         // current iteration
    int rounds;                       // remeshing rounds
  } param;

  struct {
    int    scale;                     // capacity scale factor
    size_t bucket;                    // bucket capacity
    size_t node;                      // node capacity
    size_t elem;                      // elem capacity
  } capa;

  struct {
    int*     deg;                     // nodal degree
    int*     off;                     // offset per thread for prefix sum
    char*    fixes;                   // nodal marker for topology fixes
    char*    activ;                   // nodal marker for kernel propagation
    uint8_t* tags;                    // nodal tags for geometrical features [bound|corner]
    Time     tic;                     // time point for profiling purposes
  } sync;

#ifdef DEFERRED_UPDATES
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
