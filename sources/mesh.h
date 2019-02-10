/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "header.h"
#include "timer.h"
#include "sync.h"
#include "numeric.h"
/* ------------------------------------*/
namespace trinity {

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
  int nb_nodes_;
  int nb_elems_;
  int nb_cores_;
  int nb_bucket_;

  // could use a std::set but performance suffers
  std::vector<std::vector<int>> stenc_;
  std::vector<std::vector<int>> vicin_;

  std::vector<int>     elems_;
  std::vector<double> points_;      // stride=2
  std::vector<double> tensor_;      // stride=3
  std::vector<double>  solut_;       // only on metric field calculation
  std::vector<double> qualit_;      //

  int _scale;        // capacity scale factor
  int _bucket;       // bucket capacity
  int _depth;        // max refinement/smoothing level
  int _verb;         // verbosity level
  int _iter;
  int _rounds;       // remeshing rounds
  size_t max_node;   // node max capacity
  size_t max_elem;   // elem max capacity
  int nb_activ_elem;

  int*      deg_;         // stencil degree
  int*      off_;         // offset per thread for prefix sum
  char*   fixes_;       // node marker for topology fixes
  char*   activ_;       // node marker for processFlips propagation
  uint8_t* tags_;     // node tags for geometrical features [isBoundary|isCorner]

  Time tic;

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
