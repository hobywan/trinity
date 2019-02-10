/*
 *                          'mesh.cpp'
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

#include "mesh.h"
#include "tools.h"
/* ------------------------------------ */
namespace trinity {
/* ------------------------------------ */
Mesh::Mesh(int size[2], int bucket, int depth, int verbosity, int rounds) {

  nb.nodes     = size[0];
  nb.elems     = size[1];
  nb.cores     = omp_get_max_threads();
  capa.scale  = (depth - 1) * 3;
  capa.bucket  = (size_t) bucket;
  param.depth  = depth;
  param.verb   = verbosity;
  param.iter   = 0;
  param.rounds = rounds;
  sync.deg     = nullptr;
  sync.off     = nullptr;
  sync.fixes   = nullptr;
  sync.activ   = nullptr;
  sync.tags    = nullptr;

#pragma omp parallel
  reallocMemory();
}

/* ------------------------------------ */
Mesh::~Mesh() {

  std::free(sync.activ);
  std::free(sync.deg);
  std::free(sync.off);
  std::free(sync.fixes);
  std::free(sync.tags);
}

/* ------------------------------------ */
void Mesh::reallocMemory() {
#pragma omp master
  {
    if (param.verb)
      std::printf("Mesh memory alloc ... ");
    sync.tic = timer::now();
  }

#pragma omp single
  {
    assert(nb.nodes);   // expected nb of nodes
    assert(nb.elems);   // expected nb of elems
    assert(capa.scale);
    capa.node = (size_t) nb.nodes * capa.scale;
    capa.elem = (size_t) nb.elems * capa.scale;

  #pragma omp task
    topo.stenc.resize(capa.node);
  #pragma omp task
    topo.vicin.resize(capa.node);
  #pragma omp task
    topo.elems.resize(capa.elem * 3);
  #pragma omp task
    geom.points.resize(capa.node * 2);
  #pragma omp task
    geom.tensor.resize(capa.node * 3);
  #pragma omp task
    geom.solut.resize(nb.nodes);  // only for metric calculation
  #pragma omp task
    geom.qualit.resize(capa.elem);
  #pragma omp task
    sync.deg = (int*) std::realloc(sync.deg, capa.node * sizeof(int));
  #pragma omp task
    sync.off = (int*) std::realloc(sync.off, nb.cores * sizeof(int));
  #pragma omp task
    sync.activ = (char*) std::realloc(sync.activ, capa.elem);
  #pragma omp task
    sync.fixes = (char*) std::realloc(sync.fixes, capa.node);
  #pragma omp task
    sync.tags = (uint8_t*) std::realloc(sync.tags, capa.node);
  #pragma omp taskwait
  }
  // first-touch
  doFirstTouch();

#ifdef DEFERRED_UPDATES
  initUpdates();
#endif

#pragma omp master
  {
    size_t memory_print = 0;
    memory_print += (capa.bucket * capa.node * sizeof(int)); // stenc
    memory_print += (6 * capa.node * sizeof(int));      // vicin
    memory_print += (3 * capa.elem * sizeof(int));      // elems
    memory_print += (2 * capa.node * sizeof(double));   // points
    memory_print += (3 * capa.node * sizeof(double));   // tensor
    memory_print += (1 * nb.nodes * sizeof(double));   // solut
    memory_print += (1 * capa.elem * sizeof(double));   // qualit
    memory_print += (1 * capa.node * sizeof(int));      // deg
    memory_print += (1 * nb.cores * sizeof(int));      // off
    memory_print += (1 * capa.elem);                    // activ
    memory_print += (1 * capa.node);                    // fixes
    memory_print += (1 * capa.node);                    // tags
    // 1 megabyte: 10^6, 1 mebibyte: 2^20=1048576
    if (param.verb) {
      std::printf(
        "%d MB \e[32m(%.2f s)\e[0m\n",
        (int) std::ceil(memory_print / 1e6), (float) timer::elapsed_ms(sync.tic) / 1e3
      );
    } else
      std::printf("\r= Preprocess ...  33 %% =");
    std::fflush(stdout);
  }
}

/* ------------------------------------ */
void Mesh::doFirstTouch() {
#pragma omp barrier

  assert(capa.node);
  assert(capa.elem);
  assert(nb.cores);

  const int chunk = nb.nodes / nb.cores;
  const int block = nb.elems / nb.cores;

#pragma omp for schedule(static, 1) nowait
  for (int i = 0; i < nb.cores; ++i)
    sync.off[i] = 0;
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < nb.nodes; ++i)
    geom.solut[i] = 0.;
#pragma omp for schedule(static, chunk*2) nowait
  for (int i = 0; i < capa.node * 2; ++i)
    geom.points[i] = 0.;
#pragma omp for schedule(static, chunk*3) nowait
  for (int i = 0; i < capa.node * 3; ++i)
    geom.tensor[i] = 0.;

  // (!)
#ifdef DEFERRED_UPDATES
#pragma omp for schedule(static,chunk) nowait
  for(int i=0; i < nb_nodes; ++i)
    topo.stenc[i].resize(64,-1);
#pragma omp for schedule(static,chunk) nowait
  for(int i=nb_nodes; i < capa.node; ++i)
    topo.stenc[i].reserve(64);
#else
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < capa.node; ++i)
    topo.stenc[i].resize(capa.bucket, -1);
#endif

#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < capa.node; ++i)
    sync.deg[i] = 0;
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < capa.node; ++i)
    sync.tags[i] = mask::unset;
#pragma omp for schedule(static, block*3) nowait
  for (int i = 0; i < capa.elem * 3; ++i)
    topo.elems[i] = -1;
#pragma omp for schedule(static, block) nowait
  for (int i = 0; i < capa.elem; ++i)
    sync.activ[i] = 0;
#pragma omp for schedule(static, block)
  for (int i = 0; i < capa.elem; ++i)
    geom.qualit[i] = 0.;
}

/* ------------------------------------ */
void Mesh::rebuildTopology() {

#pragma omp for schedule(static, nb.nodes/nb.cores)
  for (int i = 0; i < capa.node; ++i)
    sync.deg[i] = 0;

#pragma omp for
  for (int i = 0; i < nb.elems; ++i) {
    const int* n = getElem(i);
    if (__builtin_expect(*n < 0, 0))
      continue;

    // manually unrolled
    topo.stenc[n[0]][sync::fetchAndAdd(sync.deg + n[0], 1)] = i;
    topo.stenc[n[1]][sync::fetchAndAdd(sync.deg + n[1], 1)] = i;
    topo.stenc[n[2]][sync::fetchAndAdd(sync.deg + n[2], 1)] = i;
    assert(sync.deg[n[0]] < topo.stenc[n[0]].size());
    assert(sync.deg[n[1]] < topo.stenc[n[1]].size());
    assert(sync.deg[n[2]] < topo.stenc[n[2]].size());
  }

  std::vector<int> heap;

#pragma omp for
  for (int i = 0; i < nb.nodes; ++i) {
    heap.clear();
    heap.reserve(15);
    assert(sync.deg[i]);

    for (auto t = topo.stenc[i].begin(); t < topo.stenc[i].begin() + sync.deg[i]; ++t) {
      const int* n = getElem(*t);
      //manually unrolled
      if (n[0] == i) {
        heap.push_back(n[1]);
        heap.push_back(n[2]);
        continue;
      }
      if (n[1] == i) {
        heap.push_back(n[2]);
        heap.push_back(n[0]);
        continue;
      }
      if (n[2] == i) {
        heap.push_back(n[0]);
        heap.push_back(n[1]);
        continue;
      }
    }
    std::sort(heap.begin(), heap.end());
    heap.erase(std::unique(heap.begin(), heap.end()), heap.end());
    topo.vicin[i].swap(heap);
    std::sort(topo.stenc[i].begin(), topo.stenc[i].begin() + sync.deg[i]);
#ifdef DEFERRED_UPDATES
    topo.stenc[i].resize(deg[i]);
#endif
  }
}

/* ------------------------------------ */
bool Mesh::verifyTopology() const {

#pragma omp for
  for (int i = 0; i < nb.elems; ++i) {
    const int* n = getElem(i);
    if (__builtin_expect(*n > -1, 1)) {

      if (std::find(topo.stenc[n[0]].begin(), topo.stenc[n[0]].end(), i) == topo.stenc[n[0]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[0], i, n[0], n[1], n[2]);
        tools::display(topo.stenc[n[0]]);
      }
      if (std::find(topo.stenc[n[1]].begin(), topo.stenc[n[1]].end(), i) == topo.stenc[n[1]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[1], i, n[0], n[1], n[2]);
        tools::display(topo.stenc[n[1]]);
      }
      if (std::find(topo.stenc[n[2]].begin(), topo.stenc[n[2]].end(), i) == topo.stenc[n[2]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[2], i, n[0], n[1], n[2]);
        tools::display(topo.stenc[n[2]]);
      }
      // abort immediately if an error occured
      assert(n[0] not_eq n[1] and n[1] not_eq n[2] and n[2] not_eq n[0]);
      assert(std::find(topo.stenc[n[0]].begin(), topo.stenc[n[0]].end(), i) not_eq topo.stenc[n[0]].end());
      assert(std::find(topo.stenc[n[1]].begin(), topo.stenc[n[1]].end(), i) not_eq topo.stenc[n[1]].end());
      assert(std::find(topo.stenc[n[2]].begin(), topo.stenc[n[2]].end(), i) not_eq topo.stenc[n[2]].end());
    }
  }
  return true;
}

/* ------------------------------------ */
void Mesh::initActivElems() {

#pragma omp for
  for (int i = 0; i < nb.nodes; ++i)
    sync.activ[i] = static_cast<char>(__builtin_expect(topo.vicin[i].empty(), 0) ? 0 : 1);
}

/* ------------------------------------ */
void Mesh::extractPrimalGraph() {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb.nodes; ++i) {
    if (__builtin_expect(sync.activ[i], 1)) {
      topo.vicin[i].clear();
      topo.vicin[i].reserve(sync.deg[i] * 2);

      for (auto t = topo.stenc[i].begin(); t < topo.stenc[i].begin() + sync.deg[i]; ++t) {
        const int* n = getElem(*t);
        // manually unrolled
        if (i == n[0]) {
          topo.vicin[i].push_back(n[1]);
          topo.vicin[i].push_back(n[2]);
          continue;
        }
        if (i == n[1]) {
          topo.vicin[i].push_back(n[2]);
          topo.vicin[i].push_back(n[0]);
          continue;
        }
        if (i == n[2]) {
          topo.vicin[i].push_back(n[0]);
          topo.vicin[i].push_back(n[1]);
          continue;
        }
      }
      std::sort(topo.vicin[i].begin(), topo.vicin[i].end());
      topo.vicin[i].erase(std::unique(topo.vicin[i].begin(), topo.vicin[i].end()), topo.vicin[i].end());
    }
  }
}

/* ------------------------------------ */
void Mesh::extractDualGraph(Graph* dual) const {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb.elems; ++i) {
    const int* n = getElem(i);
    if (*n < 0)
      continue;

    auto& list = dual->at(i);
    list.clear();

    std::set_intersection(topo.stenc[n[0]].begin(), topo.stenc[n[0]].begin() + sync.deg[n[0]],
                          topo.stenc[n[1]].begin(), topo.stenc[n[1]].begin() + sync.deg[n[1]],
                          std::back_inserter(list));
    std::set_intersection(topo.stenc[n[0]].begin(), topo.stenc[n[0]].begin() + sync.deg[n[0]],
                          topo.stenc[n[2]].begin(), topo.stenc[n[2]].begin() + sync.deg[n[2]],
                          std::back_inserter(list));
    std::set_intersection(topo.stenc[n[1]].begin(), topo.stenc[n[1]].begin() + sync.deg[n[1]],
                          topo.stenc[n[2]].begin(), topo.stenc[n[2]].begin() + sync.deg[n[2]],
                          std::back_inserter(list));

    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
    std::swap(*(std::find(list.begin(), list.end(), i)), list.front());  // don't shift just swap
    assert(list[0] == i);
  }
}
/* ------------------------------------ */
const int* Mesh::getElem(int i) const { return topo.elems.data() + (i * 3); }
/* ------------------------------------ */
bool Mesh::isActiveNode(int i) const { return !topo.stenc[i].empty(); }
/* ------------------------------------ */
bool Mesh::isActiveElem(int i) const { return topo.elems[i * 3] > -1; }
/* ------------------------------------ */
bool Mesh::isBoundary(int i) const { return (sync.tags[i] & mask::bound); }
/* ------------------------------------ */
bool Mesh::isCorner(int i) const { return (sync.tags[i] & mask::corner); }
/* ------------------------------------ */
int Mesh::getCapaNode() const { return static_cast<int>(capa.node); }
/* ------------------------------------ */
int Mesh::getCapaElem() const { return static_cast<int>(capa.elem); }
/* ------------------------------------ */
Patch Mesh::getVicinity(int id, int dist) const {

  // verify step
  assert(dist > 1);
  assert(sync.deg[id] > 0);

  std::set<int> all[2];
  std::set<int> last[2];
  std::set<int> cur[2];

  // init
  last[0].insert(topo.vicin[id].begin(), topo.vicin[id].end());
  all[0].insert(topo.vicin[id].begin(), topo.vicin[id].end());
  all[1].insert(topo.stenc[id].begin(), topo.stenc[id].begin() + sync.deg[id]);

  for (int i = 1; i < dist; ++i) {
    cur[0].clear();
    cur[1].clear();
    // retrieve vicin of each vertex
    for (int k : last[0]) {
      cur[0].insert(topo.vicin[k].begin(), topo.vicin[k].end());
      cur[1].insert(topo.stenc[k].begin(), topo.stenc[k].begin() + sync.deg[k]);
    }
    for (int j = 0; j < 2; ++j) {
      last[j].clear();
      std::set_difference(cur[j].begin(), cur[j].end(),
                          all[j].begin(), all[j].end(),
                          std::inserter(last[j], last[j].begin()));
      all[j].insert(last[j].begin(), last[j].end());
    }
  }

  // remove init vertex
  all[0].erase(std::find(all[0].begin(), all[0].end(), id));

  // copy back from sets to patch
  Patch patch;
  patch.node.reserve(all[0].size());
  patch.node.insert(patch.node.begin(), all[0].begin(), all[0].end());
  patch.elem.reserve(all[1].size());
  patch.elem.insert(patch.elem.begin(), all[1].begin(), all[1].end());
  return patch;
}

/* ------------------------------------ */
int Mesh::getElemNeigh(int id, int i, int j) const {

  assert(not topo.stenc[i].empty());
  for (auto t = topo.stenc[i].begin(); t < topo.stenc[i].end() and *t > -1; ++t) {
    if (__builtin_expect(*t not_eq id, 1)) {
      const int* n = getElem(*t);
      // manually unrolled
      if ((n[0] == j and n[1] == i) or
          (n[1] == j and n[2] == i) or
          (n[2] == j and n[0] == i)) {
        return *t;
      }
    }
  }
  assert(isBoundary(i) and isBoundary(j));
  return -1;
}

/* ------------------------------------ */
void Mesh::replaceElem(int id, const int* v) {
  assert(v not_eq nullptr);
  std::memcpy(topo.elems.data() + (id * 3), v, sizeof(int) * 3);
}

/* ------------------------------------ */
void Mesh::eraseElem(int id) {
  // nb: std::memsetting with -1 is ok
  std::memset(topo.elems.data() + (id * 3), -1, sizeof(int) * 3);
}

/* ------------------------------------ */
void Mesh::updateStencil(int i, int t) {
  int k = sync::fetchAndAdd(sync.deg + i, 1);
  sync::reallocBucket(topo.stenc.data(), i, (k + 1), param.verb);
  topo.stenc[i][k] = t;
}

/* ------------------------------------ */
void Mesh::updateStencil(int i, const std::initializer_list<int>& t) {
  int j = sync::fetchAndAdd(sync.deg + i, (int) t.size());
  if (topo.stenc[i].empty()) {
    assert(topo.vicin[i].size() == 0);
    assert(topo.elems[i * 3] == -1);
  }
  //assert(stenc[i].size()>0);
  sync::reallocBucket(topo.stenc.data(), i, j + t.size(), param.verb);
  //
  if ((j + t.size()) >= topo.stenc[i].capacity()) {
    std::fprintf(stderr, "stenc[%d] not fixed, stenc.size: %lu\n", i, topo.stenc[i].size());
    std::fflush(stderr);
  }

  assert((j + t.size()) < topo.stenc[i].capacity());
  for (size_t k = 0; k < t.size(); ++k)
    topo.stenc[i][j + k] = *(t.begin() + k);
}

/* ------------------------------------ */
void Mesh::copyStencil(int i, int j, int nb_rm) {
  int chunk = sync.deg[i] - nb_rm;
  int k = sync::fetchAndAdd(sync.deg + j, chunk);
  // threads attempting to insert will spin until reallocation was done
  sync::reallocBucket(topo.stenc.data(), i, k + chunk, param.verb);

  assert((k + chunk) < topo.stenc[i].size());
  std::memcpy(topo.stenc[j].data() + k, topo.stenc[i].data(), chunk * sizeof(int));
}

/* ------------------------------------ */
#ifdef DEFERRED_UPDATES
/* ------------------------------------ */
void Mesh::deferredAppend(int tid, int i, int t){
  const int key = tools::hash(i)% (def_scale_fact*nb_cores);
  deferred[tid][key].add.push_back(i);
  deferred[tid][key].add.push_back(t);
}
/* ------------------------------------ */
void Mesh::deferredAppend(int tid, int i, const std::initializer_list<int>& t){

  const int key = tools::hash(i)% (def_scale_fact*nb_cores);
  for(size_t k=0; k < t.size(); ++k){
    auto found = std::find(stenc[i].begin(), stenc[i].end(), *(t.begin()+k));
    assert(found == stenc[i].end());
    deferred[tid][key].add.push_back(i);
    deferred[tid][key].add.push_back(*(t.begin()+k));
  }
}
/* ------------------------------------ */
void Mesh::deferredRemove(int tid, int i, int t){
  const int key = tools::hash(i)% (def_scale_fact*nb_cores);
  auto found = std::find(stenc[i].begin(), stenc[i].end(), t);
  assert(found not_eq stenc[i].end());
  deferred[tid][key].rem.push_back(i);
  deferred[tid][key].rem.push_back(t);
}
/* ------------------------------------ */
void Mesh::initUpdates(){

  int size = nb_cores * def_scale_fact;

#pragma omp single
  deferred.resize(nb_cores);

#pragma omp for
  for(int i=0; i < nb_cores; ++i){
    deferred[i].resize(size);
    for(int j=0; j < size; ++j){
      deferred[i][j].add.reserve(64);
      deferred[i][j].rem.reserve(64);
    }
  }
}
/* ------------------------------------ */
void Mesh::commitUpdates(){

  // from 'pragmatic'
#pragma omp for schedule(guided)
  for(int j=0; j < nb_cores*def_scale_fact; ++j){
    for(int i=0; i < nb_cores; ++i){
      // only one thread will update the adjacency list of a given node
      for(auto v = deferred[i][j].rem.begin(); v < deferred[i][j].rem.end(); v+=2){
        auto found = std::find(stenc[*v].begin(), stenc[*v].end(), *(v+1));
        assert(found not_eq stenc[*v].end());
        stenc[*v].erase(found);
      }
      deferred[i][j].rem.clear();

      for(auto v = deferred[i][j].add.begin(); v < deferred[i][j].add.end(); v+=2)
        stenc[*v].push_back(*(v+1));
      deferred[i][j].add.clear();
    }
  }

#pragma omp for
  for(int i=0; i < nb_nodes; ++i)
    deg[i] = (int) stenc[i].size();

  verify();
}
/* ------------------------------------ */
void Mesh::resetUpdates(){

#pragma omp for schedule(guided)
  for(int j=0; j < (nb_cores * def_scale_fact); ++j){
    for(int i=0; i < nb_cores; ++i){
      deferred[i][j].add.clear();
      deferred[i][j].rem.clear();
    }
  }
}
/* ------------------------------------ */
#endif

/* ------------------------------------ */
const int* Mesh::getElemCoord(int id, double* p) const {

  // manually unrolled
  const int* n = getElem(id);
  std::memcpy(p    , geom.points.data() + (n[0] * 2), 2 * sizeof(double));
  std::memcpy(p + 2, geom.points.data() + (n[1] * 2), 2 * sizeof(double));
  std::memcpy(p + 4, geom.points.data() + (n[2] * 2), 2 * sizeof(double));
  return n;
}

/* ------------------------------------ */
double Mesh::computeLength(int i, int j) const {

  const double* p1 = geom.points.data() + (i * 2);
  const double* p2 = geom.points.data() + (j * 2);
  const double* M1 = geom.tensor.data() + (i * 3);
  const double* M2 = geom.tensor.data() + (j * 3);

  return numeric::approxRiemannDist(p1, p2, M1, M2);
}

/* ------------------------------------ */
double Mesh::computeQuality(const int* t) const {

  assert(t[0] > -1 and t[1] > -1 and t[2] > -1);

  // zero-copy
  const double* pa = geom.points.data() + (t[0] * 2);
  const double* pb = geom.points.data() + (t[1] * 2);
  const double* pc = geom.points.data() + (t[2] * 2);
  const double* ma = geom.tensor.data() + (t[0] * 3);
  const double* mb = geom.tensor.data() + (t[1] * 3);
  const double* mc = geom.tensor.data() + (t[2] * 3);

  return numeric::computeQuality(pa, pb, pc, ma, mb, mc);
}

/* ------------------------------------ */
double Mesh::computeQuality(int id) const {
  return computeQuality(getElem(id));
}

/* ------------------------------------ */
void Mesh::computeQuality(double* q) {

  double q_min = 3.;
  double q_max = 0.;
  double q_tot = 0.;
  int count = 0;

#pragma omp single
  nb.activ_elem = 0;

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb.elems; ++i) {
    const int* t = getElem(i);
    if (*t < 0)
      continue;

    geom.qualit[i] = computeQuality(t);
    if (q_min > geom.qualit[i]) q_min = geom.qualit[i];
    if (q_max < geom.qualit[i]) q_max = geom.qualit[i];
    q_tot += geom.qualit[i];
    count++;
  }

#pragma omp critical
  {
    q[0] = std::min(q[0], q_min);
    q[1] = std::max(q[1], q_max);
    q[2] += q_tot;
    nb.activ_elem += count;
  }
#pragma omp barrier
#pragma omp single
  q[2] /= nb.activ_elem;

}

/* ------------------------------------ */
void Mesh::computeSteinerPoint(int i, int j, double* p, double* M) const {

  assert(p not_eq nullptr);
  assert(M not_eq nullptr);

  const double* p1 = geom.points.data() + (i * 2);
  const double* p2 = geom.points.data() + (j * 2);
  const double* m1 = geom.tensor.data() + (i * 3);
  const double* m2 = geom.tensor.data() + (j * 3);

  numeric::computeSteinerPoint(p1, p2, m1, m2, p, M);
}

/* ------------------------------------ */
bool Mesh::isCounterclockwise(const int* t) const {

  const double* pa = geom.points.data() + (t[0] * 2);
  const double* pb = geom.points.data() + (t[1] * 2);
  const double* pc = geom.points.data() + (t[2] * 2);

  double s[4];
  s[0] = pb[0] - pa[0];
  s[1] = pc[0] - pa[0];
  s[2] = pb[1] - pa[1];
  s[3] = pc[1] - pa[1];

  return 0.5 * (s[0] * s[3] - s[1] * s[2]) > EPSILON;
}

/* ------------------------------------ */
void Mesh::fixTagged() {
#ifdef DEFERRED_UPDATES
#pragma omp for schedule(guided)
  for(int i=0; i < nb_nodes; ++i)
    fixes[i] = 0;
  //
  commitUpdates();

#else
  std::vector<int> elem;
  elem.resize(capa.bucket);
  int k;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb.nodes; ++i) {
    if (__builtin_expect(sync.fixes[i], 0)) {
      // reset
      sync.fixes[i] = 0;
      k = 0;
      if (__builtin_expect(elem.size() < sync.deg[i], 0))
        elem.resize(sync.deg[i]);

      for (auto t = topo.stenc[i].begin(); t < topo.stenc[i].begin() + sync.deg[i]; ++t) {
        const int* n = getElem(*t);
        if (__builtin_expect(i == n[0] or i == n[1] or i == n[2], 1))
          elem[k++] = *t;
      }
      if (k >= elem.size())
        std::printf("k: %d, elem.size: %lu\n", k, elem.size());
      assert(k < elem.size());
      // remove duplicates and adjust count
      std::sort(elem.begin(), elem.begin() + k);
      sync.deg[i] = static_cast<int>(std::unique(elem.begin(), elem.begin() + k) - elem.begin());
      std::fill(elem.begin() + sync.deg[i], elem.end(), -1);
      // swap arrays O(1)
      topo.stenc[i].swap(elem);
    }
  }
#endif
}

/* ------------------------------------ */
void Mesh::fixAll() {

#pragma omp for
  for (int i = 0; i < nb.nodes; ++i)
    sync.fixes[i] = static_cast<char>(topo.stenc[i].empty() ? 0 : 1);

  this->fixTagged();
}

} // namespace trinity