/*
 *                          'mesh.cpp'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *                Copyright 2016, Hoby Rakotoarivelo.
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

#include "trinity/mesh.h"
#include "trinity/tools.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
Mesh::Mesh(int size[2], int bucket, int depth, int verbosity, int rounds) {

  nb.nodes     = size[0];
  nb.elems     = size[1];
  nb.cores     = omp_get_max_threads();
  capa.scale   = (depth - 1) * 3;
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

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() {

  std::free(sync.activ);
  std::free(sync.deg);
  std::free(sync.off);
  std::free(sync.fixes);
  std::free(sync.tags);
}

/* -------------------------------------------------------------------------- */
void Mesh::reallocMemory() {
#pragma omp master
  {
    if (param.verb)
      std::printf("Mesh memory alloc ... ");
    sync.tic = timer::now();
  }

#pragma omp single
  {
    assert(nb.nodes > 0);   // expected nb of nodes
    assert(nb.elems > 0);   // expected nb of elems
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
    geom.solut.resize((size_t) nb.nodes);  // only for metric calculation
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

#if DEFER_UPDATES
  initUpdates();
#endif

#pragma omp master
  {
    // (!) count only containers
    size_t memory_print = 0;
    memory_print += (capa.bucket * capa.node * sizeof(int)); // stenc
    memory_print += (6 * capa.node * sizeof(int));           // vicin
    memory_print += (3 * capa.elem * sizeof(int));           // elems
    memory_print += (2 * capa.node * sizeof(double));        // points
    memory_print += (3 * capa.node * sizeof(double));        // tensor
    memory_print += (1 * nb.nodes  * sizeof(double));        // solut
    memory_print += (1 * capa.elem * sizeof(double));        // qualit
    memory_print += (1 * capa.node * sizeof(int));           // deg
    memory_print += (1 * nb.cores  * sizeof(int));           // off
    memory_print += (1 * capa.elem);                         // activ
    memory_print += (1 * capa.node);                         // fixes
    memory_print += (1 * capa.node);                         // tags

    if (param.verb) {
      // 1 megabyte: 10^6, 1 mebibyte: 2^20=1048576
      auto size = (int) std::ceil(memory_print / 1e6);
      auto secs = (float) timer::elapsed_ms(sync.tic) / 1e3;
      std::printf("%d MB \e[32m(%.2f s)\e[0m\n", size, secs);
    } else {
      std::printf("\r= Preprocess ...  33 %% =");
    }
    std::fflush(stdout);
  }
}

/* -------------------------------------------------------------------------- */
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
  for (unsigned i = 0; i < capa.node * 2; ++i)
    geom.points[i] = 0.;
#pragma omp for schedule(static, chunk*3) nowait
  for (unsigned i = 0; i < capa.node * 3; ++i)
    geom.tensor[i] = 0.;

  // (!)
#if DEFER_UPDATES
#pragma omp for schedule(static,chunk) nowait
  for(int i = 0; i < nb.nodes; ++i)
    topo.stenc[i].resize(64, -1);
#pragma omp for schedule(static,chunk) nowait
  for(unsigned i = (unsigned) nb.nodes; i < capa.node; ++i)
    topo.stenc[i].reserve(64);
#else
#pragma omp for schedule(static, chunk) nowait
  for (unsigned i = 0; i < capa.node; ++i)
    topo.stenc[i].resize(capa.bucket, -1);
#endif

#pragma omp for schedule(static, chunk) nowait
  for (unsigned i = 0; i < capa.node; ++i)
    sync.deg[i] = 0;
#pragma omp for schedule(static, chunk) nowait
  for (unsigned i = 0; i < capa.node; ++i)
    sync.tags[i] = mask::unset;
#pragma omp for schedule(static, block*3) nowait
  for (unsigned i = 0; i < capa.elem * 3; ++i)
    topo.elems[i] = -1;
#pragma omp for schedule(static, block) nowait
  for (unsigned i = 0; i < capa.elem; ++i)
    sync.activ[i] = 0;
#pragma omp for schedule(static, block)
  for (unsigned i = 0; i < capa.elem; ++i)
    geom.qualit[i] = 0.;
}

/* -------------------------------------------------------------------------- */
void Mesh::rebuildTopology() {

#pragma omp for schedule(static, nb.nodes/nb.cores)
  for (int i = 0; i < capa.node; ++i)
    sync.deg[i] = 0;

#pragma omp for
  for (int i = 0; i < nb.elems; ++i) {
    auto const n = getElem(i);
    if (__builtin_expect(*n < 0, 0))
      continue;

    // manually unrolled
    topo.stenc[n[0]][sync::fetchAndAdd(sync.deg + n[0], 1)] = i;
    topo.stenc[n[1]][sync::fetchAndAdd(sync.deg + n[1], 1)] = i;
    topo.stenc[n[2]][sync::fetchAndAdd(sync.deg + n[2], 1)] = i;
    assert(sync.deg[n[0]] < (int) topo.stenc[n[0]].size());
    assert(sync.deg[n[1]] < (int) topo.stenc[n[1]].size());
    assert(sync.deg[n[2]] < (int) topo.stenc[n[2]].size());
  }

  std::vector<int> heap;

#pragma omp for
  for (int i = 0; i < nb.nodes; ++i) {
    heap.clear();
    heap.reserve(15);
    assert(sync.deg[i]);

    for (auto t = topo.stenc[i].begin();
              t < topo.stenc[i].begin() + sync.deg[i]; ++t)
    {
      auto const n = getElem(*t);
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
    #if DEFER_UPDATES
      topo.stenc[i].resize(sync.deg[i]);
    #endif
  }
}

/* -------------------------------------------------------------------------- */
bool Mesh::verifyTopology() const {

#pragma omp for
  for (int i = 0; i < nb.elems; ++i) {
    auto t = getElem(i);
    if (__builtin_expect(t[0] > -1, 1)) {

      for(int j = 0; j < 3; ++j){
        #if DEBUG
          if (std::find(
            topo.stenc[t[j]].begin(), topo.stenc[t[j]].end(), i) ==
            topo.stenc[t[j]].end()
          )
          #pragma omp critical
            {
              std::fprintf(stderr,
                "n: %d, t: %d [%d,%d,%d]", t[0], i, t[0], t[1], t[2]
              );
              tools::display(topo.stenc[t[j]]);
              std::exit(EXIT_FAILURE);
            }
        #else
          // abort immediately if an error occured
          assert(t[0] != t[1] and t[1] != t[2] and t[2] != t[0]);
          auto cur = t[j];
          auto itr = std::find(
            topo.stenc[cur].begin(), topo.stenc[cur].end(), i
          );
          assert(itr != topo.stenc[cur].end());
        #endif
      }
    }
  }
  return true;
}

/* -------------------------------------------------------------------------- */
void Mesh::initActivElems() {

#pragma omp for
  for (int i = 0; i < nb.nodes; ++i)
    sync.activ[i] = char(__builtin_expect(topo.vicin[i].empty(), 0) ? 0 : 1);
}

/* -------------------------------------------------------------------------- */
void Mesh::extractPrimalGraph() {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb.nodes; ++i) {
    if (__builtin_expect(sync.activ[i], 1)) {
      auto size = (size_t) sync.deg[i] * 2;
      topo.vicin[i].clear();
      topo.vicin[i].reserve(size);

      for (auto t = topo.stenc[i].begin();
                t < topo.stenc[i].begin() + sync.deg[i]; ++t) {
        auto n = getElem(*t);
        // manually unrolled
        if (i == n[0]) {
          topo.vicin[i].push_back(n[1]);
          topo.vicin[i].push_back(n[2]);
        }
        else if (i == n[1]) {
          topo.vicin[i].push_back(n[2]);
          topo.vicin[i].push_back(n[0]);
        }
        else if (i == n[2]) {
          topo.vicin[i].push_back(n[0]);
          topo.vicin[i].push_back(n[1]);
        }
      }
      std::sort(topo.vicin[i].begin(), topo.vicin[i].end());
      topo.vicin[i].erase(
        std::unique(topo.vicin[i].begin(), topo.vicin[i].end()),
        topo.vicin[i].end()
      );
    }
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::extractDualGraph(Graph* dual) const {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb.elems; ++i) {
    auto const n = getElem(i);
    if (*n < 0)
      continue;

    auto& list = dual->at(i);
    list.clear();

    std::set_intersection(
      topo.stenc[n[0]].begin(), topo.stenc[n[0]].begin() + sync.deg[n[0]],
      topo.stenc[n[1]].begin(), topo.stenc[n[1]].begin() + sync.deg[n[1]],
      std::back_inserter(list)
    );

    std::set_intersection(
      topo.stenc[n[0]].begin(), topo.stenc[n[0]].begin() + sync.deg[n[0]],
      topo.stenc[n[2]].begin(), topo.stenc[n[2]].begin() + sync.deg[n[2]],
      std::back_inserter(list)
    );

    std::set_intersection(
      topo.stenc[n[1]].begin(), topo.stenc[n[1]].begin() + sync.deg[n[1]],
      topo.stenc[n[2]].begin(), topo.stenc[n[2]].begin() + sync.deg[n[2]],
      std::back_inserter(list)
    );

    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
    std::swap(*(std::find(list.begin(), list.end(), i)), list.front());
    assert(list[0] == i);
  }
}
/* -------------------------------------------------------------------------- */
const int* Mesh::getElem(int i) const { return topo.elems.data() + (i * 3); }

/* -------------------------------------------------------------------------- */
bool Mesh::isActiveNode(int i) const { return not topo.stenc[i].empty(); }

/* -------------------------------------------------------------------------- */
bool Mesh::isActiveElem(int i) const { return topo.elems[i * 3] > -1; }

/* -------------------------------------------------------------------------- */
bool Mesh::isBoundary(int i) const { return (sync.tags[i] & mask::bound); }

/* -------------------------------------------------------------------------- */
bool Mesh::isCorner(int i) const { return (sync.tags[i] & mask::corner); }

/* -------------------------------------------------------------------------- */
int Mesh::getCapaNode() const { return static_cast<int>(capa.node); }

/* -------------------------------------------------------------------------- */
int Mesh::getCapaElem() const { return static_cast<int>(capa.elem); }

/* -------------------------------------------------------------------------- */
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
      std::set_difference(
        cur[j].begin(), cur[j].end(),
        all[j].begin(), all[j].end(),
        std::inserter(last[j], last[j].begin())
      );
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

/* -------------------------------------------------------------------------- */
int Mesh::getElemNeigh(int id, int i, int j) const {

  assert(not topo.stenc[i].empty());
  for (auto t = topo.stenc[i].begin(); t < topo.stenc[i].end() and *t > -1; ++t) {
    if (__builtin_expect(*t != id, 1)) {
      auto const n = getElem(*t);
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

/* -------------------------------------------------------------------------- */
void Mesh::replaceElem(int id, const int* v) {
  assert(v != nullptr);
  std::memcpy(topo.elems.data() + (id * 3), v, sizeof(int) * 3);
}

/* -------------------------------------------------------------------------- */
void Mesh::eraseElem(int id) {
  // nb: memset with -1 is ok
  std::memset(topo.elems.data() + (id * 3), -1, sizeof(int) * 3);
}

/* -------------------------------------------------------------------------- */
void Mesh::updateStencil(int i, int t) {
  int k = sync::fetchAndAdd(sync.deg + i, 1);
  sync::reallocBucket(topo.stenc.data(), i, (k + 1), param.verb);
  topo.stenc[i][k] = t;
}

/* -------------------------------------------------------------------------- */
void Mesh::updateStencil(int i, const std::initializer_list<int>& t) {
  int j = sync::fetchAndAdd(sync.deg + i, (int) t.size());
  if (topo.stenc[i].empty()) {
    assert(topo.vicin[i].empty());
    assert(topo.elems[i * 3] == -1);
  }
  //assert(stenc[i].size()>0);
  sync::reallocBucket(topo.stenc.data(), i, j + t.size(), param.verb);
  //
  if ((j + t.size()) >= topo.stenc[i].capacity()) {
    std::fprintf(stderr,
      "stenc[%d] not fixed, stenc.size: %lu\n", i, topo.stenc[i].size()
    );
    std::fflush(stderr);
  }

  assert((j + t.size()) < topo.stenc[i].capacity());
  for (size_t k = 0; k < t.size(); ++k)
    topo.stenc[i][j + k] = *(t.begin() + k);
}

/* -------------------------------------------------------------------------- */
void Mesh::copyStencil(int i, int j, int nb_rm) {
  int chunk = sync.deg[i] - nb_rm;
  int k = sync::fetchAndAdd(sync.deg + j, chunk);
  // threads attempting to insert will spin until reallocation was done
  sync::reallocBucket(topo.stenc.data(), i, k + chunk, param.verb);

  assert(size_t(k + chunk) < topo.stenc[i].size());
  auto const bytes = chunk * sizeof(int);
  std::memcpy(topo.stenc[j].data() + k, topo.stenc[i].data(), bytes);
}

/* -------------------------------------------------------------------------- */
#if DEFER_UPDATES
/* -------------------------------------------------------------------------- */
void Mesh::deferredAppend(int tid, int i, int t){
  const int key = tools::hash(i) % (def_scale_fact * nb.cores);
  deferred[tid][key].add.push_back(i);
  deferred[tid][key].add.push_back(t);
}
/* -------------------------------------------------------------------------- */
void Mesh::deferredAppend(int tid, int i, const std::initializer_list<int>& t){

  const int key = tools::hash(i) % (def_scale_fact * nb.cores);
  for(size_t k = 0; k < t.size(); ++k){
    auto found = std::find(
      topo.stenc[i].begin(), topo.stenc[i].end(), *(t.begin()+k)
    );
    assert(found == stenc[i].end());
    deferred[tid][key].add.push_back(i);
    deferred[tid][key].add.push_back(*(t.begin()+k));
  }
}
/* -------------------------------------------------------------------------- */
void Mesh::deferredRemove(int tid, int i, int t){
  const int key = tools::hash(i)% (def_scale_fact * nb.cores);
  auto found = std::find(topo.stenc[i].begin(), topo.stenc[i].end(), t);
  assert(found != stenc[i].end());
  deferred[tid][key].rem.push_back(i);
  deferred[tid][key].rem.push_back(t);
}
/* -------------------------------------------------------------------------- */
void Mesh::initUpdates(){

  int const size = nb.cores * def_scale_fact;

#pragma omp single
  deferred.resize((size_t) nb.cores);

#pragma omp for
  for (int i = 0; i < nb.cores; ++i) {
    deferred[i].resize((size_t) size);
    for (int j = 0; j < size; ++j) {
      deferred[i][j].add.reserve(64);
      deferred[i][j].rem.reserve(64);
    }
  }
}
/* -------------------------------------------------------------------------- */
void Mesh::commitUpdates(){

  // from 'pragmatic'
#pragma omp for schedule(guided)
  for (int j = 0; j < nb.cores * def_scale_fact; ++j) {
    for (int i = 0; i < nb.cores; ++i) {

      auto& defer_rem = deferred[i][j].rem;
      auto& defer_add = deferred[i][j].add;

      // only one thread will update the adjacency list of a given node
      for (auto v = defer_rem.begin(); v < defer_rem.end(); v += 2) {
        auto& stenc = topo.stenc[*v];
        auto const found = std::find(stenc.begin(), stenc.end(), *(v + 1));
        assert(found != stenc.end());
        stenc.erase(found);
      }

      for(auto v = defer_add.begin(); v < defer_add.end(); v += 2) {
        topo.stenc[*v].push_back(*(v + 1));
      }

      defer_rem.clear();
      defer_add.clear();
    }
  }

#pragma omp for
  for(int i = 0; i < nb.nodes; ++i)
    sync.deg[i] = (int) topo.stenc[i].size();

  this->verifyTopology();
}
/* -------------------------------------------------------------------------- */
void Mesh::resetUpdates(){

#pragma omp for schedule(guided)
  for(int j = 0; j < (nb.cores * def_scale_fact); ++j){
    for(int i = 0; i < nb.cores; ++i){
      deferred[i][j].add.clear();
      deferred[i][j].rem.clear();
    }
  }
}
/* -------------------------------------------------------------------------- */
#endif

/* -------------------------------------------------------------------------- */
const int* Mesh::getElemCoord(int id, double* p) const {

  // manually unrolled
  auto const n = getElem(id);
  auto const bytes = 2 * sizeof(double);
  std::memcpy(p    , geom.points.data() + (n[0] * 2), bytes);
  std::memcpy(p + 2, geom.points.data() + (n[1] * 2), bytes);
  std::memcpy(p + 4, geom.points.data() + (n[2] * 2), bytes);
  return n;
}

/* -------------------------------------------------------------------------- */
double Mesh::computeLength(int i, int j) const {

  auto const p1 = geom.points.data() + (i * 2);
  auto const p2 = geom.points.data() + (j * 2);
  auto const M1 = geom.tensor.data() + (i * 3);
  auto const M2 = geom.tensor.data() + (j * 3);

  return numeric::approxRiemannDist(p1, p2, M1, M2);
}

/* -------------------------------------------------------------------------- */
double Mesh::computeQuality(const int* t) const {

  assert(t[0] > -1 and t[1] > -1 and t[2] > -1);

  // zero-copy
  auto const pa = geom.points.data() + (t[0] * 2);
  auto const pb = geom.points.data() + (t[1] * 2);
  auto const pc = geom.points.data() + (t[2] * 2);
  auto const ma = geom.tensor.data() + (t[0] * 3);
  auto const mb = geom.tensor.data() + (t[1] * 3);
  auto const mc = geom.tensor.data() + (t[2] * 3);

  return numeric::computeQuality(pa, pb, pc, ma, mb, mc);
}

/* -------------------------------------------------------------------------- */
double Mesh::computeQuality(int id) const {
  return computeQuality(getElem(id));
}

/* -------------------------------------------------------------------------- */
void Mesh::computeQuality(double quality[3]) {

  double q_min = 3.;
  double q_max = 0.;
  double q_tot = 0.;
  int count = 0;

#pragma omp single
  nb.activ_elem = 0;

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb.elems; ++i) {
    auto const t = getElem(i);
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
    quality[0] = std::min(quality[0], q_min);
    quality[1] = std::max(quality[1], q_max);
    quality[2] += q_tot;
    nb.activ_elem += count;
  }
#pragma omp barrier
#pragma omp single
  quality[2] /= nb.activ_elem;

}

/* -------------------------------------------------------------------------- */
void Mesh::computeSteinerPoint(
  int edge_i, int edge_j, double point[2], double metric[3]
) const {

  assert(point != nullptr);
  assert(metric != nullptr);

  auto const p1 = geom.points.data() + (edge_i * 2);
  auto const p2 = geom.points.data() + (edge_j * 2);
  auto const m1 = geom.tensor.data() + (edge_i * 3);
  auto const m2 = geom.tensor.data() + (edge_j * 3);

  numeric::computeSteinerPoint(p1, p2, m1, m2, point, metric);
}

/* -------------------------------------------------------------------------- */
bool Mesh::isCounterclockwise(const int* elem) const {

  auto const point_a = geom.points.data() + (elem[0] * 2);
  auto const point_b = geom.points.data() + (elem[1] * 2);
  auto const point_c = geom.points.data() + (elem[2] * 2);

  double coef[4];
  coef[0] = point_b[0] - point_a[0];
  coef[1] = point_c[0] - point_a[0];
  coef[2] = point_b[1] - point_a[1];
  coef[3] = point_c[1] - point_a[1];

  return 0.5 * (coef[0] * coef[3] - coef[1] * coef[2]) > EPSILON;
}

/* -------------------------------------------------------------------------- */
void Mesh::fixTagged() {
#if DEFER_UPDATES
#pragma omp for schedule(guided)
  for(int i = 0; i < nb.nodes; ++i)
    sync.fixes[i] = 0;
  //
  commitUpdates();

#else
  std::vector<int> elem;
  elem.resize(capa.bucket);

#pragma omp for schedule(guided)
  for (int i = 0; i < nb.nodes; ++i) {
    if (__builtin_expect(sync.fixes[i], 0)) {
      // reset
      sync.fixes[i] = 0;
      int k = 0;
      if (__builtin_expect((int) elem.size() < sync.deg[i], 0))
        elem.resize((size_t) sync.deg[i]);

      auto const& stenc = topo.stenc[i];
      for (auto t = stenc.begin(); t < stenc.begin() + sync.deg[i]; ++t) {
        auto const n = getElem(*t);
        if (__builtin_expect(i == n[0] or i == n[1] or i == n[2], 1))
          elem[k++] = *t;
      }
      #if DEBUG
        if (k >= (int) elem.size()) {
          std::printf("k: %d, elem.size: %lu\n", k, elem.size());
        }
      #endif
      assert(k < (int) elem.size());
      // remove duplicates and adjust count
      std::sort(elem.begin(), elem.begin() + k);
      sync.deg[i] = int(
        std::unique(elem.begin(), elem.begin() + k) - elem.begin()
      );
      std::fill(elem.begin() + sync.deg[i], elem.end(), -1);
      // swap arrays O(1)
      topo.stenc[i].swap(elem);
    }
  }
#endif
}

/* -------------------------------------------------------------------------- */
void Mesh::fixAll() {

#pragma omp for
  for (int i = 0; i < nb.nodes; ++i)
    sync.fixes[i] = char(topo.stenc[i].empty() ? 0 : 1);

  this->fixTagged();
}

/* -------------------------------------------------------------------------- */
} // namespace trinity
