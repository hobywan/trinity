/* ------------------------------------ */
#include "mesh.h"
#include "tools.h"
/* ------------------------------------ */
namespace trinity {
/* ------------------------------------ */
Mesh::Mesh(int size[2], int bucket, int depth, int verbosity, int rounds)
  : nb_nodes_(size[0]),
    nb_elems_(size[1]),
    nb_cores_(omp_get_max_threads()),
    _scale  ((depth - 1) * 3),
    _bucket (bucket),
    _depth  (depth),
    _verb   (verbosity),
    _iter   (0),
    _rounds (rounds),
    deg_     (nullptr),
    off_     (nullptr),
    fixes_   (nullptr),
    activ_   (nullptr),
    tags_    (nullptr)
{
#pragma omp parallel
  reallocMemory();
}

/* ------------------------------------ */
Mesh::~Mesh() {

  std::free(activ_);
  std::free(deg_);
  std::free(off_);
  std::free(fixes_);
  std::free(tags_);
}

/* ------------------------------------ */
void Mesh::reallocMemory() {
  // cstdlib
  using std::realloc;

#pragma omp master
  {
    if (_verb)
      std::printf("Mesh memory alloc ... ");
    tic = timer::now();
  }

#pragma omp single
  {
    assert(nb_nodes_);   // expected nb of nodes
    assert(nb_elems_);   // expected nb of elems
    assert(_scale);
    max_node = (size_t) nb_nodes_ * _scale;
    max_elem = (size_t) nb_elems_ * _scale;

  #pragma omp task
    stenc_.resize(max_node);
  #pragma omp task
    vicin_.resize(max_node);
  #pragma omp task
    elems_.resize(max_elem * 3);
  #pragma omp task
    points_.resize(max_node * 2);
  #pragma omp task
    tensor_.resize(max_node * 3);
  #pragma omp task
    solut_.resize(nb_nodes_);  // only for metric calculation
  #pragma omp task
    qualit_.resize(max_elem);
  #pragma omp task
    deg_ = (int*) std::realloc(deg_, max_node * sizeof(int));
  #pragma omp task
    off_ = (int*) std::realloc(off_, nb_cores_ * sizeof(int));
  #pragma omp task
    activ_ = (char*) std::realloc(activ_, max_elem);
  #pragma omp task
    fixes_ = (char*) std::realloc(fixes_, max_node);
  #pragma omp task
    tags_ = (uint8_t*) std::realloc(tags_, max_node);
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
    memory_print += (_bucket * max_node * sizeof(int)); // stenc
    memory_print += (6 * max_node * sizeof(int));      // vicin
    memory_print += (3 * max_elem * sizeof(int));      // elems
    memory_print += (2 * max_node * sizeof(double));   // points
    memory_print += (3 * max_node * sizeof(double));   // tensor
    memory_print += (1 * nb_nodes_ * sizeof(double));   // solut
    memory_print += (1 * max_elem * sizeof(double));   // qualit
    memory_print += (1 * max_node * sizeof(int));      // deg
    memory_print += (1 * nb_cores_ * sizeof(int));      // off
    memory_print += (1 * max_elem);                    // activ
    memory_print += (1 * max_node);                    // fixes
    memory_print += (1 * max_node);                    // tags
    // 1 megabyte: 10^6, 1 mebibyte: 2^20=1048576
    if (_verb) {
      std::printf(
        "%d MB \e[32m(%.2f s)\e[0m\n",
        (int) std::ceil(memory_print / 1e6), (float) timer::elapsed_ms(tic) / 1e3
      );
    } else
      std::printf("\r= Preprocess ...  33 %% =");
    std::fflush(stdout);
  }
}

/* ------------------------------------ */
void Mesh::doFirstTouch() {
#pragma omp barrier

  assert(max_node);
  assert(max_elem);
  assert(nb_cores_);

  const int chunk = nb_nodes_ / nb_cores_;
  const int block = nb_elems_ / nb_cores_;

#pragma omp for schedule(static, 1) nowait
  for (int i = 0; i < nb_cores_; ++i)
    off_[i] = 0;
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < nb_nodes_; ++i)
    solut_[i] = 0.;
#pragma omp for schedule(static, chunk*2) nowait
  for (int i = 0; i < max_node * 2; ++i)
    points_[i] = 0.;
#pragma omp for schedule(static, chunk*3) nowait
  for (int i = 0; i < max_node * 3; ++i)
    tensor_[i] = 0.;

  // (!)
#ifdef DEFERRED_UPDATES
#pragma omp for schedule(static,chunk) nowait
  for(int i=0; i < nb_nodes; ++i)
    stenc[i].resize(64,-1);
#pragma omp for schedule(static,chunk) nowait
  for(int i=nb_nodes; i < max_node; ++i)
    stenc[i].reserve(64);
#else
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < max_node; ++i)
    stenc_[i].resize(_bucket, -1);
#endif

#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < max_node; ++i)
    deg_[i] = 0;
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < max_node; ++i)
    tags_[i] = mask::unset;
#pragma omp for schedule(static, block*3) nowait
  for (int i = 0; i < max_elem * 3; ++i)
    elems_[i] = -1;
#pragma omp for schedule(static, block) nowait
  for (int i = 0; i < max_elem; ++i)
    activ_[i] = 0;
#pragma omp for schedule(static, block)
  for (int i = 0; i < max_elem; ++i)
    qualit_[i] = 0.;
}

/* ------------------------------------ */
void Mesh::rebuildTopology() {

#pragma omp for schedule(static, nb_nodes_/nb_cores_)
  for (int i = 0; i < max_node; ++i)
    deg_[i] = 0;

#pragma omp for
  for (int i = 0; i < nb_elems_; ++i) {
    const int* n = getElem(i);
    if (__builtin_expect(*n < 0, 0))
      continue;

    // manually unrolled
    stenc_[n[0]][sync::fetchAndAdd(deg_ + n[0], 1)] = i;
    stenc_[n[1]][sync::fetchAndAdd(deg_ + n[1], 1)] = i;
    stenc_[n[2]][sync::fetchAndAdd(deg_ + n[2], 1)] = i;
    assert(deg_[n[0]] < stenc_[n[0]].size());
    assert(deg_[n[1]] < stenc_[n[1]].size());
    assert(deg_[n[2]] < stenc_[n[2]].size());
  }

  std::vector<int> heap;

#pragma omp for
  for (int i = 0; i < nb_nodes_; ++i) {
    heap.clear();
    heap.reserve(15);
    assert(deg_[i]);

    for (auto t = stenc_[i].begin(); t < stenc_[i].begin() + deg_[i]; ++t) {
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
    vicin_[i].swap(heap);
    std::sort(stenc_[i].begin(), stenc_[i].begin() + deg_[i]);
#ifdef DEFERRED_UPDATES
    stenc[i].resize(deg[i]);
#endif
  }
}

/* ------------------------------------ */
bool Mesh::verifyTopology() const {

#pragma omp for
  for (int i = 0; i < nb_elems_; ++i) {
    const int* n = getElem(i);
    if (__builtin_expect(*n > -1, 1)) {

      if (std::find(stenc_[n[0]].begin(), stenc_[n[0]].end(), i) == stenc_[n[0]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[0], i, n[0], n[1], n[2]);
        tools::display(stenc_[n[0]]);
      }
      if (std::find(stenc_[n[1]].begin(), stenc_[n[1]].end(), i) == stenc_[n[1]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[1], i, n[0], n[1], n[2]);
        tools::display(stenc_[n[1]]);
      }
      if (std::find(stenc_[n[2]].begin(), stenc_[n[2]].end(), i) == stenc_[n[2]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[2], i, n[0], n[1], n[2]);
        tools::display(stenc_[n[2]]);
      }
      // abort immediately if an error occured
      assert(n[0] not_eq n[1] and n[1] not_eq n[2] and n[2] not_eq n[0]);
      assert(std::find(stenc_[n[0]].begin(), stenc_[n[0]].end(), i) not_eq stenc_[n[0]].end());
      assert(std::find(stenc_[n[1]].begin(), stenc_[n[1]].end(), i) not_eq stenc_[n[1]].end());
      assert(std::find(stenc_[n[2]].begin(), stenc_[n[2]].end(), i) not_eq stenc_[n[2]].end());
    }
  }
  return true;
}

/* ------------------------------------ */
void Mesh::initActivElems() {

#pragma omp for
  for (int i = 0; i < nb_nodes_; ++i)
    activ_[i] = static_cast<char>(__builtin_expect(vicin_[i].empty(), 0) ? 0 : 1);
}

/* ------------------------------------ */
void Mesh::extractPrimalGraph() {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_nodes_; ++i) {
    if (__builtin_expect(activ_[i], 1)) {
      vicin_[i].clear();
      vicin_[i].reserve(deg_[i] * 2);

      for (auto t = stenc_[i].begin(); t < stenc_[i].begin() + deg_[i]; ++t) {
        const int* n = getElem(*t);
        // manually unrolled
        if (i == n[0]) {
          vicin_[i].push_back(n[1]);
          vicin_[i].push_back(n[2]);
          continue;
        }
        if (i == n[1]) {
          vicin_[i].push_back(n[2]);
          vicin_[i].push_back(n[0]);
          continue;
        }
        if (i == n[2]) {
          vicin_[i].push_back(n[0]);
          vicin_[i].push_back(n[1]);
          continue;
        }
      }
      std::sort(vicin_[i].begin(), vicin_[i].end());
      vicin_[i].erase(std::unique(vicin_[i].begin(), vicin_[i].end()), vicin_[i].end());
    }
  }
}

/* ------------------------------------ */
void Mesh::extractDualGraph(Graph* dual) const {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_elems_; ++i) {
    const int* n = getElem(i);
    if (*n < 0)
      continue;

    auto& list = dual->at(i);
    list.clear();

    std::set_intersection(stenc_[n[0]].begin(), stenc_[n[0]].begin() + deg_[n[0]],
                          stenc_[n[1]].begin(), stenc_[n[1]].begin() + deg_[n[1]],
                          std::back_inserter(list));
    std::set_intersection(stenc_[n[0]].begin(), stenc_[n[0]].begin() + deg_[n[0]],
                          stenc_[n[2]].begin(), stenc_[n[2]].begin() + deg_[n[2]],
                          std::back_inserter(list));
    std::set_intersection(stenc_[n[1]].begin(), stenc_[n[1]].begin() + deg_[n[1]],
                          stenc_[n[2]].begin(), stenc_[n[2]].begin() + deg_[n[2]],
                          std::back_inserter(list));

    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
    std::swap(*(std::find(list.begin(), list.end(), i)), list.front());  // don't shift just swap
    assert(list[0] == i);
  }
}
/* ------------------------------------ */
const int* Mesh::getElem(int i) const { return elems_.data() + (i * 3); }
/* ------------------------------------ */
bool Mesh::isActiveNode(int i) const { return !stenc_[i].empty(); }
/* ------------------------------------ */
bool Mesh::isActiveElem(int i) const { return elems_[i * 3] > -1; }
/* ------------------------------------ */
bool Mesh::isBoundary(int i) const { return (tags_[i] & mask::bound); }
/* ------------------------------------ */
bool Mesh::isCorner(int i) const { return (tags_[i] & mask::corner); }
/* ------------------------------------ */
int Mesh::getCapaNode() const { return static_cast<int>(max_node); }
/* ------------------------------------ */
int Mesh::getCapaElem() const { return static_cast<int>(max_elem); }
/* ------------------------------------ */
Patch Mesh::getVicinity(int id, int dist) const {

  // verify step
  assert(dist > 1);
  assert(deg_[id] > 0);

  std::set<int> all[2];
  std::set<int> last[2];
  std::set<int> cur[2];

  // init
  last[0].insert(vicin_[id].begin(), vicin_[id].end());
  all[0].insert(vicin_[id].begin(), vicin_[id].end());
  all[1].insert(stenc_[id].begin(), stenc_[id].begin() + deg_[id]);

  for (int i = 1; i < dist; ++i) {
    cur[0].clear();
    cur[1].clear();
    // retrieve vicin of each vertex
    for (int k : last[0]) {
      cur[0].insert(vicin_[k].begin(), vicin_[k].end());
      cur[1].insert(stenc_[k].begin(), stenc_[k].begin() + deg_[k]);
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

  assert(not stenc_[i].empty());
  for (auto t = stenc_[i].begin(); t < stenc_[i].end() and *t > -1; ++t) {
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
  std::memcpy(elems_.data() + (id * 3), v, sizeof(int) * 3);
}

/* ------------------------------------ */
void Mesh::eraseElem(int id) {
  // nb: std::memsetting with -1 is ok
  std::memset(elems_.data() + (id * 3), -1, sizeof(int) * 3);
}

/* ------------------------------------ */
void Mesh::updateStencil(int i, int t) {
  int k = sync::fetchAndAdd(deg_ + i, 1);
  sync::reallocBucket(stenc_.data(), i, (k + 1), _verb);
  stenc_[i][k] = t;
}

/* ------------------------------------ */
void Mesh::updateStencil(int i, const std::initializer_list<int>& t) {
  int j = sync::fetchAndAdd(deg_ + i, (int) t.size());
  if (stenc_[i].empty()) {
    assert(vicin_[i].size() == 0);
    assert(elems_[i * 3] == -1);
  }
  //assert(stenc[i].size()>0);
  sync::reallocBucket(stenc_.data(), i, j + t.size(), _verb);
  //
  if ((j + t.size()) >= stenc_[i].capacity()) {
    std::fprintf(stderr, "stenc[%d] not fixed, stenc.size: %lu\n", i, stenc_[i].size());
    std::fflush(stderr);
  }

  assert((j + t.size()) < stenc_[i].capacity());
  for (size_t k = 0; k < t.size(); ++k)
    stenc_[i][j + k] = *(t.begin() + k);
}

/* ------------------------------------ */
void Mesh::copyStencil(int i, int j, int nb_rm) {
  int chunk = deg_[i] - nb_rm;
  int k = sync::fetchAndAdd(deg_ + j, chunk);
  // threads attempting to insert will spin until reallocation was done
  sync::reallocBucket(stenc_.data(), i, k + chunk, _verb);

  assert((k + chunk) < stenc_[i].size());
  std::memcpy(stenc_[j].data() + k, stenc_[i].data(), chunk * sizeof(int));
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
  std::memcpy(p    , points_.data() + (n[0] * 2), 2 * sizeof(double));
  std::memcpy(p + 2, points_.data() + (n[1] * 2), 2 * sizeof(double));
  std::memcpy(p + 4, points_.data() + (n[2] * 2), 2 * sizeof(double));
  return n;
}

/* ------------------------------------ */
double Mesh::computeLength(int i, int j) const {

  const double* p1 = points_.data() + (i * 2);
  const double* p2 = points_.data() + (j * 2);
  const double* M1 = tensor_.data() + (i * 3);
  const double* M2 = tensor_.data() + (j * 3);

  return numeric::approxRiemannDist(p1, p2, M1, M2);
}

/* ------------------------------------ */
double Mesh::computeQuality(const int* t) const {

  assert(t[0] > -1 and t[1] > -1 and t[2] > -1);

  // zero-copy
  const double* pa = points_.data() + (t[0] * 2);
  const double* pb = points_.data() + (t[1] * 2);
  const double* pc = points_.data() + (t[2] * 2);
  const double* ma = tensor_.data() + (t[0] * 3);
  const double* mb = tensor_.data() + (t[1] * 3);
  const double* mc = tensor_.data() + (t[2] * 3);

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
  int nb = 0;

#pragma omp single
  nb_activ_elem = 0;

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_elems_; ++i) {
    const int* t = getElem(i);
    if (*t < 0)
      continue;

    qualit_[i] = computeQuality(t);
    if (q_min > qualit_[i]) q_min = qualit_[i];
    if (q_max < qualit_[i]) q_max = qualit_[i];
    q_tot += qualit_[i];
    nb++;
  }

#pragma omp critical
  {
    q[0] = std::min(q[0], q_min);
    q[1] = std::max(q[1], q_max);
    q[2] += q_tot;
    nb_activ_elem += nb;
  }
#pragma omp barrier
#pragma omp single
  q[2] /= nb_activ_elem;

}

/* ------------------------------------ */
void Mesh::computeSteinerPoint(int i, int j, double* p, double* M) const {

  assert(p not_eq nullptr);
  assert(M not_eq nullptr);

  const double* p1 = points_.data() + (i * 2);
  const double* p2 = points_.data() + (j * 2);
  const double* m1 = tensor_.data() + (i * 3);
  const double* m2 = tensor_.data() + (j * 3);

  numeric::computeSteinerPoint(p1, p2, m1, m2, p, M);
}

/* ------------------------------------ */
bool Mesh::isCounterclockwise(const int* t) const {

  const double* pa = points_.data() + (t[0] * 2);
  const double* pb = points_.data() + (t[1] * 2);
  const double* pc = points_.data() + (t[2] * 2);

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
  elem.resize(_bucket);
  int k;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_nodes_; ++i) {
    if (__builtin_expect(fixes_[i], 0)) {
      // reset
      fixes_[i] = 0;
      k = 0;
      if (__builtin_expect(elem.size() < deg_[i], 0))
        elem.resize(deg_[i]);

      for (auto t = stenc_[i].begin(); t < stenc_[i].begin() + deg_[i]; ++t) {
        const int* n = getElem(*t);
        if (__builtin_expect(i == n[0] or i == n[1] or i == n[2], 1))
          elem[k++] = *t;
      }
      if (k >= elem.size())
        std::printf("k: %d, elem.size: %lu\n", k, elem.size());
      assert(k < elem.size());
      // remove duplicates and adjust count
      std::sort(elem.begin(), elem.begin() + k);
      deg_[i] = static_cast<int>(std::unique(elem.begin(), elem.begin() + k) - elem.begin());
      std::fill(elem.begin() + deg_[i], elem.end(), -1);
      // swap arrays O(1)
      stenc_[i].swap(elem);
    }
  }
#endif
}

/* ------------------------------------ */
void Mesh::fixAll() {

#pragma omp for
  for (int i = 0; i < nb_nodes_; ++i)
    fixes_[i] = static_cast<char>(stenc_[i].empty() ? 0 : 1);

  this->fixTagged();
}

} // namespace trinity