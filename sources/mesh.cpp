/* ------------------------------------ */
#include "mesh.h"
#include "tools.h"
/* ------------------------------------ */
namespace trinity {
/* ------------------------------------ */
mesh_t::mesh_t(int size[2], int bucket, int depth, int verbosity, int rounds)
  : nb_nodes(size[0]),
    nb_elems(size[1]),
    nb_cores(omp_get_max_threads()),
    _scale((depth - 1) * 3),
    _bucket(bucket),
    _depth(depth),
    _verb(verbosity),
    _iter(0),
    _rounds(rounds),
    deg(nullptr),
    off(nullptr),
    fixes(nullptr),
    activ(nullptr),
    tags(nullptr) {
#pragma omp parallel
  realloc();
}

/* ------------------------------------ */
mesh_t::~mesh_t() {

  free(activ);
  free(deg);
  free(off);
  free(fixes);
  free(tags);
}

/* ------------------------------------ */
void mesh_t::realloc() {
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
    assert(nb_nodes);   // expected nb of nodes
    assert(nb_elems);   // expected nb of elems
    assert(_scale);
    max_node = (size_t) nb_nodes * _scale;
    max_elem = (size_t) nb_elems * _scale;

#pragma omp task
    stenc.resize(max_node);
#pragma omp task
    vicin.resize(max_node);
#pragma omp task
    elems.resize(max_elem * 3);
#pragma omp task
    points.resize(max_node * 2);
#pragma omp task
    tensor.resize(max_node * 3);
#pragma omp task
    solut.resize(nb_nodes);  // only for metric calculation
#pragma omp task
    qualit.resize(max_elem);
#pragma omp task
    deg = (int*) std::realloc(deg, max_node * sizeof(int));
#pragma omp task
    off = (int*) std::realloc(off, nb_cores * sizeof(int));
#pragma omp task
    activ = (char*) std::realloc(activ, max_elem);
#pragma omp task
    fixes = (char*) std::realloc(fixes, max_node);
#pragma omp task
    tags = (uint8_t*) std::realloc(tags, max_node);
#pragma omp taskwait
  }
  // first-touch
  round_robin_flush();

#ifdef DEFERRED_UPDATES
  init_updates();
#endif

#pragma omp master
  {
    size_t memory_print = 0;
    memory_print += (_bucket * max_node * sizeof(int)); // stenc
    memory_print += (6 * max_node * sizeof(int));      // vicin
    memory_print += (3 * max_elem * sizeof(int));      // elems
    memory_print += (2 * max_node * sizeof(double));   // points
    memory_print += (3 * max_node * sizeof(double));   // tensor
    memory_print += (1 * nb_nodes * sizeof(double));   // solut
    memory_print += (1 * max_elem * sizeof(double));   // qualit
    memory_print += (1 * max_node * sizeof(int));      // deg
    memory_print += (1 * nb_cores * sizeof(int));      // off
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
void mesh_t::round_robin_flush() {
#pragma omp barrier

  assert(max_node);
  assert(max_elem);
  assert(nb_cores);

  const int chunk = nb_nodes / nb_cores;
  const int block = nb_elems / nb_cores;

#pragma omp for schedule(static, 1) nowait
  for (int i = 0; i < nb_cores; ++i)
    off[i] = 0;
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < nb_nodes; ++i)
    solut[i] = 0.;
#pragma omp for schedule(static, chunk*2) nowait
  for (int i = 0; i < max_node * 2; ++i)
    points[i] = 0.;
#pragma omp for schedule(static, chunk*3) nowait
  for (int i = 0; i < max_node * 3; ++i)
    tensor[i] = 0.;

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
    stenc[i].resize(_bucket, -1);
#endif

#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < max_node; ++i)
    deg[i] = 0;
#pragma omp for schedule(static, chunk) nowait
  for (int i = 0; i < max_node; ++i)
    tags[i] = mask::unset;
#pragma omp for schedule(static, block*3) nowait
  for (int i = 0; i < max_elem * 3; ++i)
    elems[i] = -1;
#pragma omp for schedule(static, block) nowait
  for (int i = 0; i < max_elem; ++i)
    activ[i] = 0;
#pragma omp for schedule(static, block)
  for (int i = 0; i < max_elem; ++i)
    qualit[i] = 0.;
}

/* ------------------------------------ */
void mesh_t::rebuild() {

#pragma omp for schedule(static, nb_nodes/nb_cores)
  for (int i = 0; i < max_node; ++i)
    deg[i] = 0;

#pragma omp for
  for (int i = 0; i < nb_elems; ++i) {
    const int* n = get_elem(i);
    if (__builtin_expect(*n < 0, 0))
      continue;

    // manually unrolled
    stenc[n[0]][sync::fetch_and_add(deg + n[0], 1)] = i;
    stenc[n[1]][sync::fetch_and_add(deg + n[1], 1)] = i;
    stenc[n[2]][sync::fetch_and_add(deg + n[2], 1)] = i;
    assert(deg[n[0]] < stenc[n[0]].size());
    assert(deg[n[1]] < stenc[n[1]].size());
    assert(deg[n[2]] < stenc[n[2]].size());
  }

  std::vector<int> heap;

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i) {
    heap.clear();
    heap.reserve(15);
    assert(deg[i]);

    for (auto t = stenc[i].begin(); t < stenc[i].begin() + deg[i]; ++t) {
      const int* n = get_elem(*t);
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
    vicin[i].swap(heap);
    std::sort(stenc[i].begin(), stenc[i].begin() + deg[i]);
#ifdef DEFERRED_UPDATES
    stenc[i].resize(deg[i]);
#endif
  }
}

/* ------------------------------------ */
bool mesh_t::verif() const {

#pragma omp for
  for (int i = 0; i < nb_elems; ++i) {
    const int* n = get_elem(i);
    if (__builtin_expect(*n > -1, 1)) {

      if (std::find(stenc[n[0]].begin(), stenc[n[0]].end(), i) == stenc[n[0]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[0], i, n[0], n[1], n[2]);
        tools::display(stenc[n[0]]);
      }
      if (std::find(stenc[n[1]].begin(), stenc[n[1]].end(), i) == stenc[n[1]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[1], i, n[0], n[1], n[2]);
        tools::display(stenc[n[1]]);
      }
      if (std::find(stenc[n[2]].begin(), stenc[n[2]].end(), i) == stenc[n[2]].end())
#pragma omp critical
      {
        std::fprintf(stderr, "n: %d, t: %d [%d,%d,%d]", n[2], i, n[0], n[1], n[2]);
        tools::display(stenc[n[2]]);
      }
      // abort immediately if an error occured
      assert(n[0] not_eq n[1] and n[1] not_eq n[2] and n[2] not_eq n[0]);
      assert(std::find(stenc[n[0]].begin(), stenc[n[0]].end(), i) not_eq stenc[n[0]].end());
      assert(std::find(stenc[n[1]].begin(), stenc[n[1]].end(), i) not_eq stenc[n[1]].end());
      assert(std::find(stenc[n[2]].begin(), stenc[n[2]].end(), i) not_eq stenc[n[2]].end());
    }
  }
  return true;
}

/* ------------------------------------ */
void mesh_t::init_activ() {

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    activ[i] = static_cast<char>(__builtin_expect(vicin[i].empty(), 0) ? 0 : 1);
}

/* ------------------------------------ */
void mesh_t::extract_primal() {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_nodes; ++i) {
    if (__builtin_expect(activ[i], 1)) {
      vicin[i].clear();
      vicin[i].reserve(deg[i] * 2);

      for (auto t = stenc[i].begin(); t < stenc[i].begin() + deg[i]; ++t) {
        const int* n = get_elem(*t);
        // manually unrolled
        if (i == n[0]) {
          vicin[i].push_back(n[1]);
          vicin[i].push_back(n[2]);
          continue;
        }
        if (i == n[1]) {
          vicin[i].push_back(n[2]);
          vicin[i].push_back(n[0]);
          continue;
        }
        if (i == n[2]) {
          vicin[i].push_back(n[0]);
          vicin[i].push_back(n[1]);
          continue;
        }
      }
      std::sort(vicin[i].begin(), vicin[i].end());
      vicin[i].erase(std::unique(vicin[i].begin(), vicin[i].end()), vicin[i].end());
    }
  }
}

/* ------------------------------------ */
void mesh_t::extract_dual(graph_t* dual) const {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_elems; ++i) {
    const int* n = get_elem(i);
    if (*n < 0)
      continue;

    auto& list = dual->at(i);
    list.clear();

    std::set_intersection(stenc[n[0]].begin(), stenc[n[0]].begin() + deg[n[0]],
                          stenc[n[1]].begin(), stenc[n[1]].begin() + deg[n[1]],
                          std::back_inserter(list));
    std::set_intersection(stenc[n[0]].begin(), stenc[n[0]].begin() + deg[n[0]],
                          stenc[n[2]].begin(), stenc[n[2]].begin() + deg[n[2]],
                          std::back_inserter(list));
    std::set_intersection(stenc[n[1]].begin(), stenc[n[1]].begin() + deg[n[1]],
                          stenc[n[2]].begin(), stenc[n[2]].begin() + deg[n[2]],
                          std::back_inserter(list));

    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
    std::swap(*(std::find(list.begin(), list.end(), i)), list.front());  // don't shift just swap
    assert(list[0] == i);
  }
}

/* ------------------------------------ */
patch_t mesh_t::vicin_dist(int id, int dist) const {

  // verif step
  assert(dist > 1);
  assert(deg[id] > 0);

  std::set<int> all[2];
  std::set<int> last[2];
  std::set<int> cur[2];

  // init
  last[0].insert(vicin[id].begin(), vicin[id].end());
  all[0].insert(vicin[id].begin(), vicin[id].end());
  all[1].insert(stenc[id].begin(), stenc[id].begin() + deg[id]);

  for (int i = 1; i < dist; ++i) {
    cur[0].clear();
    cur[1].clear();
    // retrieve vicin of each vertex
    for (int k : last[0]) {
      cur[0].insert(vicin[k].begin(), vicin[k].end());
      cur[1].insert(stenc[k].begin(), stenc[k].begin() + deg[k]);
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
  patch_t patch;
  patch.node.reserve(all[0].size());
  patch.node.insert(patch.node.begin(), all[0].begin(), all[0].end());
  patch.elem.reserve(all[1].size());
  patch.elem.insert(patch.elem.begin(), all[1].begin(), all[1].end());
  return patch;
}

/* ------------------------------------ */
int mesh_t::elem_neigh(int id, int i, int j) const {

  assert(!stenc[i].empty());
  for (auto t = stenc[i].begin(); t < stenc[i].end() and *t > -1; ++t) {
    if (__builtin_expect(*t not_eq id, 1)) {
      const int* n = get_elem(*t);
      // manually unrolled
      if ((n[0] == j and n[1] == i) or (n[1] == j and n[2] == i) or (n[2] == j and n[0] == i))
        return *t;
    }
  }
  assert(bound(i) and bound(j));
  return -1;
}

/* ------------------------------------ */
void mesh_t::replace_elem(int id, const int* v) {
  assert(v not_eq nullptr);
  std::memcpy(elems.data() + (id * 3), v, sizeof(int) * 3);
}

/* ------------------------------------ */
void mesh_t::erase_elem(int id) {
  // nb: std::memsetting with -1 is ok
  std::memset(elems.data() + (id * 3), -1, sizeof(int) * 3);
}

/* ------------------------------------ */
void mesh_t::update_stenc(int i, int t) {
  int k = sync::fetch_and_add(deg + i, 1);
  sync::check_realloc(stenc.data(), i, (k + 1), _verb);
  stenc[i][k] = t;
}

/* ------------------------------------ */
void mesh_t::update_stenc(int i, const std::initializer_list<int>& t) {
  int j = sync::fetch_and_add(deg + i, t.size());
  if (stenc[i].size() == 0) {
    assert(vicin[i].size() == 0);
    assert(elems[i * 3] == -1);
  }
  //assert(stenc[i].size()>0);
  sync::check_realloc(stenc.data(), i, j + t.size(), _verb);
  //
  if ((j + t.size()) >= stenc[i].capacity()) {
    std::fprintf(stderr, "stenc[%d] not fixed, stenc.size: %lu\n", i, stenc[i].size());
    std::fflush(stderr);
  }

  assert((j + t.size()) < stenc[i].capacity());
  for (size_t k = 0; k < t.size(); ++k)
    stenc[i][j + k] = *(t.begin() + k);
}

/* ------------------------------------ */
void mesh_t::copy_stenc(int i, int j, int nb_rm) {
  int chunk = deg[i] - nb_rm;
  int k = sync::fetch_and_add(deg + j, chunk);
  // threads attempting to insert will spin until reallocation was done
  sync::check_realloc(stenc.data(), i, k + chunk, _verb);

  assert((k + chunk) < stenc[i].size());
  std::memcpy(stenc[j].data() + k, stenc[i].data(), chunk * sizeof(int));
}

/* ------------------------------------ */
#ifdef DEFERRED_UPDATES
/* ------------------------------------ */
void mesh_t::deferred_append(int tid, int i, int t){
  const int key = tools::hash(i)% (def_scale_fact*nb_cores);
  deferred[tid][key].add.push_back(i);
  deferred[tid][key].add.push_back(t);
}
/* ------------------------------------ */
void mesh_t::deferred_append(int tid, int i, const std::initializer_list<int>& t){

  const int key = tools::hash(i)% (def_scale_fact*nb_cores);
  for(size_t k=0; k < t.size(); ++k){
    auto found = std::find(stenc[i].begin(), stenc[i].end(), *(t.begin()+k));
    assert(found == stenc[i].end());
    deferred[tid][key].add.push_back(i);
    deferred[tid][key].add.push_back(*(t.begin()+k));
  }
}
/* ------------------------------------ */
void mesh_t::deferred_remove(int tid, int i, int t){
  const int key = tools::hash(i)% (def_scale_fact*nb_cores);
  auto found = std::find(stenc[i].begin(), stenc[i].end(), t);
  assert(found not_eq stenc[i].end());
  deferred[tid][key].rem.push_back(i);
  deferred[tid][key].rem.push_back(t);
}
/* ------------------------------------ */
void mesh_t::init_updates(){

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
void mesh_t::commit_updates(){

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
    deg[i] = stenc[i].size();

  verif();
}
/* ------------------------------------ */
void mesh_t::reset_updates(){

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
const int* mesh_t::elem_coord(int id, double* p) const {

  // manually unrolled
  const int* n = get_elem(id);
  std::memcpy(p, points.data() + (n[0] * 2), 2 * sizeof(double));
  std::memcpy(p + 2, points.data() + (n[1] * 2), 2 * sizeof(double));
  std::memcpy(p + 4, points.data() + (n[2] * 2), 2 * sizeof(double));
  return n;
}

/* ------------------------------------ */
double mesh_t::edge_length(int i, int j) const {

  const double* p1 = points.data() + (i * 2);
  const double* p2 = points.data() + (j * 2);
  const double* M1 = tensor.data() + (i * 3);
  const double* M2 = tensor.data() + (j * 3);

  return numeric::riemannian_distance(p1, p2, M1, M2);
}

/* ------------------------------------ */
double mesh_t::elem_qual(const int* t) const {

  assert(t[0] > -1 and t[1] > -1 and t[2] > -1);

  // zero-copy
  const double* pa = points.data() + (t[0] * 2);
  const double* pb = points.data() + (t[1] * 2);
  const double* pc = points.data() + (t[2] * 2);
  const double* ma = tensor.data() + (t[0] * 3);
  const double* mb = tensor.data() + (t[1] * 3);
  const double* mc = tensor.data() + (t[2] * 3);

  return numeric::quality(pa, pb, pc, ma, mb, mc);
}

/* ------------------------------------ */
double mesh_t::elem_qual(int id) const {
  return elem_qual(get_elem(id));
}

/* ------------------------------------ */
void mesh_t::compute_quality(double q[3]) {

  double q_min = 3.;
  double q_max = 0.;
  double q_tot = 0.;
  int nb = 0;

#pragma omp single
  nb_activ_elem = 0;

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_elems; ++i) {
    const int* t = get_elem(i);
    if (*t < 0)
      continue;

    qualit[i] = elem_qual(t);
    if (q_min > qualit[i]) q_min = qualit[i];
    if (q_max < qualit[i]) q_max = qualit[i];
    q_tot += qualit[i];
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
void mesh_t::calcul_steiner(int i, int j, double* P, double* M) const {

  assert(P not_eq nullptr);
  assert(M not_eq nullptr);

  const double* p1 = points.data() + (i * 2);
  const double* p2 = points.data() + (j * 2);
  const double* m1 = tensor.data() + (i * 3);
  const double* m2 = tensor.data() + (j * 3);

  numeric::steiner_point(p1, p2, m1, m2, P, M);
}

/* ------------------------------------ */
bool mesh_t::counterclockwise(const int* t) const {

  const double* pa = points.data() + (t[0] * 2);
  const double* pb = points.data() + (t[1] * 2);
  const double* pc = points.data() + (t[2] * 2);

  double s[4];
  s[0] = pb[0] - pa[0];
  s[1] = pc[0] - pa[0];
  s[2] = pb[1] - pa[1];
  s[3] = pc[1] - pa[1];

  return 0.5 * (s[0] * s[3] - s[1] * s[2]) > EPSILON;
}

/* ------------------------------------ */
void mesh_t::fix() {
#ifdef DEFERRED_UPDATES
#pragma omp for schedule(guided)
  for(int i=0; i < nb_nodes; ++i)
    fixes[i] = 0;
  //
  commit_updates();

#else
  std::vector<int> elem;
  elem.resize(_bucket);
  int k;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_nodes; ++i) {
    if (__builtin_expect(fixes[i], 0)) {
      // reset
      fixes[i] = 0;
      k = 0;
      if (__builtin_expect(elem.size() < deg[i], 0))
        elem.resize(deg[i]);

      for (auto t = stenc[i].begin(); t < stenc[i].begin() + deg[i]; ++t) {
        const int* n = get_elem(*t);
        if (__builtin_expect(i == n[0] or i == n[1] or i == n[2], 1))
          elem[k++] = *t;
      }
      if (k >= elem.size())
        std::printf("k: %d, elem.size: %lu\n", k, elem.size());
      assert(k < elem.size());
      // remove duplicates and adjust count
      std::sort(elem.begin(), elem.begin() + k);
      deg[i] = static_cast<int>(std::unique(elem.begin(), elem.begin() + k) - elem.begin());
      std::fill(elem.begin() + deg[i], elem.end(), -1);
      // swap arrays O(1)
      stenc[i].swap(elem);
    }
  }
#endif
}

/* ------------------------------------ */
void mesh_t::fix_all() {

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    fixes[i] = static_cast<char>(stenc[i].empty() ? 0 : 1);

  this->fix();
}

} // namespace trinity