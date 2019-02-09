/* ------------------------------------*/
#include "coarsening.h"
/* ------------------------------------ */
using namespace trinity;

/* ------------------------------------*/
coarse_t::coarse_t(mesh_t* input, partit_t* algo)
  : mesh(input),
    off(input->off),
    fixes(input->fixes),
    activ(input->activ),
    cores(input->nb_cores),
    nb_nodes(input->nb_nodes),
    nb_elems(input->nb_elems),
    verbose(input->_verb),
    iter(input->_iter),
    rounds(input->_rounds),
    heuris(algo),
    indep(algo->subset[0]),
    target(algo->subset[3]),
    filter(algo->subset[4]),
    nb_indep(algo->card[0]),
    depth(20) {
  primal.resize(input->max_node);
  // for profiling only
  alg = new indep_t(mesh, algo);
}

/* ------------------------------------ */
coarse_t::~coarse_t() {
  delete alg;
}

/* ------------------------------------ */
void coarse_t::run(stats_t* tot) {

  init();

  int time[] = {0, 0, 0, 0, 0, 0};
  int form[] = {0, 0};
  int stat[] = {0, 0, 0};

#pragma omp parallel
  {
    int level = 0;
    std::vector<int> heap;
    heap.reserve((size_t) nb_nodes / cores);

    preprocess();
    timer::save(tic, time);

    do {
      mesh->extract_primal();
      timer::save(tic, time);

      filtering(&heap);
      timer::save(tic, time + 1);

      if (!nb_tasks)
        break;

      extract_primal();
      timer::save(tic, time + 2);

      heuris->indep_subset(primal, nb_tasks);
      timer::save(tic, time + 3);

      kernel();
      save_stat(level, stat, form);
      timer::save(tic, time + 4);

      mesh->fix();
      timer::save(tic, time + 5);

      show_stat(level++, form);

    } while (nb_indep);

    recap(time, stat, form, tot);
    mesh->verif();
  }
}

/* ------------------------------------ */
void coarse_t::identify(int i) {

  const auto& vicin = mesh->vicin[i];
  const auto stenc_beg = mesh->stenc[i].begin();
  const auto stenc_end = mesh->stenc[i].begin() + mesh->deg[i];
  const bool bound_src = mesh->bound(i);

  // 1) filtering step
  std::multimap<double, int> bucket;

  double l[] = {0., 0., 0.};
  for (auto v = vicin.begin(); v < vicin.end(); ++v) {
    assert(i not_eq *v);
    if (__builtin_expect(bound_src and !mesh->bound(*v), 0))
      continue;

    *l = mesh->edge_length(i, *v);
    if (*l < l_min)
      bucket.insert(std::make_pair(*l, *v));
  }

  int j = -1;
  int f[] = {-1, -1, -1};
  int shared = 0;
  bool skip = true;
  bool bound_dst;

  while (skip and !bucket.empty()) {
    //
    j = bucket.begin()->second;
    bucket.erase(bucket.begin());

    skip = false;
    shared = 0;
    bound_dst = mesh->bound(j);

    // simulate collapse (i->j)
    for (auto t = stenc_beg; !skip and t < stenc_end; ++t) {

      const int* n = mesh->get_elem(*t);

      // a) skip surrounding elems of (i,*it)
      if ((i == n[0] and j == n[1]) or (i == n[1] and j == n[0]) or
          (i == n[1] and j == n[2]) or (i == n[2] and j == n[1]) or
          (i == n[2] and j == n[0]) or (i == n[0] and j == n[2])) {
        // (!) check if diagonal
        skip = (++shared > 1 and bound_src and bound_dst);   // decompose into 2 steps ?
        continue;
      }

      // b) simulate collapse and eval orientation and lengths
      // manually unrolled
      f[0] = (n[0] == i ? j : n[0]);
      f[1] = (n[1] == i ? j : n[1]);
      f[2] = (n[2] == i ? j : n[2]);
      if (not (f[0] not_eq f[1] and f[1] not_eq f[2] and f[2] not_eq f[0])) {
        std::printf("problem identify: v: %d, t: %d, f:[%d,%d,%d]\n", i, *t, f[0], f[1], f[2]);
        tools::display(mesh->stenc[i]);
        if (mesh->bound(i))
          std::printf("boundary vertex %d\n", i);
      }

      assert(f[0] not_eq f[1] and f[1] not_eq f[2] and f[2] not_eq f[0]);

      if (mesh->counterclockwise(f)) {
        l[0] = mesh->edge_length(f[0], f[1]);
        l[1] = mesh->edge_length(f[1], f[2]);
        l[2] = mesh->edge_length(f[2], f[0]);
        skip = (l[0] > l_max or l[1] > l_max or l[2] > l_max);
        continue;
      }
      // reset
      skip = true;
    } // end for each elem K in N(i)
  }
  target[i] = (skip ? -1 : j);
}

/* ------------------------------------ */
void coarse_t::collapse(int i, int j) {

#ifdef DEFERRED_UPDATES
  int tid = omp_get_thread_num();
#endif

  auto& stenc = mesh->stenc[i];
  auto& vicin = mesh->vicin[i];

  int k = 0;
  // branches can be factorized for genericity but performance suffers
  if (__builtin_expect(mesh->bound(i), 0)) {
    k++;
    int opp = -1;
    auto trash = stenc.end();
    // retrieve shell elems
    for (auto t = stenc.begin(); t < stenc.end(); ++t) {
      const int* n = mesh->get_elem(*t);
      // manually unrolled
      if ((i == n[0] and j == n[1]) or (j == n[0] and i == n[1])) {
        opp = n[2];
        trash = t;
        break;
      }
      if ((i == n[1] and j == n[2]) or (j == n[1] and i == n[2])) {
        opp = n[0];
        trash = t;
        break;
      }
      if ((i == n[2] and j == n[0]) or (j == n[2] and i == n[0])) {
        opp = n[1];
        trash = t;
        break;
      }
    }
#ifdef DEFERRED_UPDATES
    mesh->deferred_remove(tid,  j, *trash);
    mesh->deferred_remove(tid,opp, *trash);
#endif
    // erase common patch
    mesh->erase_elem(*trash);
    stenc.erase(trash);
    // mark opposite vertex as to be fixed
    sync::compare_and_swap(fixes + opp, 0, 1);
  } else {
    int opp[] = {-1, -1};
    std::vector<int>::iterator trash[] = {stenc.end(), stenc.end()}; // (!) volatile

    // retrieve shell elems
    for (auto t = stenc.begin(); k < 2 and t < stenc.end(); ++t) {
      const int* n = mesh->get_elem(*t);
      // manually unrolled
      if ((i == n[0] and j == n[1]) or (j == n[0] and i == n[1])) {
        opp[k] = n[2];
        trash[k++] = t;
        continue;
      }
      if ((i == n[1] and j == n[2]) or (j == n[1] and i == n[2])) {
        opp[k] = n[0];
        trash[k++] = t;
        continue;
      }
      if ((i == n[2] and j == n[0]) or (j == n[2] and i == n[0])) {
        opp[k] = n[1];
        trash[k++] = t;
        continue;
      }
    }
    assert(k == 2);
#ifdef DEFERRED_UPDATES
    mesh->deferred_remove(tid, j, *trash[0]);
    mesh->deferred_remove(tid, j, *trash[1]);
    mesh->deferred_remove(tid, opp[0], *trash[0]);
    mesh->deferred_remove(tid, opp[1], *trash[1]);
#endif
    // erase common patch
    mesh->erase_elem(*trash[0]);
    mesh->erase_elem(*trash[1]);
    stenc.erase(trash[0]);
    stenc.erase(trash[1] - 1);  // (!) shift: valid because trash[0] < trash[1]
    // mark opposite vertices as to be fixed
    sync::compare_and_swap(fixes + opp[0], 0, 1);
    sync::compare_and_swap(fixes + opp[1], 0, 1);
  }

#ifdef DEFERRED_UPDATES
  for(const int& t : stenc)
    mesh->deferred_append(tid, j, t);
#else
  // N[j] += (N[i]-shell)
  mesh->copy_stenc(i, j, k);
#endif

  for (auto t = stenc.begin(); t < stenc.end() and *t > -1; ++t) {
    int* n = (int*) mesh->get_elem(*t);
    // manually unrolled
    if (n[0] == i) { n[0] = j; continue; }
    if (n[1] == i) { n[1] = j; continue; }
    if (n[2] == i) { n[2] = j; continue; }
    //
    assert(n[0] not_eq n[1] and n[1] not_eq n[2] and n[2] not_eq n[0]);
  }

  // mark v[j] as to be fixed too
  sync::compare_and_swap(fixes + j, 0, 1);
  // mark neighboring nodes for propagation
  for (auto v = vicin.begin(); v < vicin.end(); ++v)
    sync::compare_and_swap(activ + (*v), 0, 1);

  vicin.clear();
  stenc.clear();
}

/* ------------------------------------ */
void coarse_t::preprocess() {
#pragma omp single
  {
    nb_tasks = 0;
    nb_indep = 0;
  }
  // first-touch
#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i) {
    if (__builtin_expect(mesh->corner(i), 0))
      activ[i] = -1;
    else
      activ[i] = static_cast<char>(mesh->stenc[i].empty() ? 0 : 1);
  }

#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    fixes[i] = 0;
#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    target[i] = -1;
#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    filter[i] = 0;
}

/* ------------------------------------ */
void coarse_t::filtering(std::vector<int>* heap) {

#pragma omp master
  nb_activ = nb_tasks = nb_indep = 0;

  int count = 0;

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_nodes; ++i) {
    if (__builtin_expect(activ[i] > 0, 1)) {
      count++;
      activ[i] = 0;
      identify(i);
      if (target[i] > -1)
        heap->push_back(i);
    }
  }

  sync::task_reduction(filter, heap, &nb_tasks, off);
  sync::fetch_and_add(&nb_activ, count);
}

/* ------------------------------------ */
void coarse_t::extract_primal() {
#pragma omp barrier

#pragma omp for
  for (int i = 0; i < nb_tasks; ++i) {
    const int& k = filter[i];
    size_t deg = mesh->vicin[k].size();
    primal[i].resize(deg + 1);
    primal[i][0] = k;
    std::memcpy(primal[i].data() + 1, mesh->vicin[k].data(), sizeof(int) * deg);
  }
}

/* ------------------------------------ */
void coarse_t::kernel() {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_indep; ++i) {
    const int& k = indep[i];
    collapse(k, target[k]);
  }
}

/* ------------------------------------ */
void coarse_t::init() {
#pragma omp master
  {
    if (verbose == 1)
      std::printf("%-18s%s", "= coarsening", "...");

    else if (verbose == 2)
      std::printf("Process coarsening ... ");

    std::fflush(stdout);
    start = round = tic = timer::now();
  }
}

/* ------------------------------------ */
void coarse_t::save_stat(int level, int* stat, int* form) {
#pragma omp master
  {
    stat[0] += nb_activ;
    stat[1] += nb_tasks;
    stat[2] += nb_indep;
    if (!level) {
      form[0] = tools::format(nb_activ);
      form[1] = tools::format(nb_tasks);
      form[2] = tools::format(nb_indep);
    }
  }
}

/* ------------------------------------ */
void coarse_t::show_stat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d filt. \e[0m(%2d %%)\e[0m, "
                  "%*d indep \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  level + 1, form[0], nb_activ, nb_activ * 100 / nb_nodes,
                  form[1], nb_tasks, nb_tasks * 100 / nb_activ,
                  form[2], nb_indep, nb_indep * 100 / nb_tasks, timer::round(round));
      std::fflush(stdout);
    }
  }
}

/* ------------------------------------ */
void coarse_t::recap(int* time, int* stat, int* form, stats_t* tot) {
#pragma omp master
  {
    int end = timer::elapsed_ms(start);
    int span = 0;

    // reduction
    tot->eval += stat[1];
    tot->task += stat[2];
    tot->elap += end;

    tot->step[0] += time[1];          // filter
    tot->step[1] += time[0] + time[2];  // vicin + primal
    tot->step[2] += time[3];          // indep
    tot->step[3] += time[4];          // kernel
    tot->step[4] += time[5];          // repair


    for (int i = 0; i < 6; ++i)
      span = std::max(span, time[i]);
    *form = tools::format(span);

    int new_nodes = nb_nodes - stat[2];
    int new_elems = nb_elems - (2 * stat[2]);

    if (!verbose)
      std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100 * (++iter) / (4 * rounds + 1)));

    else if (verbose == 1) {
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m [%.1f %%]\e[0m\n",
                  (int) std::floor(stat[1] / (end * 1e-3)), (float) end / 1e3,
                  (float) (new_nodes - nb_nodes) * 100 / nb_elems);
    } else if (verbose == 2) {
      std::printf("\n\n");
      std::printf("= nodes: %d old, %d new \e[0m(%.1f %%)\e[0m\n", nb_nodes, new_nodes,
                  (float) (new_nodes - nb_nodes) * 100 / nb_nodes);
      std::printf("= elems: %d old, %d new \e[0m(%.1f %%)\e[0m\n", nb_elems, new_elems,
                  (float) (new_elems - nb_elems) * 100 / nb_elems);
      std::printf("= rate : %d merge/sec (%d tasks) \n", (int) std::floor(stat[2] / (end * 1e-3)), stat[2]);
      std::printf("= time per step\n");
      std::printf("  %2d %% vicin  \e[32m(%*d ms)\e[0m\n", time[0] * 100 / end, *form, time[0]);
      std::printf("  %2d %% filter \e[32m(%*d ms)\e[0m\n", time[1] * 100 / end, *form, time[1]);
      std::printf("  %2d %% primal \e[32m(%*d ms)\e[0m\n", time[2] * 100 / end, *form, time[2]);
      std::printf("  %2d %% indep  \e[32m(%*d ms)\e[0m\n", time[3] * 100 / end, *form, time[3]);
      std::printf("  %2d %% kernel \e[32m(%*d ms)\e[0m\n", time[4] * 100 / end, *form, time[4]);
      std::printf("  %2d %% fixes  \e[32m(%*d ms)\e[0m\n", time[5] * 100 / end, *form, time[5]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}
