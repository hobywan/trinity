/*
 *                        'coarsening.cpp'
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
/* --------------------------------------------------------------------------- */
#include "coarsening.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
Coarse::Coarse(Mesh* input, Partit* algo)
  : mesh    (input),
    cores   (mesh->nb.cores),
    nb_nodes(mesh->nb.nodes),
    nb_elems(mesh->nb.elems),
    verbose (mesh->param.verb),
    iter    (mesh->param.iter),
    rounds  (mesh->param.rounds),
    heuris  (algo),
    nb_indep(algo->task.cardin[0])
{
  sync.off   = mesh->sync.off;
  sync.activ = mesh->sync.activ;
  sync.fixes = mesh->sync.fixes;

  task.depth  = 20;
  task.indep  = heuris->task.subset[0];
  task.target = heuris->task.subset[3];
  task.filter = heuris->task.subset[4];
  primal.resize(input->capa.node);
}

/* --------------------------------------------------------------------------- */
Coarse::~Coarse() {}

/* --------------------------------------------------------------------------- */
void Coarse::run(Stats* total) {

  initialize();

  int elap[] = {0, 0, 0, 0, 0, 0};
  int form[] = {0, 0};
  int stat[] = {0, 0, 0};

#pragma omp parallel
  {
    int level = 0;
    std::vector<int> heap;
    heap.reserve((size_t) nb_nodes / cores);

    preProcess();
    timer::save(time.tic, elap);

    do {
      mesh->extractPrimalGraph();
      timer::save(time.tic, elap);

      filterPoints(&heap);
      timer::save(time.tic, elap + 1);

      if (!nb.tasks)
        break;

      extractSubGraph();
      timer::save(time.tic, elap + 2);

      heuris->extractIndepSet(primal, nb.tasks);
      timer::save(time.tic, elap + 3);

      processPoints();
      saveStat(level, stat, form);
      timer::save(time.tic, elap + 4);

      mesh->fixTagged();
      timer::save(time.tic, elap + 5);

      showStat(level++, form);

    } while (nb_indep);

    recap(elap, stat, form, total);
    mesh->verifyTopology();
  }
}

/* --------------------------------------------------------------------------- */
void Coarse::identifyTarget(int i) {

  const auto& vicin = mesh->topo.vicin[i];
  const auto stenc_beg = mesh->topo.stenc[i].begin();
  const auto stenc_end = mesh->topo.stenc[i].begin() + mesh->sync.deg[i];
  const bool bound_src = mesh->isBoundary(i);

  // 1) filterElems step
  std::multimap<double, int> bucket;

  double l[] = {0., 0., 0.};
  for (auto v = vicin.begin(); v < vicin.end(); ++v) {
    assert(i not_eq *v);
    if (__builtin_expect(bound_src and !mesh->isBoundary(*v), 0))
      continue;

    *l = mesh->computeLength(i, *v);
    if (*l < l_min)
      bucket.insert(std::make_pair(*l, *v));
  }

  int j = -1;
  int f[] = {-1, -1, -1};
  int shared = 0;
  bool skip = true;
  bool bound_dst;

  while (skip and not bucket.empty()) {
    //
    j = bucket.begin()->second;
    bucket.erase(bucket.begin());

    skip = false;
    shared = 0;
    bound_dst = mesh->isBoundary(j);

    // simulate collapse (i->j)
    for (auto t = stenc_beg; not(skip) and t < stenc_end; ++t) {
      const int* n = mesh->getElem(*t);

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

      #ifdef DEBUG
        if (f[0] == f[1] or f[1] == f[2] or f[2] == f[0]) {
          std::fprintf(stderr, "error: identify: v: %d, t: %d, f:[%d,%d,%d]\n",
                        i, *t, f[0], f[1], f[2]);
          tools::display(mesh->topo.stenc[i]);
          if (mesh->isBoundary(i)) {
            std::fprintf(stderr, "boundary vertex %d\n", i);
          }
          std::fflush(stderr);
          std::exit(EXIT_FAILURE);
        }
      #endif
      assert(f[0] not_eq f[1] and f[1] not_eq f[2] and f[2] not_eq f[0]);

      if (mesh->isCounterclockwise(f)) {
        l[0] = mesh->computeLength(f[0], f[1]);
        l[1] = mesh->computeLength(f[1], f[2]);
        l[2] = mesh->computeLength(f[2], f[0]);
        skip = (l[0] > l_max or l[1] > l_max or l[2] > l_max);
        continue;
      }
      // reset
      skip = true;
    } // end for each elem K in N(i)
  }
  task.target[i] = (skip ? -1 : j);
}

/* --------------------------------------------------------------------------- */
void Coarse::collapseEdge(int i, int j) {

#ifdef DEFERRED_UPDATES
  int tid = omp_get_thread_num();
#endif

  auto& stenc = mesh->topo.stenc[i];
  auto& vicin = mesh->topo.vicin[i];

  int k = 0;
  // branches can be factorized for genericity but performance suffers
  if (__builtin_expect(mesh->isBoundary(i), 0)) {
    k++;
    int opp = -1;
    auto trash = stenc.end();
    // retrieve shell elems
    for (auto t = stenc.begin(); t < stenc.end(); ++t) {
      const int* n = mesh->getElem(*t);
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
    mesh->deferredRemove(tid, j, *trash);
    mesh->deferredRemove(tid, opp, *trash);
#endif
    // erase common patch
    mesh->eraseElem(*trash);
    stenc.erase(trash);
    // mark opposite vertex as to be fixed
    sync::compareAndSwap(sync.fixes + opp, 0, 1);
  } else {
    int opp[] = {-1, -1};
    std::vector<int>::iterator trash[] = {stenc.end(), stenc.end()}; // (!) volatile

    // retrieve shell elems
    for (auto t = stenc.begin(); k < 2 and t < stenc.end(); ++t) {
      const int* n = mesh->getElem(*t);
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
    mesh->deferredRemove(tid, j, *trash[0]);
    mesh->deferredRemove(tid, j, *trash[1]);
    mesh->deferredRemove(tid, opp[0], *trash[0]);
    mesh->deferredRemove(tid, opp[1], *trash[1]);
#endif
    // erase common patch
    mesh->eraseElem(*trash[0]);
    mesh->eraseElem(*trash[1]);
    stenc.erase(trash[0]);
    stenc.erase(trash[1] - 1);  // (!) shift: valid because trash[0] < trash[1]
    // mark opposite vertices as to be fixed
    sync::compareAndSwap(sync.fixes + opp[0], 0, 1);
    sync::compareAndSwap(sync.fixes + opp[1], 0, 1);
  }

#ifdef DEFERRED_UPDATES
  for(const int& t : stenc)
    mesh->deferredAppend(tid, j, t);
#else
  // N[j] += (N[i]-shell)
  mesh->copyStencil(i, j, k);
#endif

  for (auto t = stenc.begin(); t < stenc.end() and *t > -1; ++t) {
    int* n = (int*) mesh->getElem(*t);
    // manually unrolled
    if (n[0] == i) { n[0] = j; continue; }
    if (n[1] == i) { n[1] = j; continue; }
    if (n[2] == i) { n[2] = j; continue; }
    //
    assert(n[0] not_eq n[1] and n[1] not_eq n[2] and n[2] not_eq n[0]);
  }

  // mark v[j] as to be fixed too
  sync::compareAndSwap(sync.fixes + j, 0, 1);
  // mark neighboring nodes for propagation
  for (auto v = vicin.begin(); v < vicin.end(); ++v)
    sync::compareAndSwap(sync.activ + (*v), 0, 1);

  vicin.clear();
  stenc.clear();
}

/* --------------------------------------------------------------------------- */
void Coarse::preProcess() {
#pragma omp single
  {
    nb.tasks = 0;
    nb_indep = 0;
  }
  // first-touch
#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i) {
    if (__builtin_expect(mesh->isCorner(i), 0))
      sync.activ[i] = -1;
    else
      sync.activ[i] = (char) (mesh->topo.stenc[i].empty() ? 0 : 1);
  }

#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    sync.fixes[i] = 0;
#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    task.target[i] = -1;
#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    task.target[i] = 0;
}

/* --------------------------------------------------------------------------- */
void Coarse::filterPoints(std::vector<int>* heap) {

#pragma omp master
  nb.activ = nb.tasks = nb_indep = 0;

  int count = 0;

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_nodes; ++i) {
    if (__builtin_expect(sync.activ[i] > 0, 1)) {
      count++;
      sync.activ[i] = 0;
      identifyTarget(i);
      if (task.target[i] > -1)
        heap->push_back(i);
    }
  }

  sync::reduceTasks(task.filter, heap, &nb.tasks, sync.off);
  sync::fetchAndAdd(&nb.activ, count);
}

/* --------------------------------------------------------------------------- */
void Coarse::extractSubGraph() {
#pragma omp barrier

#pragma omp for
  for (int i = 0; i < nb.tasks; ++i) {
    const int& k = task.filter[i];
    size_t deg = mesh->topo.vicin[k].size();
    primal[i].resize(deg + 1);
    primal[i][0] = k;
    std::memcpy(primal[i].data() + 1, mesh->topo.vicin[k].data(), sizeof(int) * deg);
  }
}

/* --------------------------------------------------------------------------- */
void Coarse::processPoints() {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_indep; ++i) {
    const int& k = task.indep[i];
    collapseEdge(k, task.target[k]);
  }
}

/* --------------------------------------------------------------------------- */
void Coarse::initialize() {
#pragma omp master
  {
    if (verbose == 1)
      std::printf("%-18s%s", "= coarsening", "...");

    else if (verbose == 2)
      std::printf("Process coarsening ... ");

    std::fflush(stdout);
    time.start = time.iter = time.tic = timer::now();
  }
}

/* --------------------------------------------------------------------------- */
void Coarse::saveStat(int level, int* stat, int* form) {
#pragma omp master
  {
    stat[0] += nb.activ;
    stat[1] += nb.tasks;
    stat[2] += nb_indep;
    if (!level) {
      form[0] = tools::format(nb.activ);
      form[1] = tools::format(nb.tasks);
      form[2] = tools::format(nb_indep);
    }
  }
}

/* --------------------------------------------------------------------------- */
void Coarse::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d filt. \e[0m(%2d %%)\e[0m, "
                  "%*d indep \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  level + 1, form[0], nb.activ, nb.activ * 100 / nb_nodes,
                  form[1], nb.tasks, nb.tasks * 100 / nb.activ,
                  form[2], nb_indep, nb_indep * 100 / nb.tasks, timer::round(time.iter));
      std::fflush(stdout);
    }
  }
}

/* --------------------------------------------------------------------------- */
void Coarse::recap(int* elap, int* stat, int* form, Stats* total) {
#pragma omp master
  {
    int end = std::max(timer::elapsed_ms(time.start), 1);
    int span = 0;

    if (total not_eq nullptr) {
      // reduction
      total->eval += stat[1];
      total->task += stat[2];
      total->elap += end;

      total->step[0] += elap[1];            // filtering
      total->step[1] += elap[0] + elap[2];  // vicin + primal
      total->step[2] += elap[3];            // indep
      total->step[3] += elap[4];            // kernel
      total->step[4] += elap[5];            // repair
    }

    for (int i = 0; i < 6; ++i)
      span = std::max(span, elap[i]);
    *form = tools::format(span);

    int new_nodes = nb_nodes - stat[2];
    int new_elems = nb_elems - (2 * stat[2]);

    if (not verbose) {
      std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100 * (++iter) / (4 * rounds + 1)));
    }
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
      std::printf("  %2d %% vicin  \e[32m(%*d ms)\e[0m\n", elap[0] * 100 / end, *form, elap[0]);
      std::printf("  %2d %% filter \e[32m(%*d ms)\e[0m\n", elap[1] * 100 / end, *form, elap[1]);
      std::printf("  %2d %% primal \e[32m(%*d ms)\e[0m\n", elap[2] * 100 / end, *form, elap[2]);
      std::printf("  %2d %% indep  \e[32m(%*d ms)\e[0m\n", elap[3] * 100 / end, *form, elap[3]);
      std::printf("  %2d %% kernel \e[32m(%*d ms)\e[0m\n", elap[4] * 100 / end, *form, elap[4]);
      std::printf("  %2d %% fixes  \e[32m(%*d ms)\e[0m\n", elap[5] * 100 / end, *form, elap[5]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}

} // namespace trinity
