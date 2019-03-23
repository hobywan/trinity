/*
 *                        'coarsening.cpp'
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
/* -------------------------------------------------------------------------- */
#include "trinity/coarsening.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
Coarse::Coarse(Mesh* input, Partit* algo)
  : mesh    (input),
    heuris  (algo),
    cores   (mesh->nb.cores),
    nb_nodes(mesh->nb.nodes),
    nb_elems(mesh->nb.elems),
    nb_indep(algo->task.cardin[0]),
    verbose (mesh->param.verb),
    iter    (mesh->param.iter),
    rounds  (mesh->param.rounds)
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

/* -------------------------------------------------------------------------- */
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

      if (not nb.tasks)
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

/* -------------------------------------------------------------------------- */
void Coarse::identifyTarget(int source) {

  bool bound_source = mesh->isBoundary(source);
  bool bound_destin;

  // 1) filtering step
  std::multimap<double, int> bucket;

  for (auto&& neigh : mesh->topo.vicin[source]) {
    assert(source != neigh);
    bound_destin = mesh->isBoundary(neigh);
    if (__builtin_expect(bound_source and not bound_destin, 0))
      continue;

    auto len = mesh->computeLength(source, neigh);
    if (len < l_min)
      bucket.insert(std::make_pair(len, neigh));
  }

  int destin;
  int face[3];
  double length[3];
  bool skip = true;
  auto const stenc_beg = mesh->topo.stenc[source].begin();
  auto const stenc_end = stenc_beg + mesh->sync.deg[source];

  while (skip and not bucket.empty()) {
    //
    destin = bucket.begin()->second;
    bucket.erase(bucket.begin());

    skip = false;
    int shared = 0;
    bound_destin = mesh->isBoundary(destin);

    // simulate collapse (i->j)
    for (auto t = stenc_beg; not(skip) and t < stenc_end; ++t) {
      auto n = mesh->getElem(*t);

      // a) skip surrounding elems of (i,*it)
      if ((source == n[0] and destin == n[1])
       or (source == n[1] and destin == n[0])
       or (source == n[1] and destin == n[2])
       or (source == n[2] and destin == n[1])
       or (source == n[2] and destin == n[0])
       or (source == n[0] and destin == n[2])) {
        // (!) check if diagonal
        // EDIT: decompose into 2 steps ?
        skip = (++shared > 1 and bound_source and bound_destin);
        continue;
      }

      // b) simulate collapse and eval orientation and lengths
      // manually unrolled
      face[0] = (n[0] == source ? destin : n[0]);
      face[1] = (n[1] == source ? destin : n[1]);
      face[2] = (n[2] == source ? destin : n[2]);

      #ifdef DEBUG
        if (f[0] == f[1] or f[1] == f[2] or f[2] == f[0]) {
          std::fprintf(stderr,
            "error: identify: v: %d, t: %d, f:[%d,%d,%d]\n",
            i, *t, f[0], f[1], f[2]
          );
          tools::display(mesh->topo.stenc[i]);
          if (bound_src) {
            std::fprintf(stderr, "boundary vertex %d\n", i);
          }
          std::fflush(stderr);
          std::exit(EXIT_FAILURE);
        }
      #endif
      assert(face[0] != face[1] and face[1] != face[2] and face[2] != face[0]);

      if (mesh->isCounterclockwise(face)) {
        length[0] = mesh->computeLength(face[0], face[1]);
        length[1] = mesh->computeLength(face[1], face[2]);
        length[2] = mesh->computeLength(face[2], face[0]);
        skip = (length[0] > l_max or length[1] > l_max or length[2] > l_max);
        continue;
      }
      // reset
      skip = true;
    } // end for each elem K in N(i)
  }
  task.target[source] = (skip ? -1 : destin);
}

/* -------------------------------------------------------------------------- */
void Coarse::collapseEdge(int source, int destin) {

#ifdef DEFERRED_UPDATES
  int tid = omp_get_thread_num();
#endif

  auto& stenc = mesh->topo.stenc[source];
  auto& vicin = mesh->topo.vicin[source];

  int k = 0;
  // branches can be factorized for genericity but performance suffers
  if (__builtin_expect(mesh->isBoundary(source), 0)) {
    k++;
    int opp = -1;
    auto trash = stenc.end();
    // retrieve shell elems
    for (auto t = stenc.begin(); t < stenc.end(); ++t) {
      auto n = mesh->getElem(*t);
      // manually unrolled
      if ((source == n[0] and destin == n[1])
       or (destin == n[0] and source == n[1])) {
        opp = n[2];
        trash = t;
        break;
      }
      if ((source == n[1] and destin == n[2])
       or (destin == n[1] and source == n[2])) {
        opp = n[0];
        trash = t;
        break;
      }
      if ((source == n[2] and destin == n[0])
       or (destin == n[2] and source == n[0])) {
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
    // (!) volatile
    std::vector<int>::iterator trash[] = {stenc.end(), stenc.end()};

    // retrieve shell elems
    for (auto t = stenc.begin(); k < 2 and t < stenc.end(); ++t) {
      auto n = mesh->getElem(*t);
      // manually unrolled
      if ((source == n[0] and destin == n[1]) or
          (destin == n[0] and source == n[1])) {
        opp[k] = n[2];
        trash[k++] = t;
        continue;
      }
      if ((source == n[1] and destin == n[2]) or
          (destin == n[1] and source == n[2])) {
        opp[k] = n[0];
        trash[k++] = t;
        continue;
      }
      if ((source == n[2] and destin == n[0]) or
          (destin == n[2] and source == n[0])) {
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
    stenc.erase(trash[1] - 1); // (!) shift: valid because trash[0] < trash[1]
    // mark opposite vertices as to be fixed
    sync::compareAndSwap(sync.fixes + opp[0], 0, 1);
    sync::compareAndSwap(sync.fixes + opp[1], 0, 1);
  }

#ifdef DEFERRED_UPDATES
  for(const int& t : stenc)
    mesh->deferredAppend(tid, j, t);
#else
  // N[j] += (N[i]-shell)
  mesh->copyStencil(source, destin, k);
#endif

  for (auto t = stenc.begin(); t < stenc.end() and *t > -1; ++t) {
    auto n = (int*) mesh->getElem(*t);
    // manually unrolled
         if (n[0] == source) { n[0] = destin; }
    else if (n[1] == source) { n[1] = destin; }
    else if (n[2] == source) { n[2] = destin; }
    //
    assert(n[0] != n[1] and n[1] != n[2] and n[2] != n[0]);
  }

  // mark v[j] as to be fixed too
  sync::compareAndSwap(sync.fixes + destin, 0, 1);
  // mark neighboring nodes for propagation
  for (auto v = vicin.begin(); v < vicin.end(); ++v)
    sync::compareAndSwap(sync.activ + (*v), 0, 1);

  vicin.clear();
  stenc.clear();
}

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
void Coarse::extractSubGraph() {
#pragma omp barrier

#pragma omp for
  for (int i = 0; i < nb.tasks; ++i) {
    const int& k = task.filter[i];
    size_t deg = mesh->topo.vicin[k].size();
    primal[i].resize(deg + 1);
    primal[i][0] = k;
    auto const bytes = sizeof(int) * deg;
    std::memcpy(primal[i].data() + 1, mesh->topo.vicin[k].data(), bytes);
  }
}

/* -------------------------------------------------------------------------- */
void Coarse::processPoints() {

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_indep; ++i) {
    const int& k = task.indep[i];
    collapseEdge(k, task.target[k]);
  }
}

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
void Coarse::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {

      auto rounds = level + 1;
      int const percent[] = {
        nb.activ * 100 / nb_nodes,
        nb.tasks * 100 / nb.activ,
        nb_indep * 100 / nb.tasks
      };
      int const secs = timer::round(time.iter);

      std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, "
                  "%*d filt. \e[0m(%2d %%)\e[0m, "
                  "%*d indep \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  rounds,
                  form[0], nb.activ, percent[0],
                  form[1], nb.tasks, percent[1],
                  form[2], nb_indep, percent[2], secs);
      std::fflush(stdout);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Coarse::recap(int* elap, int* stat, int* form, Stats* total) {
#pragma omp master
  {
    int end = std::max(timer::elapsed_ms(time.start), 1);
    int span = 0;

    if (total != nullptr) {
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

    for (int i = 0; i < 6; ++i) {
      span = std::max(span, elap[i]);
    }

    *form = tools::format(span);

    int const old[] = { nb_nodes, nb_elems };
    int const cur[] = { nb_nodes - stat[2], nb_elems - (2 * stat[2]) };

    if (not verbose) {
      auto percent = (int) std::floor(100 * (++iter) / (4 * rounds + 1));
      std::printf("\r= Remeshing  ... %3d %% =", percent);
    }
    else if (verbose == 1) {

      auto rate  = (int) std::floor(stat[1] / (end * 1e-3));
      auto secs  = (float) end / 1E3;
      auto ratio = (float) (cur[0] - old[0]) * 100 / old[0];

      std::printf("%10d task/sec \e[32m(%4.2f s)", rate, secs);
      std::printf("\e[0m [%.1f %%]\e[0m\n", ratio);

    } else if (verbose == 2) {

      float ratio[2];
      for (int i = 0; i < 2; ++i)
        ratio[i] = (float) (cur[i] - old[i]) * 100 / old[i];

      std::printf("\n\n");
      std::printf("= nodes: %d old, %d new ", old[0], cur[0]);
      std::printf("\e[0m(%.1f %%)\e[0m\n", ratio[0]);
      std::printf("= elems: %d old, %d new ", old[1], cur[1]);
      std::printf("\e[0m(%.1f %%)\e[0m\n", ratio[1]);

      int const rate = (int) std::floor(stat[2] / (end * 1e-3));
      std::printf("= rate : %d merge/sec (%d tasks) \n", rate, stat[2]);

      int step[6];
      for (int i = 0; i < 6; ++i)
        step[i] = elap[i] * 100 / end;

      std::printf("= time per step\n");
      std::printf("  %2d %% vicin  \e[32m(%*d ms)\e[0m\n", step[0], *form, elap[0]);
      std::printf("  %2d %% filter \e[32m(%*d ms)\e[0m\n", step[1], *form, elap[1]);
      std::printf("  %2d %% primal \e[32m(%*d ms)\e[0m\n", step[2], *form, elap[2]);
      std::printf("  %2d %% indep  \e[32m(%*d ms)\e[0m\n", step[3], *form, elap[3]);
      std::printf("  %2d %% kernel \e[32m(%*d ms)\e[0m\n", step[4], *form, elap[4]);
      std::printf("  %2d %% fixes  \e[32m(%*d ms)\e[0m\n", step[5], *form, elap[5]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}

} // namespace trinity
