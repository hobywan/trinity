/*
 *                          'swapping.cpp'
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

#include <trinity.h>
#include "trinity/swapping.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
Swap::Swap(Mesh* input)
  : mesh    (input),
    cores   (mesh->nb.cores),
    nb_nodes(mesh->nb.nodes),
    nb_elems(mesh->nb.elems),
    verbose (mesh->param.verb),
    iter    (mesh->param.iter),
    rounds  (mesh->param.rounds)
{
  geom.qualit = mesh->geom.qualit.data();
  sync.off    = mesh->sync.off;
  sync.fixes  = mesh->sync.fixes;
  sync.activ  = mesh->sync.activ;

  size_t capacity = input->capa.elem;
  task.match = new int[capacity];
  task.list  = new int[capacity];
  task.depth = 20;

  dual.resize(capacity);
  heuris.initialize(capacity, task.match, sync.off);
}

/* -------------------------------------------------------------------------- */
Swap::~Swap() {
  delete[] task.match;
  delete[] task.list;
}

/* -------------------------------------------------------------------------- */
int Swap::swap(int id1, int id2, int index) {

  int j, k;
  int cur1[] = {-1, -1, -1};
  int cur2[] = {-1, -1, -1};

  const int* old1 = mesh->getElem(id1);
  const int* old2 = mesh->getElem(id2);

  __builtin_prefetch(old1, 1);
  __builtin_prefetch(old2, 1);

  // rotate (t1,t2) and find shared edge/opposite vertices
  for (int i = 0; i < 3; ++i) {
    j = (i + 1) % 3;
    k = (i + 2) % 3;

    // manually unrolled
    if (old1[i] == old2[1] and old1[j] == old2[0]) {
      cur1[0] = old1[i];
      cur1[1] = old2[2];
      cur1[2] = old1[k];
      cur2[0] = old1[k];
      cur2[1] = old2[2];
      cur2[2] = old1[j];
      break;
    }
    if (old1[i] == old2[2] and old1[j] == old2[1]) {
      cur1[0] = old1[i];
      cur1[1] = old2[0];
      cur1[2] = old1[k];
      cur2[0] = old1[k];
      cur2[1] = old2[0];
      cur2[2] = old1[j];
      break;
    }
    if (old1[i] == old2[0] and old1[j] == old2[2]) {
      cur1[0] = old1[i];
      cur1[1] = old2[1];
      cur1[2] = old1[k];
      cur2[0] = old1[k];
      cur2[1] = old2[1];
      cur2[2] = old1[j];
      break;
    }
  }
  // check inversions
  if (not mesh->isCounterclockwise(cur1) or not mesh->isCounterclockwise(cur2))
    return 0;

  // eval quality improvement
  const double q[] = { mesh->computeQuality(cur1), mesh->computeQuality(cur2) };
  const double q_old = std::min(geom.qualit[id1], geom.qualit[id2]);
  const double q_new = std::min(q[0], q[1]);

  if (q_new < q_old)
    return 0;

  // update mesh
  mesh->replaceElem(id1, cur1);
  mesh->replaceElem(id2, cur2);

  #if DEFER_UPDATES
    int tid = omp_get_thread_num();
    mesh->deferredRemove(tid, cur1[0], id2);
    mesh->deferredRemove(tid, cur2[2], id1);
    mesh->deferredAppend(tid, cur1[2], id2); //  N(s[2],k2)
    mesh->deferredAppend(tid, cur2[1], id1); //  N(s[3],k1)
  #else
    mesh->updateStencil(cur2[0], id2); //  N(s[2],k2)
    mesh->updateStencil(cur2[1], id1); //  N(s[3],k1)
  #endif

  // mark nodes as to be fixed
  sync::compareAndSwap(sync.fixes + cur1[0], 0, 1);
  sync::compareAndSwap(sync.fixes + cur1[1], 0, 1);
  sync::compareAndSwap(sync.fixes + cur1[2], 0, 1);
  sync::compareAndSwap(sync.fixes + cur2[2], 0, 1);

  // update cached qualit
  geom.qualit[id1] = q[0];
  geom.qualit[id2] = q[1];

  // propagate
  for (const int& t : dual[index]) {
    sync::compareAndSwap(sync.activ + t, 0, 1);
  }

  if (task.match[id2] > -1) {
    for (const int& t : dual[task.match[id2]]) {
      sync::compareAndSwap(sync.activ + t, 0, 1);
    }
  }

  return 2;
}

/* -------------------------------------------------------------------------- */
void Swap::cacheQuality() {

#pragma omp single
  nb.tasks = nb.commit = 0;

#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    sync.fixes[i] = 0;
#pragma omp for nowait
  for (int i = 0; i < nb_elems; ++i)
    task.match[i] = -1;
#pragma omp for
  for (int i = 0; i < nb_elems; ++i)
    sync.activ[i] = 0;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_elems; ++i) {
    const int* n = mesh->getElem(i);
    if (__builtin_expect(*n > -1, 1)) {
      geom.qualit[i] = mesh->computeQuality(i);
      sync.activ[i] = 1;
    }
  }

#pragma omp master
  time.iter = timer::now();
}

/* -------------------------------------------------------------------------- */
void Swap::filterElems(std::vector<int>* heap) {

#pragma omp master
  nb.activ = 0;

  int count = 0;
  heap->reserve((size_t) nb_elems / cores);

//#pragma omp for schedule(dynamic,chunk) nowait
#pragma omp for nowait
  for (int i = 0; i < nb_elems; ++i) {
    if (__builtin_expect(sync.activ[i], 0)) {
      sync.activ[i] = 0;
      count++;
      if (geom.qualit[i] < trinity::q_min)
        heap->push_back(i);
    }
  }
  sync::reduceTasks(task.list, heap, &nb.tasks, sync.off);
  sync::fetchAndAdd(&nb.activ, count);
}

/* -------------------------------------------------------------------------- */
void Swap::extractDualGraph() {

  const auto& stenc = mesh->topo.stenc;
  const auto* deg = mesh->sync.deg;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb.tasks; ++i) {

    const int& k = task.list[i];
    const int* n = mesh->getElem(k);

    dual[i].clear();
    dual[i].reserve(4);
    dual[i].push_back(k);

    for (auto t = stenc[*n].begin();
         dual[i].size() < 3 and t < stenc[*n].begin() + deg[*n]; ++t) {
      if (__builtin_expect(*t == k, 0))
        continue;

      // manually unrolled
      const int* v = mesh->getElem(*t);
      if ((v[0] == n[0] and v[1] == n[2]) or (v[0] == n[1] and v[1] == n[0]) or
          (v[1] == n[0] and v[2] == n[2]) or (v[1] == n[1] and v[2] == n[0]) or
          (v[2] == n[0] and v[0] == n[2]) or (v[2] == n[1] and v[0] == n[0])) {
        dual[i].push_back(*t);
      }
    }

    for (auto t = stenc[n[2]].begin(); t < stenc[n[2]].begin() + deg[n[2]]; ++t) {
      if (__builtin_expect(*t == k, 0))
        continue;

      // manually unrolled
      const int* v = mesh->getElem(*t);
      if ((v[0] == n[2] and v[1] == n[1]) or
          (v[1] == n[2] and v[2] == n[1]) or
          (v[2] == n[2] and v[0] == n[1])) {
        dual[i].push_back(*t);
        break;
      }
    }
    task.match[k] = i;
  }
}

/* -------------------------------------------------------------------------- */
void Swap::processFlips() {

  int succ = 0;

#pragma omp single
  nb.commit = 0;

#pragma omp for schedule(guided) nowait
  for (int index = 0; index < nb.tasks; ++index) {
    const int& k1 = dual[index][0];
    const int& k2 = heuris.task.matched[k1];
    if (__builtin_expect(k2 > -1, 1))
      succ += swap(k1, k2, index);
  }

  // update nb.commit
  sync::fetchAndAdd(&nb.commit, succ);
#pragma omp barrier
}

/* -------------------------------------------------------------------------- */
void Swap::initialize() {
#pragma omp master
  {
    if (verbose == 1)
      std::printf("%-18s%s", "= swapping", "...");

    else if (verbose == 2)
      std::printf("Process swapping ... ");

    std::fflush(stdout);
    time.start = time.iter = time.tic = timer::now();
  }
}

/* -------------------------------------------------------------------------- */
void Swap::saveStat(int level, int* stat, int* form) {
#pragma omp master
  {
    stat[0] += nb.activ;
    stat[1] += nb.tasks;
    stat[2] += nb.commit;
    if (!level) {
      form[0] = tools::format(nb.activ);
      form[1] = tools::format(nb.tasks);
      form[2] = tools::format(nb.commit);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Swap::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      int const round = level + 1;
      int const percent[] = {
        nb.activ  * 100 / nb_nodes,
        nb.tasks  * 100 / nb.activ,
        nb.commit * 100 / nb.tasks
      };
      int const secs = timer::round(time.iter);

      std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d"
                  " filt. \e[0m(%2d %%)\e[0m, "
                  "%*d comm. \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  round,
                  form[0], nb.activ, percent[0],
                  form[1], nb.tasks, percent[1],
                  form[2], nb.commit, percent[2], secs);
      std::fflush(stdout);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Swap::recap(int* elap, int* stat, int* form, Stats* total) {
#pragma omp master
  {
    int end = std::max(timer::elapsed_ms(time.start), 1);

    if (total != nullptr) {
      total->eval += stat[0];
      total->task += stat[1];
      total->elap += end;
      // manually unrolled
      total->step[0] += elap[0] + elap[1];  // qualit + filtering
      total->step[1] += elap[2];            // dual
      total->step[2] += elap[3];            // match
      total->step[3] += elap[4];            // kernel
      total->step[4] += elap[5];            // repair
    }

    int span = 0;
    for (int i = 0; i < 6; ++i)
      span = std::max(span, elap[i]);
    *form = tools::format(span);

    if (not verbose) {
      auto percent = (int) std::floor(100 * (++iter) / (4 * rounds + 1));
      std::printf("\r= Remeshing  ... %3d %% =", percent);

    } else if (verbose == 1) {
      auto rate = (int) std::floor(stat[1] / (end * 1e-3));
      auto secs = (float) end / 1E3;
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m\n", rate, secs);

    } else if (verbose == 2) {
      int step[6];
      auto const rate = (int) std::floor(stat[1] / (end * 1e-3));
      for (int i = 0; i < 6; ++i)
        step[i] = elap[i] * 100 / end;

      std::printf("\n\n");
      std::printf("= rate : %d flip/sec (%d tasks) \n", rate, stat[1]);

      std::printf("= time per step\n");
      std::printf("= time per step\n");
      std::printf("  %2d %% qualit \e[32m(%*d ms)\e[0m\n", step[0], *form, elap[0]);
      std::printf("  %2d %% filter \e[32m(%*d ms)\e[0m\n", step[1], *form, elap[1]);
      std::printf("  %2d %% dual   \e[32m(%*d ms)\e[0m\n", step[2], *form, elap[2]);
      std::printf("  %2d %% match  \e[32m(%*d ms)\e[0m\n", step[3], *form, elap[3]);
      std::printf("  %2d %% kernel \e[32m(%*d ms)\e[0m\n", step[4], *form, elap[4]);
      std::printf("  %2d %% fixes  \e[32m(%*d ms)\e[0m\n", step[5], *form, elap[5]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}

/* -------------------------------------------------------------------------- */
void Swap::run(Stats* total) {

  initialize();

  int elap[] = {0, 0, 0, 0, 0, 0};
  int form[] = {0, 0, 0};
  int stat[] = {0, 0, 0};

#pragma omp parallel
  {
    std::vector<int> heap;

    cacheQuality();
    timer::save(time.tic, elap);

    int level = 0;
    do {

      filterElems(&heap);
      timer::save(time.tic, elap + 1);

      if (!nb.tasks)
        break;

      extractDualGraph();
      timer::save(time.tic, elap + 2);

      heuris.computeGreedyMatching(dual, nb.tasks);
      heuris.getRatio(dual, nb.tasks, &nb.count);
      timer::save(time.tic, elap + 3);

      processFlips();
      saveStat(level, stat, form);
      timer::save(time.tic, elap + 4);

      mesh->fixTagged();
      timer::save(time.tic, elap + 5);

      showStat(level++, form);

    } while (nb.commit);

    recap(elap, stat, form, total);
    mesh->verifyTopology();
  }
}
/* -------------------------------------------------------------------------- */
} // namespace trinity
