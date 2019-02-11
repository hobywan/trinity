/*
 *                          'swapping.cpp'
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

#include "swapping.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
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
  heuris.init(capacity, task.match, sync.off);
}

/* --------------------------------------------------------------------------- */
Swap::~Swap() {
  delete[] task.match;
  delete[] task.list;
}

/* --------------------------------------------------------------------------- */
void Swap::run(Stats* tot) {

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

      heuris.computeKarpSipser(dual, nb.tasks);
      heuris.getRatio(dual, nb.tasks, &nb.count);
      timer::save(time.tic, elap + 3);

      processFlips();
      saveStat(level, stat, form);
      timer::save(time.tic, elap + 4);

      mesh->fixTagged();
      timer::save(time.tic, elap + 5);

      showStat(level++, form);

    } while (nb.commit);

    recap(elap, stat, form, tot);
    mesh->verifyTopology();
  }
}

/* --------------------------------------------------------------------------- */
int Swap::swap(int k1, int k2, int index) {

  int j, k;
  int f1[] = {-1, -1, -1};
  int f2[] = {-1, -1, -1};

  const int* t1 = mesh->getElem(k1);
  const int* t2 = mesh->getElem(k2);

  __builtin_prefetch(t1, 1);
  __builtin_prefetch(t2, 1);

  // rotate (t1,t2) and find shared edge/opposite vertices
  for (int i = 0; i < 3; ++i) {
    j = (i + 1) % 3;
    k = (i + 2) % 3;

    // manually unrolled
    if (t1[i] == t2[1] and t1[j] == t2[0]) {
      f1[0] = t1[i];
      f1[1] = t2[2];
      f1[2] = t1[k];
      f2[0] = t1[k];
      f2[1] = t2[2];
      f2[2] = t1[j];
      break;
    }
    if (t1[i] == t2[2] and t1[j] == t2[1]) {
      f1[0] = t1[i];
      f1[1] = t2[0];
      f1[2] = t1[k];
      f2[0] = t1[k];
      f2[1] = t2[0];
      f2[2] = t1[j];
      break;
    }
    if (t1[i] == t2[0] and t1[j] == t2[2]) {
      f1[0] = t1[i];
      f1[1] = t2[1];
      f1[2] = t1[k];
      f2[0] = t1[k];
      f2[1] = t2[1];
      f2[2] = t1[j];
      break;
    }
  }
  // check inversions
  if (not mesh->isCounterclockwise(f1) or not mesh->isCounterclockwise(f2))
    return 0;

  // eval computeQuality improvement
  const double q[]   = { mesh->computeQuality(f1), mesh->computeQuality(f2) };
  const double q_old = std::min(geom.qualit[k1], geom.qualit[k2]);
  const double q_new = std::min(q[0], q[1]);

  if (q_new < q_old)
    return 0;

  // update mesh
  mesh->replaceElem(k1, f1);
  mesh->replaceElem(k2, f2);

  #ifdef DEFERRED_UPDATES
    int tid = omp_get_thread_num();
    mesh->deferredRemove(tid, f1[0], k2);
    mesh->deferredRemove(tid, f2[2], k1);
    mesh->deferredAppend(tid, f1[2], k2); //  N(s[2],k2)
    mesh->deferredAppend(tid, f2[1], k1); //  N(s[3],k1)
  #else
    mesh->updateStencil(f2[0], k2); //  N(s[2],k2)
    mesh->updateStencil(f2[1], k1); //  N(s[3],k1)
  #endif

  // mark nodes as to be fixed
  sync::compareAndSwap(sync.fixes + f1[0], 0, 1);
  sync::compareAndSwap(sync.fixes + f1[1], 0, 1);
  sync::compareAndSwap(sync.fixes + f1[2], 0, 1);
  sync::compareAndSwap(sync.fixes + f2[2], 0, 1);

  // update cached qualit
  geom.qualit[k1] = q[0];
  geom.qualit[k2] = q[1];

  // propagate
  for (const int& t : dual[index]) {
    sync::compareAndSwap(sync.activ + t, 0, 1);
  }

  if (task.match[k2] > -1) {
    for (const int& t : dual[task.match[k2]]) {
      sync::compareAndSwap(sync.activ + t, 0, 1);
    }
  }

  return 2;
}

/* --------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------- */
void Swap::extractDualGraph() {

  const auto& stenc = mesh->topo.stenc;
  const int* deg = mesh->sync.deg;

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

/* --------------------------------------------------------------------------- */
void Swap::processFlips() {

  int succ = 0;

#pragma omp single
  nb.commit = 0;

#pragma omp for schedule(guided) nowait
  for (int index = 0; index < nb.tasks; ++index) {
    const int& k1 = dual[index][0];
    const int& k2 = heuris.matched[k1];
    if (__builtin_expect(k2 > -1, 1))
      succ += swap(k1, k2, index);
  }

  // update nb.commit
  sync::fetchAndAdd(&nb.commit, succ);
#pragma omp barrier
}

/* --------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------- */
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

/* --------------------------------------------------------------------------- */
void Swap::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d filt. \e[0m(%2d %%)\e[0m, "
                  "%*d comm. \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  level + 1, form[0], nb.activ, (int) (nb.activ * 100 / nb_elems),
                  form[1], nb.tasks, (int) (nb.tasks * 100 / nb.activ),
                  form[2], nb.commit, (int) (nb.commit * 100 / nb.tasks), timer::round(time.iter));
      std::fflush(stdout);
    }
  }
}

/* --------------------------------------------------------------------------- */
void Swap::recap(int* elap, int* stat, int* form, Stats* tot) {
#pragma omp master
  {
    int end = std::max(timer::elapsed_ms(time.start), 1);

    tot->eval += stat[0];
    tot->task += stat[1];
    tot->elap += end;

    // manually unrolled
    tot->step[0] += elap[0] + elap[1];  // qualit + filterElems
    tot->step[1] += elap[2];          // dual
    tot->step[2] += elap[3];          // match
    tot->step[3] += elap[4];          // processFlips
    tot->step[4] += elap[5];          // repair

    int span = 0;
    for (int i = 0; i < 6; ++i)
      span = std::max(span, elap[i]);
    *form = tools::format(span);

    if (!verbose)
      std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100 * (++iter) / (4 * rounds + 1)));

    else if (verbose == 1) {
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m\n", (int) std::floor(stat[1] / (end * 1e-3)), (float) end / 1e3);
    } else if (verbose == 2) {
      std::printf("\n\n");
      std::printf("= rate : %d flip/sec (%d tasks) \n", (int) std::floor(stat[1] / (end * 1e-3)), stat[1]);
      std::printf("= time per step\n");
      std::printf("  %2d %% qualit \e[32m(%*d ms)\e[0m\n", (int) elap[0] * 100 / end, *form, elap[0]);
      std::printf("  %2d %% filter \e[32m(%*d ms)\e[0m\n", (int) elap[1] * 100 / end, *form, elap[1]);
      std::printf("  %2d %% dual   \e[32m(%*d ms)\e[0m\n", (int) elap[2] * 100 / end, *form, elap[2]);
      std::printf("  %2d %% match  \e[32m(%*d ms)\e[0m\n", (int) elap[3] * 100 / end, *form, elap[3]);
      std::printf("  %2d %% kernel \e[32m(%*d ms)\e[0m\n", (int) elap[4] * 100 / end, *form, elap[4]);
      std::printf("  %2d %% fixes  \e[32m(%*d ms)\e[0m\n", (int) elap[5] * 100 / end, *form, elap[5]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}
} // namespace trinity