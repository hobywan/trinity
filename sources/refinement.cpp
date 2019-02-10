/*
 *                          'refinement.cpp'
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

#include "refinement.h"
/* ------------------------------------*/
namespace trinity {
/* ------------------------------------*/
Refine::Refine(Mesh* input, int level)
  :
  mesh(input),
  off(input->sync.off),
  activ(input->sync.activ),
  cores(input->nb.cores),
  steiner(input->capa.node, input->capa.bucket, 2),
  nb_nodes(input->nb.nodes),
  nb_elems(input->nb.elems),
  verbose(input->param.verb),
  iter(input->param.iter),
  rounds(input->param.rounds),
  depth(level)
{
  const auto size = input->capa.elem;
  index = new int[cores + 1];
  edges = new int[size * 3];   // 3 data/edge, 3 edges/elem
  elems = new int[size];
  pattern = new char[size];
}

/* ------------------------------------*/
Refine::~Refine() {

  delete[] index;
  delete[] edges;
  delete[] elems;
  delete[] pattern;
}

/* ------------------------------------ */
void Refine::run(Stats* tot) {

  init();

  int time[] = {0, 0, 0, 0, 0};
  int stat[] = {0, 0};
  int form[] = {0, 0, 0};

#pragma omp parallel
  {
    std::vector<int> heap[2];

    int level = 0;
    int tid = omp_get_thread_num();

    preProcess(heap);
    timer::save(tic, time);

    do {
      filterElems(heap);
      timer::save(tic, time + 1);

      if (!nb_tasks) {
        break;
      }

      computeSteinerPoints();
      saveStat(level, stat, form);
      timer::save(tic, time + 2);

      processElems(tid);
      timer::save(tic, time + 3);

#ifdef DEFERRED_UPDATES
      mesh->commitUpdates();
      timer::save(tic, time+4);
#endif
      showStat(level, form);

    } while (++level < depth);

#ifndef DEFERRED_UPDATES
    mesh->fixAll();
    timer::save(tic, time + 4);
#endif

    recap(time, stat, form, tot);
    mesh->verifyTopology();
  }
}

/* ------------------------------------*/
void Refine::cutElem(int id, int* offset) {

#ifdef DEFERRED_UPDATES
  int tid = omp_get_thread_num();
#endif

  const int* n = mesh->getElem(id);
  // the i'th edge is opposite the i'th node in the element.
  const int s[] = {
    steiner.getValue(n[1], n[2]),
    steiner.getValue(n[2], n[0]),
    steiner.getValue(n[0], n[1])
  };

  if (pattern[id] == 1) {

    const int t[] = {id, *offset};

    for (int i = 0; i < 3; ++i) {
      if (s[i] > -1) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        assert(s[i] not_eq n[i]);
        //
        const int elem0[] = {n[i], n[j], s[i]};
        const int elem1[] = {n[i], s[i], n[k]};
        assert(mesh->isCounterclockwise(elem0));
        assert(mesh->isCounterclockwise(elem1));
        //
#ifdef DEFERRED_UPDATES
        mesh->deferredAppend(tid, s[i], {t[0], t[1]});
        mesh->deferredAppend(tid, n[i], t[1]);
        mesh->deferredRemove(tid, n[k], t[0]); // replacement
        mesh->deferredAppend(tid, n[k], t[1]);
#else
        mesh->updateStencil(s[i], {t[0], t[1]});
        mesh->updateStencil(n[i], t[1]);
        mesh->updateStencil(n[k], t[1]);
#endif
        //
        mesh->replaceElem(t[0], elem0);
        mesh->replaceElem(t[1], elem1);
        // manually unrolled
        activ[t[0]] = activ[t[1]] = 1;
        break;
      }
    }
  } else if (pattern[id] == 2) {

    const int t[] = {id, *offset, *offset + 1};

    for (int i = 0; i < 3; ++i) {
      if (s[i] < 0) {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;

        double diag[2];
        diag[0] = mesh->computeLength(n[j], s[j]);
        diag[1] = mesh->computeLength(n[k], s[k]);
        // selected diagonal: v[s],n[s+1]
        int v = (diag[0] < diag[1] ? j : k);
        int w = (v == j ? k : j);

        const int elem0[] = {n[i], s[k], s[j]};
        const int elem1[] = {n[j], n[k], s[v]};
        const int elem2[] = {s[j], s[k], n[v]};
        assert(mesh->isCounterclockwise(elem0));
        assert(mesh->isCounterclockwise(elem1));
        assert(mesh->isCounterclockwise(elem2));
        //
#ifdef DEFERRED_UPDATES
        mesh->deferredAppend(tid, s[v], {t[0], t[1], t[2]});
        mesh->deferredAppend(tid, s[w], {t[0], t[2]});
        mesh->deferredRemove(tid, n[v], t[0]);
        mesh->deferredAppend(tid, n[v], {t[1], t[2]});
        mesh->deferredRemove(tid, n[w], t[0]);
        mesh->deferredAppend(tid, n[w], t[1]);
#else
        mesh->updateStencil(s[v], {t[0], t[1], t[2]});
        mesh->updateStencil(s[w], {t[0], t[2]});
        mesh->updateStencil(n[v], {t[1], t[2]});
        mesh->updateStencil(n[w], t[1]);
#endif
        //
        mesh->replaceElem(t[0], elem0);
        mesh->replaceElem(t[1], elem1);
        mesh->replaceElem(t[2], elem2);
        // manually unrolled
        activ[t[0]] = activ[t[1]] = activ[t[2]] = 1;
        break;
      }
    }
  } else {

    const int t[] = {id, *offset, *offset + 1, *offset + 2};
    //
    const int elem0[] = {n[0], s[2], s[1]};
    const int elem1[] = {n[1], s[0], s[2]};
    const int elem2[] = {n[2], s[1], s[0]};
    const int elem3[] = {s[0], s[1], s[2]};
    assert(mesh->isCounterclockwise(elem0));
    assert(mesh->isCounterclockwise(elem1));
    assert(mesh->isCounterclockwise(elem2));
    assert(mesh->isCounterclockwise(elem3));
    //
#ifdef DEFERRED_UPDATES
    mesh->deferredAppend(tid, s[0], {t[1], t[2], t[3]});
    mesh->deferredAppend(tid, s[1], {t[0], t[2], t[3]});
    mesh->deferredAppend(tid, s[2], {t[0], t[1], t[3]});
    mesh->deferredRemove(tid, n[1], t[0]);
    mesh->deferredAppend(tid, n[1], t[1]);
    mesh->deferredRemove(tid, n[2], t[0]);
    mesh->deferredAppend(tid, n[2], t[2]);
#else
    mesh->updateStencil(s[0], {t[1], t[2], t[3]});
    mesh->updateStencil(s[1], {t[0], t[2], t[3]});
    mesh->updateStencil(s[2], {t[0], t[1], t[3]});
    mesh->updateStencil(n[1], t[1]);
    mesh->updateStencil(n[2], t[2]);
#endif
    //
    mesh->replaceElem(t[0], elem0);
    mesh->replaceElem(t[1], elem1);
    mesh->replaceElem(t[2], elem2);
    mesh->replaceElem(t[3], elem3);
    // manually unrolled
    activ[t[0]] = activ[t[1]] = activ[t[2]] = activ[t[3]] = 1;
  }
  *offset += pattern[id];
}

/* ------------------------------------ */
void Refine::preProcess(std::vector<int> heap[2]) {

  heap[0].reserve(nb_elems);
  heap[1].reserve(nb_elems / cores);

#pragma omp master
  {
    old_node = nb_nodes;
    old_elem = nb_elems;
    shift = 0;
    nb_adds = 0;
    nb_split = 0;   // number of elems to be appended
    nb_eval = 0;
    nb_tasks = 0;
    nb_steiner = 0;
  }

#pragma omp for
  for (int i = 0; i < nb_elems; ++i) {
    activ[i] = static_cast<char>(mesh->isActiveElem(i) ? 1 : 0);
  }

  steiner.reset();

#pragma omp single
  round = timer::now();
}

/* ------------------------------------ */
void Refine::filterElems(std::vector<int> heap[2]) {
#pragma omp single
  {
    nb_eval = nb_adds = nb_split = nb_tasks = nb_steiner = 0;
    std::memset(index, 0, (cores + 1) * sizeof(int));
  }

  // step 1: filter elems
  int count[] = {0, 0};
  double len[] = {0, 0, 0};

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_elems; ++i) {
    pattern[i] = 0;
    // i > offset => new elems
    if (__builtin_expect(activ[i], i > shift)) {
      auto n = mesh->getElem(i);
      len[0] = mesh->computeLength(n[0], n[1]);
      len[1] = mesh->computeLength(n[1], n[2]);
      len[2] = mesh->computeLength(n[2], n[0]);

      for (int j = 0; j < 3; ++j) {
        const int k = (j + 1) % 3;
        // (!) avoid duplicated steiner-points but don't forget boundary edges
        if (len[j] > trinity::l_max) {
          ++pattern[i];
          // (!) be aware of diagonals: dont rely only on node activs
          const int opp = mesh->getElemNeigh(i, n[j], n[k]);
          if (n[j] < n[k] or opp < 0) {
            heap[0].push_back(n[j]);
            heap[0].push_back(n[k]);
            heap[0].push_back(opp);
          }
        }
      }
      if (pattern[i]) {
        count[1] += pattern[i];
        heap[1].push_back(i);
      }
      count[0]++;
    }
  }
  sync::fetchAndAdd(&nb_eval, count[0]);
  sync::fetchAndAdd(&nb_adds, count[1]);

  // reductions on steiner point & bad elems list
  sync::reduceTasks(edges, heap, &nb_split, 3);
  sync::reduceTasks(elems, heap + 1, &nb_tasks, off);
}

/* ------------------------------------ */
void Refine::computeSteinerPoints() {

  int count = 0;

  // step 2: compute steiner point and fix new cells ID
#pragma omp for nowait
  for (int i = 0; i < nb_split; ++i) {
    // implicit index
    const int id = nb_nodes + i;
    const int& v1 = edges[i * 3];
    const int& v2 = edges[i * 3 + 1];
    const int& opp = edges[i * 3 + 2];
    // 1) calculate and insert the steiner point
    double* P = mesh->geom.points.data() + (id * 2);
    double* M = mesh->geom.tensor.data() + (id * 3);
    mesh->computeSteinerPoint(v1, v2, P, M);

    count++;
    // (!) map (i,j)->id
    int v = std::min(v1, v2);
    int w = std::max(v1, v2);
    steiner.push(v, {w, id});

    if (opp < 0) {
      mesh->sync.tags[id] = mask::bound;
    }
  }

  sync::fetchAndAdd(&nb_steiner, count);
#pragma omp barrier
}

/* ------------------------------------ */
void Refine::processElems(int tid) {
#pragma omp master
  {
    index[0] = shift = nb_elems;
    nb_nodes += nb_split;
    nb_elems += nb_adds;
  }

#pragma omp for schedule(static)
  for (int i = 0; i < nb_tasks; ++i)
    index[tid + 1] += pattern[elems[i]];  // no data race

  sync::prefixSum(index, cores, 16);

  // step 3: refine bad elems
#pragma omp for schedule(static)
  for (int i = 0; i < nb_tasks; ++i)
    cutElem(elems[i], index + tid);
}

/* ------------------------------------ */
void Refine::init() {
#pragma omp master
  {
    if (verbose == 1) {
      std::printf("%-18s%s", "= refinement", "...");
    } else if (verbose == 2) {
      std::printf("Process refinement ... ");
    }

    std::fflush(stdout);
    start = round = tic = timer::now();
  }
}

/* ------------------------------------ */
void Refine::saveStat(int level, int* stat, int* form) {
#pragma omp master
  {
    cur_elem = nb_elems;
    stat[0] += nb_eval;
    stat[1] += nb_tasks;
    if (!level) {
      form[0] = tools::format(nb_eval);
      form[1] = tools::format(nb_tasks);
      form[2] = tools::format(nb_steiner);
    }
  }
}

/* ------------------------------------ */
void Refine::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d filt. \e[0m(%2d %%)\e[0m, "
                  "%*d stein \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  level + 1, form[0], nb_eval, (int) (nb_eval * 100 / cur_elem),
                  form[1], nb_tasks, (int) (nb_tasks * 100 / nb_eval),
                  form[2], nb_steiner, (int) (nb_steiner * 100 / nb_nodes), timer::round(round));
      std::fflush(stdout);
    }
  }
}

/* ------------------------------------ */
void Refine::recap(int* time, int* stat, int* form, Stats* tot) {
#pragma omp master
  {
    int end = std::max(timer::elapsed_ms(start), 1);
    int span = 0;

    tot->eval += stat[0];
    tot->task += stat[1];
    tot->elap += end;

    tot->step[0] += time[0] + time[1];  // filterElems
    tot->step[1] += time[2];          // steiner
    tot->step[2] += time[3];          // processFlips
    tot->step[3] += time[4];          // repair

    for (int i = 0; i < 5; ++i)
      span = std::max(span, time[i]);
    *form = tools::format(span);

    if (not verbose) {
      std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100 * (++iter) / (4 * rounds + 1)));
    } else if (verbose == 1) {
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m [+%.1f %%]\n",
                  (int) std::floor(stat[1] / (end * 1e-3)), (float) end / 1e3,
                  (float) (nb_elems - old_elem) * 100 / nb_elems);
    } else if (verbose == 2) {
      std::printf("\n\n");
      std::printf("= nodes: %d old, %d new \e[0m(+%.1f %%)\e[0m\n", old_node, nb_nodes,
                  (float) (nb_nodes - old_node) * 100 / nb_nodes);
      std::printf("= elems: %d old, %d new \e[0m(+%.1f %%)\e[0m\n", old_elem, nb_elems,
                  (float) (nb_elems - old_elem) * 100 / nb_elems);
      std::printf("= rate : %d split/sec (%d tasks) \n", (int) std::floor(stat[1] / (end * 1e-3)), stat[1]);
      std::printf("= time per step\n");
      std::printf("  %2d %% preproc \e[32m(%*d ms)\e[0m\n", (int) time[0] * 100 / end, *form, time[0]);
      std::printf("  %2d %% filter  \e[32m(%*d ms)\e[0m\n", (int) time[1] * 100 / end, *form, time[1]);
      std::printf("  %2d %% steiner \e[32m(%*d ms)\e[0m\n", (int) time[2] * 100 / end, *form, time[2]);
      std::printf("  %2d %% kernel  \e[32m(%*d ms)\e[0m\n", (int) time[3] * 100 / end, *form, time[3]);
      std::printf("  %2d %% fixes   \e[32m(%*d ms)\e[0m\n", (int) time[4] * 100 / end, *form, time[4]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}
} // namespace trinity
