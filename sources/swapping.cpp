/* ------------------------------------*/
#include "swapping.h"
/* ------------------------------------ */
using namespace trinity;

/* ------------------------------------ */
Swap::Swap(Mesh* input) :
  mesh    (input),
  off     (input->off_),
  fixes   (input->fixes_),
  activ   (input->activ_),
  cores   (input->nb_cores_),
  qualit  (input->qualit_.data()),
  nb_nodes(input->nb_cores_),
  nb_elems(input->nb_elems_),
  verbose (input->_verb),
  iter    (input->_iter),
  rounds  (input->_rounds),
  depth   (20)
{
  auto size = input->max_elem;
  map = new int[size];
  tasks = new int[size];
  dual.resize(size);
  heuris.init(size, map, off);
}

/* ------------------------------------ */
Swap::~Swap() {

  delete[] map;
  delete[] tasks;
}

/* ------------------------------------ */
void Swap::run(Stats* tot) {

  init();

  int time[] = {0, 0, 0, 0, 0, 0};
  int form[] = {0, 0, 0};
  int stat[] = {0, 0, 0};

#pragma omp parallel
  {
    std::vector<int> heap;

    cacheQuality();
    timer::save(tic, time);

    int level = 0;
    do {

      filterElems(&heap);
      timer::save(tic, time + 1);

      if (!nb_tasks)
        break;

      extractDualGraph();
      timer::save(tic, time + 2);

      heuris.computeKarpSipser(dual, nb_tasks);
      heuris.getRatio(dual, nb_tasks, &count);
      timer::save(tic, time + 3);

      processFlips();
      saveStat(level, stat, form);
      timer::save(tic, time + 4);

      mesh->fixTagged();
      timer::save(tic, time + 5);

      showStat(level++, form);

    } while (nb_comms);

    recap(time, stat, form, tot);
    mesh->verifyTopology();
  }
}

/* ------------------------------------ */
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
  const double q[] = { mesh->computeQuality(f1), mesh->computeQuality(f2) };
  const double q_old = std::min(qualit[k1], qualit[k2]);
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
  sync::compareAndSwap(fixes + f1[0], 0, 1);
  sync::compareAndSwap(fixes + f1[1], 0, 1);
  sync::compareAndSwap(fixes + f1[2], 0, 1);
  sync::compareAndSwap(fixes + f2[2], 0, 1);
  // update cached qualit
  qualit[k1] = q[0];
  qualit[k2] = q[1];

  // propagate
  for (const int& t : dual[index])
    sync::compareAndSwap(activ + t, 0, 1);

  if (map[k2] > -1)
    for (const int& t : dual[map[k2]])
      sync::compareAndSwap(activ + t, 0, 1);

  return 2;
}

/* ------------------------------------ */
void Swap::cacheQuality() {

#pragma omp single
  nb_tasks = nb_comms = 0;

#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    fixes[i] = 0;
#pragma omp for nowait
  for (int i = 0; i < nb_elems; ++i)
    map[i] = -1;
#pragma omp for
  for (int i = 0; i < nb_elems; ++i)
    activ[i] = 0;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_elems; ++i) {
    const int* n = mesh->getElem(i);
    if (__builtin_expect(*n > -1, 1)) {
      qualit[i] = mesh->computeQuality(i);
      activ[i] = 1;
    }
  }

#pragma omp master
  round = timer::now();
}

/* ------------------------------------ */
void Swap::filterElems(std::vector<int>* heap) {

#pragma omp master
  nb_activ = 0;

  int count = 0;
  heap->reserve(nb_elems / cores);

//#pragma omp for schedule(dynamic,chunk) nowait
#pragma omp for nowait
  for (int i = 0; i < nb_elems; ++i) {
    if (__builtin_expect(activ[i], 0)) {
      activ[i] = 0;
      count++;
      if (qualit[i] < trinity::q_min)
        heap->push_back(i);
    }
  }
  sync::reduceTasks(tasks, heap, &nb_tasks, off);
  sync::fetchAndAdd(&nb_activ, count);
}

/* ------------------------------------ */
void Swap::extractDualGraph() {

  const auto& stenc = mesh->stenc_;
  const int* deg = mesh->deg_;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_tasks; ++i) {

    const int& k = tasks[i];
    const int* n = mesh->getElem(k);

    dual[i].clear();
    dual[i].reserve(4);
    dual[i].push_back(k);

    for (auto t = stenc[*n].begin(); dual[i].size() < 3 and t < stenc[*n].begin() + deg[*n]; ++t) {
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
    map[k] = i;
  }
}

/* ------------------------------------ */
void Swap::processFlips() {

  int succ = 0;

#pragma omp single
  nb_comms = 0;

#pragma omp for schedule(guided) nowait
  for (int index = 0; index < nb_tasks; ++index) {
    const int& k1 = dual[index][0];
    const int& k2 = heuris.matched[k1];
    if (__builtin_expect(k2 > -1, 1))
      succ += swap(k1, k2, index);
  }

  // update nb_comms
  sync::fetchAndAdd(&nb_comms, succ);
#pragma omp barrier
}

/* ------------------------------------ */
void Swap::init() {
#pragma omp master
  {
    if (verbose == 1)
      std::printf("%-18s%s", "= swapping", "...");

    else if (verbose == 2)
      std::printf("Process swapping ... ");

    std::fflush(stdout);
    start = round = tic = timer::now();
  }
}

/* ------------------------------------ */
void Swap::saveStat(int level, int* stat, int* form) {
#pragma omp master
  {
    stat[0] += nb_activ;
    stat[1] += nb_tasks;
    stat[2] += nb_comms;
    if (!level) {
      form[0] = tools::format(nb_activ);
      form[1] = tools::format(nb_tasks);
      form[2] = tools::format(nb_comms);
    }
  }
}

/* ------------------------------------ */
void Swap::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d filt. \e[0m(%2d %%)\e[0m, "
                  "%*d comm. \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  level + 1, form[0], nb_activ, (int) (nb_activ * 100 / nb_elems),
                  form[1], nb_tasks, (int) (nb_tasks * 100 / nb_activ),
                  form[2], nb_comms, (int) (nb_comms * 100 / nb_tasks), timer::round(round));
      std::fflush(stdout);
    }
  }
}

/* ------------------------------------ */
void Swap::recap(int* time, int* stat, int* form, Stats* tot) {
#pragma omp master
  {
    int end = timer::elapsed_ms(start);

    tot->eval += stat[0];
    tot->task += stat[1];
    tot->elap += end;

    // manually unrolled
    tot->step[0] += time[0] + time[1];  // qualit + filterElems
    tot->step[1] += time[2];          // dual
    tot->step[2] += time[3];          // match
    tot->step[3] += time[4];          // processFlips
    tot->step[4] += time[5];          // repair

    int span = 0;
    for (int i = 0; i < 6; ++i)
      span = std::max(span, time[i]);
    *form = tools::format(span);

    if (!verbose)
      std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100 * (++iter) / (4 * rounds + 1)));

    else if (verbose == 1) {
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m\n", (int) std::floor(stat[1] / (end * 1e-3)), (float) end / 1e3);
    } else if (verbose == 2) {
      std::printf("\n\n");
      std::printf("= rate : %d flip/sec (%d tasks) \n", (int) std::floor(stat[1] / (end * 1e-3)), stat[1]);
      std::printf("= time per step\n");
      std::printf("  %2d %% qualit \e[32m(%*d ms)\e[0m\n", (int) time[0] * 100 / end, *form, time[0]);
      std::printf("  %2d %% filterElems \e[32m(%*d ms)\e[0m\n", (int) time[1] * 100 / end, *form, time[1]);
      std::printf("  %2d %% dual   \e[32m(%*d ms)\e[0m\n", (int) time[2] * 100 / end, *form, time[2]);
      std::printf("  %2d %% match  \e[32m(%*d ms)\e[0m\n", (int) time[3] * 100 / end, *form, time[3]);
      std::printf("  %2d %% processFlips \e[32m(%*d ms)\e[0m\n", (int) time[4] * 100 / end, *form, time[4]);
      std::printf("  %2d %% fixes  \e[32m(%*d ms)\e[0m\n", (int) time[5] * 100 / end, *form, time[5]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}
