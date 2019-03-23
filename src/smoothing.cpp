/*
 *                          'smoothing.cpp'
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

#include "trinity/smoothing.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
Smooth::Smooth(Mesh* input, Partit* algo, int level)
  : mesh    (input),
    heuris  (algo),
    cores   (mesh->nb.cores),
    nb_nodes(mesh->nb.nodes),
    nb_elems(mesh->nb.elems),
    verbose (mesh->param.verb),
    iter    (mesh->param.iter),
    rounds  (mesh->param.rounds)
{
  sync.activ  = mesh->sync.activ;
  geom.qualit = mesh->geom.qualit.data();
  task.depth  = level;
}

/* -------------------------------------------------------------------------- */
int Smooth::moveSmartLaplacian(int i) {

  const auto& vicin = mesh->topo.vicin[i];
  const auto& stenc = mesh->topo.stenc[i];
  const int deg = mesh->sync.deg[i];
  const int nb_neigh = (int) vicin.size();
  const int nb_cells = (int) stenc.size();

  double q_min = std::numeric_limits<double>::max();

  for (auto t = stenc.begin(); t < stenc.begin() + deg; ++t)
    q_min = std::min(q_min, geom.qualit[*t]);

  double metric[nb_neigh * 3];
  double qualit[nb_cells];
  auto point_a  = mesh->geom.points.data() + (i * 2);
  auto tensor_a = mesh->geom.tensor.data() + (i * 3);

  double length;
  double point_opt[] = {0, 0};

  // 1) compute average optimal position of v[i]
  for (int k = 0; k < nb_neigh; ++k) {

    const int& neigh = vicin[k];
    const auto point_b  = mesh->geom.points.data() + (neigh * 2);
    const auto tensor_b = mesh->geom.tensor.data() + (neigh * 3);

    // a. reduction on each local Pimal position of v[i]
    length = mesh->computeLength(i, neigh);
    assert(length);
    // unrolled for perfs
    point_opt[0] += point_a[0] + ((point_b[0] - point_a[0]) / length);
    point_opt[1] += point_a[1] + ((point_b[1] - point_a[1]) / length);

    // b. storeFile tensor locally for interpolation
    std::memcpy(metric + (k * 3), tensor_b, sizeof(double) * 3);
    assert(std::isfinite(metric[k * 3]));
  }
  point_opt[0] /= nb_neigh;
  point_opt[1] /= nb_neigh;

  // 2) interpolate tensor according to stencil
  const double point_ini[]  = {point_a[0], point_a[1]};
  const double tensor_ini[] = {tensor_a[0], tensor_a[1], tensor_a[2]};
  numeric::interpolateTensor(metric, tensor_a, nb_neigh);

#ifdef DEBUG
  std::printf("interpolated[%d]: (%.2f,%.2f,%.2f)\n", i, ma[0], ma[1], ma[2]);
#endif

  // 3) adjust coef iteratively, cf. dobrynzski
  bool valid;
  int retry = 0;
  double weight = 1;

  do {
    // reset
    valid = true;
    // update weighted coords
    point_a[0] = (1 - weight) * point_ini[0] + weight * point_opt[0];
    point_a[1] = (1 - weight) * point_ini[1] + weight * point_opt[1];

    int j = 0;
    // stencil convexity and qualit improvement test
    for (auto t = stenc.begin(); valid and t < stenc.begin() + deg; ++t, ++j) {
      const int* n = mesh->getElem(*t);
      valid = false;
      if (mesh->isCounterclockwise(n)) {
        qualit[j] = mesh->computeQuality(*t);
        valid = (qualit[j] > q_min);
      }
    }
    weight *= 0.5;
  } while (not valid and ++retry < 5);

  if (valid) {
    int j = 0;
    for (auto t = stenc.begin(); t < stenc.begin() + deg; ++t, ++j)
      geom.qualit[*t] = qualit[j];
  } else {
    std::memcpy( point_a,  point_ini, sizeof(double) * 2);
    std::memcpy(tensor_a, tensor_ini, sizeof(double) * 3);
  }

  return (int) valid;
}

/* -------------------------------------------------------------------------- */
void Smooth::preProcess() {

#pragma omp single
  nb.tasks = nb.commit = 0;

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    if (not mesh->topo.stenc[i].empty()) {
      sync.activ[i] = char(__builtin_expect(mesh->isBoundary(i), 0) ? 0 : 1);
    }

  mesh->extractPrimalGraph();

}

/* -------------------------------------------------------------------------- */
void Smooth::cacheQuality() {

#pragma omp for
  for (int i = 0; i < mesh->nb.elems; ++i)
    if (mesh->isActiveElem(i)) {
      geom.qualit[i] = mesh->computeQuality(i);
    }

#pragma omp master
  time.iter = time.tic;
}

/* -------------------------------------------------------------------------- */
void Smooth::movePoints() {

  int success = 0;
  int total = 0;

#pragma omp single
  nb.tasks = nb.commit = 0;

  for (int i = 0; i < heuris->nb.parts; ++i) {
#pragma omp for schedule(guided)
    for (int j = 0; j < heuris->task.cardin[i]; ++j) {
      const int& k = heuris->task.subset[i][j];
      if (sync.activ[k]) {
        success += moveSmartLaplacian(k);
        total++;
      }
    }
  }
  sync::fetchAndAdd(&nb.commit, success);
  sync::fetchAndAdd(&nb.tasks, total);
#pragma omp barrier
}

/* -------------------------------------------------------------------------- */
void Smooth::initialize() {
#pragma omp master
  {
    if (verbose == 1) {
      std::printf("%-18s%s", "= smoothing", "...");
    } else if (verbose == 2) {
      std::printf("Process smoothing ... ");
    }

    std::fflush(stdout);
    time.start = time.iter = time.tic = timer::now();
  }
}

/* -------------------------------------------------------------------------- */
void Smooth::saveStat(int level, int* stat, int* form) {
#pragma omp single
  {
    stat[0] += nb.tasks;
    stat[1] += nb.commit;
    if (!level) {
      form[0] = tools::format(nb.tasks);
      form[1] = tools::format(nb.commit);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Smooth::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      int const round = level + 1;
      int const percent = nb.commit * 100 / nb.tasks;
      int const secs = timer::round(time.iter);

      std::printf("\n= round %2d. %*d tasks \e[0m(100 %%)\e[0m, "
                  "%*d comm. \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  round,
                  form[0], nb.tasks,
                  form[1], nb.commit, percent, secs);
      std::fflush(stdout);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Smooth::recap(int* elap, int* stat, int* form, Stats* total) {
#pragma omp master
  {
    int end = std::max(timer::elapsed_ms(time.start), 1);

    if (total != nullptr) {
      total->eval += stat[0];
      total->task += stat[1];
      total->elap += end;
      for (int i = 0; i < 4; ++i)
        total->step[i] += elap[i];
    }

    *form = tools::format(elap[3]);
    if (not verbose) {
      auto percent = (int) std::floor(100 * (++iter) / (4 * rounds + 1));
      std::printf("\r= Remeshing  ... %3d %% =", percent);

    } else if (verbose == 1) {
      auto rate = (int) std::floor(stat[0] / (end * 1e-3));
      auto secs = (float) end / 1E3;
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m\n", rate, secs);

    } else if (verbose == 2) {
      int rate = (int) std::floor(stat[0] / (end * 1e-3));
      std::printf("\n\n");
      std::printf("= rate : %d move/sec (%d tasks) \n", rate, stat[0]);

      int step[4];
      for (int i = 0; i < 4; ++i)
        step[i] = elap[i] * 100 / end;

      std::printf("= time per step\n");
      std::printf("  %2d %% primal \e[32m(%*d ms)\e[0m\n", step[0], *form, elap[0]);
      std::printf("  %2d %% color  \e[32m(%*d ms)\e[0m\n", step[1], *form, elap[1]);
      std::printf("  %2d %% qualit \e[32m(%*d ms)\e[0m\n", step[2], *form, elap[2]);
      std::printf("  %2d %% kernel \e[32m(%*d ms)\e[0m\n", step[3], *form, elap[3]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}

/* -------------------------------------------------------------------------- */
void Smooth::run(Stats* total) {

  initialize();

  int elap[] = {0, 0, 0, 0};
  int form[] = {0, 0};
  int stat[] = {0, 0};

#pragma omp parallel
  {
    preProcess();
    timer::save(time.tic, elap);

    heuris->extractColoring(mesh);
    timer::save(time.tic, elap + 1);

    cacheQuality();
    timer::save(time.tic, elap + 2);

    for (int level = 0; level < task.depth; ++level) {
      movePoints();
      saveStat(level, stat, form);
      showStat(level, form);
    }
    timer::save(time.tic, elap + 3);
    // finalize
    recap(elap, stat, form, total);
  }
}
/* -------------------------------------------------------------------------- */
} // namespace trinity
