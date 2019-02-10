/*
 *                          'smoothing.cpp'
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

#include "smoothing.h"
/* ------------------------------------ */
namespace trinity {
/* ------------------------------------ */
Smooth::Smooth(Mesh* input, Partit* algo, int level)
  : mesh(input),
    cores(input->nb_cores_),
    activ(input->activ_),
    qualit(input->qualit_.data()),
    nb_nodes(input->nb_nodes_),
    nb_elems(input->nb_elems_),
    verbose(input->_verb),
    iter(input->_iter),
    rounds(input->_rounds),
    heuris(algo),
    depth(level)
{}

/* ------------------------------------ */
Smooth::~Smooth() {}

/* ------------------------------------ */
void Smooth::run(Stats* tot) {

  init();

  int time[] = {0, 0, 0, 0};
  int form[] = {0, 0};
  int stat[] = {0, 0};

#pragma omp parallel
  {
    preProcess();
    timer::save(tic, time);

    heuris->extractColoring(mesh);
    timer::save(tic, time + 1);

    cacheQuality();
    timer::save(tic, time + 2);

    for (int level = 0; level < depth; ++level) {
      movePoints();
      saveStat(level, stat, form);
      showStat(level, form);
    }
    timer::save(tic, time + 3);
    // finalize
    recap(time, stat, form, tot);
  }
}

/* ------------------------------------ */
int Smooth::moveLaplacian(int i) {

  const auto& vicin = mesh->vicin_[i];
  const auto& stenc = mesh->stenc_[i];
  const int deg = mesh->deg_[i];
  const int nb = vicin.size();

  double q_min = std::numeric_limits<double>::max();

  for (auto t = stenc.begin(); t < stenc.begin() + deg; ++t)
    q_min = std::min(q_min, qualit[*t]);

  double* M = new double[nb * 3];
  double* q = new double[stenc.size()];
  double* pa = mesh->points_.data() + (i * 2);
  double* ma = mesh->tensor_.data() + (i * 3);

  double len = 0.;
  double p_opt[] = {0, 0};

  // 1) compute average optimal position of v[i]
  for (int k = 0; k < nb; ++k) {

    const int& j = vicin[k];
    const double* pb = mesh->points_.data() + (j * 2);
    const double* mb = mesh->tensor_.data() + (j * 3);

    // a. reduction on each local Pimal position of v[i]
    len = mesh->computeLength(i, j);
    assert(len);
    p_opt[0] += pa[0] + ((pb[0] - pa[0]) / len);  // unrolled for perfs
    p_opt[1] += pa[1] + ((pb[1] - pa[1]) / len);

    // b. storeFile tensor locally for interpolation
    std::memcpy(M + (k * 3), mb, sizeof(double) * 3);
    assert(std::isfinite(M[k * 3]));
  }
  p_opt[0] /= nb;
  p_opt[1] /= nb;

  // 2) interpolate tensor according to stencil
  const double p_ini[] = {pa[0], pa[1]};
  const double m_ini[] = {ma[0], ma[1], ma[2]};
  numeric::interpolateTensor(M, ma, nb);

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
    pa[0] = (1 - weight) * p_ini[0] + weight * p_opt[0];
    pa[1] = (1 - weight) * p_ini[1] + weight * p_opt[1];

    int j = 0;
    // stencil convexity and qualit improvement test
    for (auto t = stenc.begin(); valid and t < stenc.begin() + deg; ++t, ++j) {
      const int* n = mesh->getElem(*t);
      valid = false;
      if (mesh->isCounterclockwise(n)) {
        q[j] = mesh->computeQuality(*t);
        valid = (q[j] > q_min);
      }
    }
    weight *= 0.5;
  } while (not valid and ++retry < 5);

  if (valid) {
    int j = 0;
    for (auto t = stenc.begin(); t < stenc.begin() + deg; ++t, ++j)
      qualit[*t] = q[j];
  } else {
    std::memcpy(pa, p_ini, sizeof(double) * 2);
    std::memcpy(ma, m_ini, sizeof(double) * 3);
  }

  delete[] M;
  delete[] q;
  return (int) valid;
}

/* ------------------------------------ */
void Smooth::preProcess() {

  int count = 0;

#pragma omp single
  nb_tasks = nb_comms = 0;

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    if (not mesh->stenc_[i].empty()) {
      activ[i] = static_cast<char>(__builtin_expect(mesh->isBoundary(i), 0) ? 0 : 1);
    }

  mesh->extractPrimalGraph();

}

/* ------------------------------------ */
void Smooth::cacheQuality() {

#pragma omp for
  for (int i = 0; i < mesh->nb_elems_; ++i)
    if (mesh->isActiveElem(i)) {
      qualit[i] = mesh->computeQuality(i);
    }

#pragma omp master
  round = tic;
}

/* ------------------------------------ */
void Smooth::movePoints() {

  int succ = 0;
  int tota = 0;

#pragma omp single
  nb_tasks = nb_comms = 0;

  for (int i = 0; i < heuris->parts; ++i) {
#pragma omp for schedule(guided)
    for (int j = 0; j < heuris->card[i]; ++j) {
      const int& k = heuris->subset[i][j];
      if (activ[k]) {
        succ += moveLaplacian(k);
        tota++;
      }
    }
  }
  sync::fetchAndAdd(&nb_comms, succ);
  sync::fetchAndAdd(&nb_tasks, tota);
#pragma omp barrier
}

/* ------------------------------------ */
void Smooth::init() {
#pragma omp master
  {
    if (verbose == 1) {
      std::printf("%-18s%s", "= smoothing", "...");
    } else if (verbose == 2) {
      std::printf("Process smoothing ... ");
    }

    std::fflush(stdout);
    start = round = tic = timer::now();
  }
}

/* ------------------------------------ */
void Smooth::saveStat(int level, int* stat, int* form) {
#pragma omp single
  {
    stat[0] += nb_tasks;
    stat[1] += nb_comms;
    if (!level) {
      form[0] = tools::format(nb_tasks);
      form[1] = tools::format(nb_comms);
    }
  }
}

/* ------------------------------------ */
void Smooth::showStat(int level, int* form) {
#pragma omp single
  {
    if (verbose == 2) {
      std::printf("\n= round %2d. %*d tasks \e[0m(100 %%)\e[0m, "
                  "%*d comm. \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
                  level + 1, form[0], nb_tasks, form[1], nb_comms,
                  (int) (nb_comms * 100 / nb_tasks), timer::round(round));
      std::fflush(stdout);
    }
  }
}

/* ------------------------------------ */
void Smooth::recap(int* time, int* stat, int* form, Stats* tot) {
#pragma omp master
  {
    int end = timer::elapsed_ms(start);

    tot->eval += stat[0];
    tot->task += stat[1];
    tot->elap += end;
    for (int i = 0; i < 4; ++i)
      tot->step[i] += time[i];

    *form = tools::format(time[3]);
    if (not verbose) {
      std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100 * (++iter) / (4 * rounds + 1)));
    } else if (verbose == 1) {
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m\n", (int) std::floor(stat[0] / (end * 1e-3)), (float) end / 1e3);
    } else if (verbose == 2) {
      std::printf("\n\n");
      std::printf("= rate : %d move/sec (%d tasks) \n", (int) std::floor(stat[0] / (end * 1e-3)), stat[0]);
      std::printf("= time per step\n");
      std::printf("  %2d %% primal \e[32m(%*d ms)\e[0m\n", (int) time[0] * 100 / end, *form, time[0]);
      std::printf("  %2d %% color  \e[32m(%*d ms)\e[0m\n", (int) time[1] * 100 / end, *form, time[1]);
      std::printf("  %2d %% qualit \e[32m(%*d ms)\e[0m\n", (int) time[2] * 100 / end, *form, time[2]);
      std::printf("  %2d %% processFlips \e[32m(%*d ms)\e[0m\n", (int) time[3] * 100 / end, *form, time[3]);
      std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}
} // namespace trinity