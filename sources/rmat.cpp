/*
 *                          'rmat.h'
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

#include "rmat.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
RMAT::RMAT() { reset(); }
/* --------------------------------------------------------------------------- */
RMAT::~RMAT() { reset(); }
/* --------------------------------------------------------------------------- */
void RMAT::reset() {
  start = timer::now();
  nb_nodes = 0;
  nb_edges = 0;
  nb_rounds = 0;
  nb_error = 0;
  deg_max = std::numeric_limits<int>::max();
  deg_avg = 0;
  ratio = 0.;
  graph.clear();
}

/* --------------------------------------------------------------------------- */
void RMAT::load(const std::string path) {

  //
  reset();

  std::ifstream file(path, std::ios::in);
  assert(file.is_open());
  assert(file.good());

  std::printf("Loading '\e[36m%s\e[0m'...", path.data());
  saveChrono();

  // init counters
  file >> nb_nodes >> nb_edges;
  assert(nb_nodes);
  assert(nb_edges);
  graph.resize(nb_nodes);

  int k, u, v;

  // first-touch
#pragma omp parallel for schedule(static)
  for (int i = 0; i < nb_nodes; ++i) {
    graph[i].reserve(200);
    graph[i].push_back(i);
  }

  while (k < nb_edges and file >> u >> v) {
    graph[u].push_back(v);
    graph[v].push_back(u);
    ++k;
  }

  file.close();

  // degree stats
#pragma omp parallel
  {
    int deg_max_local = std::numeric_limits<int>::max();
    int deg_sum_local = 0;
#pragma omp for schedule(static)
    for (int i = 0; i < nb_nodes; ++i) {
      graph[i].shrink_to_fit();
      int deg = (int) graph[i].size() - 1;
      deg_max_local = std::max(deg, deg_max_local);
      deg_sum_local += deg;
    }

#pragma omp critical
    {
      if (deg_max < deg_max_local) {
        deg_max = deg_max_local;
      }
      deg_avg += deg_sum_local;
    }
  }  // implicit barrier here

  deg_avg /= nb_nodes;

  std::printf("|V|=%d, |E|=%d, deg_max=%d, deg_avg=%d \e[32m(%d ms)\e[0m\n",
              nb_nodes, nb_edges, deg_max, deg_avg, elapsed());

}

/* --------------------------------------------------------------------------- */
void RMAT::info(const std::string testcase) {

  std::printf("|V|=%d, |E|=%d, rounds=%d, errors=%d, colors=%d, "
              "deg_max=%d, deg_avg=%d, ratio=%.3f, %s "
              "[%2d threads] \e[32m(%d ms)\e[0m\n",
              nb_nodes,
              nb_edges,
              nb_rounds,
              nb_error,
              nb_color,
              deg_max,
              deg_avg,
              ratio,
              testcase.data(),
              omp_get_num_threads(),
              elapsed());
}

/* --------------------------------------------------------------------------- */
void RMAT::saveChrono() {
  start = timer::now();
}

/* --------------------------------------------------------------------------- */
int RMAT::elapsed() {
  return timer::elapsed_ms(start);
}
} // namespace trinity
