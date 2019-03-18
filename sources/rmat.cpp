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

#include "trinity/rmat.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
RMAT::RMAT() { reset(); }

/* --------------------------------------------------------------------------- */
RMAT::~RMAT() { reset(); }

/* --------------------------------------------------------------------------- */
void RMAT::reset() {
  time.start = timer::now();
  nb.nodes   = 0;
  nb.edges   = 0;
  nb.rounds  = 0;
  nb.error   = 0;
  deg.max    = std::numeric_limits<int>::max();
  deg.avg    = 0;
  stat.ratio = 0.;
  // reset graph
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
  file >> nb.nodes >> nb.edges;
  assert(nb.nodes);
  assert(nb.edges);
  graph.resize((size_t) nb.nodes);

  int k, u, v;

  // first-touch
#pragma omp parallel for schedule(static)
  for (int i = 0; i < nb.nodes; ++i) {
    graph[i].reserve(200);
    graph[i].push_back(i);
  }

  while (k < nb.edges and file >> u >> v) {
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
    for (int i = 0; i < nb.nodes; ++i) {
      graph[i].shrink_to_fit();
      int deg = (int) graph[i].size() - 1;
      deg_max_local = std::max(deg, deg_max_local);
      deg_sum_local += deg;
    }

#pragma omp critical
    {
      if (deg.max < deg_max_local) {
        deg.max = deg_max_local;
      }
      deg.avg += deg_sum_local;
    }
  }  // implicit barrier here

  deg.avg /= nb.nodes;

  std::printf("|V|=%d, |E|=%d, deg.max=%d, deg.avg=%d \e[32m(%d ms)\e[0m\n",
              nb.nodes, nb.edges, deg.max, deg.avg, elapsed());

}

/* --------------------------------------------------------------------------- */
void RMAT::info(const std::string testcase) {

  std::printf("|V|=%d, |E|=%d, rounds=%d, errors=%d, colors=%d, "
              "deg.max=%d, deg.avg=%d, ratio=%.3f, %s "
              "[%2d threads] \e[32m(%d ms)\e[0m\n",
              nb.nodes,
              nb.edges,
              nb.rounds,
              nb.error,
              nb.color,
              deg.max,
              deg.avg,
              stat.ratio,
              testcase.data(),
              omp_get_num_threads(),
              elapsed());
}

/* --------------------------------------------------------------------------- */
void RMAT::saveChrono() {
  time.start = timer::now();
}

/* --------------------------------------------------------------------------- */
int RMAT::elapsed() {
  return timer::elapsed_ms(time.start);
}
/* --------------------------------------------------------------------------- */
} // namespace trinity
