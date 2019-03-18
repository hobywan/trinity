/*
 *                          'matching.cpp'
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

#include "trinity/matching.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
Match::Match() {
  max.capa     = 0;
  max.depth    = 0;
  param.cores  = 1;
  param.found  = false;
  task.mapping = nullptr;
  task.matched = nullptr;
  task.lists   = nullptr;
  sync.visited = nullptr;
  sync.degree  = nullptr;
  sync.off     = nullptr;
}

/* --------------------------------------------------------------------------- */
void Match::initialize(size_t capacity, int* mapping, int* index) {

  max.capa    = (int) capacity;
  max.depth   = 0;
  param.cores = omp_get_max_threads();
  param.found = false;
  sync.off    = index;
  // memalloc
  task.mapping = mapping;
  task.matched = new int[max.capa];
  sync.visited = new char[max.capa];
  sync.degree  = new char[max.capa];

#pragma omp parallel
  reset();
}

/* --------------------------------------------------------------------------- */
Match::~Match() {

  delete[] task.matched;
  delete[] sync.visited;
  delete[] sync.degree;
}

/* --------------------------------------------------------------------------- */
void Match::reset() {

#pragma omp master
  task.cardin[0] = task.cardin[1] = 0;

#pragma omp for nowait
  for (int i = 0; i < max.capa; ++i)
    sync.visited[i] = 0;

#pragma omp for nowait
  for (int i = 0; i < max.capa; ++i)
    sync.degree[i] = 0;

#pragma omp for
  for (int i = 0; i < max.capa; ++i)
    task.matched[i] = -1;
}

/* --------------------------------------------------------------------------- */
int* Match::computeGreedyMatching(const Graph& graph, int nb) {

  reset();

  std::stack<int> stack;

#pragma omp for
  for (int i = 0; i < nb; ++i) {
    const int& u = graph[i][0];
    sync.degree[u] = (char) (graph[i].size() - 1);
    assert(sync.degree[u]);
  }

#pragma omp for schedule(guided)
  for (int i = 0; i < nb; ++i)
    // find maximal set of vertex-disjoint augmenting paths via DFS
    matchAndUpdate(i, graph, &stack);

  return task.matched;
}

/* --------------------------------------------------------------------------- */
void Match::matchAndUpdate(int i, const Graph& graph, std::stack<int>* stack) {

  stack->push(i);

  int j;
  int k;
  int u;

  do {
    j = stack->top();
    u = graph[j][0];
    stack->pop();

    if (!sync::compareAndSwap(sync.visited + u, 0, 1))
      continue;

    for (auto v = graph[j].begin() + 1; v < graph[j].end(); ++v) {
      if (sync::compareAndSwap(sync.visited + (*v), 0, 1)) {
        task.matched[u] = *v;
        task.matched[*v] = MATCHED;  // avoid duplicates

        k = task.mapping[*v];
        if (k < 0)
          continue;

        // update the degree of neighbors of v
        // and recursive call to match the new vertex w of degree=1
        for (auto w = graph[k].begin() + 1; w < graph[k].end(); ++w) {
          const int& nxt = task.mapping[*w];
          if (nxt > -1 and sync::fetchAndSub(sync.degree + (*w), char(1)) == 2)
            stack->push(nxt);
        }
        break;
      }
    }
  } while (!stack->empty());
}

/* --------------------------------------------------------------------------- */
int Match::getRatio(const Graph& graph, int nb, int* count) {

  int local_matched = 0;

#pragma omp single
  *count = 0;

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb; ++i) {
    const int& u = graph[i][0];
    if (-1 < task.matched[u] or task.matched[u] == MATCHED)
      local_matched++;
  }
#pragma omp critical
  *count += local_matched;
#pragma omp barrier
#ifdef DEBUG
#pragma omp master
  std::printf("ratio task.matched: %.2f %%\n", float((*count)*100)/nb);
#endif
  return *count;
}

/* --------------------------------------------------------------------------- */
int* Match::localSearchBipartite(const Graph& graph, int nb) {

  auto tic = timer::now();

  // retrieve a greedy matching using karp-sipser heuristic
  computeGreedyMatching(graph, nb);

  //
  int* look_ahead = task.lists[0]; // reuse array
#pragma omp for
  for (int i = 0; i < max.capa; ++i)
    look_ahead[i] = 1;

  bool found;
  int level = 0;

  std::stack<int> stack;

  do {
    // (!) mandatory barrier
#pragma omp barrier
#pragma omp master
    param.found = false;
    found = false;

#pragma omp for
    for (int i = 0; i < max.capa; ++i)
      sync.visited[i] = 0;

#pragma omp for schedule(guided) nowait
    for (size_t i = 0; i < graph.size(); ++i) {
      const int& u = graph[i][0];

      if (task.matched[u] < 0)
        found = lookAheadDFS(i, graph, &stack);
    }
#pragma omp critical
    param.found |= found;
#pragma omp barrier
  } while (param.found);

  tools::showElapsed(tic, "pothen-fan maximal matching done", 2);
  return task.matched;
}

/* --------------------------------------------------------------------------- */
bool Match::lookAheadDFS(int id, const Graph& graph, std::stack<int>* stack) {

  int* look_ahead = task.lists[0];

  int j;
  int u;
  int k;

  stack->push(id);

  do {
    j = stack->top();
    u = graph[j][0];
    stack->pop();

    // look ahead step
    for (auto v = graph[j].begin() + look_ahead[u]; v < graph[j].end(); ++v) {
      look_ahead[u]++;
      if (task.matched[*v] < 0) {
        if (sync::compareAndSwap(sync.visited + u, 0, 1)) {
          task.matched[u] = *v;
          task.matched[*v] = u;
          return true;
        }
      }
    }
    // scan unmatched but sync.visited neighbors if not found
    for (auto v = graph[j].begin() + 1; v < graph[j].end(); ++v) {
      if (sync::compareAndSwap(sync.visited + u, 0, 1)) {
        const int& w = task.mapping[task.matched[*v]];
        if (w > -1)
          stack->push(w);
      }
    }
  } while (!stack->empty());
  return false;
}
/* --------------------------------------------------------------------------- */
} // namespace trinity
