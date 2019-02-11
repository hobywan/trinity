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

#include "matching.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
Match::Match()
  :
  size(0),
  depth(0),
  cores(1),
  path(false),
  mapping(nullptr),
  matched(nullptr),
  visited(nullptr),
  degree(nullptr),
  off(nullptr),
  tasks(nullptr) {}

/* --------------------------------------------------------------------------- */
void Match::init(size_t capa, int* map, int* idx) {
  size = capa;
  depth = 0;
  cores = omp_get_max_threads();
  path = false;
  off = idx;
  // memalloc
  mapping = map;
  matched = new int[size];
  visited = new char[size];
  degree = new char[size];
#pragma omp parallel
  flush();
}

/* --------------------------------------------------------------------------- */
Match::~Match() {

  delete[] matched;
  delete[] visited;
  delete[] degree;
}

/* --------------------------------------------------------------------------- */
void Match::flush() {

#pragma omp master
  card[0] = card[1] = 0;

#pragma omp for nowait
  for (int i = 0; i < size; ++i)
    visited[i] = 0;

#pragma omp for nowait
  for (int i = 0; i < size; ++i)
    degree[i] = 0;

#pragma omp for
  for (int i = 0; i < size; ++i)
    matched[i] = -1;
}

/* --------------------------------------------------------------------------- */
int* Match::computeKarpSipser(const Graph& graph, int nb) {

  flush();

  std::stack<int> stack;

#pragma omp for
  for (int i = 0; i < nb; ++i) {
    const int& u = graph[i][0];
    degree[u] = static_cast<char>(graph[i].size() - 1);
    assert(degree[u]);
  }

#pragma omp for schedule(guided)
  for (int i = 0; i < nb; ++i)
    // find maximal set of vertex-disjoint augmenting paths via DFS
    matchAndUpdate(i, graph, &stack);

  return matched;
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

    if (!sync::compareAndSwap(visited + u, 0, 1))
      continue;

    for (auto v = graph[j].begin() + 1; v < graph[j].end(); ++v) {
      if (sync::compareAndSwap(visited + (*v), 0, 1)) {
        matched[u] = *v;
        matched[*v] = MATCHED;  // avoid duplicates

        k = mapping[*v];
        if (k < 0)
          continue;

        // update the degree of neighbors of v
        // and recursive call to match the new vertex w of degree=1
        for (auto w = graph[k].begin() + 1; w < graph[k].end(); ++w) {
          const int& nxt = mapping[*w];
          if (nxt > -1 and sync::fetchAndSub(degree + (*w), char(1)) == 2)
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
    if (-1 < matched[u] or matched[u] == MATCHED)
      local_matched++;
  }
#pragma omp critical
  *count += local_matched;
#pragma omp barrier
#ifdef DEBUG
#pragma omp master
  std::printf("ratio matched: %.2f %%\n", float((*count)*100)/nb);
#endif
  return *count;
}

/* --------------------------------------------------------------------------- */
int* Match::computePothenFan(const Graph& graph, int nb) {

  auto tic = timer::now();

  // retrieve a greedy matching using karp-sipser heuristic
  computeKarpSipser(graph, nb);

  //
  int* look_ahead = tasks[0]; // reuse array
#pragma omp for
  for (int i = 0; i < size; ++i)
    look_ahead[i] = 1;

  bool found;
  int level = 0;

  std::stack<int> stack;

  do {
    // (!) mandatory barrier
#pragma omp barrier
#pragma omp master
    path = false;
    found = false;

#pragma omp for
    for (int i = 0; i < size; ++i)
      visited[i] = 0;

#pragma omp for schedule(guided) nowait
    for (size_t i = 0; i < graph.size(); ++i) {
      const int& u = graph[i][0];

      if (matched[u] < 0)
        found = DFS_lookAhead(i, graph, &stack);
    }
#pragma omp critical
    path |= found;
#pragma omp barrier
  } while (path);

  tools::showElapsed(tic, "pothen-fan maximal matching done", 2);
  return matched;
}

/* --------------------------------------------------------------------------- */
bool Match::DFS_lookAhead(int i, const Graph& graph, std::stack<int>* stack) {

  int* look_ahead = tasks[0];

  int j;
  int u;
  int k;

  stack->push(i);

  do {
    j = stack->top();
    u = graph[j][0];
    stack->pop();

    // look ahead step
    for (auto v = graph[j].begin() + look_ahead[u]; v < graph[j].end(); ++v) {
      look_ahead[u]++;
      if (matched[*v] < 0) {
        if (sync::compareAndSwap(visited + u, 0, 1)) {
          matched[u] = *v;
          matched[*v] = u;
          return true;
        }
      }
    }
    // scan unmatched but visited neighbors if not found
    for (auto v = graph[j].begin() + 1; v < graph[j].end(); ++v) {
      if (sync::compareAndSwap(visited + u, 0, 1)) {
        const int& w = mapping[matched[*v]];
        if (w > -1)
          stack->push(w);
      }
    }
  } while (!stack->empty());
  return false;
}
} // namespace trinity
