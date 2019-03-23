/*
 *                          'matching.cpp'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *                Copyright 2016, Hoby Rakotoarivelo.
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

#include "trinity/matching.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
Match::~Match() {

  delete[] task.matched;
  delete[] sync.visited;
  delete[] sync.degree;
}

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
void Match::matchAndUpdate(int vertex, const Graph& graph, std::stack<int>* stack) {

  stack->push(vertex);

  int cur;
  int nxt;
  int current;

  do {
    cur = stack->top();
    current = graph[cur][0];
    stack->pop();

    if (not sync::compareAndSwap(sync.visited + current, 0, 1))
      continue;

    for (auto neigh = graph[cur].begin() + 1; neigh < graph[cur].end(); ++neigh) {
      if (sync::compareAndSwap(sync.visited + (*neigh), 0, 1)) {
        task.matched[current] = *neigh;
        task.matched[*neigh] = MATCHED;  // avoid duplicates

        nxt = task.mapping[*neigh];
        if (nxt < 0)
          continue;

        // update the degree of neighbors of v
        // and recursive call to match the new vertex w of degree=1
        for (auto far = graph[nxt].begin() + 1; far < graph[nxt].end(); ++far) {
          const int& next = task.mapping[*far];
          if (next > -1 and sync::fetchAndSub(sync.degree + (*far), char(1)) == 2)
            stack->push(next);
        }
        break;
      }
    }
  } while (not stack->empty());
}

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------- */
bool Match::lookAheadDFS(int vertex, const Graph& graph, std::stack<int>* stack) {

  auto look_ahead = task.lists[0];

  int cur;
  int current;

  stack->push(vertex);

  do {
    cur = stack->top();
    current = graph[cur][0];
    stack->pop();

    // look ahead step
    for (auto check = graph[cur].begin() + look_ahead[current];
              check < graph[cur].end(); ++check) {
      look_ahead[current]++;
      if (task.matched[*check] < 0) {
        if (sync::compareAndSwap(sync.visited + current, 0, 1)) {
          task.matched[current] = *check;
          task.matched[*check] = current;
          return true;
        }
      }
    }
    // scan unmatched but sync.visited neighbors if not found
    for (auto neigh = graph[cur].begin() + 1; neigh < graph[cur].end(); ++neigh) {
      if (sync::compareAndSwap(sync.visited + current, 0, 1)) {
        const int& next = task.mapping[task.matched[*neigh]];
        if (next > -1)
          stack->push(next);
      }
    }
  } while (not stack->empty());

  return false;
}
/* -------------------------------------------------------------------------- */
} // namespace trinity
