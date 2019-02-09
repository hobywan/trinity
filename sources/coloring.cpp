/* ------------------------------------*/
#include "partition.h"

using namespace trinity;

/* ------------------------------------ */
void partit_t::pseudo_color(const graph_t& graph, std::vector<int>& forbidden, int i) {

  // use alias for clarity
  int* color = mapping;

  const int& v = graph[i][0];
  for (size_t j = 1; j < graph[i].size(); ++i) {
    const int& w = graph[i][j];
    forbidden[color[w]] = v;
  }
  for (int c = 1; c < max_p; ++c) {
    if (forbidden[c] not_eq v) {
      color[v] = c;
      break;
    }
  }
  assert(color[v]);
}

/* ------------------------------------ */
bool partit_t::detect_error(const graph_t& graph, std::vector<int>& conflicts, int i) {

  // use alias for clarity
  int* color = mapping;
  //
  const int& v = graph[i][0];
  for (size_t j = 1; j < graph[i].size(); ++i) {
    const int& w = graph[i][j];
    if (v < w and color[v] == color[w]) {
      conflicts.push_back(i);   // (!) index not the node value
      return true;
    }
  }
  return false;
}

/* ------------------------------------ */
void partit_t::reduce_maxcol(int nb_nodes) {

  // -- POST PROCESS
  // (!) OMP reduction doesn't work for class variables
  int nb_col = 0;
#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i) {
    const int& k = mapping[i];
    if (k > nb_col)  // k > 0 automatically
      nb_col = k;
  }

#pragma omp critical
  if (parts < nb_col)
    parts = nb_col;
#pragma omp barrier

}

/* ------------------------------------ */
void partit_t::catalyurek(RMAT* rmat) {

#pragma omp master
  rmat->save_chrono();

  // re-init containers
  flush();

  std::vector<int> forbidden;     // forbidden colors for vertex 'i'
  std::vector<int> conflicts;

  // use aliases for clarity
  int* color = mapping;
  const int& nb_nodes = rmat->nb_nodes;
  const auto& graph = rmat->graph;

  // max_p is resolved here
  forbidden.resize(max_p, std::numeric_limits<int>::max());
  conflicts.reserve(nb_nodes);

  // --------------------------------

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    pseudo_color(graph, forbidden, i);
#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_nodes; ++i)
    detect_error(graph, conflicts, i);

  sync::task_reduction(tasks[0], &conflicts, remain, off);

  // - propagation stage
  int k = 0;  // current task list index

  while (remain[k]) {
#pragma omp single
    {
      rounds++;
      defect += remain[k];
    }

#pragma omp for nowait
    for (int i = 0; i < remain[k]; ++i)
      pseudo_color(graph, forbidden, tasks[k][i]);
#pragma omp for schedule(guided) nowait
    for (int i = 0; i < remain[k]; ++i)
      detect_error(graph, conflicts, tasks[k][i]);

#pragma omp single
    remain[k ^ 1] = 0;
    // switch _tasklist
    k ^= 1;
    sync::task_reduction(tasks[k], &conflicts, remain + k, off);
  }

  reduce_maxcol(nb_nodes);

  // final step: print stats
#pragma omp master
  {
    rmat->nb_rounds = rounds;
    rmat->nb_error = defect;
    rmat->nb_color = parts;
    rmat->info("catalyurek et al.");
  }
}

/* ------------------------------------ */
void partit_t::gebremedhin(RMAT* rmat) {


#pragma omp master
  rmat->save_chrono();

  // re-init containers
  flush();

  std::vector<int> forbidden;     // forbidden colors for vertex 'i'
  std::vector<int> conflicts;

  // use aliases for clarity
  int* color = mapping;
  const int& nb_nodes = rmat->nb_nodes;
  const auto& graph = rmat->graph;

  // max_p is resolved here
  forbidden.resize(max_p, std::numeric_limits<int>::max());
  conflicts.reserve(nb_nodes);

  // ----------------------------
#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    pseudo_color(graph, forbidden, i);
#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_nodes; ++i)
    detect_error(graph, conflicts, i);


  // step 3: serial resolution of conflicts
  int nb_conflicts = conflicts.size();

#pragma omp atomic
  defect += nb_conflicts;

  for (int i = 0; i < nb_conflicts; ++i)
    pseudo_color(graph, forbidden, conflicts[i]);
#pragma omp barrier

  reduce_maxcol(nb_nodes);

  // final step: print stats
#pragma omp master
  {
    rmat->nb_rounds = rounds;
    rmat->nb_error = defect;
    rmat->nb_color = parts;
    rmat->info("gebremedhin-manne");
  }
}

/* ------------------------------------ */
void partit_t::rokos_gorman(RMAT* rmat) {

#pragma omp master
  rmat->save_chrono();

  // re-init containers
  flush();

  std::vector<int> forbidden;     // forbidden colors for vertex 'i'
  std::vector<int> conflicts;

  // use aliases for clarity
  int* color = mapping;
  const int& nb_nodes = rmat->nb_nodes;
  const auto& graph = rmat->graph;

  // max_p is resolved here
  forbidden.resize(max_p, std::numeric_limits<int>::max());
  conflicts.reserve(nb_nodes);
  // --------------------------------


#pragma omp for schedule(static)
  for (int i = 0; i < nb_nodes; ++i)
    pseudo_color(graph, forbidden, i);

#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_nodes; ++i) {
    // if defective, recolor immediately
    if (detect_error(graph, conflicts, i))
      pseudo_color(graph, forbidden, i);
  }
  sync::task_reduction(tasks[0], &conflicts, remain, off);

  int tid = omp_get_thread_num();
  int k = 0;

  // - propagation stage
  while (remain[k]) {
#pragma omp single
    {
      rounds++;
      defect += remain[k];
    }

#pragma omp for schedule(static)
    for (int i = 0; i < nb_nodes; ++i)
      pseudo_color(graph, forbidden, tasks[k][i]);

#pragma omp for schedule(guided) nowait
    for (int i = 0; i < nb_nodes; ++i) {
      // if defective, recolor immediately
      if (detect_error(graph, conflicts, tasks[k][i]))
        pseudo_color(graph, forbidden, tasks[k][i]);
    }

#pragma omp single
    remain[k ^ 1] = 0;
    // switch _tasklist
    k ^= 1;
    sync::task_reduction(tasks[k], &conflicts, remain + k, off);
  }

  reduce_maxcol(nb_nodes);

  // final step: print stats
#pragma omp master
  {
    rmat->nb_rounds = rounds;
    rmat->nb_error = defect;
    rmat->nb_color = parts;
    rmat->info("rokos-gorman");
  }

}

/* ------------------------------------ */
void partit_t::process_benchmark(int nb_rounds) {

  std::string path[] = {"data/RMAT_ER", "data/RMAT_G", "data/RMAT_B"};

  RMAT graph;
  for (const auto& i : path) {
    graph.load(i);

    // (!) don't merge loops to ease stats postprocessing
    for (int iter = 0; iter < nb_rounds; ++iter)
      catalyurek(&graph);
    for (int iter = 0; iter < nb_rounds; ++iter)
      gebremedhin(&graph);
    for (int iter = 0; iter < nb_rounds; ++iter)
      rokos_gorman(&graph);
  }
}

