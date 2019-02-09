/* ------------------------------------ */
#include "indep.h"
/* ------------------------------------ */
using namespace trinity;

/* ------------------------------------ */
indep_t::indep_t(mesh_t* input, partit_t* algo) :
  mesh(input),
  cores(input->nb_cores),
  added(input->fixes),
  weight(input->qualit.data()),
  remain(0),
  rounds(0),
  nb_nodes(0),
  subset(algo->subset[0]),
  nb_indep(algo->card[0]) {
  auto size = mesh->max_node;
  activ = new char[size];
  added = new char[size];
  subset = new int[size];
}

/* ------------------------------------ */
indep_t::~indep_t() {
  delete[] activ;
  delete[] added;
  delete[] subset;
}

/* ------------------------------------ */
void indep_t::flush() {

  assert(nb_nodes);

#pragma omp for nowait
  for (int i = 0; i < mesh->max_node; ++i)
    subset[i] = 0;

#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    activ[i] = 0;

#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    weight[i] = 0.;

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    added[i] = 0;
}

/* ------------------------------------ */
void indep_t::metivier(const graph_t& graph, int nb) {

  // init
  int tid = omp_get_thread_num();
  int nb_tasks = 0;
  int nb_added = 0;
  double max;

  trinity::random_engine engine(tid);

#pragma omp master
  nb_nodes = nb;
#pragma omp single
  rounds = remain = nb_indep = 0;

  flush();

#pragma omp for nowait
  for (int i = 0; i < nb_nodes; ++i)
    activ[graph[i][0]] = 1;

  do {

#pragma omp barrier
    nb_tasks = 0;
    nb_added = 0;


#pragma omp single
    {
      rounds++;
       std::printf("remain: %d\n", remain);
      remain = 0;
    }

    // generate weights
#pragma omp for schedule(guided)
    for (int i = 0; i < nb_nodes; ++i) {
      const int& u = graph[i][0];
      if (activ[u]) {
        weight[u] = (engine() * trinity::max_f);
        nb_tasks++;
      } else
        weight[u] = 0.;   // reset
    }

    // vertex selection
#pragma omp for schedule(guided) nowait
    for (int i = 0; i < nb_nodes; ++i) {
      const int& u = graph[i][0];
      if (activ[u]) {
        max = 0.;
        for (auto j = graph[i].begin() + 1; j < graph[i].end(); ++j)
          if (max < weight[*j])
            max = weight[*j];

        // i is a local maximum
        if (weight[u] > max) {
          added[u] = 1;  // deterministic
          activ[u] = 0;
        }
      }
    }
    sync::fetch_and_add(&remain, nb_tasks);
#pragma omp barrier

#pragma omp for schedule(guided)
    for (int i = 0; i < nb_nodes; ++i) {
      const int& u = graph[i][0];
      if (added[u])
        for (auto j = graph[i].begin() + 1; j < graph[i].end(); ++j)
          sync::compare_and_swap(activ + (*j), 1, 0);
    }
  } while (remain and rounds < 15);

  verif(graph);

  std::vector<int> heap;
  heap.reserve(nb_nodes / cores);

  // post-process
#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_nodes; ++i) {
    const int& u = graph[i][0];
    if (__builtin_expect(added[u], 0))
      heap.push_back(u);
  }
  sync::task_reduction(subset, &heap, &nb_indep, mesh->off);

#pragma omp single
   std::printf("rounds: %d, nb_indep: %d \e[33m(%d %%)\e[0m\n",
         rounds, nb_indep, (int) (nb_indep * 100 / nb_nodes));

}

/* ------------------------------------ */
void indep_t::verif(const graph_t& graph) {

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i) {
    const int& u = graph[i][0];
    if (added[u]) {
      for (auto j = graph[i].begin() + 1; j < graph[i].end(); ++j) {
        if (added[*j])
           std::printf("w[%d]: %f, w[%d]: %f\n", u, weight[u], *j, weight[*j]);
        assert(!added[*j]);
      }
    }
  }
}

