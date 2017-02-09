/* ------------------------------------*/ 
#include "partition.h"
using namespace trigen;

/* ------------------------------------ */
partit_t::partit_t() :
  size   (0),
  max_p  (0),
  cores  (0),
  off    (nullptr),
  card   (nullptr),
  mapping(nullptr),
  tasks  (nullptr),
  subset (nullptr)
{}
/* ------------------------------------ */
partit_t::partit_t(int max_size, int max_part)
{
  size  = max_size;
  max_p = std::max(6,max_part);
  cores = omp_get_max_threads();
  //
  off    = new int[cores];
  card   = new int[max_p];
  subset = new int*[max_p];
  for(int i=0; i < max_p; ++i)
    subset[i] = new int[size];

  // reuse arrays
  tasks   = subset;    // 2 tasklist maximum
  mapping = subset[2];

#pragma omp parallel
  flush(); 
}
/* ------------------------------------ */
partit_t::~partit_t(){

  for(int i=0; i < max_p; ++i)
    delete [] subset[i];
  delete [] subset;
  delete [] card;
  delete [] off;
}

/* ------------------------------------ */
void partit_t::flush(){

#pragma omp master        
  {
    parts = defect = rounds = 0;
    memset(remain, 0, sizeof(int) * 2);
    memset(  card, 0, sizeof(int) * max_p);
  } 

  for(int i=0; i < 3; ++i)
#pragma omp for
    for(int j=0; j < size; ++j)
      subset[i][j]=0;
}

/* ------------------------------------ */
// nb : nested in a parallel region, so update params accordingly
void partit_t::catalyurek(const mesh_t* mesh){

  // re-init containers
  flush();

  std::vector<int> forbidden;     // forbidden colors for vertex 'i'
  std::vector<int> conflicts;                       
 
  forbidden.resize(max_p,std::numeric_limits<int>::max());
  conflicts.reserve(mesh->nb_nodes);

  // use aliases for clarity
  int* color = mapping;
  const auto& primal = mesh->vicin;

#pragma omp for
  for(int i=0; i < mesh->nb_nodes; ++i){
    if(primal[i].empty())
      continue;

    for(const int& j : primal[i])
      forbidden[color[j]] = i;
    
    for(int c=1; c < max_p; ++c){
      if(forbidden[c] != i){
        color[i] = c;
        break;
      }
    }
    assert(color[i]);
  }

#pragma omp for schedule(guided) nowait
  for(int i=0; i < mesh->nb_nodes; ++i){
    for(const int& j : primal[i]){
      if(i < j && color[i] == color[j]){
        conflicts.push_back(i);   // (!) index not the node value
        break;
      }
    }
  }
  sync::task_reduction(tasks[0], &conflicts, remain,off);

  // - propagation stage
  int k = 0;  // current task list index
 
  while(remain[k]) {
#pragma omp single
    { rounds++; defect += remain[k]; }

#pragma omp for nowait
    for(int i=0; i < remain[k]; ++i){
      const int& v = tasks[k][i];
      if(primal[v].empty())
        continue;      
      for(const int& w : primal[v])
        forbidden[color[w]] = v;

      for(int c=1; c < max_p; ++c){
        if(forbidden[c] != v){
          color[v] = c;
          break;
        }
      }
      assert(color[v]);
    }
      
#pragma omp single
    remain[k^1]=0;

    int cur = 0; // current vertex index in tasklist

#pragma omp for schedule(guided) nowait
    for(int i=0; i < remain[k]; ++i){
      const int& v = tasks[k][i];
      for(const int& w : primal[v]){
        if(v < w && color[v] == color[w]){
          conflicts.push_back(v);   // (!) index not the node value
          break;
        }
      }
    }
    // switch _tasklist
    k ^= 1;
    sync::task_reduction(tasks[k], &conflicts, remain+k,off);
  }

  // -- POST PROCESS
  std::vector<int> *list = new std::vector<int>[max_p];  // purely local

  // (!) OMP reduction doesn't work for class variables
  int nb_col=0;

#pragma omp for nowait
  for(int i=0; i < mesh->nb_nodes; ++i){
    const int& k = mapping[i];
    if(k){
      list[k-1].push_back(i);
      nb_col = std::max(k,nb_col);
    }
  }

#pragma omp critical
 if(parts < nb_col)
   parts = nb_col;
#pragma omp barrier
  
  // b) populate 'subset'
  for(int k=0; k < parts; ++k)
    sync::task_reduction(subset[k], list+k, card+k, off);

  delete [] list;  
}

/* ------------------------------------ */
void partit_t::indep_subset(const graph_t& graph, int nb_nodes){

  if(__builtin_expect(1==nb_nodes,0)){
    subset[0][0] = graph[0][0];
    card[0] = 1;
    return;
  }

  // re-init containers
  flush();
  
  int* color = mapping;

  std::vector<int> forbidden;     // forbidden colors for vertex 'i'
  std::vector<int> heap;                       
 
  forbidden.resize(max_p,std::numeric_limits<int>::max());
  heap.reserve(nb_nodes/cores);

#pragma omp for
  for(int i=0; i < nb_nodes; ++i){
    const int& v = graph[i][0];
    for(auto w = graph[i].begin()+1; w < graph[i].end(); ++w)
      forbidden[color[*w]] = v;
    
    for(int c=1; c < max_p; ++c){
      if(forbidden[c] != v){
        color[v] = c;
        break;
      }
    }
    assert(color[v]);
  }

#pragma omp for schedule(guided) nowait
  for(int i=0; i < nb_nodes; ++i){
    const int& v = graph[i][0];
    if(color[v]==1){
      for(auto w = graph[i].begin()+1; w < graph[i].end(); ++w){
        if(v < *w && color[*w]==1){
          heap.push_back(i);   // (!) index not the node value
          break;
        }
      }
    }
  }
  sync::task_reduction(tasks[0], &heap, remain, off);

  // - propagation stage
  int k = 0;  // current task list index
 
  while(remain[k]) {
#pragma omp single
    { rounds++; defect += remain[k]; }

#pragma omp for nowait
    for(int i=0; i < remain[k]; ++i){
      const int& j = tasks[k][i];
      const int& v = graph[j][0];
      for(auto w = graph[j].begin()+1; w < graph[j].end(); ++w)    
        forbidden[color[*w]] = v;

      for(int c=1; c < max_p; ++c){
        if(forbidden[c] != v){
          color[v] = c;
          break;
        }
      }
      assert(color[v]);
    }
      
#pragma omp single
    remain[k^1]=0;

#pragma omp for schedule(guided) nowait
    for(int i=0; i < remain[k]; ++i){
      const int& j = tasks[k][i];
      const int& v = graph[j][0];
      if(color[v]==1){
        for(auto w = graph[j].begin()+1; w < graph[j].end(); ++w){ 
          if(v < *w && color[*w]==1){
            heap.push_back(j);   // (!) index not the node value
            break;
          }
        }
      }
    }
    // switch _tasklist
    k ^= 1;
    sync::task_reduction(tasks[k], &heap, remain+k, off);
  }  

  // post-process

  // (!) OMP reduction doesn't work for class variables
#pragma omp for schedule(guided) nowait
  for(int i=0; i < nb_nodes; ++i){
    const int& v = graph[i][0];
    if(color[v]==1)
      heap.push_back(v);
  }
  sync::task_reduction(subset[0], &heap, card, off);
}

