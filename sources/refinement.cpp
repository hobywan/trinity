/* ------------------------------------*/
#include "refinement.h"
/* ------------------------------------*/
using namespace trinity;
/* ------------------------------------*/
refine_t::refine_t(mesh_t* input, int level) :
  mesh    (input),
  off     (input->off),
  activ   (input->activ),
  cores   (input->nb_cores),
  steiner (input->max_node,input->_bucket,2),
  nb_nodes(input->nb_nodes),
  nb_elems(input->nb_elems),
  verbose (input->_verb),
  iter    (input->_iter),
  rounds  (input->_rounds),
  depth   (level)
{
  const int size = input->max_elem;
  index = new int[cores+1];
  edges = new int[size*3];   // 3 data/edge, 3 edges/elem
  elems = new int[size];
  pattern = new char[size];
}
/* ------------------------------------*/
refine_t::~refine_t(){

  delete [] index;
  delete [] edges;
  delete [] elems;
  delete [] pattern;
}
/* ------------------------------------ */
void refine_t::run(stats_t* tot){

  init();

  int time[] = {0,0,0,0,0};
  int stat[] = {0,0};
  int form[] = {0,0,0};

#pragma omp parallel
  {
    std::vector<int> heap[2];

    int level = 0;
    int tid = omp_get_thread_num();

    preprocess(heap);
    timer::save(tic,time);

    do {
      filtering(heap);
      timer::save(tic, time+1);

      if(!nb_tasks)
        break;

      compute_steiner();
      save_stat(level,stat,form);
      timer::save(tic, time+2);

      kernel(tid);
      timer::save(tic, time+3);

#ifdef DEFERRED_UPDATES
      mesh->commit_updates();
      timer::save(tic, time+4);
#endif
      show_stat(level,form);

    } while(++level < depth);

#ifndef DEFERRED_UPDATES
    mesh->fix_all();
    timer::save(tic, time+4);
#endif

    recap(time,stat,form,tot);
    mesh->verif();
  }
}
/* ------------------------------------*/
void refine_t::dissect(int id, int* offset){

#ifdef DEFERRED_UPDATES
  int tid = omp_get_thread_num();
#endif

  const int* n = mesh->get_elem(id);
  // the i'th edge is opposite the i'th node in the element.
  const int s[] = {
    steiner.retrieve(n[1],n[2]),
    steiner.retrieve(n[2],n[0]),
    steiner.retrieve(n[0],n[1])
  };

  if(pattern[id] == 1){

    const int t[] = {id, *offset};

    for(int i=0; i < 3; ++i){
      if(s[i] > -1){
        int j = (i+1)%3;
        int k = (i+2)%3;
        assert(s[i] not_eq n[i]);
        //
        const int elem0[] = {n[i], n[j], s[i]};
        const int elem1[] = {n[i], s[i], n[k]};
        assert(mesh->counterclockwise(elem0));
        assert(mesh->counterclockwise(elem1));
        //
#ifdef DEFERRED_UPDATES
        mesh->deferred_append(tid, s[i], {t[0],t[1]});
        mesh->deferred_append(tid, n[i], t[1]);
        mesh->deferred_remove(tid, n[k], t[0]); // replacement
        mesh->deferred_append(tid, n[k], t[1]);
#else
        mesh->update_stenc(s[i], {t[0],t[1]});
        mesh->update_stenc(n[i], t[1]);
        mesh->update_stenc(n[k], t[1]);
#endif
        //
        mesh->replace_elem(t[0], elem0);
        mesh->replace_elem(t[1], elem1);
        // manually unrolled
        activ[t[0]] = activ[t[1]] = 1;
        break;
      }
    }
  }

  else if(pattern[id] == 2){

    const int t[] = {id, *offset, *offset+1};

    for(int i=0; i < 3; ++i){
      if(s[i] < 0){
        int j = (i+1)%3;
        int k = (i+2)%3;

        double diag[2];
        diag[0] = mesh->edge_length(n[j], s[j]);
        diag[1] = mesh->edge_length(n[k], s[k]);
        // selected diagonal: v[s],n[s+1]
        int v = (diag[0] < diag[1] ? j : k);
        int w = (v == j ? k : j);

        const int elem0[] = {n[i], s[k], s[j]};
        const int elem1[] = {n[j], n[k], s[v]};
        const int elem2[] = {s[j], s[k], n[v]};
        assert(mesh->counterclockwise(elem0));
        assert(mesh->counterclockwise(elem1));
        assert(mesh->counterclockwise(elem2));
        //
#ifdef DEFERRED_UPDATES
        mesh->deferred_append(tid, s[v], {t[0],t[1],t[2]});
        mesh->deferred_append(tid, s[w], {t[0],t[2]});
        mesh->deferred_remove(tid, n[v], t[0]);
        mesh->deferred_append(tid, n[v], {t[1],t[2]});
        mesh->deferred_remove(tid, n[w], t[0]);
        mesh->deferred_append(tid, n[w], t[1]);
#else
        mesh->update_stenc(s[v], {t[0],t[1],t[2]});
        mesh->update_stenc(s[w], {t[0],t[2]});
        mesh->update_stenc(n[v], {t[1],t[2]});
        mesh->update_stenc(n[w], t[1]);
#endif
        //
        mesh->replace_elem(t[0], elem0);
        mesh->replace_elem(t[1], elem1);
        mesh->replace_elem(t[2], elem2);
        // manually unrolled
        activ[t[0]] = activ[t[1]] = activ[t[2]] = 1;
        break;
      }
    }
  }

  else {

    const int t[] = {id, *offset, *offset+1, *offset+2};
    //
    const int elem0[] = {n[0], s[2], s[1]};
    const int elem1[] = {n[1], s[0], s[2]};
    const int elem2[] = {n[2], s[1], s[0]};
    const int elem3[] = {s[0], s[1], s[2]};
    assert(mesh->counterclockwise(elem0));
    assert(mesh->counterclockwise(elem1));
    assert(mesh->counterclockwise(elem2));
    assert(mesh->counterclockwise(elem3));
    //
#ifdef DEFERRED_UPDATES
    mesh->deferred_append(tid, s[0], {t[1],t[2],t[3]});
    mesh->deferred_append(tid, s[1], {t[0],t[2],t[3]});
    mesh->deferred_append(tid, s[2], {t[0],t[1],t[3]});
    mesh->deferred_remove(tid, n[1], t[0]);
    mesh->deferred_append(tid, n[1], t[1]);
    mesh->deferred_remove(tid, n[2], t[0]);
    mesh->deferred_append(tid, n[2], t[2]);
#else
    mesh->update_stenc(s[0], {t[1],t[2],t[3]});
    mesh->update_stenc(s[1], {t[0],t[2],t[3]});
    mesh->update_stenc(s[2], {t[0],t[1],t[3]});
    mesh->update_stenc(n[1], t[1]);
    mesh->update_stenc(n[2], t[2]);
#endif
    //
    mesh->replace_elem(t[0], elem0);
    mesh->replace_elem(t[1], elem1);
    mesh->replace_elem(t[2], elem2);
    mesh->replace_elem(t[3], elem3);
    // manually unrolled
    activ[t[0]] = activ[t[1]] = activ[t[2]] = activ[t[3]] = 1;
  }
  *offset += pattern[id];
}

/* ------------------------------------ */
void refine_t::preprocess(std::vector<int> heap[2]){

  heap[0].reserve(nb_elems);
  heap[1].reserve(nb_elems/cores);

#pragma omp master
{
  old_node = nb_nodes;
  old_elem = nb_elems;
  shift     = 0;
  nb_adds   = 0;
  nb_split  = 0;   // number of elems to be appended
  nb_eval   = 0;
  nb_tasks  = 0;
  nb_stein  = 0;
}

#pragma omp for
  for(int i=0; i < nb_elems; ++i)
    activ[i] = (mesh->active_elem(i) ? 1:0);

  steiner.flush();

#pragma omp single
  round = timer::now();
}
/* ------------------------------------ */
void refine_t::filtering(std::vector<int> heap[2]){
#pragma omp single
{
  nb_eval = nb_adds = nb_split = nb_tasks = nb_stein = 0;
  std::memset(index, 0, (cores+1)*sizeof(int));
}

  // step 1: filtering
  int count[] = {0,0};
  double len[] = {0.,0.,0.};

#pragma omp for schedule(guided) nowait
  for(int i=0; i < nb_elems; ++i){
    pattern[i] = 0;
    // i > offset => new elems
    if(__builtin_expect(activ[i],i > shift)){
      const int* n = mesh->get_elem(i);
      len[0] = mesh->edge_length(n[0],n[1]);
      len[1] = mesh->edge_length(n[1],n[2]);
      len[2] = mesh->edge_length(n[2],n[0]);

      for(int j=0; j < 3; ++j){
        const int k = (j+1)%3;
        // (!) avoid duplicated steiner-points but don't forget boundary edges
        if(len[j] > trinity::l_max){
          ++pattern[i];
          // (!) be aware of diagonals: dont rely only on node activs
          const int opp = mesh->elem_neigh(i,n[j],n[k]);
          if(n[j] < n[k] or opp < 0){
            heap[0].push_back(n[j]);
            heap[0].push_back(n[k]);
            heap[0].push_back(opp);
          }
        }
      }
      if(pattern[i]){
        count[1] += pattern[i];
        heap[1].push_back(i);
      }
      count[0]++;
    }
  }
  sync::fetch_and_add(&nb_eval,count[0]);
  sync::fetch_and_add(&nb_adds,count[1]);

  // reductions on steiner point & bad elems list
  sync::task_reduction(edges, heap  , &nb_split, 3);
  sync::task_reduction(elems, heap+1, &nb_tasks, off);
}
/* ------------------------------------ */
void refine_t::compute_steiner(){

  int count = 0;

  // step 2: compute steiner point and fix new cells ID
#pragma omp for nowait
  for(int i=0; i < nb_split; ++i){
    // implicit index
    const int id = nb_nodes + i;
    const int& v1  = edges[i*3];
    const int& v2  = edges[i*3+1];
    const int& opp = edges[i*3+2];
    // 1) calculate and insert the steiner point
    double* P = mesh->points.data()+(id*2);
    double* M = mesh->tensor.data()+(id*3);
    mesh->calcul_steiner(v1, v2, P, M);

    count++;
    // (!) map (i,j)->id
    int v = std::min(v1,v2);
    int w = std::max(v1,v2);
    steiner.push(v, {w,id});

    if(opp < 0)
      mesh->tags[id] = mask::bound;
  }

  sync::fetch_and_add(&nb_stein,count);
#pragma omp barrier
}
/* ------------------------------------ */
void refine_t::kernel(int tid){
#pragma omp master
  {
    index[0] = shift = nb_elems;
    nb_nodes += nb_split;
    nb_elems += nb_adds;
  }

#pragma omp for schedule(static)
  for(int i=0; i < nb_tasks; ++i)
    index[tid+1] += pattern[elems[i]];  // no data race

  sync::prefix_sum(index,cores,16);

  // step 3: refine bad elems
#pragma omp for schedule(static)
  for(int i=0; i < nb_tasks; ++i)
    dissect(elems[i], index+tid);
}

/* ------------------------------------ */
void refine_t::init(){
#pragma omp master
  {
    if(verbose==1)
       std::printf("%-18s%s","= refinement","...");

    else if(verbose==2)
       std::printf("Process refinement ... ");

    std::fflush(stdout);
    start = round = tic = timer::now();
  }
}
/* ------------------------------------ */
void refine_t::save_stat(int level, int* stat, int* form){
#pragma omp master
  {
    cur_elem = nb_elems;
    stat[0] += nb_eval;
    stat[1] += nb_tasks;
    if(!level){
      form[0] = tools::format(nb_eval);
      form[1] = tools::format(nb_tasks);
      form[2] = tools::format(nb_stein);
    }
  }
}
/* ------------------------------------ */
void refine_t::show_stat(int level, int* form){
#pragma omp single
  {
    if(verbose==2)
    {
       std::printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d filt. \e[0m(%2d %%)\e[0m, "
             "%*d stein \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
       level+1, form[0], nb_eval , (int) (nb_eval *100/cur_elem),
                form[1], nb_tasks, (int) (nb_tasks*100/nb_eval),
                form[2], nb_stein, (int) (nb_stein*100/nb_nodes), timer::round(round));
      std::fflush(stdout);
    }
  }
}
/* ------------------------------------ */
void refine_t::recap(int* time, int* stat, int* form, stats_t* tot){
#pragma omp master
  {
    int end = timer::elapsed_ms(start);
    int span = 0;

    tot->eval += stat[0];
    tot->task += stat[1];
    tot->elap += end;

    tot->step[0] += time[0]+time[1];  // filter
    tot->step[1] += time[2];          // steiner
    tot->step[2] += time[3];          // kernel
    tot->step[3] += time[4];          // repair

    for(int i=0; i < 5; ++i)
      span = std::max(span,time[i]);
    *form = tools::format(span);

    if(!verbose)
       std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100*(++iter)/(4*rounds+1)));

    else if(verbose==1)
    {
       std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m [+%.1f %%]\n",
       (int)std::floor(stat[1]/(end*1e-3)), (float)end/1e3, (float) (nb_elems-old_elem)*100/nb_elems);
    }
    else if(verbose==2)
    {
       std::printf("\n\n");
       std::printf("= nodes: %d old, %d new \e[0m(+%.1f %%)\e[0m\n", old_node, nb_nodes,
            (float) (nb_nodes-old_node)*100/nb_nodes);
       std::printf("= elems: %d old, %d new \e[0m(+%.1f %%)\e[0m\n", old_elem, nb_elems,
            (float) (nb_elems-old_elem)*100/nb_elems);
       std::printf("= rate : %d split/sec (%d tasks) \n", (int)std::floor(stat[1]/(end*1e-3)), stat[1]);
       std::printf("= time per step\n");
       std::printf("  %2d %% preproc \e[32m(%*d ms)\e[0m\n", (int) time[0]*100/end, *form, time[0]);
       std::printf("  %2d %% filter  \e[32m(%*d ms)\e[0m\n", (int) time[1]*100/end, *form, time[1]);
       std::printf("  %2d %% steiner \e[32m(%*d ms)\e[0m\n", (int) time[2]*100/end, *form, time[2]);
       std::printf("  %2d %% kernel  \e[32m(%*d ms)\e[0m\n", (int) time[3]*100/end, *form, time[3]);
       std::printf("  %2d %% fixes   \e[32m(%*d ms)\e[0m\n", (int) time[4]*100/end, *form, time[4]);
       std::printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    std::fflush(stdout);
  }
}

