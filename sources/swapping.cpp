/* ------------------------------------*/ 
#include "swapping.h"
/* ------------------------------------ */
using namespace trigen;
/* ------------------------------------ */
swap_t::swap_t(mesh_t* input) :
  
  mesh    (input),
  off     (input->off),
  fixes   (input->fixes),
  activ   (input->activ),
  cores   (input->nb_cores),
  qualit  (input->qualit.data()),
  nb_nodes(input->nb_cores),
  nb_elems(input->nb_elems),
  verbose (input->_verb),
  iter    (input->_iter),
  rounds  (input->_rounds),
  depth   (20)
{
  const int size = input->max_elem;
  map   = new int[size];
  tasks = new int[size];
  dual.resize(size);
  heuris.init(size,map,off);
}
/* ------------------------------------ */
swap_t::~swap_t(){
  
  delete [] map;
  delete [] tasks;
}
/* ------------------------------------ */
void swap_t::run(stats_t* tot){

  init();
  
  int time[] = {0,0,0,0,0,0};
  int form[] = {0,0,0};
  int stat[] = {0,0,0};

#pragma omp parallel
  {
    std::vector<int> heap;
    
    cache_qualit();
    timer::save(tic,time);

    int level=0;
    do {
   
      filter(&heap);
      timer::save(tic,time+1);

      if(!nb_tasks)
        break;
       
      extract_dual();
      timer::save(tic,time+2);

      heuris.karp_sipser(dual,nb_tasks);
      timer::save(tic,time+3);
      
      kernel();
      save_stat(level,stat,form);
      timer::save(tic,time+4);

      mesh->fix();
      timer::save(tic,time+5);

      show_stat(level++,form);

    } while(nb_comms);        
  
    recap(time,stat,form,tot);
    mesh->verif();
  }
}
/* ------------------------------------ */
int swap_t::swap(int k1, int k2, int index){

  int j,k;
  int f1[] = {-1,-1,-1};
  int f2[] = {-1,-1,-1};

  const int* t1 = mesh->get_elem(k1);
  const int* t2 = mesh->get_elem(k2);

  __builtin_prefetch(t1,1);
  __builtin_prefetch(t2,1);
 
  // rotate (t1,t2) and find shared edge/opposite vertices
  for(int i=0; i < 3; ++i){
    j = (i+1)%3;
    k = (i+2)%3;
    
    // manually unrolled 
    if(t1[i] == t2[1] && t1[j] == t2[0]){  
      f1[0]=t1[i]; f1[1]=t2[2]; f1[2]=t1[k];
      f2[0]=t1[k]; f2[1]=t2[2]; f2[2]=t1[j];
      break;
    }
    if(t1[i] == t2[2] && t1[j] == t2[1]){
      f1[0]=t1[i]; f1[1]=t2[0]; f1[2]=t1[k];
      f2[0]=t1[k]; f2[1]=t2[0]; f2[2]=t1[j];      
      break;
    }      
    if(t1[i] == t2[0] && t1[j] == t2[2]){
      f1[0]=t1[i]; f1[1]=t2[1]; f1[2]=t1[k];
      f2[0]=t1[k]; f2[1]=t2[1]; f2[2]=t1[j];      
      break;
    }
  }
  // check inversions
  if(!mesh->counterclockwise(f1) || !mesh->counterclockwise(f2))
    return 0;
   
  // eval qualit improvement
  const double q[] = {mesh->elem_qual(f1), mesh->elem_qual(f2)};
  const double q_old = std::min(qualit[k1], qualit[k2]);
  const double q_new = std::min(q[0], q[1]); 

  if(q_new < q_old)
    return 0;

  // update mesh
  mesh->replace_elem(k1,f1);
  mesh->replace_elem(k2,f2);  
#ifdef DEFERRED_UPDATES
  int tid = omp_get_thread_num();
  mesh->deferred_remove(tid, f1[0], k2);
  mesh->deferred_remove(tid, f2[2], k1);
  mesh->deferred_append(tid, f1[2], k2); //  N(s[2],k2) 
  mesh->deferred_append(tid, f2[1], k1); //  N(s[3],k1)
#else  
  mesh->update_stenc(f2[0],k2); //  N(s[2],k2)
  mesh->update_stenc(f2[1],k1); //  N(s[3],k1)
#endif  
  // mark nodes as to be fixed
  sync::compare_and_swap(fixes+f1[0],0,1);
  sync::compare_and_swap(fixes+f1[1],0,1);
  sync::compare_and_swap(fixes+f1[2],0,1);
  sync::compare_and_swap(fixes+f2[2],0,1);
  // update cached qualit
  qualit[k1] = q[0];
  qualit[k2] = q[1];
  
  // propagate
  for(const int& t : dual[index])
    sync::compare_and_swap(activ+t,0,1);

  if(map[k2] > -1)
    for(const int& t : dual[map[k2]])
      sync::compare_and_swap(activ+t,0,1); 
      
  return 2; 
}

/* ------------------------------------ */
void swap_t::cache_qualit(){

#pragma omp single
  nb_tasks = nb_comms = 0;

#pragma omp for nowait
  for(int i=0; i < nb_nodes; ++i)
    fixes[i]=0; 
#pragma omp for nowait
  for(int i=0; i < nb_elems; ++i)
    map[i]=-1;  
#pragma omp for
  for(int i=0; i < nb_elems; ++i)
    activ[i]=0; 

#pragma omp for schedule(guided)
  for(int i=0; i < nb_elems; ++i){
    const int* n = mesh->get_elem(i);
    if(__builtin_expect(*n > -1,1)){
      qualit[i] = mesh->elem_qual(i);
      activ[i] = 1;
    }
  }

#pragma omp master
  round = timer::now();  
} 
/* ------------------------------------ */
void swap_t::filter(std::vector<int>* heap){

#pragma omp master
  nb_activ = 0;

  int count=0;
  heap->reserve(nb_elems/cores);
    
//#pragma omp for schedule(dynamic,chunk) nowait
#pragma omp for nowait
  for(int i=0; i < nb_elems; ++i){
    if(__builtin_expect(activ[i],0)){
      activ[i]=0;
      count++;
      if(qualit[i] < trigen::q_min)
        heap->push_back(i);
    }
  }
  sync::task_reduction(tasks, heap, &nb_tasks, off);
  sync::fetch_and_add(&nb_activ,count);
}

/* ------------------------------------ */
void swap_t::extract_dual(){

  const auto& stenc = mesh->stenc;
  const int* deg = mesh->deg;
  
#pragma omp for schedule(guided)
  for(int i=0; i < nb_tasks; ++i){
    
    const int& k = tasks[i];
    const int* n = mesh->get_elem(k);

    dual[i].clear();
    dual[i].reserve(4);
    dual[i].push_back(k);
   
    for(auto t = stenc[*n].begin(); dual[i].size() < 3 && t < stenc[*n].begin()+deg[*n]; ++t){
      if(__builtin_expect(*t == k,0))
        continue;
      
      // manually unrolled
      const int* v = mesh->get_elem(*t);
      if((v[0]==n[0] && v[1]==n[2]) || (v[0]==n[1] && v[1]==n[0]) ||
         (v[1]==n[0] && v[2]==n[2]) || (v[1]==n[1] && v[2]==n[0]) ||     
         (v[2]==n[0] && v[0]==n[2]) || (v[2]==n[1] && v[0]==n[0])){
         dual[i].push_back(*t);      
      }
    }
    
    for(auto t = stenc[n[2]].begin(); t < stenc[n[2]].begin()+deg[n[2]]; ++t){
      if(__builtin_expect(*t == k,0))
        continue;
     
      // manually unrolled
      const int* v = mesh->get_elem(*t);
      if((v[0]==n[2] && v[1]==n[1]) ||
         (v[1]==n[2] && v[2]==n[1]) ||
         (v[2]==n[2] && v[0]==n[1])){ 
        dual[i].push_back(*t); break;
      }    
    }
    map[k] = i;
  }
}
/* ------------------------------------ */
void swap_t::kernel(){

  int succ=0;

#pragma omp single
  nb_comms = 0;  

#pragma omp for schedule(guided) nowait
  for(int index=0; index < nb_tasks; ++index){
    const int& k1 = dual[index][0];
    const int& k2 = heuris.matched[k1];
    if(__builtin_expect(k2 > -1,1))
      succ += swap(k1,k2,index); 
  }
  
  // update nb_comms
  sync::fetch_and_add(&nb_comms,succ);
#pragma omp barrier   
}

/* ------------------------------------ */
void swap_t::init(){
#pragma omp master
  {
    if(verbose==1)
      printf("%-18s%s","= swapping","...");
     
    else if(verbose==2)
      printf("Process swapping ... ");

    fflush(stdout);
    start = round = tic = timer::now();
  }
}
/* ------------------------------------ */
void swap_t::save_stat(int level, int* stat, int* form){
#pragma omp master
  {
    stat[0] += nb_activ;
    stat[1] += nb_tasks;
    stat[2] += nb_comms;
    if(!level){
      form[0] = tools::format(nb_activ);
      form[1] = tools::format(nb_tasks);
      form[2] = tools::format(nb_comms);
    }
  }
}
/* ------------------------------------ */
void swap_t::show_stat(int level, int* form){
#pragma omp single
  {
    if(verbose==2)
    {
      printf("\n= round %2d. %*d tasks \e[0m(%2d %%)\e[0m, %*d filt. \e[0m(%2d %%)\e[0m, "
             "%*d comm. \e[0m(%2d %%) \e[32m(%d ms)\e[0m",
       level+1, form[0], nb_activ, (int) (nb_activ*100/nb_elems),
                form[1], nb_tasks, (int) (nb_tasks*100/nb_activ),
                form[2], nb_comms, (int) (nb_comms*100/nb_tasks), timer::round(round));  
      fflush(stdout);
    }
  }
}
/* ------------------------------------ */
void swap_t::recap(int* time, int* stat, int* form, stats_t* tot){
#pragma omp master
  {
    int end = timer::elapsed_ms(start);
    
    tot->eval += stat[0];
    tot->task += stat[1];
    tot->elap += end;    

    // manually unrolled
    tot->step[0] += time[0]+time[1];  // qualit + filter
    tot->step[1] += time[2];          // dual
    tot->step[2] += time[3];          // match
    tot->step[3] += time[4];          // kernel    
    tot->step[4] += time[5];          // repair
  
    int span=0.;
    for(int i=0; i < 6; ++i)
      span = std::max(span,time[i]);
    *form = tools::format(span); 

    if(!verbose)
      printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100*(++iter)/(4*rounds+1))); 
    
    else if(verbose==1)
    {
      printf("%10d task/sec \e[32m(%4.2f s)\e[0m\n", (int)std::floor(stat[1]/(end*1e-3)), (float)end/1e3);
    }
    else if(verbose==2)
    {
      printf("\n\n"); 
      printf("= rate : %d flip/sec (%d tasks) \n", (int)std::floor(stat[1]/(end*1e-3)), stat[1]);
      printf("= time per step\n");
      printf("  %2d %% qualit \e[32m(%*d ms)\e[0m\n", (int) time[0]*100/end, *form, time[0]);
      printf("  %2d %% filter \e[32m(%*d ms)\e[0m\n", (int) time[1]*100/end, *form, time[1]);
      printf("  %2d %% dual   \e[32m(%*d ms)\e[0m\n", (int) time[2]*100/end, *form, time[2]);
      printf("  %2d %% match  \e[32m(%*d ms)\e[0m\n", (int) time[3]*100/end, *form, time[3]);
      printf("  %2d %% kernel \e[32m(%*d ms)\e[0m\n", (int) time[4]*100/end, *form, time[4]);
      printf("  %2d %% fixes  \e[32m(%*d ms)\e[0m\n", (int) time[5]*100/end, *form, time[5]);
      printf("done. \e[32m(%d ms)\e[0m\n", end);
      tools::separator();
    }
    fflush(stdout);
  }
}
