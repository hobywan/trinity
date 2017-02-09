/* ------------------------------------*/ 
#pragma once
/* ------------------------------------*/ 
#include "header.h"
#include "timer.h"
#include "sync.h"
#include "numeric.h"
/* ------------------------------------*/ 
namespace trigen {
  
  class mesh_t {
   
    friend class metric_t;
    friend class refine_t;
    friend class coarse_t;
    friend class swap_t;
    friend class smooth_t;
    friend class indep_t;
  
public:

    int nb_nodes;
    int nb_elems;
    int nb_cores;  
    int nb_bucket;

    // could use a std::set but performance suffers
    std::vector<std::vector<int>> stenc;
    std::vector<std::vector<int>> vicin;
    
    std::vector<int>    elems;
    std::vector<double> points;      // stride=2
    std::vector<double> tensor;      // stride=3
    std::vector<double> solut;       // only on metric field calculation
    std::vector<double> qualit;      // 

     mesh_t(int size[2], int bucket, int depth, int verbosity, int rounds);
    ~mesh_t();

    // global routines
    void realloc();
    void round_robin_flush();
    void compress_storage();
    void reorder_cells();
    void extract_primal();
    void extract_dual(graph_t* dual) const;

    // repair topology
    void rebuild();
    void fix();
    void fix_all();
    void init_activ();
    bool verif() const;

    // I/O
    void load (const std::string& path, const std::string& solu);
    void store(const std::string& path) const;
    void store_primal(const std::string& path) const;

    // topological queries
    patch_t vicin_dist(int id, int deg) const;
    const int* elem_coord(int id, double* p) const;
    int  elem_neigh(int id, int i, int j) const;

    inline const int* get_elem(int i) const { return elems.data()+(i*3); }
    inline bool active_node(int i) const { return !stenc[i].empty(); }
    inline bool active_elem(int i) const { return elems[i*3] > -1; }
    inline bool bound(int i)       const { return (tags[i] & mask::bound); }
    inline bool corner(int i)      const { return (tags[i] & mask::corner); }
    inline auto stenc_beg(int i)   const { return stenc[i].begin(); }
    inline auto stenc_end(int i)   const { return stenc[i].begin()+deg[i]; }
    inline int node_max()          const { return max_node; }  
    inline int elem_max()          const { return max_elem; }  

    // topological updates
    void replace_elem(int id, const int* v); 
    void erase_elem(int id); 
    void update_stenc(int i, int t);
    void update_stenc(int i, const std::initializer_list<int>& t);
    void copy_stenc(int i, int j, int nb_rm);

    // deferred updates (for perfs comparison)
#ifdef DEFERRED_UPDATES
    void init_updates();
    void commit_updates();
    void reset_updates();
    void deferred_append(int tid, int i, int t);
    void deferred_append(int tid, int i, const std::initializer_list<int>& t);
    void deferred_remove(int tid, int i, int t);    
#endif

    // geometrical queries
    double edge_length(int i, int j) const;
    double elem_qual(int i) const;
    double elem_qual(const int* t) const;
    bool counterclockwise(const int* t) const;
    void calcul_steiner(int i, int j, double* p, double* m) const;
    void compute_quality(double q[3]);

private:

    int _scale;        // capacity scale factor 
    int _bucket;       // bucket capacity
    int _depth;        // max refinement/smoothing level
    int _verb;         // verbosity level
    int _iter;
    int _rounds;       // remeshing rounds
    int max_node;      // node max capacity
    int max_elem;      // elem max capacity

    int nb_activ_elem;
    
    int*  deg;         // stencil degree
    int*  off;         // offset per thread for prefix sum
    char* fixes;       // node marker for topology fixes
    char* activ;       // node marker for kernel propagation
    uint8_t* tags;     // node tags for geometrical features [bound|corner]

    time_t tic;

#ifdef DEFERRED_UPDATES
    struct updates_t {
      std::vector<int> add;
      std::vector<int> rem;
    };      
    // table for pending updates
    std::vector<std::vector<updates_t>> deferred;
    // scaling factor
    const int def_scale_fact = 32;
#endif    
        
    // internal routines
    void reduce_offset_cells(int* shift, int* map, int* N, int* range, int rank);
  };
} // namespace
