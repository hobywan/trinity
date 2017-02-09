/* ------------------------------------*/ 
#pragma once
/* ------------------------------------*/ 
#include "sync.h"
#include "tools.h"
/* ------------------------------------*/ 
namespace trigen {

  template<typename type_t>
  class hashtable {

public:
    
    /* ------------------------------------*/ 
    hashtable() : bucket(nullptr), off(nullptr), stride(0), capa(0), size(0) {}      
   /* ------------------------------------*/ 
    hashtable(size_t s, size_t buck, size_t n) :
      size(s), capa(buck), stride(n)
    {
      // memalloc
      off = new int[size];
      bucket = new type_t*[size];

#pragma omp parallel for      
      for(int i=0; i < size; ++i)
        bucket[i] = new type_t[capa];
    }      
   /* ------------------------------------*/ 
    ~hashtable(){
      for(int i=0; i < size; ++i)
        delete [] bucket[i];
      delete [] bucket;
      delete [] off;
    }     
    /* ------------------------------------*/ 
    inline int generate_key(int i, int j, int scale, int nb_cores) const{
      return tools::hash(std::min(i,j)) % (scale * nb_cores);
    }
    /* ------------------------------------*/ 
    inline int capacity() const{ return size; }
    /* ------------------------------------*/ 
    void push(int key, const std::initializer_list<type_t>& val){
      assert(val.size()==stride);
      
      int j = sync::fetch_and_add(off+key,stride);
      assert((j+stride) < capa);
      
      for(int i=0; i < stride; ++i)
        bucket[key][j+i] = *(val.begin()+i);
    }
    /* ------------------------------------*/ 
    inline int retrieve(int v1, int v2) const{
      
      const int index[] = {std::min(v1,v2),std::max(v1,v2)};
    
      for(int k=0; k < off[*index]-1; k += 2)
        if(bucket[*index][k] == *(index+1))
          return bucket[*index][k+1];
      // not found      
      return -1;
    } 
    /* ------------------------------------*/ 
    inline void flush(){
#pragma omp for      
      for(int i=0; i < size; ++i)
        memset(bucket[i], -1, capa * sizeof(int));
#pragma omp for
      for(int i=0; i < size; ++i)
        off[i] = 0;      
    }      
    
private:

    type_t ** bucket;
    int* off;
    
    int capa;
    int size;
    int stride;
             
  };  
}
