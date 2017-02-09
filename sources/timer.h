#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>

namespace trigen {
  // shortcut
  using time_t = std::chrono::high_resolution_clock::time_point;

  namespace timer {
    
    inline std::chrono::high_resolution_clock::time_point now(){
      return std::chrono::high_resolution_clock::now();
    } 
    
    inline int elapsed_ms(std::chrono::high_resolution_clock::time_point& tic){
      std::chrono::high_resolution_clock::time_point toc = std::chrono::high_resolution_clock::now();
      return std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic).count();
    }     
    
    /* -------------------------------- */
    inline void save(std::chrono::high_resolution_clock::time_point& tic, int* time){
    #pragma omp master
      {
        *time += elapsed_ms(tic);
        tic = now();
      }
    }  
    inline int round(std::chrono::high_resolution_clock::time_point& tic){
      int elap = elapsed_ms(tic);
      tic = now();
      return elap;
    }             
  }
}  
#endif
