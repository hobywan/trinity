#pragma once

#include <chrono>

namespace trinity {
  using time_t = std::chrono::high_resolution_clock::time_point;
}

namespace trinity { namespace timer {

inline trinity::time_t now() {
  return std::chrono::high_resolution_clock::now();
}

inline int elapsed_ms(trinity::time_t &tic) {
  const trinity::time_t toc = now();
  return static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count());
}

inline void save(trinity::time_t &tic, int *time) {
#pragma omp master
  {
    *time += elapsed_ms(tic);
    tic = now();
  }
}

inline int round(trinity::time_t &tic) {
  int elapsed = elapsed_ms(tic);
  tic = now();
  return elapsed;
}

}} // namespace trinity::timer
