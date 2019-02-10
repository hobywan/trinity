#pragma once

#include <chrono>

namespace trinity {
  using Time = std::chrono::high_resolution_clock::time_point;
}

namespace trinity { namespace timer {

inline trinity::Time now() {
  return std::chrono::high_resolution_clock::now();
}

inline int elapsed_ms(trinity::Time &tic) {
  const trinity::Time toc = now();
  return static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count());
}

inline void save(trinity::Time &tic, int *time) {
#pragma omp master
  {
    *time += elapsed_ms(tic);
    tic = now();
  }
}

inline int round(trinity::Time &tic) {
  int elapsed = elapsed_ms(tic);
  tic = now();
  return elapsed;
}

}} // namespace trinity::timer
