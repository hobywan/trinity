/*
 *                          'timer.h'
  *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *                Copyright 2016, Hoby Rakotoarivelo
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once
/* -------------------------------------------------------------------------- */
#include <chrono>
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
using Time = std::chrono::high_resolution_clock::time_point;
} // namespace trinity
/* -------------------------------------------------------------------------- */
namespace trinity { namespace timer {
/* -------------------------------------------------------------------------- */
inline trinity::Time now() {
  return std::chrono::high_resolution_clock::now();
}
/* -------------------------------------------------------------------------- */
inline int elapsed_ms(trinity::Time &tic) {
  const trinity::Time toc = now();
  return static_cast<int>(
    std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()
  );
}
/* -------------------------------------------------------------------------- */
inline void save(trinity::Time &tic, int *time) {
#pragma omp master
  {
    *time += elapsed_ms(tic);
    tic = now();
  }
}
/* -------------------------------------------------------------------------- */
inline int round(trinity::Time &tic) {
  int elapsed = elapsed_ms(tic);
  tic = now();
  return elapsed;
}
/* -------------------------------------------------------------------------- */
}} // namespace trinity::timer
