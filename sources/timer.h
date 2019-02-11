/*
 *                          'timer.h'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *               Copyright (c) 2016 Hoby Rakotoarivelo.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
/* --------------------------------------------------------------------------- */
#include <chrono>
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
using Time = std::chrono::high_resolution_clock::time_point;
} // namespace trinity
/* --------------------------------------------------------------------------- */
namespace trinity { namespace timer {
/* --------------------------------------------------------------------------- */
inline trinity::Time now() {
  return std::chrono::high_resolution_clock::now();
}
/* --------------------------------------------------------------------------- */
inline int elapsed_ms(trinity::Time &tic) {
  const trinity::Time toc = now();
  return static_cast<int>(
    std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()
  );
}
/* --------------------------------------------------------------------------- */
inline void save(trinity::Time &tic, int *time) {
#pragma omp master
  {
    *time += elapsed_ms(tic);
    tic = now();
  }
}
/* --------------------------------------------------------------------------- */
inline int round(trinity::Time &tic) {
  int elapsed = elapsed_ms(tic);
  tic = now();
  return elapsed;
}
/* --------------------------------------------------------------------------- */
}} // namespace trinity::timer
