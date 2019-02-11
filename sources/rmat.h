/*
 *                          'rmat.h'
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
#include "header.h"
#include "timer.h"
#include "tools.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
class RMAT {

  friend class Partit;

public:

  // rule of five
  RMAT();
  RMAT(const RMAT& other) = delete;
  RMAT& operator=(RMAT other) = delete;
  RMAT(RMAT&& other) noexcept = delete;
  RMAT& operator=(RMAT&& other) noexcept = delete;
  ~RMAT();

  // utils
  void reset();
  void load(std::string path);
  void info(std::string name);
  void saveChrono();
  int elapsed();

private:
  Graph graph;
  int   nb_nodes;
  int   nb_edges;
  int   nb_rounds;
  int   nb_error;
  int   nb_color;
  int   deg_max;
  int   deg_avg;
  float ratio;
  Time  start;
};
} // namespace trinity
