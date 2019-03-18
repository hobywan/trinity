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

private:

  // utils
  void saveChrono();
  int elapsed();

  Graph graph;

  struct {
    int nodes;
    int edges;
    int rounds;
    int error;
    int color;
  } nb;

  struct {
    int max;
    int avg;
  } deg;

  struct { double ratio; } stat;
  struct { Time start; }   time;
};
/* --------------------------------------------------------------------------- */
} // namespace trinity