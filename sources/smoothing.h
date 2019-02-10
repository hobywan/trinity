/*
 *                          'smoothing.h'
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
/* ------------------------------------ */
#include "mesh.h"
#include "hashtable.h"
#include "numeric.h"
#include "partition.h"
/* ------------------------------------ */
namespace trinity {
/* ------------------------------------ */
class Smooth {

public:

  // rule of five
  Smooth() = delete;
  Smooth(const Smooth& other) = delete;
  Smooth& operator=(Smooth other) = delete;
  Smooth(Smooth&& other) noexcept = delete;
  Smooth& operator=(Smooth&& other) noexcept = delete;
  Smooth(Mesh* input, Partit* algo, int level);
  ~Smooth();

  void run(Stats* tot);

private:

  // steps
  void preProcess();
  void cacheQuality();
  void movePoints();

  Mesh* mesh;
  Partit* heuris;
  //
  char* activ;
  double* qualit;

  // counters
  int& cores;
  int& nb_nodes;
  int& nb_elems;
  int& verbose;
  int& iter;
  int& rounds;
  int depth;
  int nb_tasks;
  int nb_comms;

  // processFlips
  int moveSmartLaplacian(int id);

  // timers and stats
  Time start;
  Time round;
  Time tic;
  //
  void init();
  void saveStat(int level, int* stat, int* form);
  void showStat(int level, int* form);
  void recap(int* time, int* stat, int* form, Stats* tot);
};
} // namespace trinity
