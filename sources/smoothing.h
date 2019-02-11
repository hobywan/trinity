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
/* --------------------------------------------------------------------------- */
#include "mesh.h"
#include "hashtable.h"
#include "numeric.h"
#include "partition.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
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

  // kernel
  int moveSmartLaplacian(int id);

  // stats
  void initialize();
  void saveStat(int level, int* stat, int* form);
  void showStat(int level, int* form);
  void recap(int* elap, int* stat, int* form, Stats* tot);

  Mesh*   mesh;
  Partit* heuris;

  struct { char* activ; }    sync;
  struct { double* qualit; } geom;
  struct { int depth; }      task;

  struct {
    int tasks;
    int commit;
  } nb;

  struct {
    Time start;
    Time iter;
    Time tic;
  } time;

  // unpackable
  int& cores;
  int& nb_nodes;
  int& nb_elems;
  int& verbose;
  int& iter;
  int& rounds;
};
/* --------------------------------------------------------------------------- */
} // namespace trinity
