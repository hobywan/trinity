/*
 *                          'coarsening.h'
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
/* ------------------------------------*/
#include "mesh.h"
#include "sync.h"
#include "numeric.h"
#include "partition.h"
/* --------------------------------------------------------------------------- */
namespace trinity {

class Coarse {

public:

  // rule of five
  Coarse() = delete;
  Coarse(const Coarse& other) = delete;
  Coarse& operator=(Coarse other) = delete;
  Coarse(Coarse&& other) noexcept = delete;
  Coarse& operator=(Coarse&& other) noexcept = delete;
  Coarse(Mesh* input, Partit* algo);
  ~Coarse();

  void run(Stats* tot);

private:
  // steps
  void preProcess();
  void filterPoints(std::vector<int>* heap);
  void extractPrimalGraph();
  void processPoints();

  //
  Mesh* mesh;
  Graph primal;
  Partit* heuris;

  // tasklists
  int* target;
  int* filter;
  int* indep;
  int* off;
  char* fixes;
  char* activ;

  // counters
  int& cores;
  int& nb_nodes;
  int& nb_elems;
  int& nb_indep;
  int& verbose;
  int& iter;
  int& rounds;
  int depth;
  int nb_activ;
  int nb_tasks;

  // kernels
  void identify(int id);
  void collapse(int i, int j);

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
}
