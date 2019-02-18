/*
 *                          'swapping.h'
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
#include "matching.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
class Swap {

public:

  // rule of five
  Swap() = delete;
  Swap(const Swap& other) = delete;
  Swap& operator=(Swap other) = delete;
  Swap(Swap&& other) noexcept = delete;
  Swap& operator=(Swap&& other) noexcept = delete;
  explicit Swap(Mesh* input);
  ~Swap();

  void run(Stats* total = nullptr);

private:
  // steps
  void cacheQuality();
  void filterElems(std::vector<int>* heap);
  void extractDualGraph();
  void processFlips();

  // kernel
  int swap(int k1, int k2, int idx);

  // stats
  void initialize();
  void saveStat(int level, int* stat, int* form);
  void showStat(int level, int* form);
  void recap(int* elap, int* stat, int* form, Stats* total);

  Mesh* mesh;
  Graph dual;
  Match heuris;

  struct {
    int* match;
    int* list;
    int  depth;
  } task;

  struct {
    int*  off;
    char* fixes;
    char* activ;
  } sync;

  struct {
    int activ;
    int tasks;
    int commit;
    int count;  // for profiling only
  } nb;

  struct { double* qualit; } geom;

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