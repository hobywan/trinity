/*
 *                          'swapping.h'
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
#include "mesh.h"
#include "hashtable.h"
#include "numeric.h"
#include "matching.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
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
  int swap(int id1, int id2, int idx);

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
/* -------------------------------------------------------------------------- */
} // namespace trinity