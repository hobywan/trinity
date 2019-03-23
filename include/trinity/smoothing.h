/*
 *                          'smoothing.h'
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
#include "partition.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
class Smooth {

public:

  // rule of five
  Smooth() = delete;
  Smooth(const Smooth& other) = delete;
  Smooth& operator=(Smooth other) = delete;
  Smooth(Smooth&& other) noexcept = delete;
  Smooth& operator=(Smooth&& other) noexcept = delete;
  Smooth(Mesh* input, Partit* algo, int level);
  ~Smooth() = default;

  void run(Stats* total = nullptr);

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
  void recap(int* elap, int* stat, int* form, Stats* total);

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
/* -------------------------------------------------------------------------- */
} // namespace trinity
