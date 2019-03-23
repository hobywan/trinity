/*
 *                        'coarsening.h'
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
#include "sync.h"
#include "numeric.h"
#include "partition.h"
/* -------------------------------------------------------------------------- */
namespace trinity {

class Coarse {

public:

  // rule of five
  Coarse() = delete;
  Coarse(const Coarse& other) = delete;
  Coarse& operator=(Coarse other) = delete;
  Coarse(Coarse&& other) noexcept = delete;
  Coarse& operator=(Coarse&& other) noexcept = delete;
  ~Coarse() = default;

  Coarse(Mesh* input, Partit* algo);

  void run(Stats* total = nullptr);

private:
  // steps
  void preProcess();
  void filterPoints(std::vector<int>* heap);
  void extractSubGraph();
  void processPoints();

  // kernels
  void identifyTarget(int source);
  void collapseEdge(int source, int destin);

  // stats
  void initialize();
  void saveStat(int level, int* stat, int* form);
  void showStat(int level, int* form);
  void recap(int* elap, int* stat, int* form, Stats* total);

  //
  Mesh*   mesh;
  Graph   primal;
  Partit* heuris;

  struct {
    int* target;
    int* filter;
    int* indep;
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
  } nb;

  struct {
    Time start;
    Time iter;
    Time tic;
  } time;

  // unpackable fields
  int& cores;
  int& nb_nodes;
  int& nb_elems;
  int& nb_indep;
  int& verbose;
  int& iter;
  int& rounds;
};
/* -------------------------------------------------------------------------- */
} // namespace trinity
