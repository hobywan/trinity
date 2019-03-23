/*
 *                          'metric.h'
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
#include "numeric.h"
#include "hessian.h"
/* -------------------------------------------------------------------------- */
namespace trinity {

class Metrics {

public:

  // rule of five
  Metrics() = delete;
  Metrics(const Metrics& other) = delete;
  Metrics& operator=(Metrics other) = delete;
  Metrics(Metrics&& other) noexcept = delete;
  Metrics& operator=(Metrics&& other) noexcept = delete;
  Metrics(Mesh* input_mesh, double target_factor, int Lp_norm, double h_min, double h_max);
  ~Metrics();

  void run(Stats* total = nullptr);

private:

  // steps
  void recoverHessianField();
  void normalizeLocally();
  void computeComplexity();
  void normalizeGlobally();

  // kernels
  void computeGradient(int index);
  void computeHessian(int index);

  // stats
  void initialize();
  void recap(Stats* total);
  void clear();

  Mesh* mesh;

  struct {
    double* gradient;
    double* solut;
    double* tensor;
    Patch*  stencil;
    double  complexity;
  } field;

  struct {
    int    chunk;
    int    target;
    int    norm;
    double h_min;
    double h_max;
    struct { double min, max; }  eigen;
    struct { double fact, exp; } scale;
  } param;

  struct { Time start; } time;

  int& nb_nodes;
  int& nb_elems;
  int& nb_cores;
  int& verbose;
  int& iter;
  int& rounds;

};
} // namespace trinity
