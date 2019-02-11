/*
 *                          'metric.h'
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
#include "numeric.h"
#include "hessian.h"
/* --------------------------------------------------------------------------- */
namespace trinity {

class Metrics {

public:

  // rule of five
  Metrics() = delete;
  Metrics(const Metrics& other) = delete;
  Metrics& operator=(Metrics other) = delete;
  Metrics(Metrics&& other) noexcept = delete;
  Metrics& operator=(Metrics&& other) noexcept = delete;
  Metrics(Mesh* input, double targ_f, int norm, double min_h, double max_h);
  ~Metrics();

  void computeTensorField(Stats* tot);
  void clear();

private:
  // steps
  void recoverHessianField();
  void normalizeLocally();
  void computeComplexity();
  void normalizeGlobally();
  // kernels
  void computeGradient(int index);
  void computeHessian(int index);

  //
  Mesh* mesh;

  int target;
  int p_norm;
  int chunk;
  int& nb_nodes;
  int& nb_elems;
  int& nb_cores;
  int& verbose;
  int& iter;
  int& rounds;

  double* nabla;
  double* solut;
  double* tens;
  Patch* stenc;

  double h_min;
  double h_max;
  double scale_fact;
  double scale_exp;
  double lambda_min;
  double lambda_max;
  double complexity;

  // timers and stats
  Time start;
  void init();
  void recap(Stats* tot);

};
} // namespace trinity
