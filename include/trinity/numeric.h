/*
 *                          'numeric.h'
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
#include "header.h"
/* -------------------------------------------------------------------------- */
namespace trinity { namespace numeric {
/* -------------------------------------------------------------------------- */
void eigenDecompose(
  const double* matrix, double* val, double* vec1, double* vec2
);
void interpolateTensor(
  const double* matrix, double* result, int n
);
void doKroneckerProduct(
  const double* vec1, const double* vec2, double* matrix
);
double computeQuality(
  const double* point_a, const double* point_b, const double* point_c,
  const double* tensor_a, const double* tensor_b, const double* tensor_c
);
void computeSteinerPoint(
  const double* point_a, const double* point_b,
  const double* tensor_a, const double* tensor_b,
  double* result_point, double* result_tensor
);
double approxRiemannDist(
  const double* point_a, const double* point_b, const double* tensor
);
double approxRiemannDist(
  const double* point_a, const double* point_b,
  const double* tensor_a, const double* tensor_b
);
void approxRiemannCircum(
  const double* point_a, const double* point_b, const double* point_c,
  const double* matrix, double* circum
);
/* -------------------------------------------------------------------------- */
}} // namespace trinity::numeric

