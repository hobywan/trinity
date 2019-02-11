/*
 *                          'numeric.h'
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
#include "header.h"
/* --------------------------------------------------------------------------- */
namespace trinity { namespace numeric {
/* --------------------------------------------------------------------------- */
void eigenDecompose(const double* m, double* val, double* vec1, double* vec2);
void interpolateTensor(const double* M, double* R, int n);
void kroneckerProduct(const double* u1, const double* u2, double* M);

double computeQuality(const double* pa, const double* pb, const double* pc,
                      const double* Ma, const double* Mb, const double* Mc);
void computeSteinerPoint(const double* pa, const double* pb,
                         const double* Ma, const double* Mb, double* p, double* M);

double approxRiemannDist(const double* pa, const double* pb, const double* M);
double approxRiemannDist(const double* pa, const double* pb,
                         const double* Ma, const double* Mb);
void approxRiemannCircum(const double* pa, const double* pb, const double* pc,
                         const double* m, double* p);
/* --------------------------------------------------------------------------- */
}} // namespace trinity::numeric

