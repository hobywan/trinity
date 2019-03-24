/*
 *                          'hessian.cpp'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *               Copyright 2016, Hoby Rakotoarivelo.
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

#include "trinity/hessian.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
void Metrics::computeGradient(int index) {

  // init
  double total_area = 0.;      // stencil area
  double point[6];             // elem coords and metric tensor
  double coef[4];              // matrix for area computation
  double deriv[6];               // finite elem basis function

  auto const size = field.stencil[index].elem.size();

  double area[size];
  double grad[size * 2];
  std::fill(grad, grad + (size * 2), 0);

  // compute element-wise gradient vector
  int i = 0;
  for (auto&& elem : field.stencil[index].elem) {
    // elem attributes
    auto const v = mesh->getElemCoord(elem, point);

    // elem area
    coef[0] = point[2] - point[0];      // t[1].x - t[0].x
    coef[1] = point[4] - point[0];      // t[2].x - t[0].x
    coef[2] = point[3] - point[1];      // t[1].y - t[0].y
    coef[3] = point[5] - point[1];      // t[2].y - t[0].y
    area[i] = 0.5 * (coef[0] * coef[3] - coef[1] * coef[2]);

    assert(area[i]);
    double const norm = 2 / area[i];
    total_area += area[i];

    // derivative of basis function
    deriv[0] = norm * (point[3] - point[5]);   // t[1].y - t[2].y
    deriv[1] = norm * (point[4] - point[2]);   // t[2].x - t[1].x
    deriv[2] = norm * (point[5] - point[1]);   // t[2].y - t[0].y
    deriv[3] = norm * (point[0] - point[4]);   // t[0].x - t[2].x
    deriv[4] = norm * (point[1] - point[3]);   // t[0].y - t[1].y
    deriv[5] = norm * (point[2] - point[0]);   // t[1].x - t[0].x

    // elem_grad = sum_i=1^3 (solut[i] . psi[i])
    // constant per elem
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 2; ++k)
        grad[i * 2 + k] += field.solut[*(v + j)] * deriv[j * 2 + k];

    ++i;
  }

  int const k = index * 2;

  // 3) recover nodal gradient vector
  field.gradient[k] = field.gradient[k + 1] = 0.;
  for (i = 0; i < (int) size; ++i)
    for (int j = 0; j < 2; ++j)
      field.gradient[k + j] += (area[i] * grad[i * 2 + j]);

  // normalize by sizeil area
  field.gradient[k] /= total_area;
  field.gradient[k + 1] /= total_area;
}

/* -------------------------------------------------------------------------- */
void Metrics::computeHessian(int index) {

  // init
  double total_area = 0;        // stencil area
  double point[6];              // elem coords and metric tensor
  double coef[4];               // matrix for area computation
  double deriv[6];              // finite elem basis function
  double hessian[4] = {0};

  auto const size = field.stencil[index].elem.size();

  double area[size];
  double delta[size * 4];
  std::fill(delta, delta + (size * 4), 0);

  // compute element-wise hessian matrices
  // nb : derivative of basis functions may be already computed on
  // gradient recovery step (if also computed by L2 projection)
  // but may not if another gradient recovery method was used
  int i = 0;
  for (auto&& elem : field.stencil[index].elem) {
    // elem attributes
    auto const v = mesh->getElemCoord(elem, point);

    // elem area
    coef[0] = point[2] - point[0];      // t[1].x - t[0].x
    coef[1] = point[4] - point[0];      // t[2].x - t[0].x
    coef[2] = point[3] - point[1];      // t[1].y - t[0].y
    coef[3] = point[5] - point[1];      // t[2].y - t[0].y
    area[i] = 0.5 * (coef[0] * coef[3] - coef[1] * coef[2]);

    // sum to size total area
    assert(area[i]);
    double const norm = 2 / area[i];
    total_area += area[i];

    // compute derivative of basis function : nabla_psi[i]
    deriv[0] = norm * (point[3] - point[5]); // (t[1].y - t[2].y)
    deriv[1] = norm * (point[4] - point[2]); // (t[2].x - t[1].x)
    deriv[2] = norm * (point[5] - point[1]); // (t[2].y - t[0].y)
    deriv[3] = norm * (point[0] - point[4]); // (t[0].x - t[2].x)
    deriv[4] = norm * (point[1] - point[3]); // (t[0].y - t[1].y)
    deriv[5] = norm * (point[2] - point[0]); // (t[1].x - t[0].x)

    // hess = sum_i=1^3 (field.gradient[j] . psi[j])
    for (int j = 0; j < 3; ++j) {
      delta[i * 4]     += (field.gradient[*(v + j) * 2]     * deriv[j * 2]);
      delta[i * 4 + 1] += (field.gradient[*(v + j) * 2]     * deriv[j * 2 + 1]);
      delta[i * 4 + 2] += (field.gradient[*(v + j) * 2 + 1] * deriv[j * 2]);
      delta[i * 4 + 3] += (field.gradient[*(v + j) * 2 + 1] * deriv[j * 2 + 1]);
    }
    ++i;
  }

  // 3) recover nodal hessian matrix
  for (i = 0; i < (int) size; ++i)
    for (int j = 0; j < 4; ++j)
      hessian[j] += area[i] * delta[i * 4 + j];
  // normalize by stencil area
  for (auto&& h_coef : hessian)
    h_coef /= total_area;

  int const k = index * 3;
  // symmetrize H : 0.5 * (H^t + H)
  field.tensor[k]     = hessian[0];
  field.tensor[k + 1] = 0.5 * (hessian[1] + hessian[2]);
  field.tensor[k + 2] = hessian[3];
}


/* ------------------------------------
void least_squares(int  k,
			             double*  solut,
			             Patch* size,
			             Mesh*  mesh,
			             double*  nabla){
  assert(solut != nullptr);
  assert(size != nullptr);
  assert(mesh  != nullptr);
  assert(nabla != nullptr);

  int nb_v = size->node.size();
  assert(nb_v>=2);

  Eigen::MatrixXf A(nb_v,2);
  Eigen::VectorXf X(2),
                  B(nb_v);

  int z=0;
  // AX = B : fill matrix A and vector B
  for(int i=0; i < nb_v; ++i){
    z = size->node[i];
    A(i,0) = mesh->points[k*2]   - mesh->points[z*2];
    A(i,1) = mesh->points[k*2+1] - mesh->points[z*2+1];
    B[i]   = solut[k] - solut[z];
  }
  // solve the normal system (A^t.A).X = (A^t.B)
  X = A.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(B);
  // copy values
  // can't direclty std::memcpy without knowing the underlying struct
  nabla[0] = X[0];
  nabla[1] = X[1];
}*/
} // namespace trinity
