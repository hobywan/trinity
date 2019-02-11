/*
 *                          'hessian.cpp'
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

#include "hessian.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
void Metrics::computeGradient(int index) {

  // init
  double size = 0.;      // stencil area
  double norm = 0.;      // normalization factor
  double p[6];           // elem coords and metric tensor
  double s[4];           // matrix for area computation
  double psi[6];         // finite elem basis function

  size_t N = field.stencil[index].elem.size();

  auto* area = new double[N];
  auto* grad = new double[2 * N];
  std::memset(grad, 0, 2 * N * sizeof(double));

  // compute element-wise gradient vector
  int i = 0;
  for (auto it = field.stencil[index].elem.begin(); it < field.stencil[index].elem.end(); ++it) {
    // elem attributes
    const int* v = mesh->getElemCoord(*it, p);

    // elem area
    s[0] = p[2] - p[0];      // t[1].x - t[0].x
    s[1] = p[4] - p[0];      // t[2].x - t[0].x
    s[2] = p[3] - p[1];      // t[1].y - t[0].y
    s[3] = p[5] - p[1];      // t[2].y - t[0].y
    area[i] = 0.5 * (s[0] * s[3] - s[1] * s[2]);

    assert(area[i]);
    norm = 2 / area[i];
    size += area[i];

    // derivative of basis function
    psi[0] = norm * (p[3] - p[5]);   // t[1].y - t[2].y
    psi[1] = norm * (p[4] - p[2]);   // t[2].x - t[1].x
    psi[2] = norm * (p[5] - p[1]);   // t[2].y - t[0].y
    psi[3] = norm * (p[0] - p[4]);   // t[0].x - t[2].x
    psi[4] = norm * (p[1] - p[3]);   // t[0].y - t[1].y
    psi[5] = norm * (p[2] - p[0]);   // t[1].x - t[0].x

    // elem_grad = sum_i=1^3 (solut[i] . psi[i])
    // constant per elem
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 2; ++k)
        grad[i * 2 + k] += field.solut[*(v + j)] * psi[j * 2 + k];

    ++i;
  }

  int k = index * 2;

  // 3) recover nodal gradient vector
  field.gradient[k] = field.gradient[k + 1] = 0.;
  for (i = 0; i < N; ++i)
    for (int j = 0; j < 2; ++j)
      field.gradient[k + j] += (area[i] * grad[i * 2 + j]);

  // normalize by sizeil area
  field.gradient[k] /= size;
  field.gradient[k + 1] /= size;

  delete[] area;
  delete[] grad;
}

/* --------------------------------------------------------------------------- */
void Metrics::computeHessian(int index) {

  // init
  double size = 0.;     // sizeil area
  double norm = 0.;      // normalization factor
  double p[6];           // elem coords and metric tensor
  double s[4];           // matrix for area computation
  double psi[6];         // finite elem basis function
  double H[4] = {0., 0., 0., 0.};

  size_t N = field.stencil[index].elem.size();

  auto* area = new double[N];
  auto* delta = new double[4 * N];
  std::memset(delta, 0, 4 * N * sizeof(double));

  // compute element-wise hessian matrices
  // nb : derivative of basis functions may be already computed on
  // gradient recovery step (if also computed by L2 projection)
  // but may not if another gradient recovery method was used
  int i = 0;
  for (auto it = field.stencil[index].elem.begin(); it < field.stencil[index].elem.end(); ++it) {
    // elem attributes
    const int* v = mesh->getElemCoord(*it, p);

    // elem area
    s[0] = p[2] - p[0];      // t[1].x - t[0].x
    s[1] = p[4] - p[0];      // t[2].x - t[0].x
    s[2] = p[3] - p[1];      // t[1].y - t[0].y
    s[3] = p[5] - p[1];      // t[2].y - t[0].y
    area[i] = 0.5 * (s[0] * s[3] - s[1] * s[2]);

    // sum to size total area
    assert(area[i]);
    norm = 2 / area[i];
    size += area[i];

    // compute derivative of basis function : nabla_psi[i]
    psi[0] = norm * (p[3] - p[5]); // (t[1].y - t[2].y)
    psi[1] = norm * (p[4] - p[2]); // (t[2].x - t[1].x)
    psi[2] = norm * (p[5] - p[1]); // (t[2].y - t[0].y)
    psi[3] = norm * (p[0] - p[4]); // (t[0].x - t[2].x)
    psi[4] = norm * (p[1] - p[3]); // (t[0].y - t[1].y)
    psi[5] = norm * (p[2] - p[0]); // (t[1].x - t[0].x)

    // hess = sum_i=1^3 (field.gradient[j] . psi[j])
    for (int j = 0; j < 3; ++j) {
      delta[i * 4]     += (field.gradient[*(v + j) * 2]     * psi[j * 2]);
      delta[i * 4 + 1] += (field.gradient[*(v + j) * 2]     * psi[j * 2 + 1]);
      delta[i * 4 + 2] += (field.gradient[*(v + j) * 2 + 1] * psi[j * 2]);
      delta[i * 4 + 3] += (field.gradient[*(v + j) * 2 + 1] * psi[j * 2 + 1]);
    }
    ++i;
  }

  // 3) recover nodal hessian matrix
  for (i = 0; i < N; ++i)
    for (int j = 0; j < 4; ++j)
      H[j] += area[i] * delta[i * 4 + j];
  // normalize by stencil area
  for (int j = 0; j < 4; ++j)
    H[j] /= size;


  int k = index * 3;
  // symmetrize H : 0.5 * (H^t + H)
  field.tensor[k] = H[0];
  field.tensor[k + 1] = 0.5 * (H[1] + H[2]);
  field.tensor[k + 2] = H[3];

  delete[] area;
  delete[] delta;
}


/* ------------------------------------
void least_squares(int  k,
			             double*  solut,
			             Patch* size,
			             Mesh*  mesh,
			             double*  nabla){
  assert(solut not_eq nullptr);
  assert(size not_eq nullptr);
  assert(mesh  not_eq nullptr);
  assert(nabla not_eq nullptr);

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