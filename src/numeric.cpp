/*
 *                          'numeric.cpp'
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

#include "trinity/numeric.h"
#include "trinity/tools.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
void numeric::eigenDecompose(const double matrix[4],
                             double val[2],
                             double vec1[2],
                             double vec2[2]) {

  if (matrix[1] * matrix[2] < EPSILON) {
    val[0]  = matrix[3];
    val[1]  = matrix[0];
    vec1[0] = 0;
    vec2[0] = 1;
    vec1[1] = 1;
    vec2[1] = 0;
    return;
  }

  std::fill( val,  val + 2, 0);
  std::fill(vec1, vec1 + 2, 0);
  std::fill(vec2, vec2 + 2, 0);

  double trace = matrix[0] + matrix[3];
  double det   = matrix[0] * matrix[3] - matrix[1] * matrix[2];
  double term  = trace * 0.5;
  //
  double sqrt_delta = sqrt(term * term - det);
  val[0] = term - sqrt_delta;
  val[1] = term + sqrt_delta;

  if (not std::isnormal(val[0]) or not std::isnormal(val[1])) {
    std::fprintf(stderr,
      "\nm:(%.2f,%.2f,%.2f,%.2f)", matrix[0], matrix[1], matrix[2], matrix[3]
    );
    std::fprintf(stderr,
      "term=%.4f, sqrt_delt=%.4f, val[0]=%.8f, val[1]=%.8f\n",
      term, sqrt_delta, val[0], val[1]
    );
    std::fprintf(stderr,
      "%.5f %.5f %.5f %.5f\n", matrix[0], matrix[1], matrix[2], matrix[3]
    );
    val[0] = val[1];
  }
  assert(std::isnormal(val[0]));
  assert(std::isnormal(val[1]));

  // (!) orthogonal
  term = std::sqrt(
    std::max(0.25 * pow(matrix[0] - matrix[3], 2) + matrix[1] * matrix[2], 0.)
  );
  if (matrix[0] - matrix[3] >= 0) {
    vec1[0] = +0.5 * (matrix[0] - matrix[3]) - term;
    vec1[1] = matrix[1];
    vec2[0] = matrix[2];
    vec2[1] = -0.5 * (matrix[0] - matrix[3]) + term;
  } else {
    vec1[0] = matrix[2];
    vec1[1] = +0.5 * (matrix[0] - matrix[3]) + term;
    vec2[0] = -0.5 * (matrix[0] - matrix[3]) - term;
    vec2[1] = matrix[1];
  }

  double norm[2];
  norm[0] = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1]);
  norm[1] = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1]);
  assert(norm[0]);
  assert(norm[1]);

//#pragma omp simd
  for (int i = 0; i < 2; ++i) {
    vec1[i] /= norm[0];
    vec2[i] /= norm[1];
  }
}

/* -------------------------------------------------------------------------- */
void numeric::interpolateTensor(const double matrix[4], double result[3], int nb) {

  assert(nb);

  double copy[4];
  double cache[4];
  double val[2];
  double vec[4];
  double log_mat[4];

  std::memset( result, 0, sizeof(double) * 3);
  std::memset(log_mat, 0, sizeof(double) * 4);

  // 1) compute weighted sum of tensors logarithm
  // sum_i=1^2 coef[i]*log(M[i]) = sum_i (P[i]. coef[i]*log(diag[i]). P[i]^-1)
  for (int i = 0; i < nb; ++i) {
    // dim=2, so offset=3
    int j = 3 * i;
    copy[0] = matrix[j];
    copy[1] = matrix[j + 1];
    copy[2] = matrix[j + 1];
    copy[3] = matrix[j + 2];

    for (double k : copy)
      assert(std::isfinite(k));

    eigenDecompose(copy, val, vec, vec + 2);

    if (not std::isnormal(val[0]) or not std::isnormal(val[1])) {
      std::fprintf(stderr,
        "(val[0],val[1]) = (%.2f,%.2f)\n", val[0], val[1]
      );
      std::fprintf(stderr,
        "(vec[0],vec[1]) = [(%.2f,%.2f),(,%.2f,%.2f)]\n",
        vec[0], vec[1], vec[2], vec[3]
      );
      std::fprintf(stderr,
        "\nM:(%.2f,%.2f,%.2f)", matrix[j], matrix[j + 1], matrix[j + 2]
      );
    }

    val[0] = std::log(val[0]);
    val[1] = std::log(val[1]);
    assert(std::isfinite(val[0]));
    assert(std::isfinite(val[1]));

    // h = P.D.P^-1
    std::memset(copy, 0, sizeof(double) * 4);
    copy[0] = val[0] * vec[0];
    copy[1] = val[0] * vec[2];
    copy[2] = val[1] * vec[1];
    copy[3] = val[1] * vec[3];

    cache[0] = vec[0] * copy[0] + vec[1] * copy[2];
    cache[1] = vec[0] * copy[1] + vec[1] * copy[3];
    cache[2] = vec[2] * copy[0] + vec[3] * copy[2];
    cache[3] = vec[2] * copy[1] + vec[3] * copy[3];

#if DEBUG
    for(int k=0; k < 4; ++k){
      if(not std::isfinite(log_mat[k])){
         std::fprintf(stderr, "\nmatrix:(%.2f,%.2f,%.2f,%.2f)",
                      matrix[0], matrix[1], matrix[2], matrix[3]);
         std::fprintf(stderr, "\ncopy:(%.2f,%.2f,%.2f)",
                      copy[j], copy[j+1], copy[j+2]);
         std::fprintf(stderr, "\ncache:(%.2f,%.2f,%.2f,%.2f)",
                      cache[0], cache[1], cache[2], cache[3]);
      }
      assert(std::isfinite(log_mat[k]));  // abort
    }
#endif
#pragma omp simd
    for (int k = 0; k < 4; ++k)
      log_mat[k] += cache[k];
  }

#if DEBUG
  for(int k=0; k < 4; ++k){
    if(not std::isfinite(log_mat[k])){
      std::fprintf(stderr, "\nmatrix:(%.2f,%.2f,%.2f,%.2f)",
                   matrix[0], matrix[1], matrix[2], matrix[3]);
      std::fprintf(stderr, "\ncache:(%.2f,%.2f,%.2f,%.2f)",
                   cache[0], cache[1], cache[2], cache[3]);
      std::fprintf(stderr, "(val[0],val[1]) = (%.2f,%.2f)\n", val[0], val[1]);
      std::fprintf(stderr, "(vec[0],vec[1]) = [(%.2f,%.2f),(,%.2f,%.2f)]\n",
                   vec[0],vec[1],vec[2],vec[3]);
    }
    assert(std::isfinite(log_mat[k]));  // abort
  }
#endif

  // normalize
  // EDIT use coef[i]*val[i] instead
#pragma omp simd
  for (int i = 0; i < 4; ++i)
    log_mat[i] /= nb;

  // exponentiate 'log_m' to get the average tensor
  // (!) 'log_m' is not symmetric
  eigenDecompose(log_mat, val, vec, vec + 2);
  val[0] = std::exp(val[0]);
  val[1] = std::exp(val[1]);

  std::memset( copy, 0, sizeof(double) * 4);
  std::memset(cache, 0, sizeof(double) * 4);
  copy[0] = val[0] * vec[0];
  copy[1] = val[0] * vec[2];
  copy[2] = val[1] * vec[1];
  copy[3] = val[1] * vec[3];

  cache[0] = vec[0] * copy[0] + vec[1] * copy[2];
  cache[1] = vec[0] * copy[1] + vec[1] * copy[3];
  cache[2] = vec[2] * copy[0] + vec[3] * copy[2];
  cache[3] = vec[2] * copy[1] + vec[3] * copy[3];

  assert(std::abs(cache[1] - cache[2]) < 1.e-3);
  result[0] = cache[0];
  result[1] = cache[1];
  result[2] = cache[3];

#if DEBUG
  if(not std::isnormal(result[0]) or
     not std::isnormal(result[1]) or
     not std::isnormal(result[2])) {
#pragma omp critical
    {
      std::fprintf(stderr, "result:(%.2f,%.2f,%.2f)\n",
                   result[0], result[1], result[2]);
      for(int i=0; i < nb; ++i) {
        std::fprintf(stderr, "(%.2f,%.2f,%.2f), ",
                     cache[i * 3], cache[i * 3 + 1], cache[i * 3 + 2]);
      }
      std::printf("\nlog_mat:(%.2f,%.2f,%.2f,%.2f)",
                  log_mat[0], log_mat[1], log_mat[2], log_mat[3]);
      std::fprintf(stderr, "\nmatrix:(%.2f,%.2f,%.2f,%.2f)",
                   matrix[0], matrix[1], matrix[2], matrix[3]);
      std::fprintf(stderr, "\ncopy:(%.2f,%.2f,%.2f)",
                   copy[0], copy[1], copy[2]);
      std::fprintf(stderr, "\ncache:(%.2f,%.2f,%.2f,%.2f)",
                   cache[0], cache[1], cache[2], cache[3]);
      std::fprintf(stderr, "(val[0],val[1]) = (%.2f,%.2f)\n", val[0], val[1]);
      std::fprintf(stderr, "(vec[0],vec[1]) = [(%.2f,%.2f),(,%.2f,%.2f)]\n",
                   vec[0], vec[1], vec[2], vec[3]);
      tools::separator();
    }
  }
#endif
  assert(std::isfinite(result[0]));
  assert(std::isfinite(result[1]));
  assert(std::isfinite(result[2]));
}

/* -------------------------------------------------------------------------- */
double numeric::computeQuality(const double point_a[2],
                               const double point_b[2],
                               const double point_c[2],
                               const double tensor_a[3],
                               const double tensor_b[3],
                               const double tensor_c[3]) {

  double coef[6];
  double vec[6];
  double jacob[3];
  double perim = 0;
  double area;
  double matrix[9];

  // average metric tensor
  // m = ((ma^-1 + mb^-1 + mc^-1)^-1)/3
  std::memcpy(matrix, tensor_a, sizeof(double) * 3);
  std::memcpy(matrix + 3, tensor_b, sizeof(double) * 3);
  std::memcpy(matrix + 6, tensor_c, sizeof(double) * 3);

  for (auto&& k : matrix)
    assert(std::isfinite(k));

  interpolateTensor(matrix, jacob, 3);
  double det = jacob[0] * jacob[2] - jacob[1] * jacob[1];
  assert(det);

  // perim. in riemannian space
  coef[0] = point_b[0] - point_a[0];
  coef[1] = point_b[1] - point_a[1];
  coef[2] = point_c[0] - point_b[0];
  coef[3] = point_c[1] - point_b[1];
  coef[4] = point_a[0] - point_c[0];
  coef[5] = point_a[1] - point_c[1];

  for (int k = 0; k < 3; ++k) {
    int i = k * 2,
        j = i + 1;
    vec[i] = jacob[0] * coef[i] + jacob[1] * coef[j];
    vec[j] = jacob[1] * coef[i] + jacob[2] * coef[j];
    perim += sqrt(coef[i] * vec[i] + coef[j] * vec[j]);
  }
  assert(perim);

  // euclidian area
  coef[0] = point_b[0] - point_a[0];
  coef[1] = point_c[0] - point_a[0];
  coef[2] = point_b[1] - point_a[1];
  coef[3] = point_c[1] - point_a[1];
  area = 0.5 * std::abs(coef[0] * coef[3] - coef[1] * coef[2]);

  // normalization
  auto const mean = std::min(perim / 3, 3 / perim);
  auto const norm = 12 * sqrt(3) * pow(mean * (2. - mean), 3);
  return area * norm * sqrt(det) / pow(perim, 2);
}

/* ---------------------------------------------------------------------------
 * riemannian triangle circumcenter w.r.t metric M
 *
 * |x12  y12| |Ox|   |d1|          |m11  m12|
 * |x13  y13| |Oy| = |d2| and Mp = |m12  m22|
 *
 * ai = 2(m11 (xi - x1) + m12 (yi - y1)
 * bi = 2(m22 (yi - y1) + m12 (xi - x1)
 * di = m11 xi² + 2 m12 (xiyi) + m22 yi²
 *   - (m11 x1² + 2 m12 (x1y1) + m22 y1²)
 * cf. Bourouchaki, Georges, Frey
 * -------------------------------------------------------------------------- */
void numeric::approxRiemannCircum(const double point_a[2],
                                  const double point_b[2],
                                  const double point_c[2],
                                  const double matrix[3],
                                  double circum[2]) {
  double coef[3];
  double x[2];
  double y[2];

  x[0] = 2 * (matrix[0] * (point_b[0] - point_a[0]) +
              matrix[1] * (point_b[1] - point_a[1]));  // x12
  x[1] = 2 * (matrix[0] * (point_c[0] - point_a[0]) +
              matrix[1] * (point_c[1] - point_a[1]));  // x13
  y[0] = 2 * (matrix[2] * (point_b[1] - point_a[1]) +
              matrix[1] * (point_b[0] - point_a[0]));  // y12
  y[1] = 2 * (matrix[2] * (point_c[1] - point_a[1]) +
              matrix[1] * (point_c[0] - point_a[0]));  // y13
  double const det = x[0] * y[1] - x[1] * y[0];

  if (std::abs(det) < EPSILON) {
    std::fprintf(stderr, "Failure on anisotropic circumcenter computation");
#if DEBUG
    std::printf("|m11  m12| : |%.5f\t%.5f|\n", matrix[0], matrix[1]);
    std::printf("|m12  m22| : |%.5f\t%.5f|\n", matrix[1], matrix[2]);
    std::printf(" ---\n");
    std::printf("|x12  y12| : |%.5f\t%.5f|\n", x[0], y[0]);
    std::printf("|x13  y13| : |%.5f\t%.5f|\n", x[1], y[1]);
    std::printf("det = %.5f\n", det);
#endif
    exit(EXIT_FAILURE);
  }

  coef[0] = matrix[0] * point_a[0] * point_a[0] +
            matrix[2] * point_a[1] * point_a[1] +
            matrix[1] * point_a[0] * point_a[1] * 2;
  coef[1] = matrix[0] * point_b[0] * point_b[0] +
            matrix[2] * point_b[1] * point_b[1] +
            matrix[1] * point_b[0] * point_b[1] * 2;
  coef[2] = matrix[0] * point_c[0] * point_c[0] +
            matrix[2] * point_c[1] * point_c[1] +
            matrix[1] * point_c[0] * point_c[1] * 2;

  coef[1] = coef[1] - coef[0];
  coef[2] = coef[2] - coef[0];

  // set coords
  circum[0] = (coef[1] * y[1] - coef[2] * y[0]) / det;
  circum[1] = (coef[2] * x[0] - coef[1] * x[1]) / det;
}

/* -------------------------------------------------------------------------- */
void numeric::doKroneckerProduct(const double vec1[2],
                                 const double vec2[2],
                                 double matrix[4]) {

  assert(vec1 != nullptr);
  assert(vec2 != nullptr);
  assert(matrix != nullptr);

  matrix[0] = 0.25 * std::pow(vec1[0] + vec2[0], 2);
  matrix[1] = 0.25 * (vec1[0] + vec2[0]) * (vec1[1] * vec2[1]);
  matrix[2] = 0.25 * (vec1[0] + vec1[1]) * (vec2[0] * vec2[1]);
  matrix[3] = 0.25 * std::pow(vec1[1] + vec2[1], 2);

}

/* -------------------------------------------------------------------------- */
void numeric::computeSteinerPoint(const double point_a[2],
                                  const double point_b[2],
                                  const double tensor_a[3],
                                  const double tensor_b[3],
                                  double result_point[2],
                                  double result_tensor[3]) {
/* midpoint of a segment parametrized by e : = a + t ab.
 * parameter 't' is given by the system :
 *   | len[a] (t * ab) = (1 - t) * len[b]
 *   | len[a] (t * ab) =   0.5 * len[bar]
 */
  double length[2] = {0};
  double metric[6] = {0};
  double weight = 0;
  double coef = 1;   // or nb-2

  length[0] = approxRiemannDist(point_a, point_b, tensor_a);
  length[1] = approxRiemannDist(point_a, point_b, tensor_b);

  std::memset(result_point, 0, sizeof(double) * 2);
  // 1) parametrize
  weight = (std::abs(length[0] - length[1]) < EPSILON
    ? 0.5
    : (sqrt(coef) * sqrt(length[0] * length[1]) - length[1]) /
                              (coef * length[0] - length[1])
  );

  // 2) coordinates
  result_point[0] = (1 - weight) * point_a[0] + (weight * point_b[0]);
  result_point[1] = (1 - weight) * point_a[1] + (weight * point_b[1]);
#if DEBUG
  if(std::abs(result_point[0] - 0) < EPSILON){
     std::printf("warning: ");
     std::printf("t=%.5f, pa.x=%.5f, pa.y=%.5f, pb.x=%.5f, pb.y=%.5f\n",
                  weight, point_a[0], point_a[1], point_b[0], point_b[1]);
  }
#endif
  assert(std::isfinite(result_point[0]));
  assert(std::isfinite(result_point[1]));

  // 3) interpolate tensor
  std::memcpy(metric    , tensor_a, sizeof(double) * 3);
  std::memcpy(metric + 3, tensor_b, sizeof(double) * 3);
  interpolateTensor(metric, result_tensor, 2);

#if DEBUG
  if(not std::isnormal(metric[0]) or
     not std::isnormal(metric[1]) or
     not std::isnormal(metric[2])) {
     std::printf("ma:(%.2f,%.2f,%.2f)\n", tensor_a[0], tensor_a[1], tensor_a[2]);
     std::printf("mb:(%.2f,%.2f,%.2f)\n", tensor_b[0], tensor_b[1], tensor_b[2]);
     std::printf("m :(%.2f,%.2f,%.2f)\n", metric[0], metric[1], metric[2]);
  }
#endif

  assert(std::isnormal(result_tensor[0]));
  //assert(std::isnormal(steiner_tensor[0]));
  assert(std::isnormal(result_tensor[2]));
}

/* -------------------------------------------------------------------------- */
// len = sqrt(u M u)
// nb : M is an upper triangular matrix
double numeric::approxRiemannDist(const double point_a[2],
                                  const double point_b[2],
                                  const double tensor[3]) {

  double u[2], v[2];
  u[0] = point_b[0] - point_a[0];
  u[1] = point_b[1] - point_a[1];
  // m * u
  v[0] = tensor[0] * u[0] + tensor[1] * u[1];
  v[1] = tensor[1] * u[0] + tensor[2] * u[1];
  return sqrt(u[0] * v[0] + u[1] * v[1]);
}

/* ------------------------------------
 * log_euclidian segment length.
 * r = l_max/l_min
 * - linear    : l_max * ln(r) / (r-1)
 * - geometric : l_max * (r-1) / r * ln(r)
 */
double numeric::approxRiemannDist(const double point_a[2],
                                  const double point_b[2],
                                  const double tensor_a[3],
                                  const double tensor_b[3]) {
  double length[4];
  double vector[2];
  double matrix[4];

  // lengths w.r.t each metric
  vector[0] = point_b[0] - point_a[0];
  vector[1] = point_b[1] - point_a[1];

  matrix[0] = tensor_a[0] * vector[0] + tensor_a[1] * vector[1];
  matrix[1] = tensor_a[1] * vector[0] + tensor_a[2] * vector[1];
  matrix[2] = tensor_b[0] * vector[0] + tensor_b[1] * vector[1];
  matrix[3] = tensor_b[1] * vector[0] + tensor_b[2] * vector[1];

  length[0] = sqrt(vector[0] * matrix[0] + vector[1] * matrix[1]);
  length[1] = sqrt(vector[0] * matrix[2] + vector[1] * matrix[3]);

  if (std::abs(length[0] - length[1]) < EPSILON) { return length[0]; }

  auto const l_min = std::min(length[0], length[1]);
  auto const l_max = std::max(length[0], length[1]);
  auto const ratio = l_max / l_min;

  // geometric interpolation
  return l_max * (ratio - 1) / (ratio * std::log(ratio));
}
/* -------------------------------------------------------------------------- */
}
