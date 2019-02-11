/*
 *                          'metric.cpp'
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

#include "metric.h"
#include "mesh.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
Metrics::Metrics(Mesh* input, double targ_f, int norm, double min_h, double max_h)
  : mesh(input),
    solut(input->geom.solut.data()),
    tens(input->geom.tensor.data()),
    nb_nodes(input->nb.nodes),
    nb_elems(input->nb.elems),
    nb_cores(input->nb.cores),
    verbose(input->param.verb),
    iter(input->param.iter),
    rounds(input->param.rounds),
    chunk(input->nb.nodes / input->nb.cores),
    target((int) std::floor(input->nb.cores * targ_f)),
    p_norm(norm),
    h_min(min_h),
    h_max(max_h),
    scale_fact(0.),
    scale_exp(norm ? -1. / (2 * norm + 2) : 1.),
    lambda_min(1. / (min_h * min_h)),
    lambda_max(1. / (max_h * max_h)),
    complexity(0.) {
  stenc = new Patch[nb_nodes];
  nabla = new double[nb_nodes * 2];
}

/* --------------------------------------------------------------------------- */
Metrics::~Metrics() {}

/* --------------------------------------------------------------------------- */
void Metrics::computeTensorField(Stats* tot) {

#pragma omp parallel
  {
    init();

    recoverHessianField();
    normalizeLocally();
    computeComplexity();
    normalizeGlobally();

    recap(tot);
  }
}

/* --------------------------------------------------------------------------- */
void Metrics::clear() {

  delete[] stenc;
  delete[] nabla;
}

/* --------------------------------------------------------------------------- */
void Metrics::recoverHessianField() {

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i) {
    stenc[i] = mesh->getVicinity(i, 2);
    computeGradient(i);
  }

#pragma omp for
  for (int i = 0; i < nb_nodes; ++i)
    computeHessian(i);

#pragma omp master
  {
    mesh->geom.solut.clear();
    mesh->geom.solut.shrink_to_fit();
  }
}

/* --------------------------------------------------------------------------- */
void Metrics::normalizeLocally() {

  // 3) local normalization
  int j, k;
  double s[4], det;
  double p[6], m[9];
  double h[4];
  double val[2], vec[4];
  double scale_loc;

#pragma omp for schedule(guided)
  for (int i = 0; i < nb_nodes; ++i) {
    k = i * 3;
    std::memset(val, 0, sizeof(double) * 2);
    std::memset(vec, 0, sizeof(double) * 4);

    m[0] = tens[k];
    m[1] = tens[k + 1];
    m[2] = m[1];
    m[3] = tens[k + 2];

    // diagonalize
    numeric::eigenDecompose(m, val, vec, vec + 2);
    val[0] = std::max(std::abs(val[0]), EPSILON);
    val[1] = std::max(std::abs(val[1]), EPSILON);

    // modif here: isotropic case
    /*double temp = std::max(val[0],val[1]);
    val[0] = val[1] = temp;*/

    // compute local scale factor w.r.t to L^p, and normalize eigenvalues
    det = val[0] * val[1];
    scale_loc = (p_norm > 0 ? pow(det, scale_exp) : 1.);
    val[0] *= scale_loc;
    val[1] *= scale_loc;

    // h = P.D.P^-1
    m[0] = val[0] * vec[0];
    m[1] = val[0] * vec[2];
    m[2] = val[1] * vec[1];
    m[3] = val[1] * vec[3];

    h[0] = vec[0] * m[0] + vec[1] * m[2];
    h[1] = vec[0] * m[1] + vec[1] * m[3];
    h[2] = vec[2] * m[0] + vec[3] * m[2];
    h[3] = vec[2] * m[1] + vec[3] * m[3];
    assert(std::abs(h[1] - h[2]) < 1.e-3);
    // storeFile back
    tens[k] = h[0];
    tens[k + 1] = h[1];
    tens[k + 2] = h[3];
  }
}

/* --------------------------------------------------------------------------- */
void Metrics::computeComplexity() {

  int j, k;
  double s[4], det;
  double p[6];
  double rho;
  double area;
  double phi = 0.;
  // 4) compute complexity
#pragma omp for schedule(guided) nowait
  for (int i = 0; i < nb_elems; ++i) {
    const int* v = mesh->getElemCoord(i, p);
    //
    s[0] = p[2] - p[0];      // t[1].x - t[0].x
    s[1] = p[4] - p[0];      // t[2].x - t[0].x
    s[2] = p[3] - p[1];      // t[1].y - t[0].y
    s[3] = p[5] - p[1];      // t[2].y - t[0].y
    area = 0.5 * (s[0] * s[3] - s[1] * s[2]);

    // local average aspect ratio
    // rho=(1/3) * sum_k=1^3 det(metric[v[k]])
    rho = 0.;
    for (j = 0; j < 3; ++j) {
      k = 3 * v[j]; // offset
      det = tens[k] * tens[k + 2] - pow(tens[k + 1], 2);
      assert(det > 0);
      rho += sqrt(det);
    }
    rho /= 3;
    phi += rho * area;
  }
#pragma omp critical
  complexity += phi;
#pragma omp barrier
}

/* --------------------------------------------------------------------------- */
void Metrics::normalizeGlobally() {

  assert(complexity);

  int j, k;
  double s[4], det;
  double p[6], m[9];
  double h[4];
  double val[2], vec[4];

#pragma omp single
  scale_fact = target / complexity;

  // 5) global normalization on eigenvalues (only)
#pragma omp for schedule(guided)
  for (int i = 0; i < nb_nodes; ++i) {

    k = i * 3;
    m[0] = tens[k];
    m[1] = tens[k + 1];
    m[2] = m[1];
    m[3] = tens[k + 2];

    numeric::eigenDecompose(m, val, vec, vec + 2);
    val[0] *= scale_fact;
    val[1] *= scale_fact;
    for (j = 0; j < 2; ++j) {
      val[j] = std::min(val[j], lambda_min);
      val[j] = std::max(val[j], lambda_max);
    }

    // modif here: isotropic case
    /*double temp = std::max(val[0],val[1]);
    val[0] = val[1] = temp;*/

    m[0] = val[0] * vec[0];
    m[1] = val[0] * vec[2];
    m[2] = val[1] * vec[1];
    m[3] = val[1] * vec[3];

    // vec * diag(val) * transp(vec)
    h[0] = vec[0] * m[0] + vec[1] * m[2];
    h[1] = vec[0] * m[1] + vec[1] * m[3];
    h[2] = vec[2] * m[0] + vec[3] * m[2];
    h[3] = vec[2] * m[1] + vec[3] * m[3];
    assert(std::abs(h[1] - h[2]) < 1.e-3);

    tens[k] = h[0];
    tens[k + 1] = h[1];
    tens[k + 2] = h[3];
  }
}

/* --------------------------------------------------------------------------- */

void Metrics::init() {
#pragma omp master
  {
    assert(p_norm < 5);
    assert(h_min > 0);
    assert(h_max > 0);

    if (not verbose)
      std::printf("\n\r= Remeshing  ... %3d %% =", 0);

    else if (verbose == 1)
      std::printf("%-18s%s", "= metric field", "...");

    else if (verbose == 2)
      std::printf("Computing multi-scale metric field ... \n");

    std::fflush(stdout);
    start = timer::now();
  }

  int rank = omp_get_thread_num();
  int off = rank * chunk;

  // first-touch
  std::memset(nabla + (off * 2), 0, (chunk * 2) * sizeof(double));
  std::memset( tens + (off * 3), 0, (chunk * 3) * sizeof(double));
}

/* --------------------------------------------------------------------------- */
void Metrics::recap(Stats* tot) {
#pragma omp master
  {
    int end = timer::elapsed_ms(start);

    tot->eval += nb_nodes;
    tot->task += nb_nodes;
    tot->elap += end;

    if (!verbose)
      std::printf("\r= Remeshing  ... %3d %% =", (int) std::floor(100 * (++iter) / (4 * rounds + 1)));

    else if (verbose == 1) {
      std::printf("%10d task/sec \e[32m(%4.2f s)\e[0m\n",
                  (int) std::floor(nb_nodes / (end * 1e-3)), (float) end / 1e3);
    } else if (verbose == 2) {
      std::printf("= norm L^%s\n", p_norm <= 0 ? "inf" : std::to_string(p_norm).data());
      std::printf("= target : %.2e\n", (float) target);
      std::printf("= complex: %.1f\n", (float) std::floor(complexity));
      std::printf("= scale_f: %.1f\n", scale_fact);
      std::printf("done. \e[32m(%d ms)\e[0m\n\n", end);
    }
    std::fflush(stdout);
  }
}
} // namespace trinity