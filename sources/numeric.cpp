/* ------------------------------------ */
#include <numeric.h>
#include <tools.h>

namespace trinity {
/* ------------------------------------ */
void numeric::eigenDecomposeTensor(const double m[4],
                                   double val[2],
                                   double vec1[2],
                                   double vec2[2]) {

  if (m[1] * m[2] < EPSILON) {
    val[0] = m[3];
    val[1] = m[0];
    vec1[0] = 0;
    vec2[0] = 1;
    vec1[1] = 1;
    vec2[1] = 0;
    return;
  }

  std::memset(val, 0, sizeof(double) * 2);
  std::memset(vec1, 0, sizeof(double) * 2);
  std::memset(vec2, 0, sizeof(double) * 2);

  double trace = m[0] + m[3];
  double det = m[0] * m[3] - m[1] * m[2];
  double term = trace * 0.5;
  //
  double sqrt_delt = sqrt(term * term - det);
  val[0] = term - sqrt_delt;
  val[1] = term + sqrt_delt;

  if (!std::isnormal(val[0]) or !std::isnormal(val[1])) {
    std::printf("\nm:(%.2f,%.2f,%.2f,%.2f)", m[0], m[1], m[2], m[3]);
    std::printf("term=%.4f, sqrt_delt=%.4f, val[0]=%.8f, val[1]=%.8f\n",
                term, sqrt_delt, val[0], val[1]);
    std::printf("%.5f %.5f %.5f %.5f\n", m[0], m[1], m[2], m[3]);
    val[0] = val[1];
  }
  assert(std::isnormal(val[0]));
  assert(std::isnormal(val[1]));

  // (!) orthogonal
  term = sqrt(std::max(0.25 * pow(m[0] - m[3], 2) + m[1] * m[2], 0.));
  if (m[0] - m[3] >= 0) {
    vec1[0] = +0.5 * (m[0] - m[3]) - term;
    vec1[1] = m[1];
    vec2[0] = m[2];
    vec2[1] = -0.5 * (m[0] - m[3]) + term;
  } else {
    vec1[0] = m[2];
    vec1[1] = +0.5 * (m[0] - m[3]) + term;
    vec2[0] = -0.5 * (m[0] - m[3]) - term;
    vec2[1] = m[1];
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

/* ------------------------------------ */
void numeric::interpolateTensor(const double* M, double R[3], int nb) {

  assert(nb);

  double m[4];
  double h[4];
  double val[2];
  double vec[4];
  double log_m[4];

  std::memset(R, 0, sizeof(double) * 3);
  std::memset(log_m, 0, sizeof(double) * 4);

  // 1) compute weighted sum of tensors logarithm
  // sum_i=1^2 coef[i]*log(M[i]) = sum_i=1^2 (P[i]. coef[i]*log(diag[i]). P[i]^-1)
  for (int i = 0; i < nb; ++i) {
    // dim=2, so offset=3
    int j = 3 * i;
    m[0] = M[j];
    m[1] = M[j + 1];
    m[2] = M[j + 1];
    m[3] = M[j + 2];

    for (double k : m)
      assert(std::isfinite(k));

    eigenDecomposeTensor(m, val, vec, vec + 2);

    if (not std::isnormal(val[0]) or not std::isnormal(val[1])) {
      std::printf("(val[0],val[1]) = (%.2f,%.2f)\n", val[0], val[1]);
      std::printf("(vec[0],vec[1]) = [(%.2f,%.2f),(,%.2f,%.2f)]\n", vec[0], vec[1], vec[2], vec[3]);
      std::printf("\nM:(%.2f,%.2f,%.2f)", M[j], M[j + 1], M[j + 2]);
    }

    val[0] = std::log(val[0]);
    val[1] = std::log(val[1]);
    assert(std::isfinite(val[0]));
    assert(std::isfinite(val[1]));

    // h = P.D.P^-1
    std::memset(m, 0, sizeof(double) * 4);
    m[0] = val[0] * vec[0];
    m[1] = val[0] * vec[2];
    m[2] = val[1] * vec[1];
    m[3] = val[1] * vec[3];

    h[0] = vec[0] * m[0] + vec[1] * m[2];
    h[1] = vec[0] * m[1] + vec[1] * m[3];
    h[2] = vec[2] * m[0] + vec[3] * m[2];
    h[3] = vec[2] * m[1] + vec[3] * m[3];

#ifdef DEBUG
    for(int k=0; k < 4; ++k){
      if(!std::isfinite(log_m[k])){
         std::printf("\nm:(%.2f,%.2f,%.2f,%.2f)", m[0], m[1], m[2], m[3]);
         std::printf("\nM:(%.2f,%.2f,%.2f)", M[j], M[j+1], M[j+2]);
         std::printf("\nh:(%.2f,%.2f,%.2f,%.2f)", h[0], h[1], h[2], h[3]);
      }
      assert(std::isfinite(log_m[k]));  // abort
    }
#endif
#pragma omp simd
    for (int k = 0; k < 4; ++k)
      log_m[k] += h[k];
  }

#ifdef DEBUG
  for(int k=0; k < 4; ++k){
    if(!std::isfinite(log_m[k])){
       std::printf("\nm:(%.2f,%.2f,%.2f,%.2f)", m[0], m[1], m[2], m[3]);
       std::printf("\nh:(%.2f,%.2f,%.2f,%.2f)", h[0], h[1], h[2], h[3]);
       std::printf("(val[0],val[1]) = (%.2f,%.2f)\n", val[0],val[1]);
       std::printf("(vec[0],vec[1]) = [(%.2f,%.2f),(,%.2f,%.2f)]\n", vec[0],vec[1],vec[2],vec[3]);
    }
    assert(std::isfinite(log_m[k]));  // abort
  }
#endif

  // normalize
  // EDIT use coef[i]*val[i] instead
#pragma omp simd
  for (int i = 0; i < 4; ++i)
    log_m[i] /= nb;

  // exponentiate 'log_m' to get the average tensor
  // (!) 'log_m' is not symmetric
  eigenDecomposeTensor(log_m, val, vec, vec + 2);
  val[0] = std::exp(val[0]);
  val[1] = std::exp(val[1]);

  std::memset(m, 0, sizeof(double) * 4);
  std::memset(h, 0, sizeof(double) * 4);
  m[0] = val[0] * vec[0];
  m[1] = val[0] * vec[2];
  m[2] = val[1] * vec[1];
  m[3] = val[1] * vec[3];

  h[0] = vec[0] * m[0] + vec[1] * m[2];
  h[1] = vec[0] * m[1] + vec[1] * m[3];
  h[2] = vec[2] * m[0] + vec[3] * m[2];
  h[3] = vec[2] * m[1] + vec[3] * m[3];

  assert(std::abs(h[1] - h[2]) < 1.e-3);
  R[0] = h[0];
  R[1] = h[1];
  R[2] = h[3];  // error here

  assert(std::abs(R[0] - h[0]) < 1.e-3);
  assert(std::abs(R[1] - h[1]) < 1.e-3);
  assert(std::abs(R[2] - h[3]) < 1.e-3);

#ifdef DEBUG
  if(!std::isnormal(R[0]) or !std::isnormal(R[1]) or !std::isnormal(R[2])){
#pragma omp critical
    {
     std::printf("R:(%d,%d,%d)\n", R[0],R[1],R[2]);
    for(int i=0; i < nb; ++i)
       std::printf("(%.2f,%.2f,%.2f), ", M[i*3], M[i*3+1], M[i*3+2]);
     std::printf("\nlog_m:(%.2f,%.2f,%.2f,%.2f)", log_m[0], log_m[1], log_m[2], log_m[3]);
     std::printf("\nm:(%.2f,%.2f,%.2f,%.2f)", m[0], m[1], m[2], m[3]);
     std::printf("\nh:(%.2f,%.2f,%.2f,%.2f)", h[0], h[1], h[2], h[3]);
     std::printf("\nR:(%.2f,%.2f,%.2f)\n", R[0], R[1], R[2]);
     std::printf("(val[0],val[1]) = (%.2f,%.2f)\n", val[0],val[1]);
     std::printf("(vec[0],vec[1]) = [(%.2f,%.2f),(,%.2f,%.2f)]\n", vec[0],vec[1],vec[2],vec[3]);
    tools::separator();
    }
  }
#endif
  assert(std::isfinite(R[0]));
  assert(std::isfinite(R[1]));
  assert(std::isfinite(R[2]));
}

/* ------------------------------------ */
double numeric::computeQuality(const double pa[2], const double pb[2], const double pc[2],
                               const double Ma[3], const double Mb[3], const double Mc[3]) {

  double s[4],
    u[6],
    v[6],
    m[3],
    det,
    perim,
    area,
    z,
    phi,
    mat[9];

  det = perim = area = z = phi = 0.;

  // average metric tensor
  // m = ((ma^-1 + mb^-1 + mc^-1)^-1)/3
  std::memcpy(mat, Ma, sizeof(double) * 3);
  std::memcpy(mat + 3, Mb, sizeof(double) * 3);
  std::memcpy(mat + 6, Mc, sizeof(double) * 3);

  for (double k : mat)
    assert(std::isfinite(k));

  interpolateTensor(mat, m, 3);
  det = m[0] * m[2] - m[1] * m[1];
  assert(det);

  // perim. in riemannian space
  u[0] = pb[0] - pa[0];
  u[1] = pb[1] - pa[1];
  u[2] = pc[0] - pb[0];
  u[3] = pc[1] - pb[1];
  u[4] = pa[0] - pc[0];
  u[5] = pa[1] - pc[1];

  for (int k = 0; k < 3; ++k) {
    int i = k * 2,
      j = i + 1;
    v[i] = m[0] * u[i] + m[1] * u[j];
    v[j] = m[1] * u[i] + m[2] * u[j];
    perim += sqrt(u[i] * v[i] + u[j] * v[j]);
  }
  assert(perim);

  // euclidian area
  s[0] = pb[0] - pa[0];
  s[1] = pc[0] - pa[0];
  s[2] = pb[1] - pa[1];
  s[3] = pc[1] - pa[1];
  area = 0.5 * std::abs(s[0] * s[3] - s[1] * s[2]);

  // normalization
  z = std::min(perim / 3., 3. / perim);
  phi = pow(z * (2. - z), 3);
  assert(phi <= 1.);
  return 12. * sqrt(3) * area * phi * sqrt(det) / pow(perim, 2);
}

/* ------------------------------------
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
 */
void numeric::approxRiemannCircum(const double* pa,
                                  const double* pb,
                                  const double* pc,
                                  const double* m,
                                  double* p) {
  double det,
    d[3],
    x[2],
    y[2];

  x[0] = 2 * (m[0] * (pb[0] - pa[0]) + m[1] * (pb[1] - pa[1]));  // x12
  x[1] = 2 * (m[0] * (pc[0] - pa[0]) + m[1] * (pc[1] - pa[1]));  // x13
  y[0] = 2 * (m[2] * (pb[1] - pa[1]) + m[1] * (pb[0] - pa[0]));  // y12
  y[1] = 2 * (m[2] * (pc[1] - pa[1]) + m[1] * (pc[0] - pa[0]));  // y13
  det = x[0] * y[1] - x[1] * y[0];

  if (std::abs(det) < EPSILON) {
    perror("aniso::riemannian_circumcenter");
#ifdef DEBUG
    std::printf("|m11  m12| : |%.5f\t%.5f|\n", m[0], m[1]);
    std::printf("|m12  m22| : |%.5f\t%.5f|\n", m[1], m[2]);
    std::printf(" ---\n");
    std::printf("|x12  y12| : |%.5f\t%.5f|\n", x[0], y[0]);
    std::printf("|x13  y13| : |%.5f\t%.5f|\n", x[1], y[1]);
    std::printf("det = %.5f\n", det);
#endif
    exit(EXIT_FAILURE);
  }

  d[0] = m[0] * pa[0] * pa[0] + m[2] * pa[1] * pa[1] + (2 * m[1] * pa[0] * pa[1]);
  d[1] = m[0] * pb[0] * pb[0] + m[2] * pb[1] * pb[1] + (2 * m[1] * pb[0] * pb[1]);
  d[2] = m[0] * pc[0] * pc[0] + m[2] * pc[1] * pc[1] + (2 * m[1] * pc[0] * pc[1]);
  d[1] = d[1] - d[0];
  d[2] = d[2] - d[0];
  // set coords
  p[0] = (d[1] * y[1] - d[2] * y[0]) / det;
  p[1] = (d[2] * x[0] - d[1] * x[1]) / det;

}

/* ------------------------------------ */
void numeric::kroneckerProduct(const double* u1, const double* u2, double* M) {

  assert(u1 not_eq nullptr);
  assert(u2 not_eq nullptr);
  assert(M not_eq nullptr);

  M[0] = 0.25 * pow(u1[0] + u2[0], 2);
  M[1] = 0.25 * (u1[0] + u2[0]) * (u1[1] * u2[1]);
  M[2] = 0.25 * (u1[0] + u1[1]) * (u2[0] * u2[1]);
  M[3] = 0.25 * pow(u1[1] + u2[1], 2);

}

/* ------------------------------------ */
void numeric::computeSteinerPoint(const double pa[2], const double pb[2],
                                  const double Ma[3], const double Mb[3],
                                  double p[2], double M[3]) {
/* midpoint of a segment parametrized by e : = a + t ab.
 * parameter 't' is given by the system :
 *   | len[a] (t * ab) = (1 - t) * len[b]
 *   | len[a] (t * ab) =   0.5 * len[bar]
 */
  double c,
    t,
    l[2],
    met[6];

  t = 0.;
  c = 1.; // or nb-2

  l[0] = approxRiemannDist(pa, pb, Ma);
  l[1] = approxRiemannDist(pa, pb, Mb);

  std::memset(p, 0, sizeof(double) * 2);
  // 1) parametrize
  t = (std::abs(l[0] - l[1]) < EPSILON ? 0.5
                                       : (sqrt(c) * sqrt(l[0] * l[1]) - l[1]) / (c * l[0] - l[1]));

  // 2) coordinates
  p[0] = (1. - t) * pa[0] + (t * pb[0]);
  p[1] = (1. - t) * pa[1] + (t * pb[1]);
#ifdef DEBUG
  if(std::abs(p[0]-0.) < EPSILON){
     std::printf("warning: ");
     std::printf("t=%.5f, pa.x=%.5f, pa.y=%.5f, pb.x=%.5f, pb.y=%.5f\n",
            t, pa[0], pa[1],pb[0], pb[1]);
  }
#endif
  assert(std::isfinite(p[0]));
  assert(std::isfinite(p[1]));

  // 3) interpolate tensor
  std::memcpy(met, Ma, sizeof(double) * 3);
  std::memcpy(met + 3, Mb, sizeof(double) * 3);
  interpolateTensor(met, M, 2);

#ifdef DEBUG
  if(!std::isnormal(m[0]) or !std::isnormal(m[1]) or !std::isnormal(m[2])){
     std::printf("ma:(%.2f,%.2f,%.2f)\n", ma[0], ma[1], ma[2]);
     std::printf("mb:(%.2f,%.2f,%.2f)\n", mb[0], mb[1], mb[2]);
     std::printf("m :(%.2f,%.2f,%.2f)\n",  m[0],  m[1],  m[2]);
  }
#endif

  assert(std::isnormal(M[0]));
//  assert(std::isnormal(m[1]));
  assert(std::isnormal(M[2]));
}

/* ------------------------------------ */
// len = sqrt(u M u)
// nb : M is an upper triangular matrix
double numeric::approxRiemannDist(const double pa[2], const double pb[2], const double m[3]) {

  double u[2],
    v[2];

  u[0] = pb[0] - pa[0];
  u[1] = pb[1] - pa[1];
  // m * u
  v[0] = m[0] * u[0] + m[1] * u[1];
  v[1] = m[1] * u[0] + m[2] * u[1];
  return sqrt(u[0] * v[0] + u[1] * v[1]);
}

/* ------------------------------------
 * log_euclidian segment length.
 * r = l_max/l_min
 * - linear    : l_max * ln(r) / (r-1)
 * - geometric : l_max * (r-1) / r * ln(r)
 */
double numeric::approxRiemannDist(const double pa[2], const double pb[2],
                                  const double ma[3], const double mb[3]) {

  double l_min,
    l_max,
    r,
    l[4],
    u[2],
    v[4];

  // lengths w.r.t each metric
  u[0] = pb[0] - pa[0];
  u[1] = pb[1] - pa[1];

  v[0] = ma[0] * u[0] + ma[1] * u[1];
  v[1] = ma[1] * u[0] + ma[2] * u[1];
  v[2] = mb[0] * u[0] + mb[1] * u[1];
  v[3] = mb[1] * u[0] + mb[2] * u[1];

  l[0] = sqrt(u[0] * v[0] + u[1] * v[1]);
  l[1] = sqrt(u[0] * v[2] + u[1] * v[3]);

  if (std::abs(l[0] - l[1]) < EPSILON)
    return l[0];

  l_min = std::min(l[0], l[1]);
  l_max = std::max(l[0], l[1]);
  r = l_max / l_min;

  // geometric interpolation
  return l_max * (r - 1) / (r * std::log(r));
}

}