#include <R.h>
#include <Rmath.h>
#include <stdio.h>

// Proportional McCOIL with strain proportions shared across loci: each sample
// has M strains with proportions w ~ Dirichlet(1), each strain carries one
// allele per locus (prior Bernoulli(P_j)), and the true allele-1 fraction at a
// locus is sum_m w_m * a_mj.

#define MAX_ENUM 10  // COI moves sum over all 2^M allele subsets up to this M

static int g_k;
static double *g_A1, *g_A2, *g_P, *g_S, *g_L, *g_sum, *g_t;
static int *g_cnt;
static char *g_a;
static double g_err, g_rho, g_pout, g_tau, g_minreads;

static double logsumexp2(double x, double y) {
  double m = x > y ? x : y;
  return m + log(exp(x - m) + exp(y - m));
}

// Observed allele-1 fraction ~ Normal(S, binomial + rho * S(1-S)), mixed with a
// uniform outlier. Minor alleles below the detection limit are zeroed
// upstream, so a cell with one allele absent is censored at that limit.
static double llobs(double a1, double a2, double S) {
  if (a1 < 0 || a2 < 0 || a1 + a2 == 0) return 0.0;
  if (S < 0) S = 0;
  if (S > 1) S = 1;
  double nn = a1 + a2, p = S * (1 - g_err) + (1 - S) * g_err;
  double sd = sqrt(p * (1 - p) / nn + g_rho * S * (1 - S));
  if (a1 == 0 || a2 == 0) {
    double tau = g_tau > g_minreads / nn ? g_tau : g_minreads / nn;
    double pabs = a1 == 0 ? p : 1 - p;
    return logsumexp2(log(1 - g_pout) + pnorm(tau, pabs, sd, 1, 1),
                      log(g_pout) + log(tau));
  }
  return logsumexp2(log(1 - g_pout) + dnorm(a1 / nn, p, sd, 1), log(g_pout));
}

#define W(i, m) w[(size_t)(i) * maxM + (m)]
#define AL(i, m, j) g_a[((size_t)(i) * maxM + (m)) * g_k + (j)]
#define D1(i, j) g_A1[(size_t)(i) * g_k + (j)]
#define D2(i, j) g_A2[(size_t)(i) * g_k + (j)]

// Log likelihood of sample i with alleles summed out (M <= MAX_ENUM). With
// draw set, alleles are sampled from their conditional and S, L refreshed.
static double marg(int i, int M, const double *wi, int draw, int maxM) {
  int ns = 1 << M, s, j, m;
  g_sum[0] = 0;
  g_cnt[0] = 0;
  for (s = 1; s < ns; s++) {
    int low = __builtin_ctz(s);
    g_sum[s] = g_sum[s & (s - 1)] + wi[low];
    g_cnt[s] = g_cnt[s & (s - 1)] + 1;
  }
  double tot = 0;
  for (j = 0; j < g_k; j++) {
    double a1 = D1(i, j), a2 = D2(i, j), lP = log(g_P[j]), lQ = log(1 - g_P[j]);
    double mx = -1e300, z = 0;
    for (s = 0; s < ns; s++) {
      g_t[s] = g_cnt[s] * lP + (M - g_cnt[s]) * lQ + llobs(a1, a2, g_sum[s]);
      if (g_t[s] > mx) mx = g_t[s];
    }
    for (s = 0; s < ns; s++) z += exp(g_t[s] - mx);
    tot += mx + log(z);
    if (draw) {
      double r = runif(0, z), cum = 0;
      for (s = 0; s < ns - 1; s++) {
        cum += exp(g_t[s] - mx);
        if (cum >= r) break;
      }
      for (m = 0; m < M; m++) AL(i, m, j) = (char)((s >> m) & 1);
      g_S[(size_t)i * g_k + j] = g_sum[s];
      g_L[(size_t)i * g_k + j] = llobs(a1, a2, g_sum[s]);
    }
  }
  return tot;
}

void McCOIL_prop_joint(int *max0, int *iterations, int *n0, int *k0,
                       double *A1, double *A2, double *err0, double *rho0,
                       double *pout0, double *tau0, double *minreads0,
                       int *M0, char **file) {
  int maxM = *max0, iter = *iterations, n = *n0, k = *k0;
  int i, j, m, t, it;
  const double eps = 0.05;  // half-width of the weight-transfer proposal

  g_k = k; g_A1 = A1; g_A2 = A2;
  g_err = *err0; g_rho = *rho0; g_pout = *pout0; g_tau = *tau0; g_minreads = *minreads0;

  int *M = (int *)R_alloc(n, sizeof(int));
  double *w = (double *)R_alloc((size_t)n * maxM, sizeof(double));
  g_a = (char *)R_alloc((size_t)n * maxM * k, sizeof(char));
  g_S = (double *)R_alloc((size_t)n * k, sizeof(double));
  g_L = (double *)R_alloc((size_t)n * k, sizeof(double));
  g_P = (double *)R_alloc(k, sizeof(double));
  g_sum = (double *)R_alloc(1 << MAX_ENUM, sizeof(double));
  g_t = (double *)R_alloc(1 << MAX_ENUM, sizeof(double));
  g_cnt = (int *)R_alloc(1 << MAX_ENUM, sizeof(int));
  double *S = g_S, *L = g_L, *P = g_P;
  double *Snew = (double *)R_alloc(k, sizeof(double));
  double *Lnew = (double *)R_alloc(k, sizeof(double));
  double *logZ = (double *)R_alloc(k, sizeof(double));
  double *lz1 = (double *)R_alloc(k, sizeof(double));
  double *Lall = (double *)R_alloc((size_t)n * k, sizeof(double));
  double *wn = (double *)R_alloc(maxM + 1, sizeof(double));
  int *acc = (int *)R_alloc(n + 2, sizeof(int));
  for (i = 0; i < n + 2; i++) acc[i] = 0;

  GetRNGstate();
  // strain 1 carries each locus's major allele; any further starting strains
  // get Dirichlet(1) weights and random alleles
  for (j = 0; j < k; j++) P[j] = 0.5;
  for (i = 0; i < n; i++) {
    M[i] = M0[i] < 1 ? 1 : (M0[i] > maxM ? maxM : M0[i]);
    double tot = 0;
    for (m = 0; m < M[i]; m++) {
      W(i, m) = M[i] == 1 ? 1.0 : exp_rand();
      tot += W(i, m);
      for (j = 0; j < k; j++)
        AL(i, m, j) = (char)(m == 0 ? D1(i, j) >= D2(i, j) : unif_rand() < 0.5);
    }
    for (m = 0; m < M[i]; m++) W(i, m) /= tot;
  }

  FILE *V0 = fopen(file[0], "w");

  for (it = 1; it <= iter; it++) {
    for (i = 0; i < n; i++) {
      // recompute from scratch to avoid floating-point drift
      for (j = 0; j < k; j++) {
        double s = 0;
        for (m = 0; m < M[i]; m++) s += W(i, m) * AL(i, m, j);
        S[i * k + j] = s;
        L[i * k + j] = llobs(D1(i, j), D2(i, j), s);
      }

      // Gibbs update of each strain's allele
      for (m = 0; m < M[i]; m++) {
        for (j = 0; j < k; j++) {
          double S0 = S[i * k + j] - W(i, m) * AL(i, m, j);
          double l0 = log(1 - P[j]) + llobs(D1(i, j), D2(i, j), S0);
          double l1 = log(P[j]) + llobs(D1(i, j), D2(i, j), S0 + W(i, m));
          int x = runif(0, 1) < 1.0 / (1.0 + exp(l0 - l1));
          AL(i, m, j) = (char)x;
          S[i * k + j] = S0 + W(i, m) * x;
          L[i * k + j] = llobs(D1(i, j), D2(i, j), S[i * k + j]);
        }
      }

      // pairwise weight transfer; the Dirichlet(1) prior is flat
      for (t = 0; M[i] >= 2 && t < M[i]; t++) {
        int m1 = (int)floor(runif(0, M[i]));
        int m2 = (int)floor(runif(0, M[i] - 1));
        if (m2 >= m1) m2++;
        double d = runif(-eps, eps), w1 = W(i, m1) + d, w2 = W(i, m2) - d;
        if (w1 <= 0 || w2 <= 0) continue;
        double diff = 0;
        for (j = 0; j < k; j++) {
          int da = AL(i, m1, j) - AL(i, m2, j);
          Snew[j] = S[i * k + j] + d * da;
          Lnew[j] = da ? llobs(D1(i, j), D2(i, j), Snew[j]) : L[i * k + j];
          diff += Lnew[j] - L[i * k + j];
        }
        if (log(runif(0, 1)) < diff) {
          W(i, m1) = w1;
          W(i, m2) = w2;
          for (j = 0; j < k; j++) {
            S[i * k + j] = Snew[j];
            L[i * k + j] = Lnew[j];
          }
        }
      }

      // birth/death/split/merge of a strain. Prior ratios use Dirichlet(1)
      // densities (M-1)!; births draw u ~ Beta(1, M) and splits v ~ U(0, 1).
      int Mc = M[i], Mn = Mc, kind = (int)floor(runif(0, 4)), rr = 0, ma = 0, ms = 0;
      double lr = -1e300, u = 0, v = 0, wnew = 0;
      for (m = 0; m < Mc; m++) wn[m] = W(i, m);
      if (kind == 0 && Mc < maxM) {  // birth
        u = rbeta(1.0, (double)Mc);
        rr = (int)floor(runif(0, Mc + 1));
        Mn = Mc + 1;
        for (m = Mc; m > rr; m--) wn[m] = wn[m - 1] * (1 - u);
        for (m = 0; m < rr; m++) wn[m] *= (1 - u);
        wn[rr] = u;
        lr = log((double)Mc) + (Mc - 1) * log(1 - u) - dbeta(u, 1.0, (double)Mc, 1);
      } else if (kind == 1 && Mc > 1) {  // death
        rr = (int)floor(runif(0, Mc));
        u = wn[rr];
        for (m = rr; m < Mc - 1; m++) wn[m] = wn[m + 1];
        Mn = Mc - 1;
        for (m = 0; m < Mn; m++) wn[m] /= (1 - u);
        lr = -log((double)Mn) - (Mn - 1) * log(1 - u) + dbeta(u, 1.0, (double)Mn, 1);
      } else if (kind == 2 && Mc < maxM) {  // split
        ms = (int)floor(runif(0, Mc));
        v = runif(0, 1);
        double wold = wn[ms];
        wnew = (1 - v) * wold;
        wn[ms] = v * wold;
        rr = (int)floor(runif(0, Mc + 1));
        for (m = Mc; m > rr; m--) wn[m] = wn[m - 1];
        wn[rr] = wnew;
        Mn = Mc + 1;
        lr = log((double)Mc) + log(wold);
      } else if (kind == 3 && Mc > 1) {  // merge
        rr = (int)floor(runif(0, Mc));
        ma = (int)floor(runif(0, Mc - 1));
        if (ma >= rr) ma++;
        wn[ma] += wn[rr];
        double wm = wn[ma];
        for (m = rr; m < Mc - 1; m++) wn[m] = wn[m + 1];
        Mn = Mc - 1;
        lr = -log((double)Mn) - log(wm);
      }
      if (lr > -1e299 && Mc <= MAX_ENUM && Mn <= MAX_ENUM) {
        // alleles summed out; redrawn from their conditional on acceptance
        double Lc = marg(i, Mc, &W(i, 0), 0, maxM);
        double Ln = marg(i, Mn, wn, 0, maxM);
        if (log(runif(0, 1)) < lr + Ln - Lc) {
          M[i] = Mn;
          for (m = 0; m < Mn; m++) W(i, m) = wn[m];
          marg(i, Mn, &W(i, 0), 1, maxM);
          acc[i]++;
        }
      } else if (lr > -1e299) {
        // too many strains to enumerate: explicit alleles, with the new
        // strain's alleles proposed from their conditional given the data
        double wr = 0;
        for (j = 0; j < k; j++) {
          double a1 = D1(i, j), a2 = D2(i, j), Sb, l0, l1;
          if (kind == 0 || kind == 2) {
            double wadd = kind == 0 ? u : wnew;
            Sb = kind == 0 ? (1 - u) * S[i * k + j] : S[i * k + j] - wnew * AL(i, ms, j);
            Snew[j] = Sb;
            l0 = log(1 - P[j]) + llobs(a1, a2, Sb);
            lz1[j] = log(P[j]) + llobs(a1, a2, Sb + wadd);
            logZ[j] = logsumexp2(l0, lz1[j]);
            lr += logZ[j] - L[i * k + j];
          } else {
            wr = W(i, rr);
            Sb = S[i * k + j] - wr * AL(i, rr, j);
            Snew[j] = kind == 1 ? Sb / (1 - wr) : Sb + wr * AL(i, ma, j);
            Lnew[j] = llobs(a1, a2, Snew[j]);
            l0 = log(1 - P[j]) + llobs(a1, a2, Sb);
            l1 = log(P[j]) + llobs(a1, a2, Sb + wr);
            lr += Lnew[j] - logsumexp2(l0, l1);
          }
        }
        if (log(runif(0, 1)) < lr) {
          if (kind == 0 || kind == 2) {
            double wadd = kind == 0 ? u : wnew;
            for (m = Mc; m > rr; m--)
              for (j = 0; j < k; j++) AL(i, m, j) = AL(i, m - 1, j);
            for (j = 0; j < k; j++) {
              int x = runif(0, 1) < exp(lz1[j] - logZ[j]);
              AL(i, rr, j) = (char)x;
              S[i * k + j] = Snew[j] + wadd * x;
              L[i * k + j] = llobs(D1(i, j), D2(i, j), S[i * k + j]);
            }
          } else {
            for (m = rr; m < Mc - 1; m++)
              for (j = 0; j < k; j++) AL(i, m, j) = AL(i, m + 1, j);
            for (j = 0; j < k; j++) {
              S[i * k + j] = Snew[j];
              L[i * k + j] = Lnew[j];
            }
          }
          M[i] = Mn;
          for (m = 0; m < Mn; m++) W(i, m) = wn[m];
          acc[i]++;
        }
      }
    }

    // conjugate update of population allele frequencies (uniform prior)
    for (j = 0; j < k; j++) {
      double ones = 0, tot = 0;
      for (i = 0; i < n; i++)
        for (m = 0; m < M[i]; m++) {
          ones += AL(i, m, j);
          tot++;
        }
      P[j] = rbeta(1 + ones, 1 + tot - ones);
      if (P[j] < 1e-12) P[j] = 1e-12;
      if (P[j] > 1 - 1e-12) P[j] = 1 - 1e-12;
    }

    // noise parameters: log-uniform prior on rho over [1e-6, 0.25], uniform
    // prior on pout over [0, 0.5]
    for (t = 0; t < 2; t++) {
      double r0 = g_rho, p0 = g_pout, cur = 0, can = 0;
      if (t == 0) g_rho = exp(log(g_rho) + rnorm(0, 0.1));
      else g_pout = plogis(qlogis(g_pout, 0, 1, 1, 0) + rnorm(0, 0.2), 0, 1, 1, 0);
      int ok = g_rho > 1e-6 && g_rho < 0.25 && g_pout > 0 && g_pout < 0.5;
      for (i = 0; ok && i < n; i++)
        for (j = 0; j < k; j++) {
          cur += L[i * k + j];
          Lall[i * k + j] = llobs(D1(i, j), D2(i, j), S[i * k + j]);
          can += Lall[i * k + j];
        }
      // Jacobian of the logit proposal for pout
      double jac = t == 1 ? log(g_pout * (1 - g_pout)) - log(p0 * (1 - p0)) : 0;
      if (ok && log(runif(0, 1)) < can - cur + jac) {
        for (i = 0; i < n * k; i++) L[i] = Lall[i];
        acc[n + t]++;
      } else {
        g_rho = r0;
        g_pout = p0;
      }
    }

    fprintf(V0, "%d", it);
    for (i = 0; i < n; i++) fprintf(V0, "\t%d", M[i]);
    for (j = 0; j < k; j++) fprintf(V0, "\t%.6f", P[j]);
    fprintf(V0, "\t%.8f\t%.6f\n", g_rho, g_pout);
  }

  // same column count as the trace so it reads as one table
  fprintf(V0, "total_acceptance");
  for (i = 0; i < n; i++) fprintf(V0, "\t%d", acc[i]);
  for (j = 0; j < k; j++) fprintf(V0, "\t0");
  fprintf(V0, "\t%d\t%d\n", acc[n], acc[n + 1]);
  fclose(V0);
  PutRNGstate();
}
