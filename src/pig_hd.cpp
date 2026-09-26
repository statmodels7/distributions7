#include <Rcpp.h>
#include "d7_par.h"
#include <cstring>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Poisson-inverse Gaussian, both parametrizations, log-likelihood
// derivatives to fourth order.
//
// TWO implementations ship side by side. The exported pig1_hd_cpp and
// pig2_hd_cpp (at the bottom) are the EXPLICIT closed-form kernels -- every
// partial written out by hand -- and are what the package methods run,
// measured 2x to 36x faster than the jet route. The bivariate-jet kernels
// below stay compiled as pig*_hd_jet_cpp: a mechanical implementation that
// shares no algebra with the explicit one, so the two are compared in the
// tests with no tolerance argument to hide behind.
//
// With c = 1 + 2 sigma mu and alpha = sqrt(c)/sigma, the half-integer order
// of the Bessel function collapses K_{y-1/2} to a finite sum and the
// prefactors cancel down to
//
//   l(y) = y log mu - (y/2) log c + 1/sigma + psi(alpha) - lgamma(y+1),
//   psi(alpha) = -alpha + log S_y(alpha),
//   S_y(alpha) = sum_{k=0}^{y-1} a_{y,k} (2 alpha)^{-k},
//   a_{y,k} = Gamma(y+k) / (Gamma(k+1) Gamma(y-k)),      S_0 = 1.
//
// The derivatives of psi in alpha come from the weighted rising-factorial
// moments of k under the (positive) terms of S, summed on the log scale;
// everything else is elementary. The composition through alpha(mu, sigma)
// is carried by a bivariate jet truncated at total order four: the value
// and the fourteen partials propagate exactly through sums, products and
// the smooth univariate functions, so no chain rule is transcribed by hand.
//
// The second parametrization keeps mu and replaces sigma by alpha itself,
// which is gamlss's PIG2 (its sigma equals this alpha); there
// sigma(mu, alpha) = (mu + sqrt(mu^2 + alpha^2)) / alpha^2 and the Bessel
// argument is a seed variable, so mu and alpha are orthogonal.

struct Jet2 {
    // partials d^{a+b} f / d mu^a d s^b for a + b <= 4, stored as v[a][b]
    double v[5][5];
    Jet2() { std::memset(v, 0, sizeof(v)); }
};

static const double BIN[5][5] = {
    {1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 2, 1, 0, 0},
    {1, 3, 3, 1, 0}, {1, 4, 6, 4, 1}
};

static Jet2 jet_const(double x) { Jet2 j; j.v[0][0] = x; return j; }

static Jet2 jet_add(const Jet2& f, const Jet2& g) {
    Jet2 out;
    for (int a = 0; a <= 4; ++a)
        for (int b = 0; a + b <= 4; ++b) out.v[a][b] = f.v[a][b] + g.v[a][b];
    return out;
}

static Jet2 jet_scale(const Jet2& f, double s) {
    Jet2 out;
    for (int a = 0; a <= 4; ++a)
        for (int b = 0; a + b <= 4; ++b) out.v[a][b] = s * f.v[a][b];
    return out;
}

static Jet2 jet_mul(const Jet2& f, const Jet2& g) {
    Jet2 out;
    for (int a = 0; a <= 4; ++a)
        for (int b = 0; a + b <= 4; ++b) {
            double s = 0.0;
            for (int p = 0; p <= a; ++p)
                for (int q = 0; q <= b; ++q)
                    s += BIN[a][p] * BIN[b][q] *
                        f.v[p][q] * g.v[a - p][b - q];
            out.v[a][b] = s;
        }
    return out;
}

// h(f) for a univariate h with derivatives h0..h4 at f's value: the Taylor
// composition h0 + h1 d + h2/2 d^2 + h3/6 d^3 + h4/24 d^4 with d = f - f0,
// exact at this truncation order because d has no constant term
static Jet2 jet_compose(const double h[5], const Jet2& f) {
    Jet2 d = f;
    d.v[0][0] = 0.0;
    Jet2 out = jet_const(h[0]);
    Jet2 p = d;
    double fac = 1.0;
    for (int k = 1; k <= 4; ++k) {
        fac *= k;
        out = jet_add(out, jet_scale(p, h[k] / fac));
        if (k < 4) p = jet_mul(p, d);
    }
    return out;
}

static Jet2 jet_log(const Jet2& f) {
    double x = f.v[0][0];
    double h[5] = { std::log(x), 1 / x, -1 / (x * x), 2 / (x * x * x),
                    -6 / (x * x * x * x) };
    return jet_compose(h, f);
}

static Jet2 jet_recip(const Jet2& f) {
    double x = f.v[0][0];
    double h[5] = { 1 / x, -1 / (x * x), 2 / (x * x * x),
                    -6 / (x * x * x * x), 24 / (x * x * x * x * x) };
    return jet_compose(h, f);
}

static Jet2 jet_sqrt(const Jet2& f) {
    double x = f.v[0][0], s = std::sqrt(x);
    double h[5] = { s, 0.5 / s, -0.25 / (s * x), 0.375 / (s * x * x),
                    -0.9375 / (s * x * x * x) };
    return jet_compose(h, f);
}

// psi(alpha) = -alpha + log S_y(alpha) and four derivatives in alpha.
// The terms of S are positive, summed by the log-sum-exp anchored at the
// largest; the derivatives are cumulants of the rising-factorial moments
// m_r = E[k (k+1) ... (k+r-1)] under the normalized terms, with
// S^(r)/S = (-1)^r m_r / alpha^r.
static void psi_derivs(double y, double alpha, int order, double h[5],
                       double* logS = nullptr, double* A1out = nullptr) {
    if (y < 0.5) {
        h[0] = -alpha; h[1] = -1.0; h[2] = h[3] = h[4] = 0.0;
        if (logS) *logS = 0.0;
        if (A1out) *A1out = 0.0;
        return;
    }
    int n = (int) y;
    double l2a = std::log(2.0 * alpha);
    double mx = -1e308;
    std::vector<double> la(n);
    for (int k = 0; k < n; ++k) {
        la[k] = R::lgammafn(y + k) - R::lgammafn(k + 1.0) -
            R::lgammafn(y - k) - k * l2a;
        if (la[k] > mx) mx = la[k];
    }
    // Only the moments the requested order reads are accumulated: the
    // value needs none of them, the score only the first.
    double W = 0, m1 = 0, m2 = 0, m3 = 0, m4 = 0;
    for (int k = 0; k < n; ++k) {
        double u = std::exp(la[k] - mx);
        W += u;
        if (order < 1) continue;
        m1 += k * u;
        if (order < 2) continue;
        m2 += k * (k + 1.0) * u;
        if (order < 3) continue;
        m3 += k * (k + 1.0) * (k + 2.0) * u;
        if (order < 4) continue;
        m4 += k * (k + 1.0) * (k + 2.0) * (k + 3.0) * u;
    }
    m1 /= W; m2 /= W; m3 /= W; m4 /= W;
    double A1 = -m1 / alpha;
    double A2 = m2 / (alpha * alpha);
    double A3 = -m3 / (alpha * alpha * alpha);
    double A4 = m4 / (alpha * alpha * alpha * alpha);
    // log S alone is reported through logS. Recovering it from h[0] by
    // adding alpha back is the cancellation this exists to avoid: h[0] is
    // of size alpha and log S is of size log(y!).
    if (logS) *logS = mx + std::log(W);
    // A1 is d log S / d alpha, of size y^2/alpha^2. Recovering it from
    // h[1] by adding one back is the same cancellation again.
    if (A1out) *A1out = A1;
    // Every entry is written whatever the order, the unread ones as zero,
    // so a caller that declares h[5] and asks for the value alone never
    // reads an uninitialized double.
    h[0] = -alpha + mx + std::log(W);
    h[1] = h[2] = h[3] = h[4] = 0.0;
    if (order < 1) return;
    h[1] = -1.0 + A1;
    if (order < 2) return;
    h[2] = A2 - A1 * A1;
    if (order < 3) return;
    h[3] = A3 - 3.0 * A1 * A2 + 2.0 * A1 * A1 * A1;
    if (order < 4) return;
    h[4] = A4 - 4.0 * A1 * A3 - 3.0 * A2 * A2 +
        12.0 * A1 * A1 * A2 - 6.0 * A1 * A1 * A1 * A1;
}

// the fourteen partials, flattened in the fixed order the R side reads
static void write_row(NumericMatrix& out, int i, const Jet2& l) {
    static const int A[15] = {0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0};
    static const int B[15] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4};
    for (int j = 0; j < 15; ++j) out(i, j) = l.v[A[j]][B[j]];
}

// [[Rcpp::export]]
NumericMatrix pig1_hd_jet_cpp(NumericVector y, NumericVector mu,
                              NumericVector sigma) {
    int n = y.size();
    NumericMatrix out(n, 15);
    for (int i = 0; i < n; ++i) {
        double yy = y[i], m0 = mu[i], s0 = sigma[i];
        Jet2 mj = jet_const(m0); mj.v[1][0] = 1.0;
        Jet2 sj = jet_const(s0); sj.v[0][1] = 1.0;
        Jet2 c = jet_add(jet_const(1.0), jet_scale(jet_mul(mj, sj), 2.0));
        Jet2 al = jet_mul(jet_sqrt(c), jet_recip(sj));
        double h[5];
        psi_derivs(yy, al.v[0][0], 4, h);
        Jet2 l = jet_add(
            jet_add(jet_scale(jet_log(mj), yy),
                    jet_scale(jet_log(c), -yy / 2.0)),
            jet_add(jet_recip(sj), jet_compose(h, al)));
        l.v[0][0] -= R::lgammafn(yy + 1.0);
        write_row(out, i, l);
    }
    return out;
}

// [[Rcpp::export]]
NumericMatrix pig2_hd_jet_cpp(NumericVector y, NumericVector mu,
                              NumericVector alpha) {
    int n = y.size();
    NumericMatrix out(n, 15);
    for (int i = 0; i < n; ++i) {
        double yy = y[i], m0 = mu[i], a0 = alpha[i];
        Jet2 mj = jet_const(m0); mj.v[1][0] = 1.0;
        Jet2 aj = jet_const(a0); aj.v[0][1] = 1.0;
        // sigma(mu, alpha) = (mu + sqrt(mu^2 + alpha^2)) / alpha^2, the
        // positive root, with no cancellation anywhere in the domain
        Jet2 root = jet_sqrt(jet_add(jet_mul(mj, mj), jet_mul(aj, aj)));
        Jet2 sj = jet_mul(jet_add(mj, root),
                          jet_recip(jet_mul(aj, aj)));
        Jet2 c = jet_add(jet_const(1.0), jet_scale(jet_mul(mj, sj), 2.0));
        double h[5];
        psi_derivs(yy, a0, 4, h);
        Jet2 l = jet_add(
            jet_add(jet_scale(jet_log(mj), yy),
                    jet_scale(jet_log(c), -yy / 2.0)),
            jet_add(jet_recip(sj), jet_compose(h, aj)));
        l.v[0][0] -= R::lgammafn(yy + 1.0);
        write_row(out, i, l);
    }
    return out;
}


// ---------------------------------------------------------------------------
// The explicit closed-form kernels: what the package methods run.
// ---------------------------------------------------------------------------

// The block's column layout, which the wrappers and the guard read:
// order k occupies PIG_LEN[k] columns starting at PIG_OFF[k].
static const int PIG_OFF[5] = {0, 1, 3, 6, 10};
static const int PIG_LEN[5] = {1, 2, 3, 4, 5};

// The partials of a bivariate function in (mu, sigma), orders 0 to 4, sit in
// fifteen slots: slot PIG_OFF[i + j] + j holds d^{i+j} / dmu^i dsigma^j.
static inline int pix(int i, int j) { return PIG_OFF[i + j] + j; }

// h(F(mu, sigma)) for a univariate h with derivatives h[0..4] at F's value
// and F's partials in F[15]: Faa di Bruno written out per component, the
// same expressions the pig2 row composes psi with.
static void pig_compose(const double h[5], const double* F, int order,
                        double* out) {
    out[0] = h[0];
    if (order < 1) return;
    double m = F[1], u = F[2];
    out[1] = h[1] * m;
    out[2] = h[1] * u;
    if (order < 2) return;
    double mm = F[3], mu = F[4], uu = F[5];
    out[3] = h[2] * m * m + h[1] * mm;
    out[4] = h[2] * m * u + h[1] * mu;
    out[5] = h[2] * u * u + h[1] * uu;
    if (order < 3) return;
    double mmm = F[6], mmu = F[7], muu = F[8], uuu = F[9];
    out[6] = h[3] * m * m * m + 3.0 * h[2] * mm * m + h[1] * mmm;
    out[7] = h[3] * m * m * u + h[2] * (mm * u + 2.0 * mu * m) + h[1] * mmu;
    out[8] = h[3] * m * u * u + h[2] * (uu * m + 2.0 * mu * u) + h[1] * muu;
    out[9] = h[3] * u * u * u + 3.0 * h[2] * uu * u + h[1] * uuu;
    if (order < 4) return;
    out[10] = h[4] * m * m * m * m + 6.0 * h[3] * mm * m * m +
        h[2] * (3.0 * mm * mm + 4.0 * mmm * m) + h[1] * F[10];
    out[11] = h[4] * m * m * m * u +
        h[3] * (3.0 * mu * m * m + 3.0 * mm * m * u) +
        h[2] * (3.0 * mm * mu + 3.0 * mmu * m + mmm * u) + h[1] * F[11];
    out[12] = h[4] * m * m * u * u +
        h[3] * (mm * u * u + uu * m * m + 4.0 * mu * m * u) +
        h[2] * (mm * uu + 2.0 * mu * mu + 2.0 * muu * m + 2.0 * mmu * u) +
        h[1] * F[12];
    out[13] = h[4] * m * u * u * u +
        h[3] * (3.0 * mu * u * u + 3.0 * uu * m * u) +
        h[2] * (3.0 * uu * mu + 3.0 * muu * u + uuu * m) + h[1] * F[13];
    out[14] = h[4] * u * u * u * u + 6.0 * h[3] * uu * u * u +
        h[2] * (3.0 * uu * uu + 4.0 * uuu * u) + h[1] * F[14];
}

// log S_y as a function of w = 1/(2 alpha), with its four derivatives in w.
// S_y(w) = sum_k a_{y,k} w^k is a polynomial with positive coefficients, so
// S^(r)(w)/S = sum_{k>=r} a_{y,k} (k)_r w^{k-r} / S is finite as w -> 0 and
// is summed with the power of w already divided out; the derivatives of
// log S are the cumulants of those ratios. In alpha the same quantities
// carry alpha^-r and the chain onto sigma cancels at a small dispersion.
static void pig_logS_w(double y, double w, int order, double L[5]) {
    L[0] = L[1] = L[2] = L[3] = L[4] = 0.0;
    if (y < 0.5) return;
    int n = (int) y;
    double lw = std::log(w);
    std::vector<double> la(n);
    double mx = -1e308;
    for (int k = 0; k < n; ++k) {
        la[k] = R::lgammafn(y + k) - R::lgammafn(k + 1.0) -
            R::lgammafn(y - k) + (k == 0 ? 0.0 : k * lw);
        if (la[k] > mx) mx = la[k];
    }
    double W = 0.0, nr[5] = {0, 0, 0, 0, 0};
    for (int k = 0; k < n; ++k) {
        W += std::exp(la[k] - mx);
        double fall = 1.0;
        for (int r = 1; r <= order && r <= k; ++r) {
            fall *= (k - r + 1);
            // a_{y,k} (k)_r w^{k-r}, the power of w taken off inside the
            // exponent so a small w never underflows before it is divided
            nr[r] += fall * std::exp(la[k] - mx - r * lw);
        }
    }
    L[0] = mx + std::log(W);
    for (int r = 1; r <= order; ++r) nr[r] /= W;
    if (order < 1) return;
    L[1] = nr[1];
    if (order < 2) return;
    L[2] = nr[2] - nr[1] * nr[1];
    if (order < 3) return;
    L[3] = nr[3] - 3.0 * nr[1] * nr[2] + 2.0 * nr[1] * nr[1] * nr[1];
    if (order < 4) return;
    L[4] = nr[4] - 4.0 * nr[1] * nr[3] - 3.0 * nr[2] * nr[2] +
        12.0 * nr[1] * nr[1] * nr[2] - 6.0 * nr[1] * nr[1] * nr[1] * nr[1];
}

// One row of the pig1 surface: the log-mass and its partials in (mu, sigma).
//
// With c = 1 + 2 sigma mu, s = sqrt(c) and w = sigma/(2 s) = 1/(2 alpha),
//
//   l(y) = y log mu - (y/2) log c + G + log S_y(w) - log y!,
//   G = 1/sigma - alpha = -2 mu/(1 + s).
//
// Every piece is smooth as sigma -> 0: G tends to -mu and w to sigma/2, and
// log S_y(w) is the log of a polynomial with positive coefficients. The
// earlier composition went through psi(alpha) = -alpha + log S, whose
// partials carry sigma^-k and cancel against those of 1/sigma: the score in
// sigma lost digits as sigma^-2, 5e-10 at 1e-4 and 4.5e-2 at 1e-8.
//
// `pre`, where given, carries log S_y and its derivatives in w from the
// recurrence the expected-information kernel runs over the support.
static inline void pig1_row(double yy, double m, double sg,
                      int order, bool only, double* v,
                      const double* pre = nullptr) {
    // The support test lives here, as it does in every other
    // compiled discrete family: a response that is negative,
    // fractional or not finite has no mass, and the kernel says so
    // rather than leaving a caller to mask afterwards. The value
    // reads -Inf, which is the log-mass there; a derivative reads
    // NaN, there being no derivative to report.
    if (!R_finite(yy) || yy < 0.0 || yy != std::floor(yy)) {
        int lo = only ? PIG_OFF[order] : 0;
        int k = only ? PIG_LEN[order] : PIG_OFF[order] + PIG_LEN[order];
        for (int j = 0; j < k; ++j)
            v[j] = (lo + j == 0) ? R_NegInf : R_NaN;
        return;
    }
    double c = 1.0 + 2.0 * sg * m;
    double s = std::sqrt(c);
    double w = 0.5 * sg / s;
    double L[5];
    if (pre) {
        for (int j = 0; j < 5; ++j) L[j] = pre[j];
    } else {
        pig_logS_w(yy, w, order, L);
    }
    // c is bilinear: c_mu = 2 sigma, c_sigma = 2 mu, c_{mu sigma} = 2
    double C[15] = {0};
    C[0] = c; C[1] = 2.0 * sg; C[2] = 2.0 * m; C[4] = 2.0;
    double ic = 1.0 / c;
    double hlog[5] = {std::log(c), ic, -ic * ic, 2.0 * ic * ic * ic,
                      -6.0 * ic * ic * ic * ic};
    double hsq[5] = {s, 0.5 / s, -0.25 / (s * c), 0.375 / (s * c * c),
                     -0.9375 / (s * c * c * c)};
    double LC[15], S[15];
    pig_compose(hlog, C, order, LC);
    pig_compose(hsq, C, order, S);
    // G = -2 mu phi(s), phi(s) = 1/(1 + s)
    double ip = 1.0 / (1.0 + s);
    double hphi[5] = {ip, -ip * ip, 2.0 * ip * ip * ip,
                      -6.0 * ip * ip * ip * ip, 24.0 * ip * ip * ip * ip * ip};
    // w = (sigma/2) psi(s), psi(s) = 1/s
    double is = 1.0 / s;
    double hpsi[5] = {is, -is * is, 2.0 * is * is * is,
                      -6.0 * is * is * is * is, 24.0 * is * is * is * is * is};
    double PH[15], PS[15], Wp[15], LS[15];
    pig_compose(hphi, S, order, PH);
    pig_compose(hpsi, S, order, PS);
    // Leibniz against the linear factors mu (for G) and sigma (for w)
    double G[15];
    for (int o = 0; o <= order; ++o) {
        for (int j = 0; j <= o; ++j) {
            int i = o - j, k = pix(i, j);
            G[k] = -2.0 * (m * PH[k] + (i > 0 ? i * PH[pix(i - 1, j)] : 0.0));
            Wp[k] = 0.5 * (sg * PS[k] + (j > 0 ? j * PS[pix(i, j - 1)] : 0.0));
        }
    }
    Wp[0] = w;
    pig_compose(L, Wp, order, LS);
    double h = -yy / 2.0;
    double im = 1.0 / m;
    double ym[5] = {yy * std::log(m), yy * im, -yy * im * im,
                    2.0 * yy * im * im * im, -6.0 * yy * im * im * im * im};
    int lo = only ? PIG_OFF[order] : 0;
    int hi = PIG_OFF[order] + PIG_LEN[order];
    for (int o = 0; o <= order; ++o) {
        for (int j = 0; j <= o; ++j) {
            int i = o - j, k = pix(i, j);
            if (k < lo || k >= hi) continue;
            double val = h * LC[k] + G[k] + LS[k];
            if (j == 0) val += ym[i];
            if (k == 0) val -= R::lgammafn(yy + 1.0);
            v[k - lo] = val;
        }
    }
}

// [[Rcpp::export]]
NumericVector pig1_pdf_cpp(NumericVector y, NumericVector mu,
                           NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector out(n);
    double* op = out.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = sigma.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[1];
        pig1_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 0, true, v);
        op[i] = v[0];
    });
    return out;
}

// [[Rcpp::export]]
List pig1_gradient_cpp(NumericVector y, NumericVector mu,
                       NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n);
    double *d0 = o0.begin(), *d1 = o1.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = sigma.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[2];
        pig1_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 1, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
    });
    return List::create(
        Named("mu") = o0,
        Named("sigma") = o1
    );
}

// [[Rcpp::export]]
List pig1_hessian_cpp(NumericVector y, NumericVector mu,
                      NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n);
    double *d0 = o0.begin(), *d1 = o1.begin(), *d2 = o2.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = sigma.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[3];
        pig1_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 2, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
        d2[i] = v[2];
    });
    return List::create(
        Named("mu_mu") = o0,
        Named("sigma_sigma") = o2,
        Named("mu_sigma") = o1
    );
}

// [[Rcpp::export]]
List pig1_deriv3_cpp(NumericVector y, NumericVector mu,
                     NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n);
    double *d0 = o0.begin(), *d1 = o1.begin(), *d2 = o2.begin(),
        *d3 = o3.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = sigma.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[4];
        pig1_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 3, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
        d2[i] = v[2];
        d3[i] = v[3];
    });
    return List::create(
        Named("mu_mu_mu") = o0,
        Named("mu_mu_sigma") = o1,
        Named("mu_sigma_sigma") = o2,
        Named("sigma_sigma_sigma") = o3
    );
}

// [[Rcpp::export]]
List pig1_deriv4_cpp(NumericVector y, NumericVector mu,
                     NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n);
    double *d0 = o0.begin(), *d1 = o1.begin(), *d2 = o2.begin(),
        *d3 = o3.begin(), *d4 = o4.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (sigma.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = sigma.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[5];
        pig1_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 4, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
        d2[i] = v[2];
        d3[i] = v[3];
        d4[i] = v[4];
    });
    return List::create(
        Named("mu_mu_mu_mu") = o0,
        Named("mu_mu_mu_sigma") = o1,
        Named("mu_mu_sigma_sigma") = o2,
        Named("mu_sigma_sigma_sigma") = o3,
        Named("sigma_sigma_sigma_sigma") = o4
    );
}

// The fifteen-column block, kept as the reference the split
// kernels are held to and as the route the jet twin compares
// against. It runs the same body at order four.
// [[Rcpp::export]]
NumericMatrix pig1_hd_cpp(NumericVector y, NumericVector mu,
                          NumericVector sigma, int threads = 1) {
    int n = y.size();
    NumericMatrix out(n, 15);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[15];
        pig1_row(y[i], mu[i], sigma[i], 4, false, v);
        for (int j = 0; j < 15; ++j) out(i, j) = v[j];
    });
    return out;
}

// One row of the pig2 surface: the log-mass and its partials in (mu, alpha).
//
// The body is the block's, cut by order: each stage computes
// the tables its own derivatives need and returns, so asking
// for the value does not pay for four orders. With `only` the
// stage's block is written at v[0]; without it the fifteen
// columns are filled in the block's own layout.
static inline void pig2_row(double yy, double m, double aa,
                      int order, bool only, double* v,
                      const double* pre = nullptr) {
    // The support test lives here, as it does in every other
    // compiled discrete family: a response that is negative,
    // fractional or not finite has no mass, and the kernel says so
    // rather than leaving a caller to mask afterwards. The value
    // reads -Inf, which is the log-mass there; a derivative reads
    // NaN, there being no derivative to report.
    if (!R_finite(yy) || yy < 0.0 || yy != std::floor(yy)) {
        int lo = only ? PIG_OFF[order] : 0;
        int k = only ? PIG_LEN[order] : PIG_OFF[order] + PIG_LEN[order];
        for (int j = 0; j < k; ++j)
            v[j] = (lo + j == 0) ? R_NegInf : R_NaN;
        return;
    }
        // r = sqrt(m^2 + a^2) and its partials, written in r_m and r_a --
        // which lie in [0, 1] -- against powers of 1/r. The direct form
        // divides by r^3, r^5 and r^7 and multiplies by a^3 and a^4, and
        // a^3 passes double.xmax at alpha = 5.6e102: there r_mma read
        // Inf - Inf and the whole Hessian came back NaN. std::hypot keeps
        // r itself finite where m^2 + a^2 does not.
        double r = std::hypot(m, aa);
        double ir = 1.0 / r, ir2 = ir * ir, ir3 = ir2 * ir;
        // sigma = (m + r) * a^{-2} by the Leibniz rule; u = m + r. The
        // powers of a are formed from ia so that a^4 cannot overflow;
        // where the true partial is below the smallest double the ladder
        // underflows to zero, which is the honest reading there.
        double ia = 1.0 / aa;
        double v0 = ia * ia;
        double v1 = -2.0 * v0 * ia, v2 = 6.0 * v0 * v0,
            v3 = -24.0 * v0 * v0 * ia,
            v4 = 120.0 * v0 * v0 * v0;
        // b = 1/sigma = a / (t + sqrt(1 + t^2)) with t = m/a, which is
        // the same number as a^2/(m + r) with neither a^2 nor m + r formed.
        double t = m * ia, wt = std::hypot(1.0, t);
        double b = aa / (t + wt);
        // F(sigma) = -y log sigma + 1/sigma has F_k = d^k F / dsigma^k
        // carrying b^k, and every Faa di Bruno term of order k has exactly
        // k factors S: the powers cancel between them. The expansion is
        // therefore written in Q_k = F_k / b^k, each linear in b, against
        // T_x = b S_x, each of size 1/a. Formed separately, F_4 = 6 y b^4
        // + 24 b^5 leaves the doubles at alpha = 4.5e61 while the product
        // it belongs to is of size a^{-3}.
        double Q1 = -yy - b;
        double Q2 = yy + 2.0 * b;
        double Q3 = -2.0 * yy - 6.0 * b;
        double Q4 = 6.0 * yy + 24.0 * b;
        double p[5], logS, A1;
        // `pre` carries what psi_derivs() would return, from the recurrence
        // the expected-information kernel runs over the support
        if (pre) {
            for (int j = 0; j < 5; ++j) p[j] = pre[j];
            logS = pre[5];
            A1 = pre[6];
        } else {
            psi_derivs(yy, aa, order, p, &logS, &A1);
        }
        // pure pieces: y log m, -y log(a sigma) and 1/sigma - alpha.
        //
        // a sigma = (m + r)/a = t + sqrt(1 + t^2), so -y log(a sigma) is
        // -y asinh(t): written as -y log a - y log sigma it is a
        // difference of two quantities of size log a.
        //
        // 1/sigma - alpha = a(a - m - r)/(m + r), and alpha - r =
        // -m^2/(alpha + r), so the pair is
        //     -m [a/(m + r)] [1 + m/(a + r)],
        // both factors bounded, tending to -m as alpha grows, which is
        // the Poisson term. Written directly it is a difference of two
        // quantities of size alpha whose value is of size mu: past
        // alpha = 1e16 the log-mass read zero, that is a probability of
        // one, and at 1e18 it read +128, one ulp of 1e18, so the mass
        // over the support summed to Inf. A fit maximizing that reward
        // drove the dispersion out of range instead of stopping.
        double im = 1.0 / m;
        double core = -m * (aa / (m + r)) * (1.0 + m / (aa + r));
        int w = 0;
        v[w + 0] = yy * std::log(m) - yy * std::asinh(t) + core + logS -
            R::lgammafn(yy + 1.0);
        if (order == 0) return;
        double r_m = m * ir, r_a = aa * ir;
        double rm2 = r_m * r_m, ra2 = r_a * r_a;
        double u = m + r, u_m = 1.0 + r_m;
        double S_m = u_m * v0;
        double S_a = r_a * v0 + u * v1;
        double T_m = b * S_m, T_a = b * S_a;
        // F(sigma(m, a)) by the same written-out Faa di Bruno
        double P_m = Q1 * T_m;
        // The score in alpha, in closed form for the same reason the value
        // is. Since r^2 - m^2 = alpha^2 the pair is 1/sigma - alpha =
        // (r - m) - alpha, and r - alpha = m^2/(r + alpha), so its
        // derivative is -m^2/(r(r + alpha)) = -r_m m/(r + alpha); the
        // other piece, -y asinh(m/alpha), gives y m/(alpha r). Both are of
        // size alpha^-2, which is what the score is. Assembled as
        // -y/alpha + psi'(alpha) + F_1 S_a the three terms are each of
        // size one: measured, the score was noise past alpha = 1e8 and
        // read exactly -1 past 1e162, where 1/alpha^2 leaves the doubles,
        // so a search saw a slope where the surface is flat.
        double dD_a = -r_m * (m / (r + aa));
        w = only ? 0 : 1;
        v[w + 0] = yy * im + P_m;
        v[w + 1] = yy * m * ir * ia + A1 + dD_a;
        if (order == 1) return;
        double r_mm = ra2 * ir, r_aa = rm2 * ir, r_ma = -r_m * r_a * ir;
        double S_mm = r_mm * v0;
        double S_ma = r_ma * v0 + u_m * v1;
        double S_aa = r_aa * v0 + 2.0 * r_a * v1 + u * v2;
        double T_mm = b * S_mm, T_ma = b * S_ma, T_aa = b * S_aa;
        double P_mm = Q2 * T_m * T_m + Q1 * T_mm;
        double P_ma = Q2 * T_m * T_a + Q1 * T_ma;
        double P_aa = Q2 * T_a * T_a + Q1 * T_aa;
        w = only ? 0 : 3;
        v[w + 0] = -yy * im * im + P_mm;
        v[w + 1] = P_ma;
        v[w + 2] = yy * ia * ia + p[2] + P_aa;
        if (order == 2) return;
        double r_mmm = -3.0 * ra2 * r_m * ir2;
        double r_mma = (2.0 * r_a - 3.0 * ra2 * r_a) * ir2;
        double r_maa = (-r_m + 3.0 * r_m * ra2) * ir2;
        double r_aaa = -3.0 * rm2 * r_a * ir2;
        double S_mmm = r_mmm * v0;
        double S_mma = r_mma * v0 + r_mm * v1;
        double S_maa = r_maa * v0 + 2.0 * r_ma * v1 + u_m * v2;
        double S_aaa = r_aaa * v0 + 3.0 * r_aa * v1 + 3.0 * r_a * v2 +
            u * v3;
        double T_mmm = b * S_mmm, T_mma = b * S_mma, T_maa = b * S_maa,
            T_aaa = b * S_aaa;
        double P_mmm = Q3 * T_m * T_m * T_m + 3.0 * Q2 * T_mm * T_m +
            Q1 * T_mmm;
        double P_mma = Q3 * T_m * T_m * T_a +
            Q2 * (T_mm * T_a + 2.0 * T_ma * T_m) + Q1 * T_mma;
        double P_maa = Q3 * T_m * T_a * T_a +
            Q2 * (T_aa * T_m + 2.0 * T_ma * T_a) + Q1 * T_maa;
        double P_aaa = Q3 * T_a * T_a * T_a + 3.0 * Q2 * T_aa * T_a +
            Q1 * T_aaa;
        w = only ? 0 : 6;
        v[w + 0] = 2.0 * yy * im * im * im + P_mmm;
        v[w + 1] = P_mma;
        v[w + 2] = P_maa;
        v[w + 3] = -2.0 * yy * ia * ia * ia + p[3] + P_aaa;
        if (order == 3) return;
        double r_mmmm = (-3.0 * ra2 + 15.0 * ra2 * rm2) * ir3;
        double r_mmma = (-6.0 * r_a * r_m + 15.0 * ra2 * r_a * r_m) * ir3;
        double r_mmaa = (2.0 - 15.0 * ra2 + 15.0 * ra2 * ra2) * ir3;
        double r_maaa = (9.0 * r_m * r_a - 15.0 * r_m * ra2 * r_a) * ir3;
        double r_aaaa = (-3.0 * rm2 + 15.0 * rm2 * ra2) * ir3;
        double S_mmmm = r_mmmm * v0;
        double S_mmma = r_mmma * v0 + r_mmm * v1;
        double S_mmaa = r_mmaa * v0 + 2.0 * r_mma * v1 + r_mm * v2;
        double S_maaa = r_maaa * v0 + 3.0 * r_maa * v1 + 3.0 * r_ma * v2 +
            u_m * v3;
        double S_aaaa = r_aaaa * v0 + 4.0 * r_aaa * v1 + 6.0 * r_aa * v2 +
            4.0 * r_a * v3 + u * v4;
        double T_mmmm = b * S_mmmm, T_mmma = b * S_mmma,
            T_mmaa = b * S_mmaa, T_maaa = b * S_maaa, T_aaaa = b * S_aaaa;
        double P_mmmm = Q4 * T_m * T_m * T_m * T_m +
            6.0 * Q3 * T_mm * T_m * T_m +
            Q2 * (3.0 * T_mm * T_mm + 4.0 * T_mmm * T_m) + Q1 * T_mmmm;
        double P_mmma = Q4 * T_m * T_m * T_m * T_a +
            Q3 * (3.0 * T_ma * T_m * T_m + 3.0 * T_mm * T_m * T_a) +
            Q2 * (3.0 * T_mm * T_ma + 3.0 * T_mma * T_m + T_mmm * T_a) +
            Q1 * T_mmma;
        double P_mmaa = Q4 * T_m * T_m * T_a * T_a +
            Q3 * (T_mm * T_a * T_a + T_aa * T_m * T_m +
                  4.0 * T_ma * T_m * T_a) +
            Q2 * (T_mm * T_aa + 2.0 * T_ma * T_ma +
                  2.0 * T_maa * T_m + 2.0 * T_mma * T_a) +
            Q1 * T_mmaa;
        double P_maaa = Q4 * T_m * T_a * T_a * T_a +
            Q3 * (3.0 * T_ma * T_a * T_a + 3.0 * T_aa * T_m * T_a) +
            Q2 * (3.0 * T_aa * T_ma + 3.0 * T_maa * T_a + T_aaa * T_m) +
            Q1 * T_maaa;
        double P_aaaa = Q4 * T_a * T_a * T_a * T_a +
            6.0 * Q3 * T_aa * T_a * T_a +
            Q2 * (3.0 * T_aa * T_aa + 4.0 * T_aaa * T_a) + Q1 * T_aaaa;
        w = only ? 0 : 10;
        v[w + 0] = -6.0 * yy * im * im * im * im + P_mmmm;
        v[w + 1] = P_mmma;
        v[w + 2] = P_mmaa;
        v[w + 3] = P_maaa;
        v[w + 4] = 6.0 * yy * ia * ia * ia * ia + p[4] + P_aaaa;
}

// [[Rcpp::export]]
NumericVector pig2_pdf_cpp(NumericVector y, NumericVector mu,
                           NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector out(n);
    double* op = out.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (alpha.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = alpha.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[1];
        pig2_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 0, true, v);
        op[i] = v[0];
    });
    return out;
}

// [[Rcpp::export]]
List pig2_gradient_cpp(NumericVector y, NumericVector mu,
                       NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n);
    double *d0 = o0.begin(), *d1 = o1.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (alpha.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = alpha.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[2];
        pig2_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 1, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
    });
    return List::create(
        Named("mu") = o0,
        Named("alpha") = o1
    );
}

// [[Rcpp::export]]
List pig2_hessian_cpp(NumericVector y, NumericVector mu,
                      NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n);
    double *d0 = o0.begin(), *d1 = o1.begin(), *d2 = o2.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (alpha.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = alpha.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[3];
        pig2_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 2, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
        d2[i] = v[2];
    });
    return List::create(
        Named("mu_mu") = o0,
        Named("alpha_alpha") = o2,
        Named("mu_alpha") = o1
    );
}

// [[Rcpp::export]]
List pig2_deriv3_cpp(NumericVector y, NumericVector mu,
                     NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n);
    double *d0 = o0.begin(), *d1 = o1.begin(), *d2 = o2.begin(),
        *d3 = o3.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (alpha.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = alpha.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[4];
        pig2_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 3, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
        d2[i] = v[2];
        d3[i] = v[3];
    });
    return List::create(
        Named("mu_mu_mu") = o0,
        Named("mu_mu_alpha") = o1,
        Named("mu_alpha_alpha") = o2,
        Named("alpha_alpha_alpha") = o3
    );
}

// [[Rcpp::export]]
List pig2_deriv4_cpp(NumericVector y, NumericVector mu,
                     NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericVector o0(n), o1(n), o2(n), o3(n), o4(n);
    double *d0 = o0.begin(), *d1 = o1.begin(), *d2 = o2.begin(),
        *d3 = o3.begin(), *d4 = o4.begin();
    bool mu_scalar = (mu.size() == 1), p2_scalar = (alpha.size() == 1);
    const double *yp = y.begin(), *mp = mu.begin(), *pp = alpha.begin();
    // the scalar branch is read INSIDE the loop, never hoisted into a
    // variable the workers write: that is the data race d7_par.h warns
    // of, and it shows only where the parameter varies by observation.
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[5];
        pig2_row(yp[i], mu_scalar ? mp[0] : mp[i],
                 p2_scalar ? pp[0] : pp[i], 4, true, v);
        d0[i] = v[0];
        d1[i] = v[1];
        d2[i] = v[2];
        d3[i] = v[3];
        d4[i] = v[4];
    });
    return List::create(
        Named("mu_mu_mu_mu") = o0,
        Named("mu_mu_mu_alpha") = o1,
        Named("mu_mu_alpha_alpha") = o2,
        Named("mu_alpha_alpha_alpha") = o3,
        Named("alpha_alpha_alpha_alpha") = o4
    );
}

// The fifteen-column block, kept as the reference the split
// kernels are held to and as the route the jet twin compares
// against. It runs the same body at order four.
// [[Rcpp::export]]
NumericMatrix pig2_hd_cpp(NumericVector y, NumericVector mu,
                          NumericVector alpha, int threads = 1) {
    int n = y.size();
    NumericMatrix out(n, 15);
    d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
        double v[15];
        pig2_row(y[i], mu[i], alpha[i], 4, false, v);
        for (int j = 0; j < 15; ++j) out(i, j) = v[j];
    });
    return out;
}

// ---------------------------------------------------------------------------
// The expected information, exactly, by one pass over the support.
// ---------------------------------------------------------------------------
//
// E[l_ab] is a sum over y of the mass times the observed component, and
// each observed component reads log S_y and the rising-factorial moments
// m_r of k under the terms of S_y. Computed afresh at each y by
// psi_derivs() those cost y terms apiece, so the sum over the support is
// quadratic in its length. They satisfy a recurrence in y instead. With
// S_y(a) = sqrt(2a/pi) e^a K_{y-1/2}(a), the Bessel recurrence
// K_{v+1} = K_{v-1} + (2v/a) K_v reads
//
//   S_{y+1} = S_{y-1} + ((2y - 1)/a) S_y,          S_0 = S_1 = 1,
//
// and differentiating it r times in a, with M_r = (-a)^r S^(r) = m_r S,
//
//   M_{r,y+1} = M_{r,y-1} + ((2y - 1)/a) sum_j C(r,j) j! M_{r-j,y},
//
// every term positive. Divided through by S_{y+1} it is a combination
// with positive weights w + v = 1 of quantities already computed, so
// nothing cancels and the recurrence is forward stable. The ratio
// S_{y+1}/S_y - 1 is carried directly, the step to log S being a log1p:
// formed as rho + t - 1 it cancels at a large a, where it is of size y/a.
struct PigSeq {
    double a;
    int y;
    double mp[5], mc[5];   // moments at y - 1 and at y, m_0 = 1
    double logS;           // log S_y
    double e;              // S_y / S_{y-1} - 1
    explicit PigSeq(double alpha) : a(alpha), y(0), logS(0.0), e(0.0) {
        for (int r = 0; r < 5; ++r) mp[r] = mc[r] = (r == 0) ? 1.0 : 0.0;
    }
    // advance from y to y + 1
    void next() {
        if (y == 0) { y = 1; e = 0.0; return; }   // S_1 = S_0 = 1
        double t = (2.0 * y - 1.0) / a;
        double d = t - e / (1.0 + e);             // S_{y+1}/S_y - 1
        double ratio = 1.0 + d;
        double w = 1.0 / (ratio * (1.0 + e));     // S_{y-1}/S_{y+1}
        double v = t / ratio;                     // (t S_y)/S_{y+1}
        double s1 = mc[1] + 1.0;
        double s2 = mc[2] + 2.0 * mc[1] + 2.0;
        double s3 = mc[3] + 3.0 * mc[2] + 6.0 * mc[1] + 6.0;
        double s4 = mc[4] + 4.0 * mc[3] + 12.0 * mc[2] + 24.0 * mc[1] + 24.0;
        double n1 = w * mp[1] + v * s1, n2 = w * mp[2] + v * s2,
            n3 = w * mp[3] + v * s3, n4 = w * mp[4] + v * s4;
        for (int r = 0; r < 5; ++r) mp[r] = mc[r];
        mc[1] = n1; mc[2] = n2; mc[3] = n3; mc[4] = n4;
        logS += std::log1p(d);
        e = d;
        ++y;
    }
    // what psi_derivs() hands a row, from the same formulas:
    // h[0..4], then log S, then A1 = d log S / d a
    void pre(double out[7]) const {
        double A1 = -mc[1] / a;
        double A2 = mc[2] / (a * a);
        double A3 = -mc[3] / (a * a * a);
        double A4 = mc[4] / (a * a * a * a);
        out[0] = -a + logS;
        out[1] = -1.0 + A1;
        out[2] = A2 - A1 * A1;
        out[3] = A3 - 3.0 * A1 * A2 + 2.0 * A1 * A1 * A1;
        out[4] = A4 - 4.0 * A1 * A3 - 3.0 * A2 * A2 +
            12.0 * A1 * A1 * A2 - 6.0 * A1 * A1 * A1 * A1;
        out[5] = logS;
        out[6] = A1;
    }
};

// The same recurrence in w = 1/(2 alpha), for pig1, whose row reads the
// derivatives of log S in w. S_{y+1} = S_{y-1} + 2(2y - 1) w S_y, and with
// N_r = S^(r)(w), r times differentiated,
//
//   N_{r,y+1} = N_{r,y-1} + 2(2y - 1) (w N_{r,y} + r N_{r-1,y}),
//
// every term positive. Divided through by S_{y+1} the weights are finite as
// w -> 0, which is the Poisson limit, so nothing there is formed from a
// power of 1/w.
struct PigSeqW {
    double w;
    int y;
    double np[5], nc[5];   // N_r / S at y - 1 and at y
    double logS, e;
    explicit PigSeqW(double ww) : w(ww), y(0), logS(0.0), e(0.0) {
        for (int r = 0; r < 5; ++r) np[r] = nc[r] = (r == 0) ? 1.0 : 0.0;
    }
    void next() {
        if (y == 0) { y = 1; e = 0.0; return; }
        double two = 2.0 * (2.0 * y - 1.0);
        double t = two * w;
        double d = t - e / (1.0 + e);
        double R = 1.0 + d;
        double Wt = 1.0 / (R * (1.0 + e)), V = t / R, U = two / R;
        double n1 = Wt * np[1] + V * nc[1] + U * nc[0];
        double n2 = Wt * np[2] + V * nc[2] + 2.0 * U * nc[1];
        double n3 = Wt * np[3] + V * nc[3] + 3.0 * U * nc[2];
        double n4 = Wt * np[4] + V * nc[4] + 4.0 * U * nc[3];
        for (int r = 0; r < 5; ++r) np[r] = nc[r];
        nc[1] = n1; nc[2] = n2; nc[3] = n3; nc[4] = n4;
        logS += std::log1p(d);
        e = d;
        ++y;
    }
    // log S and its four derivatives in w, the cumulants of N_r / S
    void pre(double out[7]) const {
        double n1 = nc[1], n2 = nc[2], n3 = nc[3], n4 = nc[4];
        out[0] = logS;
        out[1] = n1;
        out[2] = n2 - n1 * n1;
        out[3] = n3 - 3.0 * n1 * n2 + 2.0 * n1 * n1 * n1;
        out[4] = n4 - 4.0 * n1 * n3 - 3.0 * n2 * n2 +
            12.0 * n1 * n1 * n2 - 6.0 * n1 * n1 * n1 * n1;
        out[5] = out[6] = 0.0;
    }
};

// the observed components of one row, indexed by how many of their
// indices fall on the second parameter
static inline double pig_g(const double* v, int k) { return v[1 + k]; }
static inline double pig_h(const double* v, int k) { return v[3 + k]; }
static inline double pig_d3(const double* v, int k) { return v[6 + k]; }

// The pairs of hess_names(), in its order: the diagonal, then the cross.
static const int PIG_PA[3] = {0, 1, 0};
static const int PIG_PB[3] = {0, 1, 1};

// One observation, `arg` being what the recurrence runs in (alpha for pig2,
// w for pig1). order 0 fills E[3]; order 1 fills D1[6]; order 2 fills
// D1[6] and D2[9], in the orders dexpected_names() and d2expected_names()
// give. Returns false where the sum could not be completed, so the caller
// writes NA rather than a partial sum.
//
// The support is summed until the tail is negligible. Past the mode the
// ratio r_y = f(y)/f(y-1) of consecutive masses either falls towards its
// limit q = 2 sigma mu / (1 + 2 sigma mu) from above (near the Poisson
// limit, where it is about mu/y) or rises towards it from below (in the
// heavy tail), so max(r_y, q) bounds every later ratio and the tail mass
// beyond y is at most f(y) R/(1 - R) with R that maximum; `tail` carries
// q/(1 - q) = 2 sigma mu, formed without the subtraction. The summands
// grow with y, polynomially and at most as y^(4 + 2 order) relative to
// their size near the mean, and the bound is widened by that factor. A rule
// on q alone stopped the Poisson-limit sums early, and a rule on the
// accumulated mass alone does not terminate: summing 1e5 terms rounds the
// running total by more than the tolerance it is compared with.
template <class Seq, class Row>
static bool pig_expected_obs(double m, double p2, double arg, double tail,
                             int order, Row row,
                             double* E, double* D1, double* D2) {
    int ord = order + 1;
    for (int j = 0; j < 3; ++j) E[j] = 0.0;
    for (int j = 0; j < 6; ++j) D1[j] = 0.0;
    for (int j = 0; j < 9; ++j) D2[j] = 0.0;
    Seq seq(arg);
    double v[15], pre[7], cum = 0.0, fprev = 0.0;
    int pw = 4 + 2 * order;
    const long cap = 100000000L;
    for (long yy = 0; yy < cap; ++yy) {
        if (yy > 0) seq.next();
        seq.pre(pre);
        row((double) yy, m, p2, ord, v, pre);
        double f = std::exp(v[0]);
        if (!(f > 0.0)) {
            // the mass underflows below the mode, where nothing it would
            // contribute is representable; past the mean it has run out
            if ((double) yy > m) return true;
            continue;
        }
        cum += f;
        // Written through the second Bartlett identity, E[l_ab] = -E[l_a l_b],
        // and its derivatives with the measure moving. The direct forms
        // E[l_ab], E[l_abc + l_ab l_c], ... sum terms whose mean cancels
        // to a higher order in 1/alpha than the terms themselves: measured
        // at mu = 2, E[l_alpha_alpha] loses every digit by alpha = 1e5,
        // where -E[l_alpha^2] converges onto -mu^2/(2 alpha^4) past 1e9.
        for (int p = 0; p < 3; ++p) {
            int a = PIG_PA[p], b = PIG_PB[p];
            double ga = pig_g(v, a), gb = pig_g(v, b), gab = ga * gb;
            if (order == 0) { E[p] -= f * gab; continue; }
            for (int c = 0; c < 2; ++c) {
                double gc = pig_g(v, c);
                D1[2 * p + c] -= f * (pig_h(v, a + c) * gb +
                    ga * pig_h(v, b + c) + gab * gc);
            }
            if (order < 2) continue;
            for (int s = 0; s < 3; ++s) {
                int c = PIG_PA[s], d = PIG_PB[s];
                double gc = pig_g(v, c), gd = pig_g(v, d);
                double hac = pig_h(v, a + c), hbc = pig_h(v, b + c),
                    had = pig_h(v, a + d), hbd = pig_h(v, b + d);
                D2[3 * p + s] -= f * (pig_d3(v, a + c + d) * gb +
                    hac * hbd + had * hbc + ga * pig_d3(v, b + c + d) +
                    (hac * gb + ga * hbc) * gd + (had * gb + ga * hbd) * gc +
                    gab * (pig_h(v, c + d) + gc * gd));
            }
        }
        if ((double) yy > m && fprev > 0.0) {
            double r = f / fprev;
            double rt = (r < 1.0) ? r / (1.0 - r) : R_PosInf;
            double bound = f * std::max(rt, tail) *
                std::pow((yy + 1.0) / (m + 1.0), pw);
            if (bound <= 1e-17 * cum) return true;
        }
        fprev = f;
    }
    return false;
}

template <class Row>
static List pig_expected_run(NumericVector y, NumericVector mu,
                             NumericVector p2, int order, int threads,
                             const char* pname, Row row, bool pig1) {
    int n = std::max(y.size(), std::max(mu.size(), p2.size()));
    int no = (order == 0) ? 3 : ((order == 1) ? 6 : 15);
    std::vector<NumericVector> out(no);
    std::vector<double*> op(no);
    for (int j = 0; j < no; ++j) { out[j] = NumericVector(n); op[j] = out[j].begin(); }
    bool ms = (mu.size() == 1), ps = (p2.size() == 1);
    const double *mp = mu.begin(), *pp = p2.begin();
    // Every observation carries the same parameters where both are scalar,
    // which is how fit_distrib() calls: one pass serves them all. The cost of
    // a pass grows as sigma mu, about 4 s at 1e5, so paying it n times over
    // for one number was the whole cost of a fit that visits that region.
    int nrun = (ms && ps) ? std::min(n, 1) : n;
    d7::par_for(nrun, threads, d7::kMinCostly, [&](std::size_t i) {
        double m = ms ? mp[0] : mp[i], x = ps ? pp[0] : pp[i];
        double E[3], D1[6], D2[9];
        bool ok = R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0;
        if (ok) {
            if (pig1) {
                double w = 0.5 * x / std::sqrt(1.0 + 2.0 * x * m);
                ok = pig_expected_obs<PigSeqW>(m, x, w, 2.0 * x * m, order,
                                               row, E, D1, D2);
            } else {
                double t = m / x;
                double b = x / (t + std::hypot(1.0, t));   // 1/sigma
                ok = pig_expected_obs<PigSeq>(m, x, x, 2.0 * m / b, order,
                                              row, E, D1, D2);
            }
        }
        double* src = (order == 0) ? E : D1;
        for (int j = 0; j < no; ++j) {
            double val = (j < 6 || order == 0) ? src[j] : D2[j - 6];
            op[j][i] = ok ? val : NA_REAL;
        }
    });
    if (nrun < n)
        for (int j = 0; j < no; ++j)
            for (int i = 1; i < n; ++i) op[j][i] = op[j][0];
    std::string P = pname;
    std::vector<std::string> hs = {"mu_mu", P + "_" + P, "mu_" + P};
    std::vector<std::string> ps2 = {"mu", P};
    List res(no);
    CharacterVector nms(no);
    for (int j = 0; j < no; ++j) res[j] = out[j];
    if (order == 0) {
        for (int j = 0; j < 3; ++j) nms[j] = hs[j];
    } else {
        for (int p = 0; p < 3; ++p)
            for (int c = 0; c < 2; ++c) nms[2 * p + c] = hs[p] + "_" + ps2[c];
        if (order == 2)
            for (int p = 0; p < 3; ++p)
                for (int s = 0; s < 3; ++s) nms[6 + 3 * p + s] = hs[p] + "_" + hs[s];
    }
    res.names() = nms;
    return res;
}

// [[Rcpp::export]]
List pig2_expected_cpp(NumericVector y, NumericVector mu, NumericVector alpha,
                       int order, int threads = 1) {
    auto row = [](double yy, double m, double x, int ord, double* v,
                  const double* pre) {
        pig2_row(yy, m, x, ord, false, v, pre);
    };
    return pig_expected_run(y, mu, alpha, order, threads, "alpha", row, false);
}

// [[Rcpp::export]]
List pig1_expected_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                       int order, int threads = 1) {
    auto row = [](double yy, double m, double x, int ord, double* v,
                  const double* pre) {
        pig1_row(yy, m, x, ord, false, v, pre);
    };
    return pig_expected_run(y, mu, sigma, order, threads, "sigma", row, true);
}
