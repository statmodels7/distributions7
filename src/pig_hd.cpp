#include <Rcpp.h>
#include "d7_par.h"
#include <cstring>
#include <cmath>
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

// One row of the pig1 surface: the log-mass and its partials in (mu, sigma).
//
// The body is the block's, cut by order: each stage computes
// the tables its own derivatives need and returns, so asking
// for the value does not pay for four orders. With `only` the
// stage's block is written at v[0]; without it the fifteen
// columns are filled in the block's own layout.
static inline void pig1_row(double yy, double m, double sg,
                      int order, bool only, double* v) {
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
        // partials of s = sqrt(c), with c_mu = 2 sg, c_sg = 2 m, c_musg = 2
        double sc = s * c, sc2 = sc * c, sc3 = sc2 * c;
        // alpha = s / sg by the Leibniz rule against sg^{-1}
        double g1 = 1.0 / sg, g2 = g1 * g1, g3 = g2 * g1, g4 = g3 * g1,
            g5 = g4 * g1;
        double a = s * g1;
        double p[5], logS;
        psi_derivs(yy, a, order, p, &logS);
        double p1 = p[1], p2 = p[2], p3 = p[3], p4 = p[4];
        // -(y/2) log c: log of a bilinear c, partition sum written out
        double cm = 2.0 * sg, cs = 2.0 * m, cms = 2.0;
        double ic = 1.0 / c, ic2 = ic * ic, ic3 = ic2 * ic, ic4 = ic3 * ic;
        double h = -yy / 2.0;
        // the elementary pure pieces
        double im = 1.0 / m;
        int w = 0;
        // 1/sigma + psi(alpha) = (1 - s)/sigma + log S with s = sqrt(c),
        // and 1 - s = (1 - c)/(1 + s) = -2 sigma mu/(1 + s), so the pair
        // is -2 mu/(1 + s) + log S. The two terms it replaces are each of
        // size 1/sigma while their difference is of size mu: written
        // directly the value loses a digit per factor of ten in alpha and
        // is worthless past sigma = 1e-15, where the mass reads one and
        // the support no longer sums to one.
        v[w + 0] = yy * std::log(m) + h * std::log(c) -
            2.0 * m / (1.0 + s) + logS - R::lgammafn(yy + 1.0);
        if (order == 0) return;
        double s_m = sg / s,            s_s = m / s;
        double a_m = s_m * g1;
        double a_s = s_s * g1 - s * g2;
        // psi(alpha(mu, sigma)) by Faa di Bruno, written out per component
        double P_m = p1 * a_m;
        double P_s = p1 * a_s;
        double G_m = h * cm * ic;
        double G_s = h * cs * ic;
        w = only ? 0 : 1;
        v[w + 0] = yy * im + G_m + P_m;
        v[w + 1] = G_s - g2 + P_s;
        if (order == 1) return;
        double s_mm = -sg * sg / sc;
        double s_ss = -m * m / sc;
        double s_ms = 1.0 / s - sg * m / sc;
        double a_mm = s_mm * g1;
        double a_ms = s_ms * g1 - s_m * g2;
        double a_ss = s_ss * g1 - 2.0 * s_s * g2 + 2.0 * s * g3;
        double P_mm = p2 * a_m * a_m + p1 * a_mm;
        double P_ms = p2 * a_m * a_s + p1 * a_ms;
        double P_ss = p2 * a_s * a_s + p1 * a_ss;
        double G_mm = -h * cm * cm * ic2;
        double G_ms = h * (cms * ic - cm * cs * ic2);
        double G_ss = -h * cs * cs * ic2;
        w = only ? 0 : 3;
        v[w + 0] = -yy * im * im + G_mm + P_mm;
        v[w + 1] = G_ms + P_ms;
        v[w + 2] = G_ss + 2.0 * g3 + P_ss;
        if (order == 2) return;
        double s_mmm = 3.0 * sg * sg * sg / sc2;
        double s_mms = -2.0 * sg / sc + 3.0 * sg * sg * m / sc2;
        double s_mss = -2.0 * m / sc + 3.0 * m * m * sg / sc2;
        double s_sss = 3.0 * m * m * m / sc2;
        double a_mmm = s_mmm * g1;
        double a_mms = s_mms * g1 - s_mm * g2;
        double a_mss = s_mss * g1 - 2.0 * s_ms * g2 + 2.0 * s_m * g3;
        double a_sss = s_sss * g1 - 3.0 * s_ss * g2 + 6.0 * s_s * g3 -
            6.0 * s * g4;
        double P_mmm = p3 * a_m * a_m * a_m + 3.0 * p2 * a_mm * a_m +
            p1 * a_mmm;
        double P_mms = p3 * a_m * a_m * a_s +
            p2 * (a_mm * a_s + 2.0 * a_ms * a_m) + p1 * a_mms;
        double P_mss = p3 * a_m * a_s * a_s +
            p2 * (a_ss * a_m + 2.0 * a_ms * a_s) + p1 * a_mss;
        double P_sss = p3 * a_s * a_s * a_s + 3.0 * p2 * a_ss * a_s +
            p1 * a_sss;
        double G_mmm = 2.0 * h * cm * cm * cm * ic3;
        double G_mms = h * (-2.0 * cm * cms * ic2 +
                            2.0 * cm * cm * cs * ic3);
        double G_mss = h * (-2.0 * cs * cms * ic2 +
                            2.0 * cs * cs * cm * ic3);
        double G_sss = 2.0 * h * cs * cs * cs * ic3;
        w = only ? 0 : 6;
        v[w + 0] = 2.0 * yy * im * im * im + G_mmm + P_mmm;
        v[w + 1] = G_mms + P_mms;
        v[w + 2] = G_mss + P_mss;
        v[w + 3] = G_sss - 6.0 * g4 + P_sss;
        if (order == 3) return;
        double s_mmmm = -15.0 * sg * sg * sg * sg / sc3;
        double s_mmms = 9.0 * sg * sg / sc2 - 15.0 * sg * sg * sg * m / sc3;
        double s_mmss = -2.0 / sc + 12.0 * sg * m / sc2 -
            15.0 * sg * sg * m * m / sc3;
        double s_msss = 9.0 * m * m / sc2 - 15.0 * m * m * m * sg / sc3;
        double s_ssss = -15.0 * m * m * m * m / sc3;
        double a_mmmm = s_mmmm * g1;
        double a_mmms = s_mmms * g1 - s_mmm * g2;
        double a_mmss = s_mmss * g1 - 2.0 * s_mms * g2 + 2.0 * s_mm * g3;
        double a_msss = s_msss * g1 - 3.0 * s_mss * g2 + 6.0 * s_ms * g3 -
            6.0 * s_m * g4;
        double a_ssss = s_ssss * g1 - 4.0 * s_sss * g2 + 12.0 * s_ss * g3 -
            24.0 * s_s * g4 + 24.0 * s * g5;
        double P_mmmm = p4 * a_m * a_m * a_m * a_m +
            6.0 * p3 * a_mm * a_m * a_m +
            p2 * (3.0 * a_mm * a_mm + 4.0 * a_mmm * a_m) + p1 * a_mmmm;
        double P_mmms = p4 * a_m * a_m * a_m * a_s +
            p3 * (3.0 * a_ms * a_m * a_m + 3.0 * a_mm * a_m * a_s) +
            p2 * (3.0 * a_mm * a_ms + 3.0 * a_mms * a_m + a_mmm * a_s) +
            p1 * a_mmms;
        double P_mmss = p4 * a_m * a_m * a_s * a_s +
            p3 * (a_mm * a_s * a_s + a_ss * a_m * a_m +
                  4.0 * a_ms * a_m * a_s) +
            p2 * (a_mm * a_ss + 2.0 * a_ms * a_ms +
                  2.0 * a_mss * a_m + 2.0 * a_mms * a_s) +
            p1 * a_mmss;
        double P_msss = p4 * a_m * a_s * a_s * a_s +
            p3 * (3.0 * a_ms * a_s * a_s + 3.0 * a_ss * a_m * a_s) +
            p2 * (3.0 * a_ss * a_ms + 3.0 * a_mss * a_s + a_sss * a_m) +
            p1 * a_msss;
        double P_ssss = p4 * a_s * a_s * a_s * a_s +
            6.0 * p3 * a_ss * a_s * a_s +
            p2 * (3.0 * a_ss * a_ss + 4.0 * a_sss * a_s) + p1 * a_ssss;
        double G_mmmm = -6.0 * h * cm * cm * cm * cm * ic4;
        double G_mmms = h * (6.0 * cm * cm * cms * ic3 -
                             6.0 * cm * cm * cm * cs * ic4);
        double G_mmss = h * (-2.0 * cms * cms * ic2 +
                             8.0 * cm * cs * cms * ic3 -
                             6.0 * cm * cm * cs * cs * ic4);
        double G_msss = h * (6.0 * cs * cs * cms * ic3 -
                             6.0 * cs * cs * cs * cm * ic4);
        double G_ssss = -6.0 * h * cs * cs * cs * cs * ic4;
        w = only ? 0 : 10;
        v[w + 0] = -6.0 * yy * im * im * im * im + G_mmmm + P_mmmm;
        v[w + 1] = G_mmms + P_mmms;
        v[w + 2] = G_mmss + P_mmss;
        v[w + 3] = G_msss + P_msss;
        v[w + 4] = G_ssss + 24.0 * g5 + P_ssss;
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
                      int order, bool only, double* v) {
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
        psi_derivs(yy, aa, order, p, &logS, &A1);
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
