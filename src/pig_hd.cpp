#include <Rcpp.h>
#include <cstring>
using namespace Rcpp;

// Poisson-inverse Gaussian, both parametrizations, log-likelihood
// derivatives to fourth order.
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
static void psi_derivs(double y, double alpha, double h[5]) {
    if (y < 0.5) {
        h[0] = -alpha; h[1] = -1.0; h[2] = h[3] = h[4] = 0.0;
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
    double W = 0, m1 = 0, m2 = 0, m3 = 0, m4 = 0;
    for (int k = 0; k < n; ++k) {
        double u = std::exp(la[k] - mx);
        W += u;
        m1 += k * u;
        m2 += k * (k + 1.0) * u;
        m3 += k * (k + 1.0) * (k + 2.0) * u;
        m4 += k * (k + 1.0) * (k + 2.0) * (k + 3.0) * u;
    }
    m1 /= W; m2 /= W; m3 /= W; m4 /= W;
    double A1 = -m1 / alpha;
    double A2 = m2 / (alpha * alpha);
    double A3 = -m3 / (alpha * alpha * alpha);
    double A4 = m4 / (alpha * alpha * alpha * alpha);
    h[0] = -alpha + mx + std::log(W);
    h[1] = -1.0 + A1;
    h[2] = A2 - A1 * A1;
    h[3] = A3 - 3.0 * A1 * A2 + 2.0 * A1 * A1 * A1;
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
NumericMatrix pig1_hd_cpp(NumericVector y, NumericVector mu,
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
        psi_derivs(yy, al.v[0][0], h);
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
NumericMatrix pig2_hd_cpp(NumericVector y, NumericVector mu,
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
        psi_derivs(yy, a0, h);
        Jet2 l = jet_add(
            jet_add(jet_scale(jet_log(mj), yy),
                    jet_scale(jet_log(c), -yy / 2.0)),
            jet_add(jet_recip(sj), jet_compose(h, aj)));
        l.v[0][0] -= R::lgammafn(yy + 1.0);
        write_row(out, i, l);
    }
    return out;
}
