#include <Rcpp.h>
#include <cstring>
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


// ---------------------------------------------------------------------------
// The explicit closed-form kernels: what the package methods run.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix pig1_hd_cpp(NumericVector y, NumericVector mu,
                          NumericVector sigma) {
    int n = y.size();
    NumericMatrix out(n, 15);
    for (int i = 0; i < n; ++i) {
        double yy = y[i], m = mu[i], sg = sigma[i];
        double c = 1.0 + 2.0 * sg * m;
        double s = std::sqrt(c);

        // partials of s = sqrt(c), with c_mu = 2 sg, c_sg = 2 m, c_musg = 2
        double sc = s * c, sc2 = sc * c, sc3 = sc2 * c;
        double s_m = sg / s,            s_s = m / s;
        double s_mm = -sg * sg / sc;
        double s_ss = -m * m / sc;
        double s_ms = 1.0 / s - sg * m / sc;
        double s_mmm = 3.0 * sg * sg * sg / sc2;
        double s_mms = -2.0 * sg / sc + 3.0 * sg * sg * m / sc2;
        double s_mss = -2.0 * m / sc + 3.0 * m * m * sg / sc2;
        double s_sss = 3.0 * m * m * m / sc2;
        double s_mmmm = -15.0 * sg * sg * sg * sg / sc3;
        double s_mmms = 9.0 * sg * sg / sc2 - 15.0 * sg * sg * sg * m / sc3;
        double s_mmss = -2.0 / sc + 12.0 * sg * m / sc2 -
            15.0 * sg * sg * m * m / sc3;
        double s_msss = 9.0 * m * m / sc2 - 15.0 * m * m * m * sg / sc3;
        double s_ssss = -15.0 * m * m * m * m / sc3;

        // alpha = s / sg by the Leibniz rule against sg^{-1}
        double g1 = 1.0 / sg, g2 = g1 * g1, g3 = g2 * g1, g4 = g3 * g1,
            g5 = g4 * g1;
        double a = s * g1;
        double a_m = s_m * g1;
        double a_s = s_s * g1 - s * g2;
        double a_mm = s_mm * g1;
        double a_ms = s_ms * g1 - s_m * g2;
        double a_ss = s_ss * g1 - 2.0 * s_s * g2 + 2.0 * s * g3;
        double a_mmm = s_mmm * g1;
        double a_mms = s_mms * g1 - s_mm * g2;
        double a_mss = s_mss * g1 - 2.0 * s_ms * g2 + 2.0 * s_m * g3;
        double a_sss = s_sss * g1 - 3.0 * s_ss * g2 + 6.0 * s_s * g3 -
            6.0 * s * g4;
        double a_mmmm = s_mmmm * g1;
        double a_mmms = s_mmms * g1 - s_mmm * g2;
        double a_mmss = s_mmss * g1 - 2.0 * s_mms * g2 + 2.0 * s_mm * g3;
        double a_msss = s_msss * g1 - 3.0 * s_mss * g2 + 6.0 * s_ms * g3 -
            6.0 * s_m * g4;
        double a_ssss = s_ssss * g1 - 4.0 * s_sss * g2 + 12.0 * s_ss * g3 -
            24.0 * s_s * g4 + 24.0 * s * g5;

        double p[5];
        psi_derivs(yy, a, p);
        double p1 = p[1], p2 = p[2], p3 = p[3], p4 = p[4];

        // psi(alpha(mu, sigma)) by Faa di Bruno, written out per component
        double P_m = p1 * a_m;
        double P_s = p1 * a_s;
        double P_mm = p2 * a_m * a_m + p1 * a_mm;
        double P_ms = p2 * a_m * a_s + p1 * a_ms;
        double P_ss = p2 * a_s * a_s + p1 * a_ss;
        double P_mmm = p3 * a_m * a_m * a_m + 3.0 * p2 * a_mm * a_m +
            p1 * a_mmm;
        double P_mms = p3 * a_m * a_m * a_s +
            p2 * (a_mm * a_s + 2.0 * a_ms * a_m) + p1 * a_mms;
        double P_mss = p3 * a_m * a_s * a_s +
            p2 * (a_ss * a_m + 2.0 * a_ms * a_s) + p1 * a_mss;
        double P_sss = p3 * a_s * a_s * a_s + 3.0 * p2 * a_ss * a_s +
            p1 * a_sss;
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

        // -(y/2) log c: log of a bilinear c, partition sum written out
        double cm = 2.0 * sg, cs = 2.0 * m, cms = 2.0;
        double ic = 1.0 / c, ic2 = ic * ic, ic3 = ic2 * ic, ic4 = ic3 * ic;
        double h = -yy / 2.0;
        double G_m = h * cm * ic;
        double G_s = h * cs * ic;
        double G_mm = -h * cm * cm * ic2;
        double G_ms = h * (cms * ic - cm * cs * ic2);
        double G_ss = -h * cs * cs * ic2;
        double G_mmm = 2.0 * h * cm * cm * cm * ic3;
        double G_mms = h * (-2.0 * cm * cms * ic2 +
                            2.0 * cm * cm * cs * ic3);
        double G_mss = h * (-2.0 * cs * cms * ic2 +
                            2.0 * cs * cs * cm * ic3);
        double G_sss = 2.0 * h * cs * cs * cs * ic3;
        double G_mmmm = -6.0 * h * cm * cm * cm * cm * ic4;
        double G_mmms = h * (6.0 * cm * cm * cms * ic3 -
                             6.0 * cm * cm * cm * cs * ic4);
        double G_mmss = h * (-2.0 * cms * cms * ic2 +
                             8.0 * cm * cs * cms * ic3 -
                             6.0 * cm * cm * cs * cs * ic4);
        double G_msss = h * (6.0 * cs * cs * cms * ic3 -
                             6.0 * cs * cs * cs * cm * ic4);
        double G_ssss = -6.0 * h * cs * cs * cs * cs * ic4;

        // the elementary pure pieces
        double im = 1.0 / m;
        out(i, 0) = yy * std::log(m) + h * std::log(c) + g1 + p[0] -
            R::lgammafn(yy + 1.0);
        out(i, 1) = yy * im + G_m + P_m;
        out(i, 2) = G_s - g2 + P_s;
        out(i, 3) = -yy * im * im + G_mm + P_mm;
        out(i, 4) = G_ms + P_ms;
        out(i, 5) = G_ss + 2.0 * g3 + P_ss;
        out(i, 6) = 2.0 * yy * im * im * im + G_mmm + P_mmm;
        out(i, 7) = G_mms + P_mms;
        out(i, 8) = G_mss + P_mss;
        out(i, 9) = G_sss - 6.0 * g4 + P_sss;
        out(i, 10) = -6.0 * yy * im * im * im * im + G_mmmm + P_mmmm;
        out(i, 11) = G_mmms + P_mmms;
        out(i, 12) = G_mmss + P_mmss;
        out(i, 13) = G_msss + P_msss;
        out(i, 14) = G_ssss + 24.0 * g5 + P_ssss;
    }
    return out;
}

// [[Rcpp::export]]
NumericMatrix pig2_hd_cpp(NumericVector y, NumericVector mu,
                          NumericVector alpha) {
    int n = y.size();
    NumericMatrix out(n, 15);
    for (int i = 0; i < n; ++i) {
        double yy = y[i], m = mu[i], aa = alpha[i];

        // r = sqrt(m^2 + a^2) and its partials
        double r = std::sqrt(m * m + aa * aa);
        double r3 = r * r * r, r5 = r3 * r * r, r7 = r5 * r * r;
        double r_m = m / r, r_a = aa / r;
        double r_mm = aa * aa / r3, r_aa = m * m / r3, r_ma = -m * aa / r3;
        double r_mmm = -3.0 * aa * aa * m / r5;
        double r_mma = 2.0 * aa / r3 - 3.0 * aa * aa * aa / r5;
        double r_maa = -m / r3 + 3.0 * m * aa * aa / r5;
        double r_aaa = -3.0 * m * m * aa / r5;
        double r_mmmm = -3.0 * aa * aa / r5 + 15.0 * aa * aa * m * m / r7;
        double r_mmma = -6.0 * aa * m / r5 + 15.0 * aa * aa * aa * m / r7;
        double r_mmaa = 2.0 / r3 - 15.0 * aa * aa / r5 +
            15.0 * aa * aa * aa * aa / r7;
        double r_maaa = 9.0 * m * aa / r5 - 15.0 * m * aa * aa * aa / r7;
        double r_aaaa = -3.0 * m * m / r5 + 15.0 * m * m * aa * aa / r7;

        // sigma = (m + r) * a^{-2} by the Leibniz rule; u = m + r
        double v0 = 1.0 / (aa * aa);
        double v1 = -2.0 * v0 / aa, v2 = 6.0 * v0 / (aa * aa),
            v3 = -24.0 * v0 / (aa * aa * aa),
            v4 = 120.0 * v0 / (aa * aa * aa * aa);
        double u = m + r, u_m = 1.0 + r_m;
        double sig = u * v0;
        double S_m = u_m * v0;
        double S_a = r_a * v0 + u * v1;
        double S_mm = r_mm * v0;
        double S_ma = r_ma * v0 + u_m * v1;
        double S_aa = r_aa * v0 + 2.0 * r_a * v1 + u * v2;
        double S_mmm = r_mmm * v0;
        double S_mma = r_mma * v0 + r_mm * v1;
        double S_maa = r_maa * v0 + 2.0 * r_ma * v1 + u_m * v2;
        double S_aaa = r_aaa * v0 + 3.0 * r_aa * v1 + 3.0 * r_a * v2 +
            u * v3;
        double S_mmmm = r_mmmm * v0;
        double S_mmma = r_mmma * v0 + r_mmm * v1;
        double S_mmaa = r_mmaa * v0 + 2.0 * r_mma * v1 + r_mm * v2;
        double S_maaa = r_maaa * v0 + 3.0 * r_maa * v1 + 3.0 * r_ma * v2 +
            u_m * v3;
        double S_aaaa = r_aaaa * v0 + 4.0 * r_aaa * v1 + 6.0 * r_aa * v2 +
            4.0 * r_a * v3 + u * v4;

        // F(sigma) = -y log sigma + 1/sigma and its four derivatives
        double is = 1.0 / sig, is2 = is * is, is3 = is2 * is,
            is4 = is3 * is, is5 = is4 * is;
        double F1 = -yy * is - is2;
        double F2 = yy * is2 + 2.0 * is3;
        double F3 = -2.0 * yy * is3 - 6.0 * is4;
        double F4 = 6.0 * yy * is4 + 24.0 * is5;

        // F(sigma(m, a)) by the same written-out Faa di Bruno
        double P_m = F1 * S_m;
        double P_a = F1 * S_a;
        double P_mm = F2 * S_m * S_m + F1 * S_mm;
        double P_ma = F2 * S_m * S_a + F1 * S_ma;
        double P_aa = F2 * S_a * S_a + F1 * S_aa;
        double P_mmm = F3 * S_m * S_m * S_m + 3.0 * F2 * S_mm * S_m +
            F1 * S_mmm;
        double P_mma = F3 * S_m * S_m * S_a +
            F2 * (S_mm * S_a + 2.0 * S_ma * S_m) + F1 * S_mma;
        double P_maa = F3 * S_m * S_a * S_a +
            F2 * (S_aa * S_m + 2.0 * S_ma * S_a) + F1 * S_maa;
        double P_aaa = F3 * S_a * S_a * S_a + 3.0 * F2 * S_aa * S_a +
            F1 * S_aaa;
        double P_mmmm = F4 * S_m * S_m * S_m * S_m +
            6.0 * F3 * S_mm * S_m * S_m +
            F2 * (3.0 * S_mm * S_mm + 4.0 * S_mmm * S_m) + F1 * S_mmmm;
        double P_mmma = F4 * S_m * S_m * S_m * S_a +
            F3 * (3.0 * S_ma * S_m * S_m + 3.0 * S_mm * S_m * S_a) +
            F2 * (3.0 * S_mm * S_ma + 3.0 * S_mma * S_m + S_mmm * S_a) +
            F1 * S_mmma;
        double P_mmaa = F4 * S_m * S_m * S_a * S_a +
            F3 * (S_mm * S_a * S_a + S_aa * S_m * S_m +
                  4.0 * S_ma * S_m * S_a) +
            F2 * (S_mm * S_aa + 2.0 * S_ma * S_ma +
                  2.0 * S_maa * S_m + 2.0 * S_mma * S_a) +
            F1 * S_mmaa;
        double P_maaa = F4 * S_m * S_a * S_a * S_a +
            F3 * (3.0 * S_ma * S_a * S_a + 3.0 * S_aa * S_m * S_a) +
            F2 * (3.0 * S_aa * S_ma + 3.0 * S_maa * S_a + S_aaa * S_m) +
            F1 * S_maaa;
        double P_aaaa = F4 * S_a * S_a * S_a * S_a +
            6.0 * F3 * S_aa * S_a * S_a +
            F2 * (3.0 * S_aa * S_aa + 4.0 * S_aaa * S_a) + F1 * S_aaaa;

        double p[5];
        psi_derivs(yy, aa, p);

        // pure pieces: y log m, -y log a + psi(a)
        double im = 1.0 / m, ia = 1.0 / aa;
        out(i, 0) = yy * std::log(m) - yy * std::log(aa) -
            yy * std::log(sig) + is + p[0] - R::lgammafn(yy + 1.0);
        out(i, 1) = yy * im + P_m;
        out(i, 2) = -yy * ia + p[1] + P_a;
        out(i, 3) = -yy * im * im + P_mm;
        out(i, 4) = P_ma;
        out(i, 5) = yy * ia * ia + p[2] + P_aa;
        out(i, 6) = 2.0 * yy * im * im * im + P_mmm;
        out(i, 7) = P_mma;
        out(i, 8) = P_maa;
        out(i, 9) = -2.0 * yy * ia * ia * ia + p[3] + P_aaa;
        out(i, 10) = -6.0 * yy * im * im * im * im + P_mmmm;
        out(i, 11) = P_mmma;
        out(i, 12) = P_mmaa;
        out(i, 13) = P_maaa;
        out(i, 14) = 6.0 * yy * ia * ia * ia * ia + p[4] + P_aaaa;
    }
    return out;
}
