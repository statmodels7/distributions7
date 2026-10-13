#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// The location-scale log-tail route of ls_log_tail_derivs(), compiled: log P
// = g(z), z = sg (x - mu) w(t), with the derivatives g^(k) of the base written
// out and the partials of z in (mu, t) put together by Faa di Bruno over the
// set partitions of the multi-index. The gaussian takes -z^2/2 out of the sum
// right of the location, as the R route does.

namespace {

// g^(1..4) of the upper tail of the standard base at z; base 0 gaussian,
// 1 laplace, 2 logistic, 3 gumbel, 4 gumbel lower tail (log F = -e^-z)
void ls_g(int base, double z, double* g) {
    switch (base) {
    case 1: {
        if (z > 0) { g[0] = -1; g[1] = g[2] = g[3] = 0; return; }
        const double u = std::exp(z) / 2, v = 1 - u;
        g[0] = -u / v;
        g[1] = -u / (v * v);
        g[2] = -u * (1 + u) / (v * v * v);
        g[3] = -u * (1 + 4 * u + u * u) / (v * v * v * v);
        return;
    }
    case 2: {
        const double p = R::plogis(z, 0, 1, 1, 0), s = R::plogis(-z, 0, 1, 1, 0);
        const double ps = p * s;
        g[0] = -p; g[1] = -ps; g[2] = -ps * (s - p); g[3] = -ps * (1 - 6 * ps);
        return;
    }
    case 3: {
        const double w = std::exp(-z);
        double t[4];
        if (w < 0.5) {
            // theta^k beta(w) = sum_n n^k B_n / n! w^n, n = 0..30
            static double c[4][31];
            static bool init = false;
            if (!init) {
                static const double Bn[31] = {
                    1, -1.0/2, 1.0/6, 0, -1.0/30, 0, 1.0/42, 0, -1.0/30, 0, 5.0/66, 0,
                    -691.0/2730, 0, 7.0/6, 0, -3617.0/510, 0, 43867.0/798, 0,
                    -174611.0/330, 0, 854513.0/138, 0, -236364091.0/2730, 0,
                    8553103.0/6, 0, -23749461029.0/870, 0, 8615841276005.0/14322};
                double f = 1.0;
                for (int n = 0; n <= 30; ++n) {
                    if (n > 0) f *= n;
                    double p = 1.0;
                    for (int k = 0; k < 4; ++k) {
                        c[k][n] = p * Bn[n] / f;
                        p *= n;
                    }
                }
                init = true;
            }
            for (int k = 0; k < 4; ++k) {
                double acc = 0.0;
                for (int n = 30; n >= 0; --n) acc = acc * w + c[k][n];
                t[k] = acc;
            }
        } else {
            const double b = w / std::expm1(w), E = 1 / (-std::expm1(-w));
            const double t1 = b * (1 - w * E);
            const double tE = -w * E * (E - 1);
            const double twE = w * E + w * tE;
            const double t2 = t1 * (1 - w * E) - b * twE;
            const double ttE = -(w * E * (E - 1) + w * tE * (E - 1) + w * E * tE);
            const double ttwE = twE + w * tE + w * ttE;
            t[0] = b; t[1] = t1; t[2] = t2;
            t[3] = t2 * (1 - w * E) - 2 * t1 * twE - b * ttwE;
        }
        g[0] = -t[0]; g[1] = t[1]; g[2] = -t[2]; g[3] = t[3];
        return;
    }
    case 4: {
        const double e = std::exp(-z);
        g[0] = e; g[1] = -e; g[2] = e; g[3] = -e;
        return;
    }
    default: {
        if (z < 9) {
            const double m = std::exp(R::dnorm(z, 0, 1, 1) - R::pnorm(z, 0, 1, 0, 1));
            const double r = m - z, m1 = m * r, m2 = m1 * r + m * (m1 - 1);
            const double m3 = m2 * r + 2 * m1 * (m1 - 1) + m * m2;
            g[0] = -m; g[1] = -m1; g[2] = -m2; g[3] = -m3;
            return;
        }
        // log P(t), t = z^-2, P = sum_n (-1)^n (2n - 1)!! t^n, n <= 40
        const double t = 1 / (z * z);
        double P[5] = {0, 0, 0, 0, 0};
        static double a[41];
        static bool init = false;
        if (!init) {
            a[0] = 1;
            for (int n = 1; n <= 40; ++n) a[n] = -a[n - 1] * (2 * n - 1);
            init = true;
        }
        for (int j = 0; j <= 4; ++j) {
            double acc = 0.0;
            for (int k = 40; k >= j; --k) {
                double ff = 1.0;
                for (int i = 0; i < j; ++i) ff *= (k - i);
                acc = acc * t + ff * a[k];
            }
            P[j] = acc;
        }
        const double L1 = P[1] / P[0];
        const double L2 = P[2] / P[0] - L1 * L1;
        const double L3 = P[3] / P[0] - 3 * P[2] * P[1] / (P[0] * P[0]) + 2 * L1 * L1 * L1;
        const double L4 = P[4] / P[0] - 4 * P[3] * P[1] / (P[0] * P[0]) -
            3 * (P[2] / P[0]) * (P[2] / P[0]) +
            12 * P[2] * P[1] * P[1] / (P[0] * P[0] * P[0]) - 6 * L1 * L1 * L1 * L1;
        const double z2 = z * z, z3 = z2 * z, z4 = z3 * z;
        const double t1 = -2 / z3, t2 = 6 / z4, t3 = -24 / (z4 * z), t4 = 120 / (z4 * z2);
        g[0] = -z - 1 / z + L1 * t1;
        g[1] = -1 + 1 / z2 + L2 * t1 * t1 + L1 * t2;
        g[2] = -2 / z3 + L3 * t1 * t1 * t1 + 3 * L2 * t1 * t2 + L1 * t3;
        g[3] = 6 / z4 + L4 * t1 * t1 * t1 * t1 + 6 * L3 * t1 * t1 * t2 +
            L2 * (3 * t2 * t2 + 4 * t1 * t3) + L1 * t4;
        return;
    }
    }
}

// w(t) and w(t)^2 with their derivatives to order 4: 0 inv, 1 rate,
// 2 invsqrt, 3 sqrt
void ls_w(int wid, double s, double* W, double* W2) {
    switch (wid) {
    case 1:
        W[0] = s; W[1] = 1; W[2] = W[3] = W[4] = 0;
        W2[0] = s * s; W2[1] = 2 * s; W2[2] = 2; W2[3] = W2[4] = 0;
        return;
    case 2: {
        const double r = 1 / std::sqrt(s);
        W[0] = r; W[1] = -0.5 * r / s; W[2] = 0.75 * r / (s * s);
        W[3] = -1.875 * r / (s * s * s); W[4] = 6.5625 * r / (s * s * s * s);
        W2[0] = 1 / s; W2[1] = -1 / (s * s); W2[2] = 2 / (s * s * s);
        W2[3] = -6 / (s * s * s * s); W2[4] = 24 / (s * s * s * s * s);
        return;
    }
    case 3: {
        const double r = std::sqrt(s);
        W[0] = r; W[1] = 0.5 / r; W[2] = -0.25 / (r * s); W[3] = 0.375 / (r * s * s);
        W[4] = -0.9375 / (r * s * s * s);
        W2[0] = s; W2[1] = 1; W2[2] = W2[3] = W2[4] = 0;
        return;
    }
    default:
        W[0] = 1 / s; W[1] = -1 / (s * s); W[2] = 2 / (s * s * s);
        W[3] = -6 / (s * s * s * s); W[4] = 24 / (s * s * s * s * s);
        W2[0] = 1 / (s * s); W2[1] = -2 / (s * s * s); W2[2] = 6 / (s * s * s * s);
        W2[3] = -24 / (s * s * s * s * s); W2[4] = 120 / (s * s * s * s * s * s);
        return;
    }
}

// the set partitions of {0..o-1} as lists of block bitmasks
void partitions(int o, std::vector<std::vector<int>>& out) {
    out.clear();
    std::vector<int> a(o, 0);
    // restricted growth strings
    while (true) {
        int nb = 0;
        for (int i = 0; i < o; ++i) nb = std::max(nb, a[i] + 1);
        std::vector<int> bl(nb, 0);
        for (int i = 0; i < o; ++i) bl[a[i]] |= 1 << i;
        out.push_back(bl);
        int i = o - 1;
        while (i > 0) {
            int mx = 0;
            for (int k = 0; k < i; ++k) mx = std::max(mx, a[k]);
            if (a[i] <= mx) { ++a[i]; break; }
            a[i] = 0;
            --i;
        }
        if (i == 0) break;
    }
}

} // namespace

// The components (j location indices, o - j of the second parameter), j from o
// down to 0, of the derivatives of order o of log P, one column each.
// [[Rcpp::export]]
NumericMatrix ls_log_tail_cpp(NumericVector x, NumericVector mu, NumericVector t,
                              int base, int wid, bool flip, bool split, int o) {
    const int n = x.size();
    NumericMatrix out(n, o + 1);
    std::vector<std::vector<int>> parts;
    partitions(o, parts);
    const double sg = flip ? -1.0 : 1.0;
    for (int i = 0; i < n; ++i) {
        double W[5], W2[5], g[4];
        ls_w(wid, t[i % t.size()], W, W2);
        const double d = sg * (x[i] - mu[i % mu.size()]);
        const double z = d * W[0];
        if (!std::isfinite(z)) { for (int c = 0; c <= o; ++c) out(i, c) = NA_REAL; continue; }
        ls_g(base, z, g);
        const bool pos = split && z > 0;
        if (pos) { g[0] += z; g[1] += 1; }
        for (int c = 0; c <= o; ++c) {
            const int j = o - c;          // positions 0..j-1 name the location
            const int locmask = (1 << j) - 1;
            double v = 0.0;
            for (const auto& bl : parts) {
                double term = g[bl.size() - 1];
                for (int b : bl) {
                    const int jb = __builtin_popcount(b & locmask);
                    const int kb = __builtin_popcount(b) - jb;
                    double zd;
                    if (jb >= 2) zd = 0.0;
                    else if (jb == 1) zd = -sg * W[kb];
                    else zd = d * W[kb];
                    term *= zd;
                    if (term == 0.0) break;
                }
                v += term;
            }
            if (pos && j <= 2) {
                const double Dd = (j == 0) ? d * d : (j == 1 ? -2 * sg * d : 2.0);
                v -= 0.5 * Dd * W2[o - j];
            }
            out(i, c) = v;
        }
    }
    return out;
}
