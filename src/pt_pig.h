#ifndef D7_PT_PIG_H
#define D7_PT_PIG_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "pig_asym.h"
#include "pt_sqrt.h"

// The Poisson-inverse Gaussian in both parametrizations: what the kernels of
// pig_hd.cpp and the scalar registry share. The special functions are
// log S_y and its derivatives, from psi_derivs() and pig_logS_w() at one
// response or from the recurrences PigSeq and PigSeqW over the support;
// pig1_slot() and the pig2 component functions below take them as
// arguments and compute one partial each, with the arithmetic of the row
// kernels pig1_row() and pig2_row(), slot for slot.

namespace d7 {

// psi(alpha) = -alpha + log S_y(alpha) and four derivatives in alpha.
// The terms of S are positive, summed by the log-sum-exp anchored at the
// largest; the derivatives are cumulants of the rising-factorial moments
// m_r = E[k (k+1) ... (k+r-1)] under the normalized terms, with
// S^(r)/S = (-1)^r m_r / alpha^r.
inline void psi_derivs(double y, double alpha, int order, double h[5],
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

// The block's column layout, which the wrappers and the guard read:
// order k occupies PIG_LEN[k] columns starting at PIG_OFF[k].
const int PIG_OFF[5] = {0, 1, 3, 6, 10};
const int PIG_LEN[5] = {1, 2, 3, 4, 5};

// The partials of a bivariate function in (mu, sigma), orders 0 to 4, sit in
// fifteen slots: slot PIG_OFF[i + j] + j holds d^{i+j} / dmu^i dsigma^j.
inline int pix(int i, int j) { return PIG_OFF[i + j] + j; }

// h(F(mu, sigma)) for a univariate h with derivatives h[0..4] at F's value
// and F's partials in F[15]: Faa di Bruno written out per component, the
// same expressions the pig2 row composes psi with.
inline void pig_compose(const double h[5], const double* F, int order,
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
inline void pig_logS_w(double y, double w, int order, double L[5]) {
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


// The threshold on tau = (1 + mu) x, x = sigma (pig1) or 1/alpha (pig2),
// below which the series of pig_asym.h replaces the sums. Measured against
// exact sums (stabilita/pig_cmp.R, pig*_dexpected_ref.py): at tau = 0.05 the
// order-20 series is within 2e-15 of every component, where the sums lose up
// to 5e-11; above 0.1 the series converges too slowly.
const double PIG_ASYM_CUT = 0.06;

// slot k of pig_compose(): one partial of h(F), from the partials of F in
// the down-set of k, with the expressions pig_compose() writes
inline double pig_compose_slot(const double h[5], const double* F, int k) {
    switch (k) {
    case 0: return h[0];
    case 1: return h[1] * F[1];
    case 2: return h[1] * F[2];
    case 3: { double m = F[1]; return h[2] * m * m + h[1] * F[3]; }
    case 4: { double m = F[1], u = F[2]; return h[2] * m * u + h[1] * F[4]; }
    case 5: { double u = F[2]; return h[2] * u * u + h[1] * F[5]; }
    case 6: {
        double m = F[1], mm = F[3];
        return h[3] * m * m * m + 3.0 * h[2] * mm * m + h[1] * F[6];
    }
    case 7: {
        double m = F[1], u = F[2], mm = F[3], mu = F[4];
        return h[3] * m * m * u + h[2] * (mm * u + 2.0 * mu * m) + h[1] * F[7];
    }
    case 8: {
        double m = F[1], u = F[2], mu = F[4], uu = F[5];
        return h[3] * m * u * u + h[2] * (uu * m + 2.0 * mu * u) + h[1] * F[8];
    }
    case 9: {
        double u = F[2], uu = F[5];
        return h[3] * u * u * u + 3.0 * h[2] * uu * u + h[1] * F[9];
    }
    case 10: {
        double m = F[1], mm = F[3], mmm = F[6];
        return h[4] * m * m * m * m + 6.0 * h[3] * mm * m * m +
            h[2] * (3.0 * mm * mm + 4.0 * mmm * m) + h[1] * F[10];
    }
    case 11: {
        double m = F[1], u = F[2], mm = F[3], mu = F[4], mmm = F[6],
            mmu = F[7];
        return h[4] * m * m * m * u +
            h[3] * (3.0 * mu * m * m + 3.0 * mm * m * u) +
            h[2] * (3.0 * mm * mu + 3.0 * mmu * m + mmm * u) + h[1] * F[11];
    }
    case 12: {
        double m = F[1], u = F[2], mm = F[3], mu = F[4], uu = F[5],
            mmu = F[7], muu = F[8];
        return h[4] * m * m * u * u +
            h[3] * (mm * u * u + uu * m * m + 4.0 * mu * m * u) +
            h[2] * (mm * uu + 2.0 * mu * mu + 2.0 * muu * m + 2.0 * mmu * u) +
            h[1] * F[12];
    }
    case 13: {
        double m = F[1], u = F[2], mu = F[4], uu = F[5], muu = F[8],
            uuu = F[9];
        return h[4] * m * u * u * u +
            h[3] * (3.0 * mu * u * u + 3.0 * uu * m * u) +
            h[2] * (3.0 * uu * mu + 3.0 * muu * u + uuu * m) + h[1] * F[13];
    }
    default: {
        double u = F[2], uu = F[5], uuu = F[9];
        return h[4] * u * u * u * u + 6.0 * h[3] * uu * u * u +
            h[2] * (3.0 * uu * uu + 4.0 * uuu * u) + h[1] * F[14];
    }
    }
}

inline bool pig_in_support(double yy) {
    return R_finite(yy) && yy >= 0.0 && yy == std::floor(yy);
}

// ---- pig1 (mu, sigma) -------------------------------------------------------

// w = 1/(2 alpha), at which log S is differentiated; pig1_row()'s expression
inline double pig1_w(double m, double sg) {
    double c = 1.0 + 2.0 * sg * m;
    double s = d7::sqrt_cr(c);
    return 0.5 * sg / s;
}

// The partial d^(i+j) l / dmu^i dsigma^j at one response, L being log S_y and
// its derivatives in w (pig_logS_w() or PigSeqW). It runs pig1_row()'s body
// on the slots of the down-set of (i, j) alone.
inline double pig1_slot(int i, int j, double yy, double m, double sg,
                        const double* L) {
    const int k = pix(i, j);
    if (!pig_in_support(yy)) return (k == 0) ? R_NegInf : R_NaN;
    double c = 1.0 + 2.0 * sg * m;
    double s = d7::sqrt_cr(c);
    double w = 0.5 * sg / s;
    double C[15] = {0};
    C[0] = c; C[1] = 2.0 * sg; C[2] = 2.0 * m; C[4] = 2.0;
    double ic = 1.0 / c;
    double hlog[5] = {std::log(c), ic, -ic * ic, 2.0 * ic * ic * ic,
                      -6.0 * ic * ic * ic * ic};
    double hsq[5] = {s, 0.5 / s, -0.25 / (s * c), 0.375 / (s * c * c),
                     -0.9375 / (s * c * c * c)};
    double ip = 1.0 / (1.0 + s);
    double hphi[5] = {ip, -ip * ip, 2.0 * ip * ip * ip,
                      -6.0 * ip * ip * ip * ip, 24.0 * ip * ip * ip * ip * ip};
    double is = 1.0 / s;
    double hpsi[5] = {is, -is * is, 2.0 * is * is * is,
                      -6.0 * is * is * is * is, 24.0 * is * is * is * is * is};
    double S[15] = {0}, PS[15] = {0}, Wp[15] = {0};
    for (int o = 0; o <= i + j; ++o) {
        for (int jj = 0; jj <= o; ++jj) {
            int ii = o - jj;
            if (ii > i || jj > j) continue;
            int kk = pix(ii, jj);
            S[kk] = pig_compose_slot(hsq, C, kk);
            PS[kk] = pig_compose_slot(hpsi, S, kk);
            Wp[kk] = 0.5 * (sg * PS[kk] + (jj > 0 ? jj * PS[pix(ii, jj - 1)] : 0.0));
        }
    }
    Wp[0] = w;
    double PHk = pig_compose_slot(hphi, S, k);
    double PHi = (i > 0) ? pig_compose_slot(hphi, S, pix(i - 1, j)) : 0.0;
    double G = -2.0 * (m * PHk + (i > 0 ? i * PHi : 0.0));
    double LC = pig_compose_slot(hlog, C, k);
    double LS = pig_compose_slot(L, Wp, k);
    double h = -yy / 2.0;
    double val = h * LC + G + LS;
    if (j == 0) {
        double im = 1.0 / m;
        double ym;
        switch (i) {
        case 0: ym = yy * std::log(m); break;
        case 1: ym = yy * im; break;
        case 2: ym = -yy * im * im; break;
        case 3: ym = 2.0 * yy * im * im * im; break;
        default: ym = -6.0 * yy * im * im * im * im;
        }
        val += ym;
    }
    if (k == 0) val -= R::lgammafn(yy + 1.0);
    return val;
}

inline double pig1_value(double yy, double m, double sg, const double* L) {
    return pig1_slot(0, 0, yy, m, sg, L);
}
inline double pig1_score_mu(double yy, double m, double sg, const double* L) {
    return pig1_slot(1, 0, yy, m, sg, L);
}
inline double pig1_score_sigma(double yy, double m, double sg, const double* L) {
    return pig1_slot(0, 1, yy, m, sg, L);
}
inline double pig1_hess_mu_mu(double yy, double m, double sg, const double* L) {
    return pig1_slot(2, 0, yy, m, sg, L);
}
inline double pig1_hess_sigma_sigma(double yy, double m, double sg,
                                    const double* L) {
    return pig1_slot(0, 2, yy, m, sg, L);
}

// ---- pig2 (mu, alpha) -------------------------------------------------------

// the per-observation pieces pig2_row() forms before its first partial
struct Pig2Base {
    double r, ir, ia, v0, v1, v2, t, b, Q1, Q2, im, r_m, r_a, rm2, ra2, u,
        u_m, T_m, T_a;
};

inline Pig2Base pig2_base(double yy, double m, double aa) {
    Pig2Base B;
    B.r = std::hypot(m, aa);
    B.ir = 1.0 / B.r;
    B.ia = 1.0 / aa;
    B.v0 = B.ia * B.ia;
    B.v1 = -2.0 * B.v0 * B.ia;
    B.v2 = 6.0 * B.v0 * B.v0;
    B.t = m * B.ia;
    double wt = std::hypot(1.0, B.t);
    B.b = aa / (B.t + wt);
    B.Q1 = -yy - B.b;
    B.Q2 = yy + 2.0 * B.b;
    B.im = 1.0 / m;
    B.r_m = m * B.ir;
    B.r_a = aa * B.ir;
    B.rm2 = B.r_m * B.r_m;
    B.ra2 = B.r_a * B.r_a;
    B.u = m + B.r;
    B.u_m = 1.0 + B.r_m;
    double S_m = B.u_m * B.v0;
    double S_a = B.r_a * B.v0 + B.u * B.v1;
    B.T_m = B.b * S_m;
    B.T_a = B.b * S_a;
    return B;
}

// logS the log of S_y; A1 = d log S / d alpha; p2 the second derivative of
// psi(alpha) = -alpha + log S, from psi_derivs() or PigSeq
inline double pig2_value(double yy, double m, double aa, double logS) {
    if (!pig_in_support(yy)) return R_NegInf;
    double r = std::hypot(m, aa);
    double ia = 1.0 / aa;
    double t = m * ia;
    double core = -m * (aa / (m + r)) * (1.0 + m / (aa + r));
    return yy * std::log(m) - yy * std::asinh(t) + core + logS -
        R::lgammafn(yy + 1.0);
}

inline double pig2_score_mu(double yy, double m, const Pig2Base& B) {
    if (!pig_in_support(yy)) return R_NaN;
    double P_m = B.Q1 * B.T_m;
    return yy * B.im + P_m;
}

inline double pig2_score_alpha(double yy, double m, double aa,
                               const Pig2Base& B, double A1) {
    if (!pig_in_support(yy)) return R_NaN;
    double dD_a = -B.r_m * (m / (B.r + aa));
    return yy * m * B.ir * B.ia + A1 + dD_a;
}

inline double pig2_hess_mu_mu(double yy, const Pig2Base& B) {
    if (!pig_in_support(yy)) return R_NaN;
    double r_mm = B.ra2 * B.ir;
    double S_mm = r_mm * B.v0;
    double T_mm = B.b * S_mm;
    double P_mm = B.Q2 * B.T_m * B.T_m + B.Q1 * T_mm;
    return -yy * B.im * B.im + P_mm;
}

inline double pig2_hess_mu_alpha(double yy, const Pig2Base& B) {
    if (!pig_in_support(yy)) return R_NaN;
    double r_ma = -B.r_m * B.r_a * B.ir;
    double S_ma = r_ma * B.v0 + B.u_m * B.v1;
    double T_ma = B.b * S_ma;
    return B.Q2 * B.T_m * B.T_a + B.Q1 * T_ma;
}

inline double pig2_hess_alpha_alpha(double yy, const Pig2Base& B, double p2) {
    if (!pig_in_support(yy)) return R_NaN;
    double r_aa = B.rm2 * B.ir;
    double S_aa = r_aa * B.v0 + 2.0 * B.r_a * B.v1 + B.u * B.v2;
    double T_aa = B.b * S_aa;
    double P_aa = B.Q2 * B.T_a * B.T_a + B.Q1 * T_aa;
    return yy * B.ia * B.ia + p2 + P_aa;
}

// ---- the routers of the scalar registry ------------------------------------

inline void pig1_score_curv(int k, double y, const double* th, double* out) {
    const double m = th[0], sg = th[1];
    double L[5] = {0, 0, 0, 0, 0};
    if (pig_in_support(y)) pig_logS_w(y, pig1_w(m, sg), 2, L);
    if (k == 0) {
        out[0] = pig1_score_mu(y, m, sg, L);
        out[1] = pig1_hess_mu_mu(y, m, sg, L);
    } else {
        out[0] = pig1_score_sigma(y, m, sg, L);
        out[1] = pig1_hess_sigma_sigma(y, m, sg, L);
    }
}

inline void pig2_score_curv(int k, double y, const double* th, double* out) {
    const double m = th[0], aa = th[1];
    double p[5] = {0, 0, 0, 0, 0}, logS = 0.0, A1 = 0.0;
    if (pig_in_support(y)) psi_derivs(y, aa, 2, p, &logS, &A1);
    const Pig2Base B = pig2_base(y, m, aa);
    if (k == 0) {
        out[0] = pig2_score_mu(y, m, B);
        out[1] = pig2_hess_mu_mu(y, B);
    } else {
        out[0] = pig2_score_alpha(y, m, aa, B, A1);
        out[1] = pig2_hess_alpha_alpha(y, B, p[2]);
    }
}

// The (k, k) expected second derivative and its derivative in parameter k,
// by pig_expected_obs()'s pass over the support, the information and its
// derivative in one pass: the information stops at its own rule (pw = 4)
// and the derivative runs on to its rule (pw = 6), each with the sums the
// kernels pig*_expected_cpp and pig*_dexpected1_cpp form. Through the
// second Bartlett identity, E[l_kk] = -E[l_k^2] and
// d_k E[l_kk] = -E[2 l_kk l_k + l_k^3]. Returns in ok[0], ok[1] whether each
// sum was completed.
template <class Seq, class Comp>
inline void pig_diag_obs(double m, double arg, double tail, Comp comp,
                         double* E, double* D, bool* ok) {
    *E = 0.0;
    *D = 0.0;
    ok[0] = ok[1] = false;
    Seq seq(arg);
    double pre[7], cum = 0.0, fprev = 0.0;
    bool e_done = false;
    const long cap = 100000000L;
    for (long yy = 0; yy < cap; ++yy) {
        if (yy > 0) seq.next();
        seq.pre(pre);
        double val, g, h;
        comp((double) yy, pre, &val, &g, &h);
        double f = std::exp(val);
        if (!(f > 0.0)) {
            if ((double) yy > m) {
                ok[0] = ok[1] = true;
                return;
            }
            continue;
        }
        cum += f;
        double gab = g * g;
        if (!e_done) *E -= f * gab;
        *D -= f * (h * g + g * h + gab * g);
        if ((double) yy > m && fprev > 0.0) {
            double r = f / fprev;
            double rt = (r < 1.0) ? r / (1.0 - r) : R_PosInf;
            double base = f * std::max(rt, tail);
            if (!e_done && base * std::pow((yy + 1.0) / (m + 1.0), 4) <=
                1e-17 * cum) {
                e_done = true;
                ok[0] = true;
            }
            if (base * std::pow((yy + 1.0) / (m + 1.0), 6) <= 1e-17 * cum) {
                ok[1] = true;
                return;
            }
        }
        fprev = f;
    }
}

inline void pig1_info_dinfo(int k, double y, const double* th, double* out) {
    (void) y;
    const double m = th[0], x = th[1];
    double E = 0.0, D = 0.0;
    bool ok[2] = {false, false};
    if (R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0) {
        if ((1.0 + m) * x <= PIG_ASYM_CUT) {
            E = (k == 0) ? pig1_asym0_0(m, x) : pig1_asym0_1(m, x);
            D = (k == 0) ? pig1_asym1_0(m, x) : pig1_asym1_3(m, x);
            ok[0] = ok[1] = true;
        } else {
            double w = 0.5 * x / d7::sqrt_cr(1.0 + 2.0 * x * m);
            pig_diag_obs<PigSeqW>(m, w, 2.0 * x * m,
                [&](double yy, const double* pre, double* val, double* g,
                    double* h) {
                    *val = pig1_value(yy, m, x, pre);
                    if (k == 0) {
                        *g = pig1_score_mu(yy, m, x, pre);
                        *h = pig1_hess_mu_mu(yy, m, x, pre);
                    } else {
                        *g = pig1_score_sigma(yy, m, x, pre);
                        *h = pig1_hess_sigma_sigma(yy, m, x, pre);
                    }
                }, &E, &D, ok);
        }
    }
    out[0] = ok[0] ? E : NA_REAL;
    out[1] = ok[1] ? D : NA_REAL;
}

inline void pig2_info_dinfo(int k, double y, const double* th, double* out) {
    (void) y;
    const double m = th[0], x = th[1];
    double E = 0.0, D = 0.0;
    bool ok[2] = {false, false};
    if (R_finite(m) && R_finite(x) && m > 0.0 && x > 0.0) {
        const double xs = 1.0 / x;
        if ((1.0 + m) * xs <= PIG_ASYM_CUT) {
            E = (k == 0) ? pig2_asym0_0(m, xs) : pig2_asym0_1(m, xs);
            D = (k == 0) ? pig2_asym1_0(m, xs) : pig2_asym1_3(m, xs);
            ok[0] = ok[1] = true;
        } else {
            double t = m / x;
            double b = x / (t + std::hypot(1.0, t));
            pig_diag_obs<PigSeq>(m, x, 2.0 * m / b,
                [&](double yy, const double* pre, double* val, double* g,
                    double* h) {
                    *val = pig2_value(yy, m, x, pre[5]);
                    const Pig2Base B = pig2_base(yy, m, x);
                    if (k == 0) {
                        *g = pig2_score_mu(yy, m, B);
                        *h = pig2_hess_mu_mu(yy, B);
                    } else {
                        *g = pig2_score_alpha(yy, m, x, B, pre[6]);
                        *h = pig2_hess_alpha_alpha(yy, B, pre[2]);
                    }
                }, &E, &D, ok);
        }
    }
    out[0] = ok[0] ? E : NA_REAL;
    out[1] = ok[1] ? D : NA_REAL;
}

} // namespace d7

#endif
