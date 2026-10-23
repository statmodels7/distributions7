#ifndef D7_PT_NEGBIN1_ASYM_H
#define D7_PT_NEGBIN1_ASYM_H

#include <cmath>

// THE NEGATIVE BINOMIAL WITH A LINEAR VARIANCE AT LARGE THETA: the quantities
// that nb1_sums0, nb1_sums1 and nb1_sums2 (pt_negbin1.h) return, A = -E[l_mu,mu]
// and G = E[l_theta,theta] with their derivatives in (mu, theta), computed
// without a sum over the support. The support sums need about theta terms at
// a ratio theta/(1 + theta), so a single row at theta = 1e8 costs seconds;
// the form below costs a few hundred operations at any theta.
//
// DERIVATION. Write r = mu/theta, delta = 1/(1 + theta) and L = log(1 +
// theta) = -log delta. With u_j = mu + theta j, A = E[a_Y] for a_k = sum_{j<k}
// 1/u_j^2, and since E[a_Y] = sum_j P(Y > j)/u_j^2, the same sum gives
//
//   W = E[V_Y] = mu delta - mu^2 A,   G = W/theta^2,
//
// because sum_j P(Y > j) = E[Y] = mu. Everything therefore follows from A.
// By the recurrence psi'(r + k) - psi'(r) = -a_k theta^2,
//
//   theta^2 A = psi'(r) - E[psi'(Y + r)],
//
// and psi'(r + k) = int_0^inf t e^(-(r+k)t)/(1 - e^(-t)) dt together with the
// probability generating function E[e^(-tY)] = (1 + theta (1 - e^(-t)))^(-r)
// give, after x = e^t,
//
//   F(r, L) := theta^2 A
//            = int_1^inf log x/(x - 1) [x^(-r) - delta^r (x - 1 + delta)^(-r)] dx.
//
// Expanding (delta + y)^(-r) in powers of delta/y for y > delta and
// log(1 + y)/y in powers of y for y < delta, and integrating each family
// termwise by the Mellin transforms
//
//   int_0^inf y^(s-1) log(1 + y) dy = pi/(s sin(pi s))   (-1 < s < 0),
//   int_0^inf y^m (delta + y)^(-r) dy = delta^(m+1-r) m! Gamma(r-m-1)/Gamma(r)
//
// (continued analytically), gives for r < 1 the convergent series in delta
//
//   F = B0(r) + (1 - e^(-rL)) q(r)/r^2
//       - q(r) e^(-rL) sum_{n>=1} a_n(r) delta^n
//       + sum_{m>=0} b_m(r) delta^(m+1),
//
//   q(r)   = pi r/sin(pi r) = 1 + 2 sum_{j>=1} eta(2j) r^(2j),
//   a_n(r) = (r+1)_(n-1)/(n! (r + n)),
//   b_m(r) = m!/((m + 1) prod_{i=1}^{m+1} (i - r)),
//   B0(r)  = psi'(1 + r) - (q(r) - 1)/r^2
//          = sum_{k>=1} c_k r^k,  c_k = (-1)^k (k+1) zeta(k+2)
//                                       - [k even] 2 eta(k+2),
//
// where eta is the Dirichlet eta function. The factor e^(-rL) = delta^r is
// the whole dependence on theta that is not a power of delta, so the
// expansion has no other logarithm. B0 is psi'(r) minus the first term of the
// delta series with its 1/r^2 poles removed, written as a power series so
// that no difference of nearly equal quantities is formed at small r;
// 1 - e^(-rL) is -expm1(-rL) for the same reason. The series was checked
// against the integral above at 28 digits and against the support sums
// (nb_ref.json) at theta = 50.
//
// DERIVATIVES. Every term is a product of a function of r and a function of
// L, so the derivatives in (r, L) up to order two are termwise, with the
// logarithmic derivatives of a_n and b_m kept as running sums,
//
//   a_n'/a_n = sum_{i=1}^{n-1} 1/(r + i) - 1/(r + n),
//   b_m'/b_m = sum_{i=1}^{m+1} 1/(i - r),
//
// and are carried to (mu, theta) through r = mu/theta and L = log1p(theta).
// The second derivative in mu of G is the one that cancels: W = mu delta -
// r^2 F, and r^2 F has a leading part of order rL whose second derivative in
// r is of order L^2, while the same quantity formed from A_mu,mu, A_mu and A
// subtracts three terms of order L/r. W is therefore built from
// H = r^2 F = q (1 - e^(-rL)) + r^2 (F - (1 - e^(-rL)) q/r^2), the first
// part with its own derivatives and the second without the 1/r^2 factor.
//
// RANGE. The series in delta converges for every theta > 0, and the series in
// r (B0 and q) for r < 1. The truncation is the one of the tables below
// (40 coefficients of B0, 18 of q) and a loop on n and on m that stops once a
// term times (n + L + 1)^2 is below 1e-18 of the leading part. Against 50
// digit values on a grid of mu from 1e-6 to 1e3 and theta from 2 to 1e8 the
// worst componentwise relative error of the 12 outputs of order two was 6e-15
// (misure_gas_congiunto/c/fase1/nb1asym/), and it stays below 7e-15 up to
// r = 0.3, where the truncation of the r series begins to show; the worst
// error is 7e-14 at r = 0.35, 1e-11 at r = 0.4 and 4e-9 at r = 0.5 (more
// terms of the tables extend the range, and r = 1 is a pole). nb1_asym_ok()
// therefore asks for
// theta >= 2 (the loops need up to 38 terms there, 9 at theta = 100 and 3 at
// 1e8) and r <= 0.3. Beyond theta = 1e50 or below r = 1e-70 the powers of
// 1/r and 1/theta leave the range of a double and the function reports false.
//
// The outputs have the layouts of nb1_sums0/1/2: A and G as the value, then
// (value, d_mu, d_theta), then (value, d_mu, d_theta, d_mu,mu, d_mu,theta,
// d_theta,theta).

// the largest r = mu/theta and the smallest theta at which the expansion
// meets 1e-13 (measured, see above)
inline double nb1_asym_rmax() { return 0.3; }
inline double nb1_asym_thmin() { return 2.0; }

// the relative size of the first omitted term of each loop
inline double nb1_asym_tol() { return 1e-18; }

inline bool nb1_asym_ok(double mu, double th) {
    return th >= nb1_asym_thmin() && th <= 1e50 &&
           mu <= nb1_asym_rmax() * th && mu >= 1e-70 * th;
}

// c_k, k = 1..40: the coefficients of B0(r) = sum c_k r^k
static const double NB1A_CB[40] = {
    -2.4041138063191885, 1.3529040421389227, -4.14771102057348,
    3.1156131273273755, -6.050095664291537, 5.036075489680314,
    -8.016067142608657, 7.01087216095382, -10.004941886041195,
    9.003191581798673, -12.001472560170942, 11.000917885063803,
    -14.000428235308298, 13.000259705461318, -16.000122195162206,
    15.000072518246284, -18.000034347828898, 17.00002003205548,
    -20.000009538659736, 19.000005484234094, -22.00000262283837,
    21.000001490190563, -24.000000715284084, 23.000000402340408,
    -26.000000193718506, 25.000000108034513, -28.000000052154473,
    27.00000002887113, -30.00000001396989, 29.000000007683425,
    -32.0000000037253, 31.00000000203727, -34.00000000098953,
    33.00000000053842, -36.000000000261934, 35.00000000014188,
    -38.00000000006912, 37.00000000003729, -40.00000000001819,
    39.00000000000978};

// 2 eta(2j), j = 1..18: q(r) = 1 + sum_j NB1A_CQ[j-1] r^(2j)
static const double NB1A_CQ[18] = {
    1.6449340668482264, 1.8940656589944918, 1.9711021825948702,
    1.9924660037052957, 1.998079015196543, 1.9995153702877164,
    1.9998783406919594, 1.9999695284298122, 1.9999923757392202,
    1.9999980932231631, 1.9999995232264616, 1.9999998807977848,
    1.999999970198464, 1.9999999925495069, 1.9999999981373646,
    1.9999999995343398, 1.999999999883585, 1.9999999999708962};

// A and G at one (mu, theta)
inline void nb1_asym_sums0(double mu, double th, double* A, double* G) {
    const double r = mu / th, r2 = r * r, L = std::log1p(th);
    const double dl = 1.0 / (1.0 + th);
    // B0 and q, and M0 - 1/r^2 = (q - 1)/r^2
    double B = 0.0, pw = r;
    for (int k = 0; k < 40; ++k) { B += NB1A_CB[k] * pw; pw *= r; }
    double q = 1.0, C0 = 0.0, p = 1.0;
    for (int j = 0; j < 18; ++j) {
        C0 += NB1A_CQ[j] * p;
        q += NB1A_CQ[j] * p * r2;
        p *= r2;
    }
    const double E = std::exp(-r * L), omE = -std::expm1(-r * L);
    const double g = (1.0 / r2 + C0) * omE, hq = q * omE;
    const double lim = nb1_asym_tol() * g;
    // S = sum_n a_n delta^n e^(-rL)
    double S = 0.0, u = 1.0, dn = 1.0;
    for (int n = 1; n <= 80; ++n) {
        const double nd = (double) n;
        dn *= dl;
        if (n > 1) u *= (r + nd - 1.0) / nd;
        const double t = u / (r + nd) * E * dn;
        S += t;
        if (t * (nd + L + 1.0) * (nd + L + 1.0) <= lim) break;
    }
    // T4 = sum_m b_m delta^(m+1)
    double T4 = 0.0, w = 1.0 / (1.0 - r), dm = 1.0;
    for (int m = 0; m < 80; ++m) {
        const double md = (double) m;
        if (m > 0) w *= md / (md + 1.0 - r);
        dm *= dl;
        const double b = w / (md + 1.0) * dm;
        T4 += b;
        if (b * (md + L + 1.0) * (md + L + 1.0) <= lim) break;
    }
    const double R = B - q * S + T4;
    const double F = g + R, H = hq + r2 * R;
    const double i2 = 1.0 / (th * th);
    *A = F * i2;
    *G = (mu * dl - H) * i2;
}

// the same with the first derivatives, each as (value, d_mu, d_theta)
inline void nb1_asym_sums1(double mu, double th, double* A, double* G) {
    const double r = mu / th, r2 = r * r, L = std::log1p(th);
    const double dl = 1.0 / (1.0 + th), it = 1.0 / th;
    double B = 0.0, B1 = 0.0, pw = 1.0;
    for (int k = 1; k <= 40; ++k) {
        const double c = NB1A_CB[k - 1];
        B += c * pw * r;
        B1 += c * (double) k * pw;
        pw *= r;
    }
    double q = 1.0, q1 = 0.0, C0 = 0.0, C01 = 0.0, p = 1.0;
    for (int j = 1; j <= 18; ++j) {
        const double c = NB1A_CQ[j - 1], jd = (double) j;
        C0 += c * p;
        q += c * p * r2;
        q1 += c * 2.0 * jd * p * r;
        if (j >= 2) C01 += c * (2.0 * jd - 2.0) * p / r;
        p *= r2;
    }
    const double M0 = 1.0 / r2 + C0, M01 = -2.0 / (r2 * r) + C01;
    const double E = std::exp(-r * L), omE = -std::expm1(-r * L);
    const double g = M0 * omE, gr = M01 * omE + M0 * L * E, gL = M0 * r * E;
    const double hq = q * omE, hq_r = q1 * omE + q * L * E, hq_L = q * r * E;
    const double lim = nb1_asym_tol() * g;
    double S = 0.0, Sr = 0.0, SL = 0.0, u = 1.0, dn = 1.0, ls = 0.0;
    for (int n = 1; n <= 80; ++n) {
        const double nd = (double) n;
        dn *= dl;
        if (n > 1) {
            u *= (r + nd - 1.0) / nd;
            ls += 1.0 / (r + nd - 1.0);
        }
        const double rn = r + nd;
        const double t = u / rn * E * dn;
        const double lam = ls - 1.0 / rn;
        S += t;
        SL -= rn * t;
        Sr += (lam - L) * t;
        if (t * (nd + L + 1.0) * (nd + L + 1.0) <= lim) break;
    }
    double T4 = 0.0, T4r = 0.0, T4L = 0.0;
    double w = 1.0 / (1.0 - r), dm = 1.0, kap = 1.0 / (1.0 - r);
    for (int m = 0; m < 80; ++m) {
        const double md = (double) m;
        if (m > 0) {
            w *= md / (md + 1.0 - r);
            kap += 1.0 / (md + 1.0 - r);
        }
        dm *= dl;
        const double b = w / (md + 1.0) * dm;
        T4 += b;
        T4r += kap * b;
        T4L -= (md + 1.0) * b;
        if (b * (md + L + 1.0) * (md + L + 1.0) <= lim) break;
    }
    // R = F - (1 - e^(-rL)) q/r^2, the part without the 1/r^2 factor
    const double R = B - q * S + T4;
    const double Rr = B1 - (q1 * S + q * Sr) + T4r;
    const double RL = -q * SL + T4L;
    const double F = g + R, Fr = gr + Rr, FL = gL + RL;
    // H = r^2 F
    const double H = hq + r2 * R;
    const double Hr = hq_r + 2.0 * r * R + r2 * Rr;
    const double HL = hq_L + r2 * RL;
    const double i2 = it * it, i3 = i2 * it;
    A[0] = F * i2;
    A[1] = Fr * it * i2;
    A[2] = (-Fr * mu * i2 + FL * dl) * i2 - 2.0 * F * i3;
    const double W0 = mu * dl - H;
    const double W1 = dl - Hr * it;
    const double W2 = -mu * dl * dl - (-Hr * mu * i2 + HL * dl);
    G[0] = W0 * i2;
    G[1] = W1 * i2;
    G[2] = W2 * i2 - 2.0 * W0 * i3;
}

// the same with the first and second derivatives, each as (value, d_mu,
// d_theta, d_mu,mu, d_mu,theta, d_theta,theta)
inline void nb1_asym_sums2(double mu, double th, double* A, double* G) {
    const double r = mu / th, r2 = r * r, L = std::log1p(th);
    const double dl = 1.0 / (1.0 + th), it = 1.0 / th;
    double B = 0.0, B1 = 0.0, B2 = 0.0, pw = 1.0;
    for (int k = 1; k <= 40; ++k) {
        const double c = NB1A_CB[k - 1], kd = (double) k;
        B += c * pw * r;
        B1 += c * kd * pw;
        if (k >= 2) B2 += c * kd * (kd - 1.0) * pw / r;
        pw *= r;
    }
    double q = 1.0, q1 = 0.0, q2 = 0.0, C0 = 0.0, C01 = 0.0, C02 = 0.0, p = 1.0;
    for (int j = 1; j <= 18; ++j) {
        const double c = NB1A_CQ[j - 1], jd = (double) j;
        C0 += c * p;
        q += c * p * r2;
        q1 += c * 2.0 * jd * p * r;
        q2 += c * 2.0 * jd * (2.0 * jd - 1.0) * p;
        if (j >= 2) {
            C01 += c * (2.0 * jd - 2.0) * p / r;
            C02 += c * (2.0 * jd - 2.0) * (2.0 * jd - 3.0) * p / r2;
        }
        p *= r2;
    }
    const double M0 = 1.0 / r2 + C0, M01 = -2.0 / (r2 * r) + C01,
                 M02 = 6.0 / (r2 * r2) + C02;
    const double E = std::exp(-r * L), omE = -std::expm1(-r * L);
    // g = M0 (1 - E) and hq = q (1 - E), with their derivatives in (r, L)
    const double g = M0 * omE, gr = M01 * omE + M0 * L * E, gL = M0 * r * E;
    const double grr = M02 * omE + 2.0 * M01 * L * E - M0 * L * L * E;
    const double grL = E * (r * M01 + M0 * (1.0 - r * L));
    const double gLL = -M0 * r2 * E;
    const double hq = q * omE, hq_r = q1 * omE + q * L * E, hq_L = q * r * E;
    const double hq_rr = q2 * omE + 2.0 * q1 * L * E - q * L * L * E;
    const double hq_rL = E * (r * q1 + q * (1.0 - r * L));
    const double hq_LL = -q * r2 * E;
    const double lim = nb1_asym_tol() * g;
    double S = 0.0, Sr = 0.0, SL = 0.0, Srr = 0.0, SrL = 0.0, SLL = 0.0;
    double u = 1.0, dn = 1.0, ls = 0.0, ls2 = 0.0;
    for (int n = 1; n <= 80; ++n) {
        const double nd = (double) n;
        dn *= dl;
        if (n > 1) {
            u *= (r + nd - 1.0) / nd;
            ls += 1.0 / (r + nd - 1.0);
            ls2 += 1.0 / ((r + nd - 1.0) * (r + nd - 1.0));
        }
        const double rn = r + nd;
        const double t = u / rn * E * dn;
        const double lam = ls - 1.0 / rn, lamp = -ls2 + 1.0 / (rn * rn);
        const double lm = lam - L;
        S += t;
        SL -= rn * t;
        SLL += rn * rn * t;
        Sr += lm * t;
        Srr += (lm * lm + lamp) * t;
        SrL -= (1.0 + rn * lm) * t;
        if (t * (nd + L + 1.0) * (nd + L + 1.0) <= lim) break;
    }
    double T4 = 0.0, T4r = 0.0, T4L = 0.0, T4rr = 0.0, T4rL = 0.0, T4LL = 0.0;
    double w = 1.0 / (1.0 - r), dm = 1.0, kap = 1.0 / (1.0 - r);
    double kap2 = kap * kap;
    for (int m = 0; m < 80; ++m) {
        const double md = (double) m;
        if (m > 0) {
            w *= md / (md + 1.0 - r);
            kap += 1.0 / (md + 1.0 - r);
            kap2 += 1.0 / ((md + 1.0 - r) * (md + 1.0 - r));
        }
        dm *= dl;
        const double b = w / (md + 1.0) * dm;
        T4 += b;
        T4L -= (md + 1.0) * b;
        T4LL += (md + 1.0) * (md + 1.0) * b;
        T4r += kap * b;
        T4rr += (kap * kap + kap2) * b;
        T4rL -= (md + 1.0) * kap * b;
        if (b * (md + L + 1.0) * (md + L + 1.0) <= lim) break;
    }
    // R = F - (1 - e^(-rL)) q/r^2 and its derivatives in (r, L)
    const double R = B - q * S + T4;
    const double Rr = B1 - (q1 * S + q * Sr) + T4r;
    const double RL = -q * SL + T4L;
    const double Rrr = B2 - (q2 * S + 2.0 * q1 * Sr + q * Srr) + T4rr;
    const double RrL = -(q1 * SL + q * SrL) + T4rL;
    const double RLL = -q * SLL + T4LL;
    const double F = g + R, Fr = gr + Rr, FL = gL + RL;
    const double Frr = grr + Rrr, FrL = grL + RrL, FLL = gLL + RLL;
    // H = r^2 F = q (1 - E) + r^2 R
    const double H0 = hq + r2 * R;
    const double Hr = hq_r + 2.0 * r * R + r2 * Rr;
    const double HL = hq_L + r2 * RL;
    const double Hrr = hq_rr + 2.0 * R + 4.0 * r * Rr + r2 * Rrr;
    const double HrL = hq_rL + 2.0 * r * RL + r2 * RrL;
    const double HLL = hq_LL + r2 * RLL;
    // (r, L) -> (mu, theta): r_mu = 1/theta, r_theta = -mu/theta^2, L_theta = delta
    const double i2 = it * it, i3 = i2 * it, i4 = i2 * i2;
    const double Fm = Fr * it, Ft = -Fr * mu * i2 + FL * dl;
    const double Fmm = Frr * i2, Fmt = (-Frr * mu * i2 + FrL * dl) * it - Fr * i2;
    const double Ftt = Frr * mu * mu * i4 - 2.0 * FrL * mu * i2 * dl
                       + FLL * dl * dl + 2.0 * Fr * mu * i3 - FL * dl * dl;
    A[0] = F * i2;
    A[1] = Fm * i2;
    A[2] = Ft * i2 - 2.0 * F * i3;
    A[3] = Fmm * i2;
    A[4] = Fmt * i2 - 2.0 * Fm * i3;
    A[5] = Ftt * i2 - 4.0 * Ft * i3 + 6.0 * F * i4;
    const double Hm = Hr * it, Ht = -Hr * mu * i2 + HL * dl;
    const double Hmm = Hrr * i2, Hmt = (-Hrr * mu * i2 + HrL * dl) * it - Hr * i2;
    const double Htt = Hrr * mu * mu * i4 - 2.0 * HrL * mu * i2 * dl
                       + HLL * dl * dl + 2.0 * Hr * mu * i3 - HL * dl * dl;
    // W = mu delta - H and G = W/theta^2
    const double W0 = mu * dl - H0, W1 = dl - Hm, W2 = -mu * dl * dl - Ht;
    const double W3 = -Hmm, W4 = -dl * dl - Hmt, W5 = 2.0 * mu * dl * dl * dl - Htt;
    G[0] = W0 * i2;
    G[1] = W1 * i2;
    G[2] = W2 * i2 - 2.0 * W0 * i3;
    G[3] = W3 * i2;
    G[4] = W4 * i2 - 2.0 * W1 * i3;
    G[5] = W5 * i2 - 4.0 * W2 * i3 + 6.0 * W0 * i4;
}

#endif
