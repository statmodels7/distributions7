// The third and fourth derivatives of the expected information of the Student
// t in its location mu, scale sigma and degrees of freedom nu (student_t1), on
// the parameter scale, in closed form.
//
// The expected Hessian is E_ab = s^p G_ab(nu) with
//   E_mm = -(nu+1)/(s^2 (nu+3)),   E_ss = -2 nu/(s^2 (nu+3)),
//   E_sn = 2/(s (nu+1)(nu+3)),     E_nn = V2(nu),  E_m,s = E_m,nu = 0,
// and it does not depend on mu. The sigma derivatives multiply by the falling
// factorial of p and lower the power of s by one per order; the nu derivatives
// act on the shape function G^(k)(nu). The shape functions Gmm, Gss, Gsn are
// rational, differentiated symbolically in factored form (misure_d34e/agents/
// student_t1/gen_t1.py). V2 carries the trigamma difference:
//   V2(nu) = (psi'(a) - psi'(b))/4 + (nu+5)/(2 nu (nu+1)(nu+3)),
//   a = (nu+1)/2, b = nu/2, with d^k/dnu^k psi'(a) = 2^-k psi^(k+1)(a).
// Above nu = 20 the trigamma difference cancels, and V2^(k) comes from the
// series of V2 in h = 1/nu (src/pt_student_t1.h), differentiated term by term:
// d/dnu h^m = -m h^(m+1). Each order is its own kernel. The components are in
// the order of d3expected_names() and d4expected_names() for student_t1.

#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {

inline double t1e_gmm_0(double nu) {
    return -(nu + 1)/(nu + 3);
}

inline double t1e_gmm_1(double nu) {
    return -2/std::pow(nu + 3, 2);
}

inline double t1e_gmm_2(double nu) {
    return 4/std::pow(nu + 3, 3);
}

inline double t1e_gmm_3(double nu) {
    return -12/std::pow(nu + 3, 4);
}

inline double t1e_gmm_4(double nu) {
    return 48/std::pow(nu + 3, 5);
}

inline double t1e_gss_0(double nu) {
    return -2*nu/(nu + 3);
}

inline double t1e_gss_1(double nu) {
    return -6/std::pow(nu + 3, 2);
}

inline double t1e_gss_2(double nu) {
    return 12/std::pow(nu + 3, 3);
}

inline double t1e_gss_3(double nu) {
    return -36/std::pow(nu + 3, 4);
}

inline double t1e_gss_4(double nu) {
    return 144/std::pow(nu + 3, 5);
}

inline double t1e_gsn_0(double nu) {
    return 2/((nu + 1)*(nu + 3));
}

inline double t1e_gsn_1(double nu) {
    return -4*(nu + 2)/(std::pow(nu + 1, 2)*std::pow(nu + 3, 2));
}

inline double t1e_gsn_2(double nu) {
    return 4*(3*std::pow(nu, 2) + 12*nu + 13)/(std::pow(nu + 1, 3)*std::pow(nu + 3, 3));
}

inline double t1e_gsn_3(double nu) {
    return -48*(nu + 2)*(std::pow(nu, 2) + 4*nu + 5)/(std::pow(nu + 1, 4)*std::pow(nu + 3, 4));
}

inline double t1e_gsn_4(double nu) {
    return 48*(5*std::pow(nu, 4) + 40*std::pow(nu, 3) + 130*std::pow(nu, 2) + 200*nu + 121)/(std::pow(nu + 1, 5)*std::pow(nu + 3, 5));
}

static inline double t1e_V2_3(double nu) {
    if (nu >= 20.0) {
        // h^(4 + k) times the term-wise derivative of the nu >= 20 series
        const double h = 1.0 / nu;
        static const double e[38] = {420.0, -2730.0, 13272.0, -59976.0, 261720.0, -1089990.0, 4328940.0, -16753308.0, 64479324.0, -244632570.0, 892817520.0, -3174431760.0, 11708703216.0, -44414523918.0, 147219778980.0, -400274972580.0, 1789882649940.0, -12070681047810.0, 21171754870728.0, 286232206857000.0, 244772264934600.0, -2.631880679196015e+16, 2775717484515612.0, 2.0880678920325102e+18, 3.095992578902286e+16, -1.977339498791409e+20, 3.4040667172978995e+17, 2.114445178570045e+22, 3.6961576033631027e+18, -2.549580821895678e+24, 3.969169244520639e+19, 3.44089193782807e+26, 4.220661144551615e+20, -5.166776397929389e+28, 4.448804990203053e+21, 8.585463721823933e+30, 4.6523333156657846e+22, -1.5710201622578877e+33};
        double acc = 0.0;
        for (int j = 38 - 1; j >= 0; --j) acc = acc * h + e[j];
        return std::pow(h, 7) * acc;
    }
    const double a = 0.5 * (nu + 1.0), b = 0.5 * nu;
    return 1.0/32.0 * (R::psigamma(a, 4.0) - R::psigamma(b, 4.0))
        + -3*(4*std::pow(nu, 7) + 74*std::pow(nu, 6) + 452*std::pow(nu, 5) + 1315*std::pow(nu, 4) + 2000*std::pow(nu, 3) + 1620*std::pow(nu, 2) + 720*nu + 135)/(std::pow(nu, 4)*std::pow(nu + 1, 4)*std::pow(nu + 3, 4));
}

static inline double t1e_V2_4(double nu) {
    if (nu >= 20.0) {
        // h^(4 + k) times the term-wise derivative of the nu >= 20 series
        const double h = 1.0 / nu;
        static const double e[38] = {-2940.0, 21840.0, -119448.0, 599760.0, -2878920.0, 13079880.0, -56276220.0, 234546312.0, -967189860.0, 3914121120.0, -15177897840.0, 57139771680.0, -222465361104.0, 888290478360.0, -3091615358580.0, 8806049396760.0, -41167300948620.0, 289696345147440.0, -529293871768200.0, -7442037378282000.0, -6608851153234200.0, 7.369265901748842e+17, -8.049580705095275e+16, -6.26420367609753e+19, -9.597576994597087e+17, 6.327486396132509e+21, -1.1233420167083067e+19, -7.189113607138152e+23, -1.293655161177086e+20, 9.178490958824441e+25, -1.4685926204726365e+21, -1.3075389363746668e+28, -1.6460578463751298e+22, 2.0667105591717558e+30, -1.8240100459832517e+23, -3.605894763166052e+32, -2.0005033257362872e+24, 6.912488713934706e+34};
        double acc = 0.0;
        for (int j = 38 - 1; j >= 0; --j) acc = acc * h + e[j];
        return std::pow(h, 8) * acc;
    }
    const double a = 0.5 * (nu + 1.0), b = 0.5 * nu;
    return 1.0/64.0 * (R::psigamma(a, 5.0) - R::psigamma(b, 5.0))
        + 12*(5*std::pow(nu, 9) + 115*std::pow(nu, 8) + 930*std::pow(nu, 7) + 3875*std::pow(nu, 6) + 9421*std::pow(nu, 5) + 14050*std::pow(nu, 4) + 13200*std::pow(nu, 3) + 7875*std::pow(nu, 2) + 2700*nu + 405)/(std::pow(nu, 5)*std::pow(nu + 1, 5)*std::pow(nu + 3, 5));
}

}  // namespace

// [[Rcpp::export]]
List student_t1_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
             NumericVector nu, int threads = 1) {
    (void) mu;
    d34e::Par S(sigma), V(nu);
    return d34e::run(y.size(), threads, d7::kMinCostly, 60, [&](std::size_t i, d34e::Row o) {
    const double s = S[i], v = V[i];
    const double si = 1.0 / s;
    const double V2_3 = t1e_V2_3(v);
    const double gmm_0 = t1e_gmm_0(v);
    const double gmm_1 = t1e_gmm_1(v);
    const double gmm_2 = t1e_gmm_2(v);
    const double gmm_3 = t1e_gmm_3(v);
    const double gsn_0 = t1e_gsn_0(v);
    const double gsn_1 = t1e_gsn_1(v);
    const double gsn_2 = t1e_gsn_2(v);
    const double gsn_3 = t1e_gsn_3(v);
    const double gss_0 = t1e_gss_0(v);
    const double gss_1 = t1e_gss_1(v);
    const double gss_2 = t1e_gss_2(v);
    const double gss_3 = t1e_gss_3(v);
    (void) si;
    o[0] = 0.0;  // mu_mu_mu_mu_mu
    o[1] = 0.0;  // mu_mu_mu_mu_sigma
    o[2] = 0.0;  // mu_mu_mu_mu_nu
    o[3] = 0.0;  // mu_mu_mu_sigma_sigma
    o[4] = 0.0;  // mu_mu_mu_sigma_nu
    o[5] = 0.0;  // mu_mu_mu_nu_nu
    o[6] = -24.0*si*si*si*si*si*gmm_0;  // mu_mu_sigma_sigma_sigma
    o[7] = 6.0*si*si*si*si*gmm_1;  // mu_mu_sigma_sigma_nu
    o[8] = -2.0*si*si*si*gmm_2;  // mu_mu_sigma_nu_nu
    o[9] = 1.0*si*si*gmm_3;  // mu_mu_nu_nu_nu
    o[10] = 0.0;  // sigma_sigma_mu_mu_mu
    o[11] = 0.0;  // sigma_sigma_mu_mu_sigma
    o[12] = 0.0;  // sigma_sigma_mu_mu_nu
    o[13] = 0.0;  // sigma_sigma_mu_sigma_sigma
    o[14] = 0.0;  // sigma_sigma_mu_sigma_nu
    o[15] = 0.0;  // sigma_sigma_mu_nu_nu
    o[16] = -24.0*si*si*si*si*si*gss_0;  // sigma_sigma_sigma_sigma_sigma
    o[17] = 6.0*si*si*si*si*gss_1;  // sigma_sigma_sigma_sigma_nu
    o[18] = -2.0*si*si*si*gss_2;  // sigma_sigma_sigma_nu_nu
    o[19] = 1.0*si*si*gss_3;  // sigma_sigma_nu_nu_nu
    o[20] = 0.0;  // nu_nu_mu_mu_mu
    o[21] = 0.0;  // nu_nu_mu_mu_sigma
    o[22] = 0.0;  // nu_nu_mu_mu_nu
    o[23] = 0.0;  // nu_nu_mu_sigma_sigma
    o[24] = 0.0;  // nu_nu_mu_sigma_nu
    o[25] = 0.0;  // nu_nu_mu_nu_nu
    o[26] = 0.0;  // nu_nu_sigma_sigma_sigma
    o[27] = 0.0;  // nu_nu_sigma_sigma_nu
    o[28] = 0.0;  // nu_nu_sigma_nu_nu
    o[29] = V2_3;  // nu_nu_nu_nu_nu
    o[30] = 0.0;  // mu_sigma_mu_mu_mu
    o[31] = 0.0;  // mu_sigma_mu_mu_sigma
    o[32] = 0.0;  // mu_sigma_mu_mu_nu
    o[33] = 0.0;  // mu_sigma_mu_sigma_sigma
    o[34] = 0.0;  // mu_sigma_mu_sigma_nu
    o[35] = 0.0;  // mu_sigma_mu_nu_nu
    o[36] = 0.0;  // mu_sigma_sigma_sigma_sigma
    o[37] = 0.0;  // mu_sigma_sigma_sigma_nu
    o[38] = 0.0;  // mu_sigma_sigma_nu_nu
    o[39] = 0.0;  // mu_sigma_nu_nu_nu
    o[40] = 0.0;  // mu_nu_mu_mu_mu
    o[41] = 0.0;  // mu_nu_mu_mu_sigma
    o[42] = 0.0;  // mu_nu_mu_mu_nu
    o[43] = 0.0;  // mu_nu_mu_sigma_sigma
    o[44] = 0.0;  // mu_nu_mu_sigma_nu
    o[45] = 0.0;  // mu_nu_mu_nu_nu
    o[46] = 0.0;  // mu_nu_sigma_sigma_sigma
    o[47] = 0.0;  // mu_nu_sigma_sigma_nu
    o[48] = 0.0;  // mu_nu_sigma_nu_nu
    o[49] = 0.0;  // mu_nu_nu_nu_nu
    o[50] = 0.0;  // sigma_nu_mu_mu_mu
    o[51] = 0.0;  // sigma_nu_mu_mu_sigma
    o[52] = 0.0;  // sigma_nu_mu_mu_nu
    o[53] = 0.0;  // sigma_nu_mu_sigma_sigma
    o[54] = 0.0;  // sigma_nu_mu_sigma_nu
    o[55] = 0.0;  // sigma_nu_mu_nu_nu
    o[56] = -6.0*si*si*si*si*gsn_0;  // sigma_nu_sigma_sigma_sigma
    o[57] = 2.0*si*si*si*gsn_1;  // sigma_nu_sigma_sigma_nu
    o[58] = -1.0*si*si*gsn_2;  // sigma_nu_sigma_nu_nu
    o[59] = 1.0*si*gsn_3;  // sigma_nu_nu_nu_nu
    });
}

// [[Rcpp::export]]
List student_t1_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
             NumericVector nu, int threads = 1) {
    (void) mu;
    d34e::Par S(sigma), V(nu);
    return d34e::run(y.size(), threads, d7::kMinCostly, 90, [&](std::size_t i, d34e::Row o) {
    const double s = S[i], v = V[i];
    const double si = 1.0 / s;
    const double V2_4 = t1e_V2_4(v);
    const double gmm_0 = t1e_gmm_0(v);
    const double gmm_1 = t1e_gmm_1(v);
    const double gmm_2 = t1e_gmm_2(v);
    const double gmm_3 = t1e_gmm_3(v);
    const double gmm_4 = t1e_gmm_4(v);
    const double gsn_0 = t1e_gsn_0(v);
    const double gsn_1 = t1e_gsn_1(v);
    const double gsn_2 = t1e_gsn_2(v);
    const double gsn_3 = t1e_gsn_3(v);
    const double gsn_4 = t1e_gsn_4(v);
    const double gss_0 = t1e_gss_0(v);
    const double gss_1 = t1e_gss_1(v);
    const double gss_2 = t1e_gss_2(v);
    const double gss_3 = t1e_gss_3(v);
    const double gss_4 = t1e_gss_4(v);
    (void) si;
    o[0] = 0.0;  // mu_mu_mu_mu_mu_mu
    o[1] = 0.0;  // mu_mu_mu_mu_mu_sigma
    o[2] = 0.0;  // mu_mu_mu_mu_mu_nu
    o[3] = 0.0;  // mu_mu_mu_mu_sigma_sigma
    o[4] = 0.0;  // mu_mu_mu_mu_sigma_nu
    o[5] = 0.0;  // mu_mu_mu_mu_nu_nu
    o[6] = 0.0;  // mu_mu_mu_sigma_sigma_sigma
    o[7] = 0.0;  // mu_mu_mu_sigma_sigma_nu
    o[8] = 0.0;  // mu_mu_mu_sigma_nu_nu
    o[9] = 0.0;  // mu_mu_mu_nu_nu_nu
    o[10] = 120.0*si*si*si*si*si*si*gmm_0;  // mu_mu_sigma_sigma_sigma_sigma
    o[11] = -24.0*si*si*si*si*si*gmm_1;  // mu_mu_sigma_sigma_sigma_nu
    o[12] = 6.0*si*si*si*si*gmm_2;  // mu_mu_sigma_sigma_nu_nu
    o[13] = -2.0*si*si*si*gmm_3;  // mu_mu_sigma_nu_nu_nu
    o[14] = 1.0*si*si*gmm_4;  // mu_mu_nu_nu_nu_nu
    o[15] = 0.0;  // sigma_sigma_mu_mu_mu_mu
    o[16] = 0.0;  // sigma_sigma_mu_mu_mu_sigma
    o[17] = 0.0;  // sigma_sigma_mu_mu_mu_nu
    o[18] = 0.0;  // sigma_sigma_mu_mu_sigma_sigma
    o[19] = 0.0;  // sigma_sigma_mu_mu_sigma_nu
    o[20] = 0.0;  // sigma_sigma_mu_mu_nu_nu
    o[21] = 0.0;  // sigma_sigma_mu_sigma_sigma_sigma
    o[22] = 0.0;  // sigma_sigma_mu_sigma_sigma_nu
    o[23] = 0.0;  // sigma_sigma_mu_sigma_nu_nu
    o[24] = 0.0;  // sigma_sigma_mu_nu_nu_nu
    o[25] = 120.0*si*si*si*si*si*si*gss_0;  // sigma_sigma_sigma_sigma_sigma_sigma
    o[26] = -24.0*si*si*si*si*si*gss_1;  // sigma_sigma_sigma_sigma_sigma_nu
    o[27] = 6.0*si*si*si*si*gss_2;  // sigma_sigma_sigma_sigma_nu_nu
    o[28] = -2.0*si*si*si*gss_3;  // sigma_sigma_sigma_nu_nu_nu
    o[29] = 1.0*si*si*gss_4;  // sigma_sigma_nu_nu_nu_nu
    o[30] = 0.0;  // nu_nu_mu_mu_mu_mu
    o[31] = 0.0;  // nu_nu_mu_mu_mu_sigma
    o[32] = 0.0;  // nu_nu_mu_mu_mu_nu
    o[33] = 0.0;  // nu_nu_mu_mu_sigma_sigma
    o[34] = 0.0;  // nu_nu_mu_mu_sigma_nu
    o[35] = 0.0;  // nu_nu_mu_mu_nu_nu
    o[36] = 0.0;  // nu_nu_mu_sigma_sigma_sigma
    o[37] = 0.0;  // nu_nu_mu_sigma_sigma_nu
    o[38] = 0.0;  // nu_nu_mu_sigma_nu_nu
    o[39] = 0.0;  // nu_nu_mu_nu_nu_nu
    o[40] = 0.0;  // nu_nu_sigma_sigma_sigma_sigma
    o[41] = 0.0;  // nu_nu_sigma_sigma_sigma_nu
    o[42] = 0.0;  // nu_nu_sigma_sigma_nu_nu
    o[43] = 0.0;  // nu_nu_sigma_nu_nu_nu
    o[44] = V2_4;  // nu_nu_nu_nu_nu_nu
    o[45] = 0.0;  // mu_sigma_mu_mu_mu_mu
    o[46] = 0.0;  // mu_sigma_mu_mu_mu_sigma
    o[47] = 0.0;  // mu_sigma_mu_mu_mu_nu
    o[48] = 0.0;  // mu_sigma_mu_mu_sigma_sigma
    o[49] = 0.0;  // mu_sigma_mu_mu_sigma_nu
    o[50] = 0.0;  // mu_sigma_mu_mu_nu_nu
    o[51] = 0.0;  // mu_sigma_mu_sigma_sigma_sigma
    o[52] = 0.0;  // mu_sigma_mu_sigma_sigma_nu
    o[53] = 0.0;  // mu_sigma_mu_sigma_nu_nu
    o[54] = 0.0;  // mu_sigma_mu_nu_nu_nu
    o[55] = 0.0;  // mu_sigma_sigma_sigma_sigma_sigma
    o[56] = 0.0;  // mu_sigma_sigma_sigma_sigma_nu
    o[57] = 0.0;  // mu_sigma_sigma_sigma_nu_nu
    o[58] = 0.0;  // mu_sigma_sigma_nu_nu_nu
    o[59] = 0.0;  // mu_sigma_nu_nu_nu_nu
    o[60] = 0.0;  // mu_nu_mu_mu_mu_mu
    o[61] = 0.0;  // mu_nu_mu_mu_mu_sigma
    o[62] = 0.0;  // mu_nu_mu_mu_mu_nu
    o[63] = 0.0;  // mu_nu_mu_mu_sigma_sigma
    o[64] = 0.0;  // mu_nu_mu_mu_sigma_nu
    o[65] = 0.0;  // mu_nu_mu_mu_nu_nu
    o[66] = 0.0;  // mu_nu_mu_sigma_sigma_sigma
    o[67] = 0.0;  // mu_nu_mu_sigma_sigma_nu
    o[68] = 0.0;  // mu_nu_mu_sigma_nu_nu
    o[69] = 0.0;  // mu_nu_mu_nu_nu_nu
    o[70] = 0.0;  // mu_nu_sigma_sigma_sigma_sigma
    o[71] = 0.0;  // mu_nu_sigma_sigma_sigma_nu
    o[72] = 0.0;  // mu_nu_sigma_sigma_nu_nu
    o[73] = 0.0;  // mu_nu_sigma_nu_nu_nu
    o[74] = 0.0;  // mu_nu_nu_nu_nu_nu
    o[75] = 0.0;  // sigma_nu_mu_mu_mu_mu
    o[76] = 0.0;  // sigma_nu_mu_mu_mu_sigma
    o[77] = 0.0;  // sigma_nu_mu_mu_mu_nu
    o[78] = 0.0;  // sigma_nu_mu_mu_sigma_sigma
    o[79] = 0.0;  // sigma_nu_mu_mu_sigma_nu
    o[80] = 0.0;  // sigma_nu_mu_mu_nu_nu
    o[81] = 0.0;  // sigma_nu_mu_sigma_sigma_sigma
    o[82] = 0.0;  // sigma_nu_mu_sigma_sigma_nu
    o[83] = 0.0;  // sigma_nu_mu_sigma_nu_nu
    o[84] = 0.0;  // sigma_nu_mu_nu_nu_nu
    o[85] = 24.0*si*si*si*si*si*gsn_0;  // sigma_nu_sigma_sigma_sigma_sigma
    o[86] = -6.0*si*si*si*si*gsn_1;  // sigma_nu_sigma_sigma_sigma_nu
    o[87] = 2.0*si*si*si*gsn_2;  // sigma_nu_sigma_sigma_nu_nu
    o[88] = -1.0*si*si*gsn_3;  // sigma_nu_sigma_nu_nu_nu
    o[89] = 1.0*si*gsn_4;  // sigma_nu_nu_nu_nu_nu
    });
}
