// The third and fourth derivatives of the expected information of the Student
// t in its location mu, standard deviation sigma and degrees of freedom nu
// (student_t2), on the parameter scale, in closed form.
//
// The expected Hessian is E_ab = s^p G_ab(nu) with
//   E_mm = -nu(nu+1)/(s^2 (nu-2)(nu+3)),   E_ss = -2 nu/(s^2 (nu+3)),
//   E_sn = -6/(s (nu-2)(nu+1)(nu+3)),      E_nn = V2(nu),  E_m,s = E_m,nu = 0,
// and it does not depend on mu. The derivative in sigma multiplies by the
// falling factorial of p and lowers the power of s by one per order, so the
// derivative in nu acts on the shape function alone: G^(k)(nu). The shape
// functions are rational (partial factored forms, differentiated symbolically
// with sympy in misure_d34e/agents/student_t2/derive.py), except V2, which
// carries the trigamma difference:
//   V2(nu) = (psi'(a) - psi'(b))/4 + (nu+4)(nu-3)/(2 (nu-2)^2 (nu+1)(nu+3)),
//   a = (nu+1)/2, b = nu/2, with d^k/dnu^k psi'(a) = 2^-k psi^(k+1)(a).
// Above nu = 20 the trigamma difference cancels, and V2^(k) is taken from the
// series of V2 in h = 1/nu (src/pt_student_t2.h), differentiated term by term:
// d/dnu h^m = -m h^(m+1), so V2^(k) = h^(4+k) sum_j (-1)^k prod_i(4+j+i) c_j h^j.
// Each order is its own kernel. The components are in the order of
// d3expected_names() and d4expected_names() for student_t2.

#include <Rcpp.h>
#include <cmath>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {

inline double t2e_gmm_0(double nu) {
    return -nu*(nu + 1)/((nu - 2)*(nu + 3));
}

inline double t2e_gmm_1(double nu) {
    return 6*(2*nu + 1)/(std::pow(nu - 2, 2)*std::pow(nu + 3, 2));
}

inline double t2e_gmm_2(double nu) {
    return -12*(3*std::pow(nu, 2) + 3*nu + 7)/(std::pow(nu - 2, 3)*std::pow(nu + 3, 3));
}

inline double t2e_gmm_3(double nu) {
    return 36*(2*nu + 1)*(2*std::pow(nu, 2) + 2*nu + 13)/(std::pow(nu - 2, 4)*std::pow(nu + 3, 4));
}

inline double t2e_gmm_4(double nu) {
    return -720*(std::pow(nu, 4) + 2*std::pow(nu, 3) + 14*std::pow(nu, 2) + 13*nu + 11)/(std::pow(nu - 2, 5)*std::pow(nu + 3, 5));
}

inline double t2e_gss_0(double nu) {
    return -2*nu/(nu + 3);
}

inline double t2e_gss_1(double nu) {
    return -6/std::pow(nu + 3, 2);
}

inline double t2e_gss_2(double nu) {
    return 12/std::pow(nu + 3, 3);
}

inline double t2e_gss_3(double nu) {
    return -36/std::pow(nu + 3, 4);
}

inline double t2e_gss_4(double nu) {
    return 144/std::pow(nu + 3, 5);
}

inline double t2e_gsn_0(double nu) {
    return -6/((nu - 2)*(nu + 1)*(nu + 3));
}

inline double t2e_gsn_1(double nu) {
    return 6*(3*std::pow(nu, 2) + 4*nu - 5)/(std::pow(nu - 2, 2)*std::pow(nu + 1, 2)*std::pow(nu + 3, 2));
}

inline double t2e_gsn_2(double nu) {
    return -12*(6*std::pow(nu, 4) + 16*std::pow(nu, 3) - 3*std::pow(nu, 2) - 12*nu + 37)/(std::pow(nu - 2, 3)*std::pow(nu + 1, 3)*std::pow(nu + 3, 3));
}

inline double t2e_gsn_3(double nu) {
    return 36*(10*std::pow(nu, 6) + 40*std::pow(nu, 5) + 35*std::pow(nu, 4) + 8*std::pow(nu, 3) + 172*std::pow(nu, 2) + 176*nu - 209)/(std::pow(nu - 2, 4)*std::pow(nu + 1, 4)*std::pow(nu + 3, 4));
}

inline double t2e_gsn_4(double nu) {
    return -144*(15*std::pow(nu, 8) + 80*std::pow(nu, 7) + 155*std::pow(nu, 6) + 228*std::pow(nu, 5) + 750*std::pow(nu, 4) + 1200*std::pow(nu, 3) - 405*std::pow(nu, 2) - 980*nu + 1309)/(std::pow(nu - 2, 5)*std::pow(nu + 1, 5)*std::pow(nu + 3, 5));
}

static inline double t2e_V2_3(double nu) {
    if (nu >= 20.0) {
        // h^(4 + k) times the term-wise derivative of the nu >= 20 series
        const double h = 1.0 / nu;
        static const double e[39] = {180.0, -630.0, 7224.0, -10584.0, 133560.0, -264330.0, 1981980.0, -4669236.0, 26910156.0, -79516710.0, 349673520.0, -1017572400.0, 4411930032.0, -17213178402.0, 54278173620.0, -67309009740.0, 652238012580.0, -8099381332350.0, 7668352960488.0, 332536644747000.0, 88383229607400.0, -2.578947252364665e+16, 1000726883331084.0, 2.0940156209964078e+18, 1.115350737713478e+16, -1.9766811963028708e+20, 1.2258836880716784e+17, 2.1145170820433506e+22, 1.3308346824055757e+18, -2.5495730594710894e+24, 1.4290123071809343e+19, 3.440892767176004e+26, 1.5194941286931376e+20, -5.166776310138384e+28, 1.601597714832547e+21, 8.585463731040045e+30, 1.6748537268037963e+22, -1.5710201622578883e+33, -1.7911176357866373e+18};
        double acc = 0.0;
        for (int j = 39 - 1; j >= 0; --j) acc = acc * h + e[j];
        return std::pow(h, 7) * acc;
    }
    const double a = 0.5 * (nu + 1.0), b = 0.5 * nu;
    return 1.0/32.0 * (R::psigamma(a, 4.0) - R::psigamma(b, 4.0))
        + -3*(4*std::pow(nu, 8) + 34*std::pow(nu, 7) - 48*std::pow(nu, 6) - 785*std::pow(nu, 5) - 1814*std::pow(nu, 4) - 1716*std::pow(nu, 3) - 2664*std::pow(nu, 2) - 3141*nu + 1506)/(std::pow(nu - 2, 5)*std::pow(nu + 1, 4)*std::pow(nu + 3, 4));
}

static inline double t2e_V2_4(double nu) {
    if (nu >= 20.0) {
        // h^(4 + k) times the term-wise derivative of the nu >= 20 series
        const double h = 1.0 / nu;
        static const double e[39] = {-1260.0, 5040.0, -65016.0, 105840.0, -1469160.0, 3171960.0, -25765740.0, 65369304.0, -403652340.0, 1272267360.0, -5944449840.0, 18316303200.0, -83826670608.0, 344263568040.0, -1139841646020.0, 1480798214280.0, -15001474289340.0, 194385151976400.0, -191708824012200.0, -8645952763422000.0, -2386347199399800.0, 7.221052306621062e+17, -2.9021079616601436e+16, -6.2820468629892235e+19, -3.457587286911782e+17, 6.325379828169186e+21, -4.045416170636539e+18, -7.189358078947393e+23, -4.657921388419515e+19, 9.178463014095921e+25, -5.2873455365694574e+20, -1.3075392515268816e+28, -5.926027101903236e+21, 2.0667105240553535e+30, -6.566550630813443e+22, -3.605894767036819e+32, -7.201871025256323e+23, 6.912488713934709e+34, 8.060029361039868e+19};
        double acc = 0.0;
        for (int j = 39 - 1; j >= 0; --j) acc = acc * h + e[j];
        return std::pow(h, 8) * acc;
    }
    const double a = 0.5 * (nu + 1.0), b = 0.5 * nu;
    return 1.0/64.0 * (R::psigamma(a, 5.0) - R::psigamma(b, 5.0))
        + 12*(5*std::pow(nu, 10) + 55*std::pow(nu, 9) - 10*std::pow(nu, 8) - 1465*std::pow(nu, 7) - 5843*std::pow(nu, 6) - 11734*std::pow(nu, 5) - 20580*std::pow(nu, 4) - 30105*std::pow(nu, 3) - 12300*std::pow(nu, 2) + 8961*nu - 11112)/(std::pow(nu - 2, 6)*std::pow(nu + 1, 5)*std::pow(nu + 3, 5));
}

}  // namespace

// [[Rcpp::export]]
List student_t2_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
             NumericVector nu, int threads = 1) {
    (void) mu;
    d34e::Par S(sigma), V(nu);
    return d34e::run(y.size(), threads, d7::kMinCostly, 60, [&](std::size_t i, d34e::Row o) {
    const double s = S[i], v = V[i];
    const double si = 1.0 / s;
    const double V2_3 = t2e_V2_3(v);
    const double gmm_0 = t2e_gmm_0(v);
    const double gmm_1 = t2e_gmm_1(v);
    const double gmm_2 = t2e_gmm_2(v);
    const double gmm_3 = t2e_gmm_3(v);
    const double gsn_0 = t2e_gsn_0(v);
    const double gsn_1 = t2e_gsn_1(v);
    const double gsn_2 = t2e_gsn_2(v);
    const double gsn_3 = t2e_gsn_3(v);
    const double gss_0 = t2e_gss_0(v);
    const double gss_1 = t2e_gss_1(v);
    const double gss_2 = t2e_gss_2(v);
    const double gss_3 = t2e_gss_3(v);
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
List student_t2_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
             NumericVector nu, int threads = 1) {
    (void) mu;
    d34e::Par S(sigma), V(nu);
    return d34e::run(y.size(), threads, d7::kMinCostly, 90, [&](std::size_t i, d34e::Row o) {
    const double s = S[i], v = V[i];
    const double si = 1.0 / s;
    const double V2_4 = t2e_V2_4(v);
    const double gmm_0 = t2e_gmm_0(v);
    const double gmm_1 = t2e_gmm_1(v);
    const double gmm_2 = t2e_gmm_2(v);
    const double gmm_3 = t2e_gmm_3(v);
    const double gmm_4 = t2e_gmm_4(v);
    const double gsn_0 = t2e_gsn_0(v);
    const double gsn_1 = t2e_gsn_1(v);
    const double gsn_2 = t2e_gsn_2(v);
    const double gsn_3 = t2e_gsn_3(v);
    const double gsn_4 = t2e_gsn_4(v);
    const double gss_0 = t2e_gss_0(v);
    const double gss_1 = t2e_gss_1(v);
    const double gss_2 = t2e_gss_2(v);
    const double gss_3 = t2e_gss_3(v);
    const double gss_4 = t2e_gss_4(v);
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
