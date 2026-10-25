// The third and fourth derivatives of the expected information of the
// beta-binomial in its mean proportion and dispersion (mu, sigma), on the
// parameter scale, for a fixed size N. The expected information is
//   E_ab = sum_{j=0}^{N} f_j l_ab(j),
// with f_j the mass and l_ab the second derivative of the log-mass in
// (mu, sigma), and the support is finite, so the derivatives of order T are
// exact sums, taken under the sum:
//   d^T E_ab / dmu^p dsigma^q = sum_j f_j sum_{i<=p, k<=q} C(p,i) C(q,k)
//                                 l_ab^(i,k)(j) Q^(p-i,q-k)(j),
// where l^(i,k) is the derivative of the log-mass with i differentiations in
// mu and k in sigma, and Q^(p,q) = f^(p,q) / f is the derivative of the mass
// relative to the mass. The form used is E[l_ab], whose highest derivative of
// the log-mass is of order T + 2 (5 at order 3, 6 at order 4).
//
// Q is built by the Leibniz recurrence
//   Q^(p+1,q) = sum_{i<=p, k<=q} C(p,i) C(q,k) Q^(p-i,q-k) l^(i+1,k)
// (and the same with the extra differentiation in sigma), from the
// derivatives of the log-mass of the mass without its first factors. The
// mass is a product of linear factors,
//   f_j = C(N, j) prod_{i<j} u_i prod_{i<N-j} v_i / prod_{i<N} w_i,
// and the factors u_0 = mu (present for j >= 1) and v_0 = 1 - mu (present
// for j < N) are the ones that are small when mu is near 0 or 1. The
// recurrence run on the full log-mass forms terms of order mu^-p that cancel
// to order mu^-1 against each other when one of them is small (1e-5 at
// mu = 0.001, sigma = 10, against 1e-14 here); so those two factors are taken
// out: Q~ is the recurrence on the remaining factors (their derivatives of
// the log-mass are the sums below with the term i = 0 left out), and the two
// factors enter by the Leibniz rule with Q_0^(p) = p! times the coefficient
// of dmu^p in (1 + dmu / mu)(1 - dmu / (1 - mu)), which has degree two at
// most. The full derivatives l^(i,k) that multiply Q are the sums with the
// term i = 0 included.
//
// The derivatives of the log-mass are the direct sums of pt_betabinom.h
// (bb1_direct): with u_i = mu + i sigma, v_i = 1 - mu + i sigma and
// w_i = 1 + i sigma,
//   l = lchoose(N, j) + sum_{i<j} log u_i + sum_{i<N-j} log v_i
//         - sum_{i<N} log w_i,
// and a derivative with a differentiations in mu and b in sigma, of order
// k = a + b, is
//   c_k [ sum_{i<j} i^b / u_i^k + (-1)^a sum_{i<N-j} i^b / v_i^k
//         - [a = 0] sum_{i<N} i^b / w_i^k ],   c_k = (-1)^(k-1) (k-1)!.
// Each sum has terms of one sign and the cumulative sums over i are formed
// once per row, so every derivative at every count is read in O(1); the
// mass is formed in the log scale from the cumulative sums of the logarithms
// and from the recurrence of lchoose. No shape, no polygamma and no map to
// the shapes enters, so nothing of order sigma^-(k+1) is formed and
// nothing cancels as sigma goes to zero. The route is used at every sigma
// (the order-2 kernels take it below bb1_direct_cut() and the shapes above).
//
// Each order is its own kernel: the order-3 kernel returns the twelve
// components of d3expected_names(), the order-4 kernel the fifteen of
// d4expected_names(), in that order: the pairs (mu, mu), (sigma, sigma),
// (mu, sigma), each with the derivatives mu^p sigma^(T-p) for p = T, ..., 0.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "d7_par.h"
#include "d34e_run.h"
using namespace Rcpp;

namespace {

// c_k = (-1)^(k-1) (k-1)!
const double kC34[7] = {0.0, 1.0, -1.0, 2.0, -6.0, 24.0, -120.0};
const double kB34[5][5] = {{1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 2, 1, 0, 0},
                           {1, 3, 3, 1, 0}, {1, 4, 6, 4, 1}};

// One row of order T (3 or 4): the 3 (T + 1) components.
template <int T>
inline void bb1_d34_row(double mu, double s, int N, d34e::Row o) {
    constexpr int K = T + 2;
    constexpr int NI = K * (K + 3) / 2;   // pairs (k, b), 1 <= k <= K, 0 <= b <= k
    constexpr int NC = 3 * (T + 1);
    const double nu = 1.0 - mu;
    // U[i * NI + idx(k, b)] = sum_{m<i} m^b / u_m^k, V likewise in v, and the
    // cumulative logarithms
    std::vector<double> U((N + 1) * NI, 0.0), V((N + 1) * NI, 0.0);
    std::vector<double> LU(N + 1, 0.0), LV(N + 1, 0.0);
    double W[K + 1];
    for (int k = 0; k <= K; ++k) W[k] = 0.0;
    double LW = 0.0;
    // powers of 1 / mu and 1 / (1 - mu), for the terms i = 0
    double imu[K + 1], inu[K + 1];
    imu[0] = inu[0] = 1.0;
    for (int k = 1; k <= K; ++k) { imu[k] = imu[k - 1] / mu; inu[k] = inu[k - 1] / nu; }
    for (int i = 0; i < N; ++i) {
        const double di = static_cast<double>(i);
        const double u = mu + di * s, v = nu + di * s, w = 1.0 + di * s;
        LU[i + 1] = LU[i] + std::log(u);
        LV[i + 1] = LV[i] + std::log(v);
        LW += std::log(w);
        const double iu = 1.0 / u, iv = 1.0 / v, iw = 1.0 / w;
        double pu = 1.0, pv = 1.0, pw = 1.0;
        const double* Ui = &U[static_cast<std::size_t>(i) * NI];
        const double* Vi = &V[static_cast<std::size_t>(i) * NI];
        double* Un = &U[static_cast<std::size_t>(i + 1) * NI];
        double* Vn = &V[static_cast<std::size_t>(i + 1) * NI];
        // the term i = 0 is left out of U and V (and adds nothing to W)
        for (int k = 1; k <= K; ++k) {
            pu *= iu; pv *= iv; pw *= iw;
            const int base = k * (k + 1) / 2 - 1;
            double pb = 1.0;
            for (int b = 0; b <= k; ++b) {
                Un[base + b] = Ui[base + b] + (i > 0 ? pb * pu : 0.0);
                Vn[base + b] = Vi[base + b] + (i > 0 ? pb * pv : 0.0);
                if (b == k) W[k] += pb * pw;
                pb *= di;
            }
        }
    }
    const int pa[3] = {2, 0, 1}, pbb[3] = {0, 2, 1};   // (mu mu), (sigma sigma), (mu sigma)
    double acc[NC];
    for (int c = 0; c < NC; ++c) acc[c] = 0.0;
    double lc = 0.0;   // log choose(N, j)
    for (int j = 0; j <= N; ++j) {
        // Lt[a][b]: a derivatives in mu and b in sigma of the log-mass at count j
        // without the factors u_0 and v_0; L is the full one
        double Lt[K + 1][K + 1], L[K + 1][K + 1];
        for (int a = 0; a <= K; ++a)
            for (int b = 0; b <= K; ++b) Lt[a][b] = 0.0;
        const double* Uj = &U[static_cast<std::size_t>(j) * NI];
        const double* Vj = &V[static_cast<std::size_t>(N - j) * NI];
        for (int k = 1; k <= K; ++k) {
            const int base = k * (k + 1) / 2 - 1;
            for (int b = 0; b <= k; ++b) {
                const int a = k - b;
                double val = Uj[base + b] + ((a & 1) ? -Vj[base + b] : Vj[base + b]);
                if (a == 0) val -= W[k];
                Lt[a][b] = kC34[k] * val;
            }
        }
        for (int a = 0; a <= K; ++a)
            for (int b = 0; b <= K; ++b) L[a][b] = Lt[a][b];
        for (int a = 1; a <= K; ++a) {
            L[a][0] += kC34[a] * ((j >= 1 ? imu[a] : 0.0) +
                                  ((a & 1) ? -1.0 : 1.0) * (j < N ? inu[a] : 0.0));
        }
        const double p = std::exp(lc + LU[j] + LV[N - j] - LW);
        // Qt[p][q]: the derivatives of the mass without u_0 and v_0, relative to it
        double Qt[T + 1][T + 1], Q[T + 1][T + 1];
        Qt[0][0] = 1.0;
        for (int d = 1; d <= T; ++d) {
            for (int pp = 0; pp <= d; ++pp) {
                const int q = d - pp;
                const int pm = pp > 0 ? pp - 1 : 0;
                const int qm = pp > 0 ? q : q - 1;
                const int da = pp > 0 ? 1 : 0, db = pp > 0 ? 0 : 1;
                double sm = 0.0;
                for (int i = 0; i <= pm; ++i)
                    for (int k = 0; k <= qm; ++k)
                        sm += kB34[pm][i] * kB34[qm][k] * Qt[pm - i][qm - k] *
                              Lt[i + da][k + db];
                Qt[pp][q] = sm;
            }
        }
        // the two exact factors: h[i] = i! times the coefficient of dmu^i
        double h[3] = {1.0, 0.0, 0.0};
        if (j >= 1 && j < N) {
            h[1] = 1.0 / mu - 1.0 / nu;
            h[2] = -2.0 / (mu * nu);
        } else if (j >= 1) {
            h[1] = 1.0 / mu;
        } else if (j < N) {
            h[1] = -1.0 / nu;
        }
        for (int pp = 0; pp <= T; ++pp) {
            for (int q = 0; pp + q <= T; ++q) {
                double sm = Qt[pp][q];
                if (pp >= 1) sm += pp * h[1] * Qt[pp - 1][q];
                if (pp >= 2) sm += kB34[pp][2] * h[2] * Qt[pp - 2][q];
                Q[pp][q] = sm;
            }
        }
        int w = 0;
        for (int r = 0; r < 3; ++r) {
            const int n0 = pa[r], n1 = pbb[r];
            for (int pp = T; pp >= 0; --pp) {
                const int q = T - pp;
                double sm = 0.0;
                for (int i = 0; i <= pp; ++i)
                    for (int k = 0; k <= q; ++k)
                        sm += kB34[pp][i] * kB34[q][k] * L[n0 + i][n1 + k] *
                              Q[pp - i][q - k];
                acc[w++] += p * sm;
            }
        }
        if (j < N) lc += std::log(static_cast<double>(N - j) / (j + 1.0));
    }
    for (int c = 0; c < NC; ++c) o[c] = acc[c];
}

}  // namespace

// [[Rcpp::export]]
List betabinom1_dexpected3_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                    double size, int threads = 1) {
    d34e::Par Mv(mu), Sv(sigma);
    const int N = static_cast<int>(size);
    return d34e::run(y.size(), threads, d7::kMinCostly, 12,
                     [&](std::size_t i, d34e::Row o) {
        bb1_d34_row<3>(Mv[i], Sv[i], N, o);
    });
}

// [[Rcpp::export]]
List betabinom1_dexpected4_cpp(NumericVector y, NumericVector mu, NumericVector sigma,
                    double size, int threads = 1) {
    d34e::Par Mv(mu), Sv(sigma);
    const int N = static_cast<int>(size);
    return d34e::run(y.size(), threads, d7::kMinCostly, 15,
                     [&](std::size_t i, d34e::Row o) {
        bb1_d34_row<4>(Mv[i], Sv[i], N, o);
    });
}
