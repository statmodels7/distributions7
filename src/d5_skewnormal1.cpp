#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "d7_par.h"
#include "pt_logphi.h"
using namespace Rcpp;

// The fifth derivatives of the skew normal log-density in (mu, sigma, alpha)
// (direct parametrization, as skewnormal1_deriv4 in skewnormal_hd.cpp):
//   l = log 2 - log sigma - z^2/2 + log phi-constant + log Phi(t),
//   z = (y - mu)/sigma, t = alpha z.
// Every component is a polynomial in z, alpha, sigma and the derivatives
// G1..G5 of log Phi at t, which logphi_derivs() supplies with the cumulant
// route below t = -1 (pt_logphi.h). The loop body is written out per
// component and derived in misure_gas_congiunto/c/fase4/D5/skewnormal1/derive.py.

// [[Rcpp::export]]
List skewnormal1_deriv5_cpp(NumericVector y, NumericVector mu,
                            NumericVector sigma, NumericVector alpha,
                            int threads = 1) {
  const int n = std::max(y.size(), std::max(mu.size(), std::max(sigma.size(), alpha.size())));
  const int n_y = y.size(), n_mu = mu.size(), n_s = sigma.size(), n_a = alpha.size();
  NumericVector o_mu_mu_mu_mu_mu(n), o_mu_mu_mu_mu_sigma(n), o_mu_mu_mu_mu_alpha(n),
    o_mu_mu_mu_sigma_sigma(n), o_mu_mu_mu_sigma_alpha(n), o_mu_mu_mu_alpha_alpha(n),
    o_mu_mu_sigma_sigma_sigma(n), o_mu_mu_sigma_sigma_alpha(n), o_mu_mu_sigma_alpha_alpha(n),
    o_mu_mu_alpha_alpha_alpha(n), o_mu_sigma_sigma_sigma_sigma(n),
    o_mu_sigma_sigma_sigma_alpha(n), o_mu_sigma_sigma_alpha_alpha(n),
    o_mu_sigma_alpha_alpha_alpha(n), o_mu_alpha_alpha_alpha_alpha(n),
    o_sigma_sigma_sigma_sigma_sigma(n), o_sigma_sigma_sigma_sigma_alpha(n),
    o_sigma_sigma_sigma_alpha_alpha(n), o_sigma_sigma_alpha_alpha_alpha(n),
    o_sigma_alpha_alpha_alpha_alpha(n), o_alpha_alpha_alpha_alpha_alpha(n);
  d7::par_for(n, threads, d7::kMinCostly, [&](std::size_t i) {
    const double yy = y[i % n_y];
    const double mm = mu[i % n_mu];
    const double ss = sigma[i % n_s];
    const double aa = alpha[i % n_a];
    const double tt = aa*((yy - mm)/ss);
    double G[6];
    d7::logphi_derivs(tt, 5, G);
    const double G1 = G[1];
    const double G2 = G[2];
    const double G3 = G[3];
    const double G4 = G[4];
    const double G5 = G[5];
        const double t0 = std::pow(aa, 5);
        const double t1 = std::pow(ss, -5);
        const double t2 = G5*t1;
        const double t3 = std::pow(aa, 4);
        const double t4 = G4*t1;
        const double t5 = 4*t4;
        const double t6 = std::pow(ss, -6);
        const double t7 = -mm + yy;
        const double t8 = G5*t7;
        const double t9 = std::pow(aa, 3);
        const double t10 = std::pow(ss, -4);
        const double t11 = G4*t10;
        const double t12 = t11*t9;
        const double t13 = t3*t7;
        const double t14 = t13*t2;
        const double t15 = G3*t9;
        const double t16 = 12*t1;
        const double t17 = G4*t6;
        const double t18 = 5*t17;
        const double t19 = 3*t12 + t14;
        const double t20 = std::pow(ss, -2);
        const double t21 = t20*t7;
        const double t22 = aa*t21;
        const double t23 = std::pow(aa, 2);
        const double t24 = G3*t23;
        const double t25 = t10*t24;
        const double t26 = t7*t9;
        const double t27 = t26*t5;
        const double t28 = 1.0/ss;
        const double t29 = t28*t7;
        const double t30 = std::pow(ss, -3);
        const double t31 = G3*t30;
        const double t32 = 6*aa;
        const double t33 = 3*t7;
        const double t34 = t11*t23;
        const double t35 = t33*t34;
        const double t36 = t23*t30;
        const double t37 = 3*G4;
        const double t38 = 24*t1;
        const double t39 = G2*t23;
        const double t40 = 20*t6;
        const double t41 = 6*t10;
        const double t42 = t24*t41;
        const double t43 = 2*G3;
        const double t44 = -t12*t7 - t36*t43;
        const double t45 = std::pow(t7, 2);
        const double t46 = G5*t45;
        const double t47 = 6*t26*t4 + t3*t46*t6 + t42;
        const double t48 = G2*aa;
        const double t49 = t24*t7;
        const double t50 = t16*t49;
        const double t51 = 4*t31;
        const double t52 = aa*t51;
        const double t53 = -t35 - t52;
        const double t54 = G2*t30;
        const double t55 = t41*t7;
        const double t56 = aa*t55;
        const double t57 = G3*t56;
        const double t58 = t45*t9;
        const double t59 = t30*t7;
        const double t60 = t43*t59;
        const double t61 = t20*t43;
        const double t62 = 2*t59;
        const double t63 = aa*t62;
        const double t64 = G1*t38;
        const double t65 = 60*t6;
        const double t66 = G2*t20;
        const double t67 = aa*t66 + t24*t59;
        const double t68 = t41*t48;
        const double t69 = 2*t54;
        const double t70 = aa*t69;
        const double t71 = -t25*t33 - t70;
        const double t72 = t1*t49;
        const double t73 = -4*t25*t7 - t4*t58 - t70;
        const double t74 = std::pow(ss, -7);
        const double t75 = std::pow(t7, 3);
        const double t76 = G5*t75;
        const double t77 = 9*t17;
        const double t78 = t3*t74*t76 + t58*t77 + t68 + 18*t72;
        const double t79 = t38*t48*t7;
        const double t80 = aa*t60 + t66;
        const double t81 = aa*t10;
        const double t82 = G3*t7*t81;
        const double t83 = t23*t4*t45;
        const double t84 = -t69 - 8*t82 - 3*t83;
        const double t85 = G2*t55;
        const double t86 = t1*t45;
        const double t87 = aa*t86;
        const double t88 = 2*t21;
        const double t89 = t6*t9;
        const double t90 = t10*t45;
        const double t91 = t43*t90;
        const double t92 = 2*aa;
        const double t93 = t11*t45;
        const double t94 = t23*t75;
        const double t95 = t11*t75;
        const double t96 = aa*t6;
        const double t97 = t25*t45 + t7*t70;
        const double t98 = t24*t45;
        const double t99 = t68*t7;
        const double t100 = t24*t86;
        const double t101 = -4*t100 - t99;
        const double t102 = t6*t98;
        const double t103 = G4*t75;
        const double t104 = t74*t9;
        const double t105 = t103*t104;
        const double t106 = -6*t100 - t103*t89 - t99;
        const double t107 = std::pow(t7, 4);
        const double t108 = G5*t107;
        const double t109 = 36*t102 + 12*t105 + t79 + t108*t3/std::pow(ss, 8);
        const double t110 = aa*t91 + t69*t7;
        const double t111 = 8*t86;
        const double t112 = G3*t45;
        const double t113 = -aa*t112*t16 - t37*t6*t94 - t85;
        const double t114 = t43*t75;
        const double t115 = t4*t75;
        o_mu_mu_mu_mu_mu[i] = -t0*t2;
        o_mu_mu_mu_mu_sigma[i] = -t0*t6*t8 - t3*t5;
        o_mu_mu_mu_mu_alpha[i] = 4*t12 + t14;
        o_mu_mu_mu_sigma_sigma[i] = -t13*t18 - t15*t16 - t19*t22;
        o_mu_mu_mu_sigma_alpha[i] = t19*t29 + 9*t25 + t27;
        o_mu_mu_mu_alpha_alpha[i] = t28*t7*(-t10*t8*t9 - t36*t37) - t31*t32 - t35;
        o_mu_mu_sigma_sigma_sigma[i] = 2*aa*t30*t44*t7 - t15*t40*t7 - t22*t47 - t22*(t27 + t42) - t38*t39 + t38;
        o_mu_mu_sigma_sigma_alpha[i] = 12*t10*t48 - t21*t44 - t22*t53 + t29*t47 + t50;
        o_mu_mu_sigma_alpha_alpha[i] = t28*t53*t7 + t28*t7*(-t2*t58 - 5*t34*t7 - t52) - 4*t54 - t57;
        o_mu_mu_alpha_alpha_alpha[i] = t29*(G4*t63 + t61) + t29*(4*G4*aa*t59 + t10*t23*t46 + t61) + t60;
        o_mu_sigma_sigma_sigma_sigma[i] = 4*aa*t30*t7*t71 + 2*aa*t30*t7*t73 - aa*t64 - t22*t78 - t22*(t50 + t68) - t22*(t18*t58 + t68 + 16*t72) - t39*t65*t7 - t56*t67 - t65*(2*mm - 2*yy);
        o_mu_sigma_sigma_sigma_alpha[i] = G1*t41 - t21*t71 - t21*t73 - t22*t84 - t22*(-t57 - t69) + t29*t78 + t62*t67 + t63*t80 + t79;
        o_mu_sigma_sigma_alpha_alpha[i] = t28*t7*t84 + t28*t7*(-t69 - t76*t89 - 10*t82 - 7*t83) - t43*t87 - t80*t88 - t85;
        o_mu_sigma_alpha_alpha_alpha[i] = t29*(t31*t33 + t92*t93) + t29*(5*aa*t93 + t2*t94 + t51*t7) + t91;
        o_mu_alpha_alpha_alpha_alpha[i] = t28*t7*(-t30*t37*t45 - t76*t81) - t95;
        o_sigma_sigma_sigma_sigma_sigma[i] = -120*G1*t7*t96 + 4*aa*t101*t30*t7 + 2*aa*t106*t30*t7 - t109*t22 - t22*(t40*t98 + t79) - t22*(30*t102 + 6*t105 + t79) - t38 - 120*t39*t45*t74 + 360*t45*t74 - t56*t97;
        o_sigma_sigma_sigma_sigma_alpha[i] = -t101*t21 - t106*t21 + t109*t29 + t110*t63 - t113*t22 - t22*(-G3*aa*t111 - t85) + 40*t45*t48*t6 + t62*t97 + t64*t7;
        o_sigma_sigma_sigma_alpha_alpha[i] = -G2*t111 - t110*t88 + t113*t28*t7 - t114*t96 + t28*t7*(-18*G3*t87 - t104*t108 - t77*t94 - t85);
        o_sigma_sigma_alpha_alpha_alpha[i] = t1*t114 + t29*(4*G3*t90 + t115*t92) + t29*(t108*t23*t6 + t112*t41 + t115*t32);
        o_sigma_alpha_alpha_alpha_alpha[i] = -t107*t4 + t28*t7*(-aa*t107*t2 - 3*t95);
        o_alpha_alpha_alpha_alpha_alpha[i] = t2*std::pow(t7, 5);
  });
  return List::create(Named("mu_mu_mu_mu_mu") = o_mu_mu_mu_mu_mu,
                      Named("mu_mu_mu_mu_sigma") = o_mu_mu_mu_mu_sigma,
                      Named("mu_mu_mu_mu_alpha") = o_mu_mu_mu_mu_alpha,
                      Named("mu_mu_mu_sigma_sigma") = o_mu_mu_mu_sigma_sigma,
                      Named("mu_mu_mu_sigma_alpha") = o_mu_mu_mu_sigma_alpha,
                      Named("mu_mu_mu_alpha_alpha") = o_mu_mu_mu_alpha_alpha,
                      Named("mu_mu_sigma_sigma_sigma") = o_mu_mu_sigma_sigma_sigma,
                      Named("mu_mu_sigma_sigma_alpha") = o_mu_mu_sigma_sigma_alpha,
                      Named("mu_mu_sigma_alpha_alpha") = o_mu_mu_sigma_alpha_alpha,
                      Named("mu_mu_alpha_alpha_alpha") = o_mu_mu_alpha_alpha_alpha,
                      Named("mu_sigma_sigma_sigma_sigma") = o_mu_sigma_sigma_sigma_sigma,
                      Named("mu_sigma_sigma_sigma_alpha") = o_mu_sigma_sigma_sigma_alpha,
                      Named("mu_sigma_sigma_alpha_alpha") = o_mu_sigma_sigma_alpha_alpha,
                      Named("mu_sigma_alpha_alpha_alpha") = o_mu_sigma_alpha_alpha_alpha,
                      Named("mu_alpha_alpha_alpha_alpha") = o_mu_alpha_alpha_alpha_alpha,
                      Named("sigma_sigma_sigma_sigma_sigma") = o_sigma_sigma_sigma_sigma_sigma,
                      Named("sigma_sigma_sigma_sigma_alpha") = o_sigma_sigma_sigma_sigma_alpha,
                      Named("sigma_sigma_sigma_alpha_alpha") = o_sigma_sigma_sigma_alpha_alpha,
                      Named("sigma_sigma_alpha_alpha_alpha") = o_sigma_sigma_alpha_alpha_alpha,
                      Named("sigma_alpha_alpha_alpha_alpha") = o_sigma_alpha_alpha_alpha_alpha,
                      Named("alpha_alpha_alpha_alpha_alpha") = o_alpha_alpha_alpha_alpha_alpha);
}
