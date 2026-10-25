// The support walk shared by the third and fourth derivatives of the
// Poisson-inverse Gaussian's expected information (src/d34e_pig1.cpp,
// src/d34e_pig2.cpp). The per-observation rows it reads are those of
// pt_pig_rows.h, the same code as pig_hd.cpp and d5_pig*.cpp.
#ifndef DISTRIBUTIONS7_PT_PIG_D34E_H
#define DISTRIBUTIONS7_PT_PIG_D34E_H

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "d7_par.h"
#include "pig_asym.h"
#include "pt_pig.h"
#include "pt_pig_d5.h"
#include "pt_sqrt.h"
#include "pt_pig_rows.h"

namespace d34e_pig {
using d7::psi_derivs;
using d7::PIG_OFF;
using d7::PIG_LEN;
using d7::pix;
using d7::pig_compose;
using d7::pig_logS_w;
using d7::PigSeq;
using d7::PigSeqW;
using d7::PIG_ASYM_CUT;
using d7pigrow::pig1_row;
using d7pigrow::pig2_row;
using d7pigrow::pig1_d5_row;
using d7pigrow::pig2_d5_row;










// ---------------------------------------------------------------------------
// The support walk of the third and fourth derivatives of the expected
// information (d34e_pig1.cpp, d34e_pig2.cpp). It is pig_expected_obs() of
// pig_hd.cpp with the order-3 or order-4 accumulation in place of the order
// 0-2 one, and with the same stopping rule. One pass over the support serves
// every component: the row is read once per response, the accumulator adds
// f times each component's polynomial in the derivatives of log f.
//
// The row supplies v[0..14], the partials of log f of orders 0-4 in the slots
// of pix(i, j), and, for K = 4, v[15..20], the fifth-order partials in the
// order of the d5 row (v[15 + j] is the partial with j indices on the second
// parameter). The fifth-order row is computed only past the first test that
// f is positive, as the mass underflows below the mode.
//
// Stopping: the tail beyond y is bounded by f max{r/(1-r), tail} times the
// growth (y+1)^(4+2K)/(mu+1)^(4+2K) of the summands; the walk ends when the
// bound falls below 1e-17 of the mass accumulated (the rule of pig_hd.cpp,
// with the exponent 4 + 2K, 10 for K = 3 and 12 for K = 4).
// ---------------------------------------------------------------------------
template <int K, class Seq, class Row4, class Row5, class Acc>
static bool d34_walk(double m, double arg, double tail, Row4 row4, Row5 row5,
                     Acc acc, double* out) {
    const int nk = (K == 3) ? 12 : 15;
    const int pw = 4 + 2 * K;
    Seq seq(arg);
    double v[21], pre[7], cum = 0.0, fprev = 0.0;
    for (int j = 0; j < nk; ++j) out[j] = 0.0;
    const long cap = 100000000L;
    for (long yy = 0; yy < cap; ++yy) {
        if (yy > 0) seq.next();
        seq.pre(pre);
        row4((double) yy, pre, v);
        const double f = std::exp(v[0]);
        if (!(f > 0.0)) {
            // the mass underflows below the mode, where nothing it would
            // contribute is representable; past the mean it has run out
            if ((double) yy > m) return true;
            continue;
        }
        if (K == 4) row5((double) yy, v);
        cum += f;
        acc(f, v, out);
        if ((double) yy > m && fprev > 0.0) {
            const double r = f / fprev;
            const double rt = (r < 1.0) ? r / (1.0 - r) : R_PosInf;
            const double bound = f * std::max(rt, tail) *
                std::pow(((double) yy + 1.0) / (m + 1.0), pw);
            if (bound <= 1e-17 * cum) return true;
        }
        fprev = f;
    }
    return false;
}

}  // namespace d34e_pig

#endif
