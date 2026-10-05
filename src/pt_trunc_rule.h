#ifndef D7_PT_TRUNC_RULE_H
#define D7_PT_TRUNC_RULE_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "pt_loc_scale.h"

// The quadrature over the retained interval [a, b] of a truncated
// continuous family, for its expected information. The rule does not look
// at the integrand: it reads the interval, a center c and a scale s of the
// parent, and the points where the parent's density has a kink, so that R
// (trunc_route_parts) and the scalar registry evaluate their integrands at
// the same nodes with the same weights.
//
// The interval is cut at c and at each kink inside it. A finite piece is
// cut again into panels whose widths double away from the end nearer the
// center, s, 2s, 4s, ..., so that the panel next to the mass has the
// parent's scale and a long piece costs a number of panels logarithmic in
// its length; each panel takes a tanh-sinh rule (step 1/16, |t| <= 3.2),
// whose nodes cluster at both ends of the panel, where a kink or an
// integrable singularity of the density sits. A piece unbounded on one side
// takes the exp-sinh half of loc_scale_rule() scaled by s from its finite
// end, beyond the doubling panels that cover the first 2^6 - 1 scales.

namespace d7 {

const double kTsStep = 1.0 / 16.0;
const double kTsHalf = 3.2;
const int kDoublings = 6;

// the tanh-sinh nodes of [u, v], appended; the distance of a node to its
// nearer end is computed directly, as fold_rule() does
inline void ts_panel(double u, double v, std::vector<double>& y,
                     std::vector<double>& w, double h = kTsStep) {
    const double L = v - u;
    if (!(L > 0)) return;
    const double hp = M_PI / 2.0;
    const int K = (int) std::floor(kTsHalf / h);
    for (int i = -K; i <= K; ++i) {
        const double t = i * h;
        const double q = hp * std::sinh(t);
        const double e = std::exp(-2.0 * std::fabs(q));
        const double near = L * e / (1.0 + e);
        const double yy = (q < 0) ? u + near : v - near;
        const double ww = h * hp * std::cosh(t) * L * 2.0 * e /
            ((1.0 + e) * (1.0 + e));
        if (!(yy > u) || !(yy < v)) continue;
        y.push_back(yy);
        w.push_back(ww);
    }
}

// a finite piece [u, v] in doubling panels from the end `from` (u or v)
inline void ts_doubling(double u, double v, bool from_u, double s,
                        std::vector<double>& y, std::vector<double>& w,
                        double h = kTsStep) {
    double width = s, at = from_u ? u : v;
    for (;;) {
        const double next = from_u ? at + width : at - width;
        // a panel that would not advance (a zero or non-finite scale, or a
        // width below the spacing of doubles at `at`) closes the piece
        const bool stuck = !(width > 0) || !R_FINITE(width) || next == at;
        const bool last = stuck || (from_u ? !(next < v) : !(next > u));
        if (last) {
            if (from_u) ts_panel(at, v, y, w, h); else ts_panel(u, at, y, w, h);
            return;
        }
        if (from_u) ts_panel(at, next, y, w, h); else ts_panel(next, at, y, w, h);
        at = next;
        width *= 2.0;
    }
}

// a piece unbounded on one side: doubling panels over the first
// 2^kDoublings - 1 scales, then the exp-sinh half rule
inline void es_tail(double u, bool right, double s, std::vector<double>& y,
                    std::vector<double>& w, double h = kTsStep,
                    int nd = kDoublings) {
    const double span = s * (std::ldexp(1.0, nd) - 1.0);
    const double end = right ? u + span : u - span;
    if (nd > 0) {
        if (right) ts_doubling(u, end, true, s, y, w, h);
        else ts_doubling(end, u, false, s, y, w, h);
    }
    const LocScaleRule R = loc_scale_rule();
    for (int j = R.n / 2; j < R.n; ++j) {
        y.push_back(right ? end + s * R.x[j] : end - s * R.x[j]);
        w.push_back(s * R.w[j]);
    }
}

// the rule over [a, b] (either end may be infinite), center c, scale s, and
// the kinks kk[0..nk-1]
inline void trunc_rule(double a, double b, double c, double s,
                       const double* kk, int nk, std::vector<double>& y,
                       std::vector<double>& w, double h = kTsStep,
                       int nd = kDoublings) {
    y.clear();
    w.clear();
    std::vector<double> cuts;
    cuts.push_back(a);
    if (c > a && c < b) cuts.push_back(c);
    for (int i = 0; i < nk; ++i)
        if (kk[i] > a && kk[i] < b && kk[i] != c) cuts.push_back(kk[i]);
    cuts.push_back(b);
    std::sort(cuts.begin(), cuts.end());
    for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
        const double u = cuts[i], v = cuts[i + 1];
        if (!R_FINITE(u) && !R_FINITE(v)) {
            es_tail(c, false, s, y, w, h, nd);
            es_tail(c, true, s, y, w, h, nd);
        } else if (!R_FINITE(v)) {
            es_tail(u, true, s, y, w, h, nd);
        } else if (!R_FINITE(u)) {
            es_tail(v, false, s, y, w, h, nd);
        } else {
            // doubling away from the end nearer the center
            const bool from_u = std::fabs(u - c) <= std::fabs(v - c);
            ts_doubling(u, v, from_u, s, y, w, h);
        }
    }
}

} // namespace d7

#endif
