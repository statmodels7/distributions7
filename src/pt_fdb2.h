#ifndef D7_PT_FDB2_H
#define D7_PT_FDB2_H

#include <algorithm>
#include <map>
#include <vector>

// Faa di Bruno for an outer function of one variable composed with an inner
// function of two variables (x1, x2), to order five. The derivative of
// f(g(x1, x2)) of order i in x1 and j in x2 is
//   sum over the set partitions of the i + j slots of
//   f^(number of blocks)(g) * prod over the blocks of d^(p_b, q_b) g,
// where p_b and q_b count the x1 and x2 slots in block b. Partitions with the
// same blocks are merged into one term with a multiplicity, so each (i, j)
// carries a short list of terms built once.

namespace d7 {

struct Fdb2Term {
    int nb;
    int p[5];
    int q[5];
    double mult;
};

inline std::vector<Fdb2Term> fdb2_build_terms(int i, int j) {
    const int k = i + j;
    std::map<std::vector<int>, Fdb2Term> acc;
    int rg[5] = {0, 0, 0, 0, 0};
    while (true) {
        int nb = 0;
        for (int s = 0; s < k; ++s) if (rg[s] + 1 > nb) nb = rg[s] + 1;
        std::vector<int> key;
        for (int b = 0; b < nb; ++b) {
            int p = 0, q = 0;
            for (int s = 0; s < k; ++s) {
                if (rg[s] == b) { if (s < i) ++p; else ++q; }
            }
            key.push_back(p * 8 + q);
        }
        std::sort(key.begin(), key.end());
        auto it = acc.find(key);
        if (it == acc.end()) {
            Fdb2Term t{nb, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, 1.0};
            for (int b = 0; b < nb; ++b) { t.p[b] = key[b] / 8; t.q[b] = key[b] % 8; }
            acc[key] = t;
        } else {
            it->second.mult += 1.0;
        }
        // the next restricted growth string
        int t2 = k - 1;
        for (; t2 >= 1; --t2) {
            int mx = 0;
            for (int s = 0; s < t2; ++s) if (rg[s] > mx) mx = rg[s];
            if (rg[t2] <= mx) break;
        }
        if (t2 < 1) break;
        ++rg[t2];
        for (int s = t2 + 1; s < k; ++s) rg[s] = 0;
    }
    std::vector<Fdb2Term> out;
    for (auto& kv : acc) out.push_back(kv.second);
    return out;
}

// The terms of order (i, j), 1 <= i + j <= 5. The table is built on the
// first call; call it once before a parallel region.
inline const std::vector<Fdb2Term>& fdb2_terms(int i, int j) {
    static const std::vector<std::vector<std::vector<Fdb2Term>>> tab = [] {
        std::vector<std::vector<std::vector<Fdb2Term>>> t(
            6, std::vector<std::vector<Fdb2Term>>(6));
        for (int a = 0; a <= 5; ++a)
            for (int b = 0; a + b <= 5; ++b)
                if (a + b >= 1) t[a][b] = fdb2_build_terms(a, b);
        return t;
    }();
    return tab[i][j];
}

// d^(i, j) f(g) from the outer derivatives fd[r] = f^(r)(g), r = 1..5, and
// the inner table dg[p][q] = d^(p, q) g.
inline double fdb2_compose(int i, int j, const double* fd,
                           const double (*dg)[6]) {
    double s = 0.0;
    for (const Fdb2Term& t : fdb2_terms(i, j)) {
        double pr = t.mult * fd[t.nb];
        for (int b = 0; b < t.nb; ++b) pr *= dg[t.p[b]][t.q[b]];
        s += pr;
    }
    return s;
}

} // namespace d7

#endif
