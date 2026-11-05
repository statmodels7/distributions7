#ifndef D7_PT_PARTS_H
#define D7_PT_PARTS_H
#include <algorithm>
#include <array>
#include <vector>
// The parts route of the registry (piano_motore.md, D7): a family's parts
// function fills a buffer of doubles at one observation with the quantities
// that cost (the special functions), up to the order a caller will read, and
// each order is then its own function over the buffer, doing arithmetic only.
// buf[0] holds the level the buffer was filled to; an order above it writes
// NaN. The buffer holds at most kPartsMax doubles.
//
// The orders three and four are written as full P^k arrays, column-major,
// out[a + P b + P^2 c (+ P^3 d)], from the unique components a family
// computes in the order of its vector kernels: the index tuples sorted
// increasingly, in lexicographic order (for P = 2 and k = 3: (1,1,1),
// (1,1,2), (1,2,2), (2,2,2)).
namespace d7 {

const int kPartsMax = 64;

// rank of a sorted tuple among the sorted tuples of length k over P values,
// in lexicographic order
inline int sorted_tuple_rank(int P, int k, const int* t) {
    int rank = 0;
    // count the tuples before t by enumeration; called only to build tables
    std::vector<int> u(k, 0);
    while (true) {
        if (std::equal(u.begin(), u.end(), t)) return rank;
        ++rank;
        int j = k - 1;
        while (j >= 0 && u[j] == P - 1) --j;
        if (j < 0) return -1;
        ++u[j];
        for (int i = j + 1; i < k; ++i) u[i] = u[j];
    }
}

typedef std::array<std::array<std::vector<int>, 5>, 5> SymTables;

inline SymTables sym_tables_build() {
    SymTables tables;
    for (int p = 1; p <= 4; ++p) {
        for (int kk = 1; kk <= 4; ++kk) {
            int total = 1;
            for (int i = 0; i < kk; ++i) total *= p;
            std::vector<int>& tab = tables[p][kk];
            tab.resize(total);
            std::vector<int> t(kk);
            for (int j = 0; j < total; ++j) {
                int r = j;
                for (int i = 0; i < kk; ++i) { t[i] = r % p; r /= p; }
                std::sort(t.begin(), t.end());
                tab[j] = sorted_tuple_rank(p, kk, t.data());
            }
        }
    }
    return tables;
}

// the table from a full index to the rank of its sorted tuple, for P <= 4 and
// k <= 4; a function-local static, so it is built once and safely from any
// thread
inline const std::vector<int>& sym_table(int P, int k) {
    static const SymTables tables = sym_tables_build();
    return tables[P][k];
}

// the full P^k array from the unique components
inline void sym_fill(int P, int k, const double* uniq, double* out) {
    const std::vector<int>& tab = sym_table(P, k);
    for (std::size_t j = 0; j < tab.size(); ++j) out[j] = uniq[tab[j]];
}

} // namespace d7
#endif
