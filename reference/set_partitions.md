# All Set Partitions of a Finite Index Set

Every way of splitting `{1, ..., n}` into disjoint non-empty blocks, as
a list of lists of integer vectors. There are \\B_n\\ of them, the Bell
number.

## Usage

``` r
set_partitions(n)
```

## Arguments

- n:

  A positive integer, the size of the index set.

## Value

A list of partitions. Each partition is a list of integer vectors, the
blocks.

## Details

Built by the usual recursion: take each partition of `{1, ..., n-1}` and
either add `n` to one of its existing blocks or open a new block for it.

This is what makes the Bartlett identities work at arbitrary order. The
identity of order \\k\\ says that summing, over all partitions \\\pi\\
of the index set, the expectation of \\\prod\_{B \in \pi} \ell_B\\ gives
zero; the single-block partition is the term wanted, so it follows from
all the others. Enumerating the partitions is therefore the whole
content of
[`expected_by_bartlett`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md).

## See also

[`expected_by_bartlett`](https://statmodels7.github.io/distributions7/reference/expected_by_bartlett.md)
