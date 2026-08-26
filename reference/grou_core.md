# The Generalized Ratio-of-Uniforms Sampler

The sampler itself, on a bare log-density over an interval.

## Usage

``` r
grou_core(lp, b, n, r)
```

## Arguments

- lp:

  A function giving the log-density.

- b:

  A length-2 numeric vector, the support.

- n:

  The number of draws wanted.

- r:

  The ratio-of-uniforms tuning parameter.

## Value

A numeric vector of draws, or `NULL` if no bounding box is found.

## Details

Kept separate from the `distrib` object so that it can also be run on a
reparametrized density, which is how the divergence transforms reuse it.

Two devices make it safe. The kernel is **recentered at the mode**,
without which a density located at \\\mu = 1000\\ produces a degenerate
bounding box; and it is **normalized to a maximum of one**. With those
it rejects far less often than expected: bimodal densities, a Student t
with half a degree of freedom and a Pareto with infinite mean are all
fine. The only genuine rejection is a density that diverges at an edge,
which is handled by transforming it away.

## See also

[`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md),
[`find_lp_anchor()`](https://statmodels7.github.io/distributions7/reference/find_lp_anchor.md)
