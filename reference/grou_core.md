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

A numeric vector of draws. An error is signaled if the log-density at
the mode is not finite, if no bounding rectangle can be built, or if no
draw is accepted in 1000 rounds.

## Details

Kept separate from the `distrib` object so that it can also be run on a
reparametrized density, which is how the divergence transforms reuse it.

Two devices make it safe. The kernel is **recentered at the mode**,
without which a density located at \\\mu = 1000\\ produces a degenerate
bounding box; and it is **normalized to a maximum of one**. With those
it accepts a Student t with half a degree of freedom and a Pareto with
infinite mean. A density that diverges at an edge is handled by
transforming it away. The kernel must be unimodal: with two
well-separated modes the search for the bounding rectangle stops at the
valley next to the first mode.

## See also

[`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md),
[`find_lp_anchor()`](https://statmodels7.github.io/distributions7/reference/find_lp_anchor.md)
