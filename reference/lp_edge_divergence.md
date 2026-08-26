# Detect and Measure a Divergence at the Edges of the Support

For each finite edge, the exponent \\\alpha\\ of a divergence \\f(y)
\sim \lvert y - a \rvert^{\alpha - 1}\\, or `NA` when the density stays
bounded there.

## Usage

``` r
lp_edge_divergence(lp, b)
```

## Arguments

- lp:

  A function giving the log-density.

- b:

  A length-2 numeric vector, the support.

## Value

A numeric vector of length 2, named `lower` and `upper`.

## Details

A density that diverges at an edge is the one case the ratio-of-uniforms
sampler cannot handle directly, so it has to be detected and then
removed by a change of variable, which needs the exponent.

Detecting and measuring are the same operation. Walking towards the edge
in decades lifts the log-density by \\(1 - \alpha)\log 10\\ per step
when it diverges, and by an amount that dies away when it does not. So
the probe establishing *whether* the density diverges also reports *how
fast*, to about four decimals, with no search, reducing the cost per
draw by several orders of magnitude for strongly divergent shapes.

## See also

[`grou_two_sided()`](https://statmodels7.github.io/distributions7/reference/grou_two_sided.md),
[`rng_grou()`](https://statmodels7.github.io/distributions7/reference/rng_grou.md)
