# Ratio-of-Uniforms for a Density Diverging at Both Edges

Samples a density unbounded at *both* ends of a finite support, by
mapping the interval through a transform behaving like a different power
at each end.

## Usage

``` r
grou_two_sided(lp, b, div, n, r)
```

## Arguments

- lp:

  A function giving the log-density.

- b:

  A length-2 numeric vector, the support.

- div:

  The exponents at the two edges, from
  [`lp_edge_divergence`](https://statmodels7.github.io/distributions7/reference/lp_edge_divergence.md).

- n:

  The number of draws wanted.

- r:

  The ratio-of-uniforms tuning parameter.

## Value

A numeric vector of draws, or `NULL` if the sampler gives up.

## Details

A single power cannot repair a two-sided divergence, since raising to a
power flattens one edge while steepening the other. This map does:
\$\$T(u) = \frac{u^p}{u^p + (1-u)^q}, \qquad Y = a + (b-a) T(U).\$\$ It
increases monotonically from 0 to 1, behaves like \\u^p\\ on the left
and like \\1 - (1-u)^q\\ on the right, and has a closed-form derivative.

The density of \\U\\ therefore carries the exponents \\p\alpha - 1\\ and
\\q\beta - 1\\ at the two ends and is bounded as soon as \\p\alpha \>
1\\ and \\q\beta \> 1\\ – indeed it vanishes there, turning the original
U-shaped density into a single-peaked one, which is exactly what the
sampler wants. Both exponents come from the same probe that detected the
divergence, so nothing has to be searched for.

## See also

[`grou_core`](https://statmodels7.github.io/distributions7/reference/grou_core.md),
[`lp_edge_divergence`](https://statmodels7.github.io/distributions7/reference/lp_edge_divergence.md)
