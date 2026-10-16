# Log-Mass of the Beta-Binomial

Returns \\\log\binom{n}{y} + \log B(y+\alpha, n-y+\beta) - \log
B(\alpha, \beta)\\ by whichever of two routes is accurate at the given
shapes. The choice depends only on the shapes, so it is deterministic
and costs one comparison.

## Usage

``` r
betabinom_log_mass(y, a, b, n)
```

## Arguments

- y:

  A numeric vector of counts, already known to lie on the support. The
  caller tests that; this function does not.

- a, b:

  The two shapes, each a numeric vector of length 1 or of the length of
  `y`, strictly positive.

- n:

  The size, a single positive integer.

## Value

A numeric vector of log-probabilities, of the recycled length of the
inputs.

## Why there are two routes

The two beta functions are each of magnitude
\\(\alpha+\beta)\log(\alpha+\beta)\\ while their difference is of order
one, so forming the mass as that difference carries an absolute error of
\\\varepsilon\\ times the larger magnitude. That route is used only
while the error stays below `1e-8`.

Beyond it the shifts \\y\\, \\n-y\\ and \\n\\ are integers, so each
log-gamma difference is an exact sum of logarithms,
\$\$\log\Gamma(\alpha+y) - \log\Gamma(\alpha) =
\sum\_{j=0}^{y-1}\log(\alpha+j),\$\$ and the mass follows from three
such sums, each of which forms only quantities no larger than
\\n\log(\alpha+\beta)\\. The sums also give the binomial limit correctly
as the shapes grow at a fixed ratio, where the beta-function route loses
one digit per factor of ten in the concentration.

The cost is \\O(n)\\ in the size, which is why the route is taken only
where it is needed.

## See also

[`distrib_pdf.BetaBinom2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BetaBinom2Distrib.md),
which calls this after testing the support, and
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
for the family.
