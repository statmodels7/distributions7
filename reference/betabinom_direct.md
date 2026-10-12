# The Beta-Binomial Derivatives From the Direct Sums

Returns the derivatives of one order of the beta-binomial log-mass from
the finite sums over the integer shifts that the compiled kernel
`betabinom_direct_cpp` forms, keyed as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
keys them. In the shapes \\(\alpha, \beta)\\ a pure derivative of order
\\K\\ in \\\alpha\\ is \\c_K\[\sum\_{i\<y}(\alpha+i)^{-K} -
\sum\_{i\<n}(\alpha+\beta+i)^{-K}\]\\ with \\c_K = (-1)^{K-1}(K-1)!\\,
the pure one in \\\beta\\ the same with \\n - y\\, and a mixed one
\\-c_K\sum\_{i\<n}(\alpha+\beta+i)^{-K}\\. In the mean and dispersion
\\(\mu, \sigma)\\ the sums run over \\\mu + i\sigma\\, \\1 - \mu +
i\sigma\\ and \\1 + i\sigma\\, with the power \\i^b\\ for \\b\\
derivatives in \\\sigma\\. Every sum has terms of one sign, so the
values keep their accuracy at any concentration; the cost is \\O(n)\\ an
observation.

## Usage

``` r
betabinom_direct(y, p1, p2, n, order, params, shapes)
```

## Arguments

- y:

  A numeric vector of counts in \\\\0, \dots, n\\\\.

- p1, p2:

  The two parameters, `mu` and `sigma` or `alpha` and `beta`, each a
  numeric vector of length 1 or of the length of `y`.

- n:

  The size, a single positive integer.

- order:

  The derivative order, an integer from 1 to 4.

- params:

  The names of the two parameters.

- shapes:

  Logical of length 1: `TRUE` for the shapes, `FALSE` for the mean and
  dispersion.

## Value

A named list of numeric vectors, one per distinct multi-index of the
order, each of the length of the recycled inputs.

## See also

[`bb_direct_cut()`](https://statmodels7.github.io/distributions7/reference/bb_direct_cut.md)
for where the families take this route;
[`betabinom2_derivs()`](https://statmodels7.github.io/distributions7/reference/betabinom2_derivs.md)
and
[`betabinom1_components()`](https://statmodels7.github.io/distributions7/reference/betabinom1_components.md),
which call it.
