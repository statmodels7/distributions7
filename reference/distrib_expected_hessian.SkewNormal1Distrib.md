# Skew Normal Expected Hessian and Its Derivatives

Returns the expected information, and through
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
its first and second derivatives in the parameters.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- y:

  A numeric vector of observations. Its length sets the length of each
  returned component; the values themselves are not read.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of the length of `y`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html).

- approx, nsim:

  Ignored.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  The thread count passed to the family's kernels.

## Value

A named list of numeric vectors of length `length(y)`: for
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
the components keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
for the two derivatives those keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
and
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

No component has an elementary form. In Azzalini's notation, with \\b =
\sqrt{2/\pi}\\ and \\a_k = \mathbb{E}\\\left\[Z^k \\\phi(\alpha
Z)/\Phi(\alpha Z)\\^2\right\]\\ for \\Z\\ a standard skew normal of
shape \\\alpha\\, \$\$-\mathbb{E}\[\ell\_{\mu\mu}\] = \frac{1 + \alpha^2
a_0}{\sigma^2},\quad -\mathbb{E}\[\ell\_{\sigma\sigma}\] = \frac{2 +
\alpha^2 a_2}{\sigma^2},\quad -\mathbb{E}\[\ell\_{\alpha\alpha}\] =
a_2,\$\$ and each \\a_k\\ is an integral that has no closed form. The
family is a location-scale family, so every component equals its value
at \\\mu = 0\\, \\\sigma = 1\\ times \\\sigma^{-k}\\, with \\k\\ the
number of indices on \\\mu\\ or \\\sigma\\. The value there depends on
\\\alpha\\ alone and is one integral over \\z\\ of the analytic observed
derivatives against the density, taken once per distinct \\\alpha\\ by
the exp-sinh rule of
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md).
It reproduces the three expressions above, computed from the \\a_k\\ by
an independent quadrature, to the last printed digit.

At \\\alpha = 0\\ the information has rank 2 and not 3: its eigenvalues
at \\\mu = 0\\, \\\sigma = 1\\ are 2, 1.637 and \\-2.6\times10^{-27}\\,
and the smallest grows like \\\alpha^4\\, reading \\4.4\times10^{-10}\\
at \\\alpha = 0.01\\ and \\1.9\times10^{-3}\\ at \\\alpha = 0.5\\. See
[`distrib_hessian.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal1Distrib.md)
for why.

`approx` and `nsim` are accepted for the generic's sake and ignored.

## See also

[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md)
for the construction,
[`distrib_hessian.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal1Distrib.md)
for the quantity this is the expectation of, and
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
a <- 2
e <- distrib_expected_hessian(d, 0, list(mu = 0, sigma = 1, alpha = a))

# Azzalini's form, with the a_k computed by an independent quadrature.
ak <- sapply(c(0, 2), function(k) integrate(function(z)
  2 * z^k * exp(dnorm(z, log = TRUE) + 2 * dnorm(a * z, log = TRUE) -
                pnorm(a * z, log.p = TRUE)), -Inf, Inf, rel.tol = 1e-12)$value)
rbind(azzalini = c(1 + a^2 * ak[1], 2 + a^2 * ak[2], ak[2]),
      package = -c(e$mu_mu, e$sigma_sigma, e$alpha_alpha))
#>              [,1]     [,2]       [,3]
#> azzalini 2.260226 2.298811 0.07470277
#> package  2.260226 2.298811 0.07470277
```
