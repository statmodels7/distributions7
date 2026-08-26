# Skew Normal Distribution Function

Computes the skew normal distribution function through Owen's T
function, with \\z = (q-\mu)/\sigma\\: \$\$F(q; \mu, \sigma, \alpha) =
\Phi(z) - 2\\T(z, \alpha).\$\$ The identity is Azzalini's. Each
evaluation costs one bounded one-dimensional quadrature, through
[`numericals7::owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.html),
where the base class would integrate the density over a semi-infinite
range; the bounded integrand is both cheaper and more accurate. A
thousand quantiles cost about 20 milliseconds.

The result is clamped to \\\[0, 1\]\\ before the tail and the logarithm
are applied, so rounding in the quadrature cannot return a probability
outside the unit interval.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- q:

  A numeric vector of quantiles, anywhere on the real line.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of the length of `q`. `sigma` must be strictly
  positive.

- lower.tail:

  Logical of length 1. When `TRUE`, the default, the value is \\P(Y \le
  q)\\; when `FALSE` it is \\P(Y \> q)\\, formed as the complement.

- log.p:

  Logical of length 1. When `TRUE` the logarithm of the probability is
  returned. Defaults to `FALSE`. The logarithm is taken after the clamp,
  so a quantile far into the light tail returns `-Inf`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of probabilities in \\\[0, 1\]\\, or their logarithms
with `log.p = TRUE`, of length
`max(length(q), length(mu), length(sigma), length(alpha))`.

## Notation

\\\Phi\\ is the standard Gaussian distribution function and \\T(h, a) =
(2\pi)^{-1}\int_0^a e^{-h^2(1+x^2)/2}(1+x^2)^{-1}dx\\ is Owen's T.
\\\alpha\\ is the shape.

## See also

[`numericals7::owen_t()`](https://statmodels7.github.io/numericals7/reference/owen_t.html)
for the quadrature,
[`distrib_pdf.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.SkewNormal1Distrib.md)
for the density it integrates,
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
for its parameter derivatives, and
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
q <- c(-2, -0.5, 0.5, 2)
th <- list(mu = 0, sigma = 1, alpha = 3)

# Owen's identity against a direct quadrature of the density.
rbind(owen = distrib_cdf(d, q, th),
      quadrature = vapply(q, function(u)
        integrate(function(v) distrib_pdf(d, v, th), -Inf, u)$value, 0))
#>                    [,1]        [,2]      [,3]      [,4]
#> owen       5.089200e-12 0.006369453 0.3892944 0.9544997
#> quadrature 5.089126e-12 0.006369453 0.3892944 0.9544997

# At shape zero, T(z, 0) = 0 and the identity is the Gaussian's.
all.equal(distrib_cdf(d, q, list(mu = 0, sigma = 1, alpha = 0)), pnorm(q))
#> [1] TRUE

# The two tails sum to one.
distrib_cdf(d, q, th) + distrib_cdf(d, q, th, lower.tail = FALSE)
#> [1] 1 1 1 1
```
