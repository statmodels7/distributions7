# Dirichlet Density

Computes the Dirichlet density \$\$f(y) =
\dfrac{\Gamma(\phi)}{\prod\_{j=1}^{p} \Gamma(\alpha_j)} \prod\_{j=1}^{p}
y_j^{\alpha_j - 1}, \qquad \alpha = \phi\mu,\$\$ one value per row of
`y`, with respect to the \\(p-1)\\-dimensional measure the simplex
carries.

The support is tested **before** the logarithm is taken. A negative
coordinate would make [`base::log()`](https://rdrr.io/r/base/Log.html)
warn, and a density evaluated off its support is 0, so a row that is not
on the simplex returns 0 without a numerical complaint.

## Arguments

- distrib:

  A `DirichletDistrib` object, from
  [`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md).

- y:

  A numeric matrix with one row per observation and \\p\\ columns, each
  row strictly positive and summing to one within `1e-8`. A single
  observation may be given as a plain vector and is read as one row. A
  row failing either test gives a density of 0, or `-Inf` with
  `log = TRUE`.

- theta:

  A named list of parameters on the parameter scale: the mean's free
  values followed by `phi`, each of length 1. A parameter may not vary
  by observation here, the mean being one point of the simplex for the
  whole sample.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities of length `nrow(y)`, one per observation.

## Notation

\\\mu\\ is the mean vector, a point of the simplex; \\\phi \> 0\\ the
concentration; \\\alpha = \phi\mu\\ the shapes; and \\p\\ the number of
coordinates.

## See also

[`distrib_gradient.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.DirichletDistrib.md)
for the derivatives of the log-density,
[`distrib_rng.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.DirichletDistrib.md)
for draws,
[`mv_marginal.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.DirichletDistrib.md)
for a coordinate's beta marginal, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- dirichlet_distrib(3)
th <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 12)
mu <- mv_location(d, th)

# The density at four draws, against the formula written out.
set.seed(1)
Y <- distrib_rng(d, 4, th)
al <- 12 * mu
all.equal(distrib_pdf(d, Y, th),
          as.numeric(exp(lgamma(12) - sum(lgamma(al))) *
                       apply(Y, 1, function(r) prod(r^(al - 1)))))
#> [1] TRUE

# A row that does not sum to one, or has a negative coordinate, is off the
# support and returns zero rather than warning.
distrib_pdf(d, rbind(c(0.5, 0.5, 0.5), c(-0.1, 0.6, 0.5)), th)
#> [1] 0 0

# A single observation may be given as a vector.
distrib_pdf(d, c(0.4, 0.25, 0.35), th)
#> [1] 9.289035
```
