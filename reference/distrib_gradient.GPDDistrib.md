# Generalized Pareto Score

Computes the two first derivatives of the log-density in closed form, in
a compiled kernel. With \\z = y/\sigma\\, \\t = 1 + \xi z\\ and \\u =
z/t\\, \$\$\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(\xi+1)u -
1}{\sigma}, \qquad \dfrac{\partial \ell}{\partial \xi} = \dfrac{\log
t}{\xi^2} - \left(1 + \dfrac{1}{\xi}\right)u.\$\$

The shape component is written this way only away from zero. Both of its
terms blow up as \\\xi \to 0\\ and their difference has the finite limit
\\z^2/2 - z\\, so the kernel evaluates it through the analytic function
\\\Lambda(u) = \log(1+u)/u\\ instead, whose derivatives come from a
recursion above \\\|u\| = 1/2\\ and from a Taylor series below it.
Measured at \\z = 1\\, the component reads \\-0.4967\\, \\-0.49997\\,
\\-0.5000000\\ at \\\xi = 10^{-2}, 10^{-4}, 10^{-8}\\ against the limit
\\-1/2\\.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation to
  the link scale is applied in the generic's body, so this method always
  returns the parameter scale.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list of two numeric vectors, `sigma` and `xi`, each of the
length of the recycled inputs.

## Notation

\\\ell\\ is the log-density of one observation, \\\sigma \> 0\\ the
scale, \\\xi\\ the shape, \\z = y/\sigma\\ and \\t = 1 + \xi z\\.

## See also

[`distrib_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GPDDistrib.md)
for the second derivatives,
[`distrib_expected_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GPDDistrib.md)
for Smith's closed form, and
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
y <- c(0.2, 1, 4)
th <- list(sigma = 1.5, xi = 0.3)
g <- distrib_gradient(d, y, th)

# Against a central difference of the log-density.
eps <- 1e-6
rbind(analytic = g$xi,
      numeric = (distrib_pdf(d, y, list(sigma = 1.5, xi = 0.3 + eps),
                             log = TRUE) -
                 distrib_pdf(d, y, list(sigma = 1.5, xi = 0.3 - eps),
                             log = TRUE)) / (2 * eps))
#>                [,1]       [,2]      [,3]
#> analytic -0.1197699 -0.3816123 0.1112099
#> numeric  -0.1197699 -0.3816123 0.1112099

# The shape component has a finite limit at zero that its written form
# does not: both of its terms diverge and the difference does not.
vapply(c(1e-2, 1e-4, 1e-8, 1e-14),
       function(x) distrib_gradient(d, 1, list(sigma = 1, xi = x))$xi, 0)
#> [1] -0.4966915 -0.4999667 -0.5000000 -0.5000000

# ...and the limit is z^2/2 - z, which at z = 1 is -0.5.
distrib_gradient(d, 1, list(sigma = 1, xi = 0))$xi
#> [1] -0.5
```
