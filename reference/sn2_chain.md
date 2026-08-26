# Derivatives of the Skew Normal in Its Centered Parametrization

Carries
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)'s
derivatives into the centered coordinates through the partition sum of
[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md),
at any order from one to four, observed or expected. The parent supplies
the derivatives in \\(\xi, \omega, \alpha)\\ and
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
supplies the map's partial derivatives; the partition sum assembles
them.

## Usage

``` r
sn2_chain(distrib, y, theta, order, expected = FALSE)
```

## Arguments

- distrib:

  A
  [SkewNormal2Distrib](https://statmodels7.github.io/distributions7/reference/SkewNormal2Distrib.md)
  object.

- y:

  A numeric vector of observations.

- theta:

  A list with `mu`, `sigma` and `gamma1`. It is aligned here, so it may
  be given in any order and by name.

- order:

  A single integer, 1, 2, 3 or 4: the derivative order.

- expected:

  Logical of length 1. When `TRUE` the parent's **expected** derivatives
  are carried instead of the observed ones. Defaults to `FALSE`.

## Value

A named list of numeric vectors, one per distinct component of the
requested order, named in the centered parameters. At order 2 the caller
subsets it by
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
to fix the ordering.

## The point that is excluded

The map to the direct parametrization runs through \\c =
\sqrt\[3\]{2\gamma_1/(4-\pi)}\\, whose derivative grows like
\\\gamma_1^{-2/3}\\. At zero skewness the map is not differentiable and
the chain rule is asked for a quantity that does not exist.

The first derivatives of the log-density survive the limit: the map's
divergent factor cancels and they approach a finite value from both
sides. The second ones do not. Measured at \\y = 0.5\\, \\\mu = 0\\,
\\\sigma = 1\\, the score in \\\gamma_1\\ runs \\-0.2152, -0.2257,
-0.2284, -0.2290\\ at \\\gamma_1 = 10^{-2}, 10^{-4}, 10^{-6}, 10^{-8}\\,
while \\\partial^2\ell/\partial\gamma_1^2\\ runs \\0.29, 11.5, 253,
5451\\ over the same values, a factor of 4.642 per decade against
\\10^{2/3} = 4.6416\\.

Zero skewness is therefore rejected here, where the map is used and the
reason can be named, with a message that points at
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
whose derivatives at \\\alpha = 0\\ are ordinary numbers.

## Where the cancellation runs out of digits

The **expected** information stays finite as \\\gamma_1 \to 0\\ and
tends to \\1/6\\ in its own component, but it is computed as a
difference of terms of size \\\gamma_1^{-2/3}\\. Measured, it holds to
seven figures down to \\\gamma_1 = 10^{-8}\\ (0.16666782 against
\\1/6\\), loses three by \\10^{-10}\\ and is **negative at**
\\10^{-12}\\, which no information can be. That is a limit of
double-precision arithmetic on a parameter value no fit visits, and it
is why the near-symmetric case is better handled by the direct
parametrization.

## Errors

Signals an error when any element of `gamma1` is exactly zero, naming
the cube root as the cause and
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
as the alternative.

## See also

[`chain_derivatives()`](https://statmodels7.github.io/distributions7/reference/chain_derivatives.md)
for the partition sum,
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
for the map's derivatives, and
[`distrib_gradient.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewNormal2Distrib.md)
for the method that calls this.

## Examples

``` r
d <- skewnormal2_distrib()
y <- c(-1, 0.3, 1.7)
th <- list(mu = 0, sigma = 1, gamma1 = 0.5)

# Order one, against the method that wraps it.
all.equal(distributions7:::sn2_chain(d, y, th, 1L),
          distrib_gradient(d, y, th))
#> [1] TRUE

# The score in the skewness stays of order one as the map's Jacobian
# diverges; the curvature does not.
t(vapply(10^-c(2, 4, 6, 8), function(g) {
  p <- list(mu = 0, sigma = 1, gamma1 = g)
  c(gamma1 = g,
    score = distrib_gradient(d, 0.5, p)$gamma1,
    curvature = distrib_hessian(d, 0.5, p)$gamma1_gamma1)
}, numeric(3)))
#>      gamma1      score   curvature
#> [1,]  1e-02 -0.2152070    0.288818
#> [2,]  1e-04 -0.2256687   11.513887
#> [3,]  1e-06 -0.2284080  252.620070
#> [4,]  1e-08 -0.2290031 5451.052094

# Zero skewness is rejected rather than approximated.
tryCatch(distrib_gradient(d, 0, list(mu = 0, sigma = 1, gamma1 = 0)),
         error = function(e) "rejected, as documented")
#> [1] "rejected, as documented"
```
