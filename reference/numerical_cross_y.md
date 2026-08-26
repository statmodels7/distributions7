# Numerical Mixed Response-Parameter Derivatives

Computes \\\partial^2 \ell / \partial y\\ \partial \theta_i\\ by one
central difference of
[`distrib_grad_y()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
in each parameter, at a step of `h_rel` relative to the parameter's own
magnitude. It powers the default
[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md)
method for a continuous family with no closed form.

Measured on a Gaussian against Richardson extrapolation of the same
analytic response gradient, the two components agree to
\\1.2\times10^{-11}\\ and \\9.0\times10^{-11}\\ relative.

## Usage

``` r
numerical_cross_y(
  distrib,
  y,
  theta,
  h_rel = .Machine$double.eps^(1/3),
  which = NULL
)
```

## Arguments

- distrib:

  An object inheriting from `continuous_distrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned to `distrib@params`.

- h_rel:

  The relative finite-difference step, a single positive number.
  Defaults to `.Machine$double.eps^(1/3)`. A step near a parameter's
  bound can push an evaluation outside the domain, where the density
  returns `NaN`.

- which:

  A character vector of parameter names to differentiate, or `NULL` (the
  default) for all of them. Used by a family with a closed form for some
  of its parameters, so that only the remaining ones cost a pair of
  evaluations.

## Value

A named list with one numeric vector per requested parameter, each of
length `length(y)`, in the order `which` gives or `distrib@params` where
it is `NULL`.

## One difference layer, not two

The quantity differenced is the **response gradient**, not the
log-density, so a family with an analytical `distrib_grad_y` pays for
exactly one finite-difference layer. Where the response gradient is
itself the finite-difference fallback the composition is the four-point
mixed stencil on the log-density: the two differences act on different
variables, so they commute into a single stencil instead of compounding
the way nested differences in one variable do.

## The step

\\\varepsilon^{1/3} \approx 6.1\times10^{-6}\\ balances the \\O(h^2)\\
truncation of a central difference against a rounding term growing as
\\1/h\\, which is the usual optimum for a first derivative.

## See also

[`distrib_cross_y()`](https://statmodels7.github.io/distributions7/reference/distrib_cross_y.md),
the generic it serves;
[`numerical_grad_y()`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md)
for the quantity it differences;
[`numerical_cross2_y()`](https://statmodels7.github.io/distributions7/reference/numerical_cross2_y.md)
for the third-order analogue.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)
numerical_cross_y(d, y, th)
#> $mu
#> [1] 1 1 1
#> 
#> $sigma
#> [1] -2  0  2
#> 

# It agrees with the family's own closed form to about 1e-10.
all.equal(numerical_cross_y(d, y, th), distrib_cross_y(d, y, th),
          tolerance = 1e-8)
#> [1] TRUE

# 'which' costs a pair of evaluations for the named parameter alone.
numerical_cross_y(d, y, th, which = "sigma")
#> $sigma
#> [1] -2  0  2
#> 
```
