# The Step a Skew t Differences the Degrees of Freedom With

Returns the finite-difference step used for the derivatives in \\\nu\\,
`pmax(1e-3 * abs(nu), 1e-6)`: relative to \\\nu\\ itself, so that the
same number of significant digits is differenced at every scale, and
floored so that it stays a number for a degrees of freedom near zero.

## Usage

``` r
skewt_nu_step(nu)
```

## Arguments

- nu:

  A numeric vector of degrees of freedom.

## Value

A numeric vector of steps, of the length of `nu`.

## Details

The relative step \\10^{-3}\\ is measured rather than assumed. Swept
over \\\nu\\ from 2 to 30 and sample sizes from 500 to 4000, it is where
the truncation error of
[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)'s
five-point stencil has fallen to the level of the rounding error and the
two are balanced. A smaller step is dominated by rounding, which the
stencil amplifies by \\18/(12h)\\, and a larger one by truncation, which
grows as \\h^4\\.

The choice shows up in the stopping rule of a fit.
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
takes its rule from the method it is given, and
[`optimizers7::crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.html)
tests the score **per observation** at a default tolerance of
\\10^{-6}\\. Measured on samples of 500 to 4000 the run converges at
between \\3\times10^{-10}\\ and \\9\times10^{-9}\\ per observation, so
this step leaves two orders of room. A three-point stencil would not:
its bias does not cancel over the sum, and a run would spend its whole
budget reporting failure at the maximum.

## See also

[`fd5_first()`](https://statmodels7.github.io/distributions7/reference/fd5_first.md)
and its siblings, which consume the step, and
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md),
the first method that needs it.

## Examples

``` r
distributions7:::skewt_nu_step(c(2, 6, 30, 1e-5))
#> [1] 2e-03 6e-03 3e-02 1e-06

# Relative above the floor, absolute below it.
nu <- c(1e-4, 1e-3, 1e-2, 1, 100)
rbind(nu = nu, step = distributions7:::skewt_nu_step(nu))
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> nu   1e-04 1e-03 1e-02 1.000 100.0
#> step 1e-06 1e-06 1e-05 0.001   0.1
```
