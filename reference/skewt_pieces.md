# The Pieces a Skew t Evaluates From

Assembles the standardized variable, the argument of the tilting
distribution function and the seven scalar functions \\A\\, \\A'\\,
\\E\\, \\B\\, \\B'\\, \\Q\\ and \\Q'\\, of which the closed-form
derivatives of the log-density in \\(\mu, \sigma, \alpha)\\ are
combinations.
[`distrib_grad_y.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.SkewTDistrib.md)
and
[`distrib_hess_y.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_y.SkewTDistrib.md)
evaluate their formulas from the result; the derivatives in the
parameters come from the compiled kernels.

## Usage

``` r
skewt_pieces(y, mu, sigma, alpha, nu)
```

## Arguments

- y:

  A numeric vector of observations.

- mu, sigma, alpha, nu:

  The four parameters, numeric vectors of length 1 or of the length of
  `y`. The arguments are not validated; `sigma` and `nu` must be
  strictly positive.

## Value

A named list of ten numeric vectors, each of the length implied by the
arguments on which it depends: `z` the standardized variable, `w` the
tilting argument, `c` the factor relating them, `a` and `da` for \\A\\
and \\A'\\, `e` for \\E\\, `b` and `db` for \\B\\ and \\B'\\, and `q`
and `dq` for \\Q\\ and \\Q'\\.

## Details

With \\z = (y-\mu)/\sigma\\, \\m = \nu + 1\\, \\s = \nu + z^2\\ and \\c
= \sqrt{m/s}\\, the tilting argument is \\w = \alpha z c\\. The
functions returned are \$\$A = \dfrac{\partial}{\partial z}\log
t\_\nu(z) = -\dfrac{m z}{s}, \qquad A' = -\dfrac{m(\nu - z^2)}{s^2},\$\$
\$\$E = \dfrac{\nu\sqrt{m}}{s^{3/2}}, \qquad B = \dfrac{\partial
w}{\partial z} = \alpha E, \qquad B' =
-\dfrac{3\alpha\nu\sqrt{m}\\z}{s^{5/2}},\$\$ and \\Q = t_m(w)/T_m(w)\\
with \\Q' = Q\\-(m+1)w/(m + w^2) - Q\\\\, the last from differentiating
the quotient.

\\Q\\ is formed as `exp(dt(log = TRUE) - pt(log.p = TRUE))`, because
\\t_m(w)\\ and \\T_m(w)\\ underflow together in the far left tail while
their ratio stays finite. The direct quotient is `NaN` there, in
particular when the degrees of freedom are large and the \\t\\ tail is
close to the Gaussian one. Where \\\|z\| \> 10^{20}\\ the functions are
written in \\\sqrt{s}\\, \\z/\sqrt{s}\\ and \\\nu/s\\, so that no power
of \\s\\ overflows.

This function does not differentiate in \\\nu\\. The derivative kernels
evaluate the same quantities in bounded variables; see
[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md).

## Notation

\\t\_\nu\\ and \\T\_\nu\\ are the standard Student \\t\\ density and
distribution function on \\\nu\\ degrees of freedom, \\\mu\\ the
location, \\\sigma\\ the scale, \\\alpha\\ the shape and \\\nu\\ the
degrees of freedom.

## See also

[`distrib_gradient.SkewTDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.SkewTDistrib.md),
which writes the score in these terms, and
[`skewt_distrib()`](https://statmodels7.github.io/distributions7/reference/skewt_distrib.md)
for the family.

## Examples

``` r
p <- distributions7:::skewt_pieces(c(-1.5, 0.4, 2.1), 0, 1, 3, 6)
names(p)
#>  [1] "z"  "w"  "c"  "a"  "da" "e"  "b"  "db" "q"  "dq"

# The score in the location is -(A + QB)/sigma.
d <- skewt_distrib()
all.equal(-(p$a + p$q * p$b) / 1,
          distrib_gradient(d, c(-1.5, 0.4, 2.1),
                           list(mu = 0, sigma = 1, alpha = 3, nu = 6))$mu)
#> [1] TRUE

# Q survives a tail where the direct quotient does not.
c(log_route = exp(dt(-60, df = 2000, log = TRUE) -
                  pt(-60, df = 2000, log.p = TRUE)),
  direct = dt(-60, df = 2000) / pt(-60, df = 2000))
#> log_route    direct 
#>  21.43451       NaN 
```
