# The Pieces a Skew t Evaluates From

Assembles the standardized variable, the argument of the tilting
distribution function and the six scalar functions every closed-form
derivative of the log-density is a combination of. Every method in this
family calls it once and then writes its own formula in terms of the
result.

## Usage

``` r
skewt_pieces(y, mu, sigma, alpha, nu)
```

## Arguments

- y:

  A numeric vector of observations.

- mu, sigma, alpha, nu:

  The four parameters, numeric vectors of length 1 or of the length of
  `y`. Nothing is validated: `sigma` and `nu` must be strictly positive.

## Value

A named list of ten numeric vectors, each of the length of the recycled
inputs: `z` the standardized variable, `w` the tilting argument, `c` the
factor relating them, `a` and `da` for \\A\\ and \\A'\\, `e` for \\E\\,
`b` and `db` for \\B\\ and \\B'\\, and `q` and `dq` for \\Q\\ and
\\Q'\\.

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

\\Q\\ is formed as `exp(dt(log = TRUE) - pt(log.p = TRUE))` because both
factors underflow together in the far left tail while the ratio stays
finite. It matters as the degrees of freedom grow and the \\t\\ tail
approaches the Gaussian's: measured at \\w = -60\\ with \\m = 2000\\,
the log route returns 21.4345 and `dt(w, m)/pt(w, m)` returns `NaN`.

Nothing here involves \\\nu\\ by differentiation; the components in
\\\nu\\ are obtained separately, by a stencil.

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
