# Derivatives of log1p(u)/u

Returns \\\Lambda(u) = \log(1+u)/u\\ and its first four derivatives, one
vector per order. The function is analytic at the origin, with
\\\Lambda(0) = 1\\, and it is the device that removes every division by
the generalized Pareto's shape from that family's survival function.

## Usage

``` r
gpd_lambda_derivs(u)
```

## Arguments

- u:

  A numeric vector, greater than \\-1\\. Values at or below \\-1\\ are
  outside the support and are masked out by the caller before they reach
  here.

## Value

A list of five numeric vectors the length of `u`, orders 0 to 4. At \\u
= 0\\ they are 1, \\-1/2\\, \\2/3\\, \\-3/2\\ and \\24/5\\.

## Two routes, and where they change over

Differentiating \\u\Lambda = \log(1+u)\\ gives the recursion
\$\$u\\\Lambda^{(r)} + r\\\Lambda^{(r-1)} =
\frac{(-1)^{r-1}(r-1)!}{(1+u)^{r}},\$\$ which is exact away from the
origin and useless at it: it divides by \\u\\ and subtracts two nearly
equal quantities. Measured against the Taylor series, its fourth
derivative is wrong by a factor of \\10^{39}\\ at \\u = 10^{-14}\\, by
1.7 at \\10^{-4}\\ and by \\3\times10^{-8}\\ at \\10^{-2}\\.

Below \\\|u\| = 1/2\\ the series \$\$\Lambda^{(r)}(u) = \sum\_{m \ge r}
(-1)^{m}\frac{m!}{(m-r)!} \frac{u^{m-r}}{m+1}\$\$ is used instead, where
the two agree to \\10^{-16}\\. Its truncation is set so that the \\m^4\\
weight at the switch point stays under the rounding.

## Why the expression is arranged this way

Differentiating \\L = -\log(1+\xi q/\sigma)/\xi\\ directly gives terms
in \\\xi^{-1-m}\\ that cancel only in the limit, which needs a guard and
is fragile whatever the guard. Writing \\L = -(q/\sigma)\Lambda(u)\\
puts the whole removable singularity inside one univariate function, and
the exponential limit \\\xi \to 0\\ becomes an ordinary point of the
formula.

## Notation

\\u = \xi q/\sigma\\ with \\\xi\\ the shape and \\\sigma\\ the scale.

## See also

[`gpd_surv_pieces()`](https://statmodels7.github.io/distributions7/reference/gpd_surv_pieces.md),
the one consumer;
[`distrib_grad_cdf.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.GPDDistrib.md)
for the family.

## Examples

``` r
# The limits at the origin, reached through the series branch.
vapply(distributions7:::gpd_lambda_derivs(1e-14), function(v) v[1],
       numeric(1))
#> [1]  1.0000000 -0.5000000  0.6666667 -1.5000000  4.8000000
```
