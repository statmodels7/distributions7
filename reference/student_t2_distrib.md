# Student t Distribution in the Standard Deviation

Creates a Student t distribution object whose second parameter is the
standard deviation rather than the scale.

## Usage

``` r
student_t2_distrib(
  link_mu = identity_link(),
  link_sigma = log_link(),
  link_nu = bounded_link(lwr = 2)
)
```

## Arguments

- link_mu:

  Link function for the location. Defaults to the identity.

- link_sigma:

  Link function for the standard deviation. Defaults to the log.

- link_nu:

  Link function for the degrees of freedom. Defaults to a link bounded
  below at two.

## Value

An S7 object of class `StudentT2Distrib`, inheriting from
`continuous_distrib`, with `params` `c("mu", "sigma", "nu")` and
`link_params` the three links given here.

## Details

The scale of
[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
is not the standard deviation: the two differ by \\\sqrt{\nu/(\nu-2)}\\.
Here the scale is \$\$s_0 = \sigma\sqrt{\dfrac{\nu-2}{\nu}},\$\$ which
exists only for \\\nu \> 2\\, so the degrees of freedom are bounded
below at two. Every derivative in \\(\mu, \sigma, \nu)\\ to order five,
the expected information and its derivatives, and the derivatives in the
response are closed forms, each order in its own compiled kernel (see
[`distrib_gradient.StudentT2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT2Distrib.md)).
This is `TF2` in gamlss.

## The distribution

\$\$f(y) = \frac{1}{s_0}\\t\_{\nu}\\\left(\frac{y-\mu}{s_0}\right),
\qquad s_0 = \sigma\sqrt{\frac{\nu-2}{\nu}}\$\$ on \\y \in \mathbb{R}\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2}\$\$

## References

Rigby, R. A. and Stasinopoulos, D. M. (2005). Generalized additive
models for location, scale and shape. *Journal of the Royal Statistical
Society, Series C* 54, 507-554.

## See also

[`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md);
[StudentT2Distrib](https://statmodels7.github.io/distributions7/reference/StudentT2Distrib.md)
for the class.

## Examples

``` r
d <- student_t2_distrib()
theta <- list(mu = 0, sigma = 2, nu = 8)
variance(d, theta)
#> [1] 4
```
