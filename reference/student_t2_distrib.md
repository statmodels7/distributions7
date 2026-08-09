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

A reparametrized distribution object.

## Details

The scale of
[`student_t1_distrib`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md)
is not the standard deviation: the two differ by \\\sqrt{\nu/(\nu-2)}\\.
Here the map is \$\$\sigma\_{\text{scale}} =
\sigma\sqrt{\dfrac{\nu-2}{\nu}},\$\$ which exists only for \\\nu \> 2\\,
and the constructor bounds \\\nu\\ there rather than letting the map
return a complex number several frames down. This is `TF2` in gamlss.

The restriction is the point rather than a limitation: a family
parametrized by a standard deviation is a family whose standard
deviation exists.

## The distribution

\$\$f(y) = \frac{1}{s_0}\\t\_{\nu}\\\left(\frac{y-\mu}{s_0}\right),
\qquad s_0 = \sigma\sqrt{\frac{\nu-2}{\nu}}\$\$ on \\y \in \mathbb{R}\\.

\$\$\mathbb{E}\[Y\] = \mu, \qquad \operatorname{Var}(Y) = \sigma^{2}\$\$

## See also

[`student_t1_distrib`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md),
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)

## Examples

``` r
d <- student_t2_distrib()
theta <- list(mu = 0, sigma = 2, nu = 8)
variance(d, theta)
#> [1] 4
```
