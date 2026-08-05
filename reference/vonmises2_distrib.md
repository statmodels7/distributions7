# von Mises Distribution in the Mean Resultant Length

Creates a von Mises distribution object parametrized by its mean
direction and its **mean resultant length** \\\rho = A(\kappa)\\, which
lives in \\(0, 1)\\.

## Usage

``` r
vonmises2_distrib(
  link_mu = bounded_link(lwr = -pi, upr = pi),
  link_rho = logit_link()
)
```

## Arguments

- link_mu:

  Link function for the mean direction. Defaults to a link bounded to
  \\(-\pi, \pi)\\.

- link_rho:

  Link function for the resultant length. Defaults to the logit, the
  natural link onto \\(0, 1)\\.

## Value

An S7 object of class
[`VonMises2Distrib`](https://statmodels7.github.io/distributions7/reference/VonMises2Distrib.md).

## Details

The concentration \\\kappa\\ of
[`vonmises1_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)
is unbounded and hard to read; the resultant length is bounded, is the
quantity circular statistics reports, and is one minus the circular
variance. The two are related by \\\rho = I_1(\kappa)/I_0(\kappa)\\, a
strictly increasing bijection.

That bijection has no closed inverse, which is why this is a family of
its own rather than a
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
of the other: \\\kappa\\ is obtained by root finding on \\\log\kappa\\,
and its four derivatives come from the inverse function rule applied to
\\A' \dots A''''\\, which the Bessel recurrences give from the same two
evaluations \\A\\ already needs.

The map touches the second parameter only, so the chain rule is the
one-variable one and the derivatives are exact. The expected information
is closed form and the two parameters are orthogonal, as they are in the
concentration parametrization.

**The moments are not the parameters.**
[`mean`](https://rdrr.io/r/base/mean.html) returns the ordinary
expectation of \\Y\\ on \\\[-\pi, \pi)\\, which differs from \\\mu\\
whenever \\\mu \ne 0\\; see
[`vonmises1_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

## See also

[`vonmises1_distrib`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md)

## Examples

``` r
d <- vonmises2_distrib()
theta <- list(mu = 0.5, rho = 0.7)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.07974218 0.40483719 0.40483719

# rho is bounded, which is what makes it readable
d@params_bounds$rho
#> [1] 0 1
```
