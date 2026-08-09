# Beta Distribution in Its Shapes

Creates a beta distribution object in its canonical parametrization, the
two shapes \\\alpha\\ and \\\beta\\.

## Usage

``` r
beta2_distrib(link_alpha = log_link(), link_beta = log_link())
```

## Arguments

- link_alpha:

  Link function for \\\alpha\\. Defaults to the log.

- link_beta:

  Link function for \\\beta\\. Defaults to the log.

## Value

An S7 object of class
[`Beta2Distrib`](https://statmodels7.github.io/distributions7/reference/Beta2Distrib.md).

## Details

The same law as
[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md),
which carries the mean and a precision: \\\alpha = \mu\varphi\\ and
\\\beta = (1-\mu)\varphi\\. The mean parametrization is the one a
regression wants; this one is the one the family is usually written in
and the one a conjugate analysis produces, the beta being conjugate for
a binomial probability.

The data enter the log-density only through \\\log y\\ and
\\\log(1-y)\\, both of which are linear in the parameters. Every
derivative beyond the first is therefore free of the data, so the
observed and the expected ones coincide at orders two, three and four,
and Fisher scoring and Newton's method take the same step on the
parameter scale.

## The distribution

\$\$f(y) = \frac{y^{\alpha-1}(1-y)^{\beta-1}}{B(\alpha, \beta)}\$\$ on
\\y \in (0, 1)\\.

\$\$\mathbb{E}\[Y\] = \frac{\alpha}{\alpha+\beta}, \qquad
\operatorname{Var}(Y) =
\frac{\alpha\beta}{(\alpha+\beta)^{2}(\alpha+\beta+1)}\$\$

## See also

[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md)

## Examples

``` r
d <- beta2_distrib()
theta <- list(alpha = 2, beta = 5)
distrib_pdf(d, c(0.1, 0.3, 0.7), theta)
#> [1] 1.9683 2.1609 0.1701
c(mean = mean(d, theta), variance = variance(d, theta))
#>      mean  variance 
#> 0.2857143 0.0255102 

# the same law as beta1 with mu = alpha/(alpha+beta), phi = alpha+beta
distrib_pdf(beta1_distrib(), 0.3, list(mu = 2 / 7, phi = 7))
#> [1] 2.1609
```
