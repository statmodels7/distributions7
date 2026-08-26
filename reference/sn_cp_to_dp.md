# From the Centered Parameters to the Direct Ones

Maps \\(\mu, \sigma, \gamma_1)\\, the mean, the standard deviation and
the skewness, to \\(\xi, \omega, \alpha)\\, the location, the scale and
the shape that
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)
takes. Every probability function of the centered family calls it and
then delegates to the direct one.

## Usage

``` r
sn_cp_to_dp(mu, sigma, gamma1, s)
```

## Arguments

- mu, sigma, gamma1:

  The centered parameters: the mean, the standard deviation and the
  skewness, each a numeric vector. `sigma` must be positive and `gamma1`
  must lie strictly inside \\(-0.9952717, 0.9952717)\\; nothing is
  validated here.

- s:

  The sign of `gamma1`, \\\pm 1\\, taken by the caller from its plain
  value.

## Value

A named list with `mu`, `sigma` and `alpha`, the direct parameters, each
of the length of the recycled inputs. The names are the parent's, so the
result can be passed to
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md)'s
methods as they stand; `mu` there is the location \\\xi\\ and `sigma`
the scale \\\omega\\.

## Details

With \\b = \sqrt{2/\pi}\\, \$\$c = \mathrm{sign}(\gamma_1)
\left(\dfrac{2\|\gamma_1\|}{4-\pi}\right)^{1/3}, \qquad \mu_z =
\dfrac{c}{\sqrt{1+c^2}}, \qquad \delta = \dfrac{\mu_z}{b}, \qquad \alpha
= \dfrac{\delta}{\sqrt{1-\delta^2}},\$\$ and then \\\omega =
\sigma/\sqrt{1-\mu_z^2}\\ and \\\xi = \mu - \omega\mu_z\\.

The caller supplies the sign, so the body reads
\\s\\(2s\gamma_1/(4-\pi))^{1/3}\\ with no
[`abs()`](https://rdrr.io/r/base/MathFun.html) in it. Away from zero the
sign is locally constant, so the expression is exact and differentiable
as written, and the derivative tables of
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
differentiate it directly.

## See also

[`sn2_theta()`](https://statmodels7.github.io/distributions7/reference/sn2_theta.md),
which supplies the sign and calls this;
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
for the derivative tables of the same map; and
[`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
for the family.

## Examples

``` r
# The mean and the standard deviation are not the location and the scale.
distributions7:::sn_cp_to_dp(0, 1, 0.5, 1)
#> $mu
#> [1] -1.052209
#> 
#> $sigma
#> [1] 1.451601
#> 
#> $alpha
#> [1] 2.173758
#> 

# At zero skewness the map is the identity on the first two, and the shape
# is zero: the Gaussian sits at the same point in both parametrizations.
distributions7:::sn_cp_to_dp(3, 2, 0, 1)
#> $mu
#> [1] 3
#> 
#> $sigma
#> [1] 2
#> 
#> $alpha
#> [1] 0
#> 

# Round trip: the direct parameters reproduce the centered moments.
dp <- distributions7:::sn_cp_to_dp(0, 1, 0.5, 1)
d1 <- skewnormal1_distrib()
c(mean = mean(d1, dp), sd = sqrt(variance(d1, dp)), skew = skewness(d1, dp))
#> mean   sd skew 
#>  0.0  1.0  0.5 
```
