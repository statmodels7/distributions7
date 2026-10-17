# The Truncation Constant and Lower Tail

Returns the tail below the interval, \\F(L^-)\\, and the retained mass
\\Z = P(L \le Y \le U)\\. Every method of a truncated distribution
divides by \\Z\\, and
[`trunc_cdf()`](https://statmodels7.github.io/distributions7/reference/trunc_cdf.md)
and
[`trunc_quantile()`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md)
shift by \\F(L^-)\\ or by the upper tail \\S(U) = P(Y \> U)\\ as well.

## Usage

``` r
trunc_constants(distrib, theta, tails = FALSE)
```

## Arguments

- distrib:

  A truncated distribution object, of either class.

- theta:

  A named list of the parent's parameters. A component may vary by
  observation.

- tails:

  Logical, default `FALSE`. When `TRUE` the result also carries `Fu`,
  \\F(U)\\, and the two survival values `Sl`, \\P(Y \ge L)\\, and `Su`,
  \\P(Y \> U)\\, which
  [`trunc_cdf()`](https://statmodels7.github.io/distributions7/reference/trunc_cdf.md)
  and
  [`trunc_quantile()`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md)
  read.

## Value

A named list with `Fl`, the tail below the interval, and `Z`, the
retained mass, and with `Fu`, `Sl` and `Su` when `tails = TRUE`. Each is
a numeric vector following the length of `theta`, or the single value
`0` or `1` for an infinite endpoint.

## Details

Both endpoints are INCLUDED in the truncated support, so any mass
sitting exactly on the lower one is added back: \\F(L^-) = F(L) - P(Y =
L)\\. That correction belongs to the ATOM case, not to the discrete
case, so it goes through
[`parent_mass_at()`](https://statmodels7.github.io/distributions7/reference/parent_mass_at.md).
An infinite endpoint contributes \\0\\ or \\1\\ without a call on the
parent.

Both quantities are vectorized in \\\theta\\, so a parameter varying by
observation gives one constant per observation.

Where \\F(L^-)\\ exceeds one half the interval lies in the upper tail,
and \\F(U) - F(L^-)\\ keeps no digits there. \\Z\\ is then computed from
the parent's survival function as \\S(L^-) - S(U)\\, with \\S(L^-) = P(Y
\ge L) = P(Y \> L) + P(Y = L)\\. An error that names the interval and
the computed mass is signaled when that mass is not finite or not
positive, which happens only when it underflows.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, which both belong to the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`parent_mass_at()`](https://statmodels7.github.io/distributions7/reference/parent_mass_at.md)
for the endpoint correction,
[`trunc_mass_derivs()`](https://statmodels7.github.io/distributions7/reference/trunc_mass_derivs.md)
for the derivatives of \\Z\\, and
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md).

## Examples

``` r
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)

cs <- distributions7:::trunc_constants(tn, theta)
unlist(cs)
#>        Fl         Z 
#> 0.1393302 0.7823795 
c(Z = cs$Z, direct = pnorm(2, 0.3, 1.2) - pnorm(-1, 0.3, 1.2))
#>         Z    direct 
#> 0.7823795 0.7823795 

# A discrete parent: the mass at the lower endpoint is added back, so the
# zero-truncated Poisson keeps everything except dpois(0, mu).
ztp <- truncated(poisson_distrib(), lower = 1)
cz <- distributions7:::trunc_constants(ztp, list(mu = 2))
c(Z = cz$Z, direct = 1 - dpois(0, 2), Fl = cz$Fl)
#>         Z    direct        Fl 
#> 0.8646647 0.8646647 0.1353353 

# The function is vectorized in theta.
distributions7:::trunc_constants(tn, list(mu = c(0, 0.5), sigma = 1.2))$Z
#> [1] 0.7498813 0.7887005

# An interval carrying no mass is reported, not returned.
far <- truncated(gaussian1_distrib(), lower = 100, upper = 200)
try(distributions7:::trunc_constants(far, theta))
#> Error : The truncation interval [100, 200] carries no probability under these parameter values (computed mass 0). A truncated distribution is not defined there.
```
