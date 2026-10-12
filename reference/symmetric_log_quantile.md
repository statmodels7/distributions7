# The Quantile of a Symmetric Family on the Log Scale

Returns the quantile of a law symmetric about its first parameter,
solved on the logarithm of the smaller tail: the root of \\\log F(q) =
\log p\\ is found below the location by
[`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) and reflected
about it when the upper tail is the smaller one, so that neither \\1 -
p\\ nor \\\exp(\log p)\\ is formed. The two pseudo-Huber families use
it.

## Usage

``` r
symmetric_log_quantile(distrib, p, theta, lower.tail, log.p, spread)
```

## Arguments

- distrib:

  A distribution object symmetric about its first parameter.

- p:

  A numeric vector of probabilities, or their logarithms.

- theta:

  A named list of parameter values.

- lower.tail, log.p:

  Logical of length 1.

- spread:

  A function of the parameter list that sets the width of the first
  bracket, which is doubled until it contains the root.

## Value

A numeric vector of quantiles.

## Examples

``` r
d <- pseudohuber_distrib()
distributions7:::symmetric_log_quantile(d, 1e-30, list(mu = 0, sigma = 1, nu = 2),
                                        FALSE, FALSE, function(th) th[[2]])
#> [1] 69.18132
```
