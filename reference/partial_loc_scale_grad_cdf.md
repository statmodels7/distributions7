# CDF Gradient When Only Some Parameters Are Location-Scale

For a family that is location-scale in \\(\mu, \sigma)\\ but carries a
further shape parameter: the two location-scale directions in closed
form, the shape direction by finite differences.

## Usage

``` r
partial_loc_scale_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters.

- lower.tail:

  Logical; whether the lower tail is wanted.

- log:

  Logical; whether derivatives of the log probability are wanted.

## Value

A named list of gradient component vectors.

## Details

The Student t and the pseudo-Huber are the two. Their shape direction is
a derivative of a hypergeometric-type integral with no elementary form,
so it is differenced; the other two are exact.

For the pseudo-Huber this is an **accuracy** gain rather than merely a
speed one. Its cdf is itself a quadrature, so differencing it is good to
only about `1e-6`, whereas \\\partial F/\partial \mu = -f\\ is exact.
