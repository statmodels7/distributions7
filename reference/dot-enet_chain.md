# The Elastic Net's Derivatives of a Given Order

Assembles the log-density's derivatives in \\(\mu, a, c)\\ and carries
them onto \\(\mu, \lambda, \alpha)\\.

## Usage

``` r
.enet_chain(y, theta, order)
```

## Arguments

- y:

  A numeric vector of observations.

- theta:

  A list containing `mu`, `lambda` and `alpha`.

- order:

  The derivative order, 1 to 4.

## Value

A named list of components.
