# Probability Density Function

Evaluates the probability density function (PDF) or probability mass
function (PMF).

## Usage

``` r
distrib_pdf(distrib, y, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- y:

  A numeric vector of observations.

- theta:

  A named list (or named numeric vector) of distribution parameters. If
  unnamed, parameters are taken in the order of `distrib@params`.

- ...:

  Additional arguments passed to the specific method (e.g., `log`).
