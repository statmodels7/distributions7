# The Jets of the Centred-to-Direct Map

The map run on jets, giving every partial derivative of the direct
parameters with respect to the centred ones, to fourth order and
exactly.

## Usage

``` r
sn2_jets(theta, n)
```

## Arguments

- theta:

  A list with `mu`, `sigma` and `gamma1`.

- n:

  The number of observations.

## Value

A list as
[`reparam_jets`](https://statmodels7.github.io/distributions7/reference/reparam_jets.md)
returns.

## Details

The sign of \\\gamma_1\\ is read off the plain value before the jets are
seeded, so the cube root is taken of a positive quantity and the result
negated. This is the step
[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
cannot take, a jet having no sign of its own.

## See also

[`skewnormal2_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md)
