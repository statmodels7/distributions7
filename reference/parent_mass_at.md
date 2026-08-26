# Probability the Parent Puts on a Single Point

Returns \\P(Y = x)\\ under the parent: the mass function for a discrete
parent, the atom's probability for a mixed one, and zero for an ordinary
continuous parent.
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
needs it because the lower endpoint is included in the truncated
support, so the tail below the interval is \\F(L^-) = F(L) - P(Y = L)\\.

## Usage

``` r
parent_mass_at(distrib, x, theta)
```

## Arguments

- distrib:

  A truncated distribution object, of either class. The branch is taken
  on it, and the mass is asked of `distrib@parent_distrib`.

- x:

  A single number, the point to evaluate at.

- theta:

  A named list of the parent's parameters.

## Value

A numeric vector of probabilities, of length one unless a parameter in
`theta` varies by observation.

## Details

The question is asked of the OBJECT, through
[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md),
never of its class. Branching on `discrete_distrib` looks right and is
wrong: the distribution function of `zero_adjusted(gamma2_distrib())`
already contains the point mass at zero, so \\F(0) \ne F(0^-)\\ although
the object is a `continuous_distrib`. Truncating such a parent from
above, with the atom retained, would then drop exactly that mass out of
\\Z\\ while every formula still read correctly, and the density would
integrate to something other than one.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`distrib_atoms()`](https://statmodels7.github.io/distributions7/reference/distrib_atoms.md)
for the declaration it reads, and
[`trunc_constants()`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
for its one caller.

## Examples

``` r
# A discrete parent: the mass function.
ztp <- truncated(poisson_distrib(), lower = 1)
c(reported = distributions7:::parent_mass_at(ztp, 1, list(mu = 2)),
  dpois = dpois(1, 2))
#>  reported     dpois 
#> 0.2706706 0.2706706 

# A mixed parent: the atom, although the object is a continuous_distrib.
tz <- truncated(zero_adjusted(gamma2_distrib()), upper = 5)
th <- list(mu = 2, sigma2 = 1, za = 0.3)
distributions7:::parent_mass_at(tz, 0, th)
#> [1] 0.3

# An ordinary continuous parent: zero, at every point.
tn <- truncated(gaussian1_distrib(), lower = -1, upper = 2)
theta <- list(mu = 0.3, sigma = 1.2)
distributions7:::parent_mass_at(tn, 0, theta)
#> [1] 0
```
