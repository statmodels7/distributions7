# Validate a Multivariate Distribution

The battery
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
runs on a multivariate distribution, where the one-dimensional checks
have no counterpart.

## Usage

``` r
check_distrib_mv(distrib, theta, n, nsim, tol)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters, already aligned.

- n:

  The number of observations drawn for the derivative checks.

- nsim:

  The Monte Carlo sample size.

- tol:

  The relative tolerance.

## Value

A list of one-row data frames, as
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
builds.

## Details

Five of the univariate checks do not survive the move to \\p\\
dimensions: the distribution function is an integral over an orthant,
the quantile function inverts an ordering that does not exist, and the
two checks built on them go with them. Running them anyway and reporting
the rejections as failures is the mistake a validator makes when it does
not know about a case, and it is worse than not checking, because a user
validating their own distribution cannot tell a real defect from it.

What replaces them are checks that do generalize:

1.  **the density is positive and finite** on a sample;

2.  **the density integrates to one**. For a family that enumerates its
    support through
    [`mv_support`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
    this is an exact sum over that support; otherwise it is importance
    sampling from the proposal
    [`mv_reference_draw`](https://statmodels7.github.io/distributions7/reference/mv_reference_draw.md)
    supplies, which by default is a gaussian with the same mean and an
    inflated covariance. The proposal is deliberately not the
    distribution itself, which would make the ratio identically one and
    the check vacuous;

3.  **the score has mean zero**, the first Bartlett identity, under the
    distribution's own generator – so it is also a check that the
    generator and the density describe the same law;

4.  **gradient and Hessian** against finite differences of the summed
    log-density;

5.  **the expected information** against the Monte Carlo average of the
    observed Hessian, and against the variance of the score, which is
    the second Bartlett identity;

6.  **the generator** against the first two moments;

7.  **the response derivatives** against finite differences in \\y\\.

The last of these is emitted only when it applies, as the univariate
battery already omits the checks that a discrete family has no
counterpart for. A family with an enumerable support is discrete and has
no derivative in the response, and the multivariate base class rejects
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
by design, so a family that has not registered one has made a choice
rather than left a gap.

## See also

[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
