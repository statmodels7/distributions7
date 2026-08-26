# Differencing the Expected Information Once

The default route behind
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md):
a central difference in each parameter of the family's own expected
information.

## Usage

``` r
numerical_dexpected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  approx = "bartlett",
  nsim = 10000,
  h_rel = numericals7::fd_step(1, 1L)
)
```

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- scale:

  Either `"parameter"` or `"link"`.

- approx, nsim:

  Passed to
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md).

- h_rel:

  The relative step, a cube root of machine epsilon by default, the
  value at which a central difference balances truncation against
  rounding.

## Value

A named list keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md).

## Details

The step is
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)',
which shrinks near a finite boundary so that both evaluation points stay
strictly inside the parameter's open domain. On the link scale the
domain is the whole line and no clamp is needed, so the step is the
plain relative one.

**It refuses rather than approximating an approximation.** Where the
expected information is obtained by quadrature or by simulation, this
would be a difference of a difference, which the package forbids
everywhere, and it would cost 2p of the dearest call the family has.

## See also

[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md),
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
