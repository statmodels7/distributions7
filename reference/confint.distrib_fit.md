# Confidence Intervals for a Maximum-Likelihood Fit

Returns Wald intervals. They are built symmetrically on the link scale,
where the parameters are unconstrained, and mapped through \\g^{-1}\\
when the parameter scale is requested, so that a limit can never leave
the parameter's domain. The two ends are sorted after mapping, because a
link need not be increasing.

## Arguments

- object:

  A
  [`distrib_fit`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- parm:

  Parameters to report, given by name or position. Defaults to all.

- level:

  Confidence level. Defaults to the level the fit was computed at, and
  any other value is obtained from the stored estimates and standard
  errors without refitting.

- scale:

  Either `"parameter"` (default) or `"link"`.

- ...:

  Unused.

## Value

A two-column matrix of confidence limits, one row per parameter.
