# Check Parameter Values Against Their Domains

Verifies that every supplied parameter value lies strictly inside the
distribution's `params_bounds`, and is finite. Domains are treated as
**open** intervals: for instance a Gaussian requires \\\sigma \> 0\\ and
a Bernoulli requires \\0 \< \mu \< 1\\, since the log-likelihood and its
derivatives are not defined at the boundary.

This is called automatically by every generic (through the internal
`align_theta()`), so passing an out-of-domain value raises an
informative error instead of silently producing `NaN`. It is exported so
that it can also be used directly, e.g. when writing an optimizer.

## Usage

``` r
check_theta_bounds(distrib, theta)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  A list of parameter values, ordered as `distrib@params`.

## Value

Invisibly `NULL`. Raises an error listing every offending parameter, the
offending value(s) and the expected domain.

## Examples

``` r
d <- gaussian_distrib()
check_theta_bounds(d, list(mu = 0, sigma = 1))
if (FALSE) { # \dontrun{
check_theta_bounds(d, list(mu = 0, sigma = -1)) # error: sigma outside (0, Inf)
} # }
```
