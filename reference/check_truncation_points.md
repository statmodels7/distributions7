# Validate the Truncation Endpoints

Checks that the interval is well formed and that truncating this parent
there leaves a model worth fitting.
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
calls it before building anything, so every rejection names the endpoint
and the support it was compared against.

## Usage

``` r
check_truncation_points(distrib, lower, upper, is_disc)
```

## Arguments

- distrib:

  The parent distribution, before wrapping. Its `bounds` and
  `distrib_name` are read.

- lower:

  The lower endpoint, or `NULL` for none.

- upper:

  The upper endpoint, or `NULL` for none.

- is_disc:

  A single logical saying whether the parent is discrete, which the
  caller has already determined. It governs the whole-number check
  alone.

## Value

Invisibly `NULL`. Every failure raises, with the endpoint and the
parent's support named.

## Details

Four conditions are tested, in order: each endpoint is a single
non-missing number; a finite endpoint on a discrete parent is a whole
number, a fractional one being ambiguous between the support points
either side; `lower` is strictly below `upper`; and neither endpoint
sits outside the parent's support, where truncation would either remove
nothing, so that the result is the parent itself, or remove everything.

Two further conditions are left to
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
which has the resolved endpoints: truncating zero away from a zero
wrapper, which cancels that wrapper's parameter out of the likelihood,
and a discrete interval leaving too few support points to identify the
parameters.

## Notation

\\L\\ and \\U\\ are the truncation endpoints, both included in the
support; \\Z(\theta) = P(L \le Y \le U)\\ is the retained mass; \\f\\
and \\F\\ are the parent's density and distribution function; \\s_i\\
and \\H\_{ij}\\ are the parent's score and observed Hessian; and
\\\mathbb{E}\_T\\ is expectation under the truncated law.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
its one caller.

## Examples

``` r
# A point at or below the support removes nothing: the gamma lives on
# (0, Inf), so truncating below zero would return the gamma itself.
try(distributions7:::check_truncation_points(
  gamma2_distrib(), lower = -2, upper = NULL, is_disc = FALSE))
#> Error : truncated() was given lower = -2, which is at or below the lower bound of the
#>   support of 'gamma2' (0). Truncating there removes no probability mass and the
#>   result would be the parent distribution. Choose a point strictly inside the
#>   support, or omit 'lower'.

# A fractional endpoint on a discrete parent is ambiguous.
try(distributions7:::check_truncation_points(
  poisson_distrib(), lower = 1.5, upper = NULL, is_disc = TRUE))
#> Error : 'lower' = 1.5 is not a point of the support. A discrete distribution is supported on
#>   the integers, so a non-integer truncation point is ambiguous: use 1 or 2.

# A well-formed interval returns nothing and raises nothing.
distributions7:::check_truncation_points(
  gaussian1_distrib(), lower = -1, upper = 2, is_disc = FALSE)
```
