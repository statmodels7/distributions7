# Check Parameter Domains, Taking the Properties as Arguments

Provides the body of
[`check_theta_bounds()`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md).
It differs only in that it receives the distribution's properties
instead of reaching for them.

## Usage

``` r
check_bounds_fast(params, param_bounds, theta, distrib_name)
```

## Arguments

- params:

  A character vector of parameter names.

- param_bounds:

  A named list of length-2 domain vectors.

- theta:

  A list of parameter values, ordered as `params`.

- distrib_name:

  The distribution's name, used in the message.

## Value

Invisibly `NULL`; raises an error if any value is outside its open
domain or not finite.

## Details

[`align_theta()`](https://statmodels7.github.io/distributions7/reference/align_theta.md)
calls this on every generic invocation, so its fixed cost is the
package's fixed cost. Each read of an S7 property has a cost of its own,
so the caller reads the properties that it needs once and passes them
on.

Every offending parameter is collected before anything is raised, so a
call with two bad values reports both rather than one at a time. At most
three distinct offending values are shown.

## See also

[`check_theta_bounds()`](https://statmodels7.github.io/distributions7/reference/check_theta_bounds.md)
