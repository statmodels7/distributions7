# Check That an Argument Is an Expectation Method

Returns `method` invisibly when it is `NULL` (if `null_ok`) or inherits
from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
and signals an error otherwise.

## Usage

``` r
check_method(method, arg = "method", null_ok = TRUE)
```

## Arguments

- method:

  The value to check.

- arg:

  The argument name used in the message.

- null_ok:

  Whether `NULL` is accepted.

## Value

`method`, invisibly.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
