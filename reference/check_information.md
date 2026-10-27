# Check That an Argument Is an Information Object

Returns `x` invisibly when it is `NULL`, an
[`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
or an
[`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
object and signals an error otherwise.

## Usage

``` r
check_information(x, arg = "information")
```

## Arguments

- x:

  The value to check.

- arg:

  The argument name used in the message.

## Value

`x`, invisibly.

## See also

[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md)
