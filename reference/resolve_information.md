# Resolve the Information Argument

Returns `information` when it is an object and
[`default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.md)
of the family when it is `NULL`.

## Usage

``` r
resolve_information(distrib, information)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- information:

  `NULL`, an
  [`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
  or an
  [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
  object.

## Value

An
[Information](https://statmodels7.github.io/distributions7/reference/Information.md)
object.

## See also

[information_objects](https://statmodels7.github.io/distributions7/reference/information_objects.md),
[`default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.md)
