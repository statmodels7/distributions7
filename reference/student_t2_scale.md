# The Scale of a Student t with a Given Standard Deviation

Returns \\s_0 = \sigma\sqrt{1 - 2/\nu}\\, the scale of the Student t
whose standard deviation is \\\sigma\\.

## Usage

``` r
student_t2_scale(theta)
```

## Arguments

- theta:

  A list with the location, the standard deviation and the degrees of
  freedom, in that order.

## Value

A numeric vector of scales.

## See also

[`distrib_pdf.StudentT2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.StudentT2Distrib.md)
