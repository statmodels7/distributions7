# Reject the Zero Skewness Where a Derivative Diverges

Signals an error when any skewness is exactly zero. The log-density's
series in \\r\\ has nonzero \\r^4\\ and \\r^5\\ terms, so every observed
derivative of order two or more in \\\gamma_1\\, and the derivatives of
the expected information in \\\gamma_1\\, grow like a negative power of
\\\gamma_1\\ and are infinite at zero.

## Usage

``` r
sn2_reject_zero(theta, what)
```

## Arguments

- theta:

  An aligned parameter list.

- what:

  A short description of the quantity, for the message.

## Value

`NULL`, invisibly, when no skewness is zero.
