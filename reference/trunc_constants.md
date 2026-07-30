# The Truncation Constant and Lower Tail

Returns \\F(\ell^-)\\ and the normalising constant \\Z = F(u) -
F(\ell^-)\\, vectorised in \\\theta\\.

## Usage

``` r
trunc_constants(distrib, theta)
```

## Arguments

- distrib:

  A truncated distribution object.

- theta:

  A named list of parameters.

## Value

A list with `Fl` and `Z`.

## Details

Both endpoints are **included** in the truncated support, so any mass
sitting exactly on the lower one has to be added back: \\F(\ell^-) =
F(\ell) - P(Y = \ell)\\. That correction is the *atom* case, not the
discrete case – see
[`parent_mass_at`](https://statmodels7.github.io/distributions7/reference/parent_mass_at.md).

An interval carrying no probability under the given parameters is
reported rather than returned, since the truncated law is not defined
there.
