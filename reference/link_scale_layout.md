# The Index Layout of the Link-Scale Assembly

Which parameter-scale components enter each link-scale component of a
given order, with the multiplicities and the lookup keys, for one vector
of parameter names.

## Usage

``` r
link_scale_layout(params, order)
```

## Arguments

- params:

  The parameter names, in the family's own order.

- order:

  The derivative order.

## Value

A list with one entry per component of that order, each carrying `name`,
the multi-index's distinct entries `uniq`, their multiplicities `mult`,
and `combos`, one entry per term of the nested sum with its exponents,
its order and its lookup key.

## Details

The multi-index list, the multiplicities and the keys of the
parameter-scale components depend only on the parameter names and the
order, not on the values.
[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
reads them from this cache instead of rebuilding them on every call,
which matters because a score-driven filter reaches it once per
observation and iteration.

The combinations are enumerated by decoding a counter in mixed radix
rather than with `expand.grid`, and the enumeration is paid once per
family and order.

The cache is keyed by the parameter names and the order, so it holds one
entry per family per order actually used. It lives in the function's own
enclosure rather than in the namespace: it is an implementation detail
with no other reader, and a package-level object would need a help topic
of its own saying so.
