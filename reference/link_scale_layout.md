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

[`to_link_scale()`](https://statmodels7.github.io/distributions7/reference/to_link_scale.md)
used to rebuild this on every call: the multi-index list, a `unique` and
a `tabulate` per component, and a `sort` and a `paste` per combination
to spell the key of the parameter-scale component to look up. None of it
depends on the values, only on the names and the order, and a profile of
a fitted score-driven model put `paste`, `sort` and `unique` among the
leaders of its self time – a filter reaches this once per observation
per iteration, so the names were being respelled millions of times per
fit.

The combinations are enumerated by decoding a counter in mixed radix
rather than with `expand.grid`, which is the same device the loop used
before and is now paid once.

The cache is keyed by the parameter names and the order, so it holds one
entry per family per order actually used. It lives in the function's own
enclosure rather than in the namespace: it is an implementation detail
with no other reader, and a package-level object would need a help topic
of its own saying so.
