# Project Overview

<br/>![Splider logo](doc/diagrams/logo_square.svg)

## Purpose

Cubic splines are parametric functions, whose parameters are the domain intervals and knot values.
Spline libraries generally allow instantiating a spline and calling it on an argument, e.g.:

```cpp
std::vector<double> u = ... // Knot abscissae
std::vector<double> v = ... // Knot values
std::vector<double> x = ... // Spline arguments
Spline spline(u, v);
for (const auto& xi : x) {
  double yi = spline(xi);
  ...
}
```

This is generally adequate, but in some cases yields huge recomputation where caching could be used.
For example, assume we need to call `spline` on the same `x` for many variations of `v`.
In this case, many spline coefficients remain valid while changing `v`.

## Approach

Splider separates the spline classes in components (subintervals, knots and arguments), which each hold their cache.
For example the subintervals hold not only their bounds but also associated byproducts like their lengths.
This allows recomputing only what has changed.

With Splider, the above example writes:

```cpp
const auto build = Spline::builder(u); // Compute domain-related coefficients
auto spline = build.spline(v); // Compute knot-related coefficients
auto y = spline(x); // Compute argument-related coefficients
```

In general splines are initialized with `u` and `v` and applied to `x`.
Yet, sometimes, one wish to initialize them with `u` and `x` and apply them to `v`.
This is known as a "cospline" in Splider:

```cpp
const auto build = Spline::builder(u); // Compute domain-related coefficients
auto cospline = build.cospline(x); // Compute arguments-related coefficients
auto y = cospline(v); // Compute knots-related coefficients
```

This is especially efficient when several functions must be interpolated with the same `x`:

```cpp
const auto build = Spline::builder(u); // Compute domain-related coefficients
auto cospline = build.cospline(x);
auto y0 = cospline(v0);
auto y1 = cospline(v1);
auto y2 = cospline(v2);
```

Moreover, Splider is compatible with any value type of a ring (i.e. with `+` and `*` operators), e.g. `std::complex`.

Extrapolation to N-dimensional splines is under development.

Splider relies on [Linx](https://github.com/kabasset/Linx) for the data structures and basic operations.

## Status

Prototype is being validated!

