# Splider

## Purpose

Cubic splines are parametric functions, whose parameters are the domain intervals and knot values.
Spline libraries generally allow instantiating a spline and calling it on an argument, e.g.:

```cpp
std::vector<double> u = ... // Knot abscissae
std::vector<double> v = ... // Knot values
Spline spline(u, v);
double x = ...
double y = spline(x);
```

This is generally sufficient, but in some cases yields huge recomputation where caching could be used.
For example, assume we need to call `spline` on the same `x` for many variations of `v`.
Worse, `x` can be a vector instead of a mere scalar.
In these cases, many spline coefficients remain valid while changing `v`.

## Approach

Splider separates the spline classes in components (intervals, knots and arguments), which each hold their cache.
For example `SplineIntervals` holds not only the knot positions but also associated byproducts like the spacing between positions.
This allows recomputing only what has changed.

For resampling several functions using spline interpolation, one can simply do:

```cpp
SplineIntervals u = ... // Compute domain-related coefficients
SplineArguments x = ... // Compute arguments-related coefficients
std::vector<std::vector<double>> vs = ...
for (const auto& v : vs) {
  auto y = Spline(u, v)(x); // Compute only what is left
  ...
}
```

Moreover, Splider is compatible with any value type of a ring (i.e. with `+` and `*` operators), e.g. `std::complex`.

## Status

Implementation is ongoing...
