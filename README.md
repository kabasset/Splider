# Project Overview

![Splider logo](doc/diagrams/logo_square.svg)

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

This is generally sufficient, but in some cases yields huge recomputation where caching could be used.
For example, assume we need to call `spline` on the same `x` for many variations of `v`.
In this case, many spline coefficients remain valid while changing `v`.

## Approach

Splider separates the spline classes in components (intervals, knots and arguments), which each hold their cache.
For example `Splider::SplineIntervals` holds not only the knot positions but also associated byproducts like the spacing between positions.
This allows recomputing only what has changed.

With Splider, the above example writes:

```cpp
Splider::SplineIntervals domain(u); // Compute domain-related coefficients
Splider::Spline<double> spline(domain, v); // Compute knot-related coefficients
std::vector<double> y = spline(x); // Compute argument-related coefficients
```

For resampling a function using spline interpolation, `x` is provided to Splider prior to `v`:

```cpp
Splider::SplineIntervals domain(u); // Compute domain-related coefficients
Splider::SplineResampler resampler(domain, x); // Compute arguments-related coefficients
std::vector<double> y = resampler(v); // Compute knots-related coefficients
```

This is especially efficient when several functions must be resampled:

```cpp
Splider::SplineIntervals domain(u);
Splider::SplineResampler resampler(domain, x);
std::vector<double> y0 = resampler(v0);
std::vector<double> y1 = resampler(v1);
std::vector<double> y2 = resampler(v2);
```

Moreover, Splider is compatible with any value type of a ring (i.e. with `+` and `*` operators), e.g. `std::complex`.

2D interpolation is also provided as `Splider::BiSplineResampler`.
It relies on [Linx](https://github.com/kabasset/Linx) for the data structures.

## Status

Prototype is being validated!

