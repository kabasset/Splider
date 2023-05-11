/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Splider/Spline.h"

#include <boost/test/unit_test.hpp>
#include <complex>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Spline_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(index_test) {
  const Splider::SplineIntervals u {1, 2, 3, 4};
  for (std::size_t i = 0; i < u.size() - 1; ++i) {
    BOOST_TEST(u.index(u[i]) == i);
  }
  BOOST_TEST(u.index(u[u.size() - 1]) == u.size() - 2);
  for (std::size_t i = 0; i < u.size() - 1; ++i) {
    BOOST_TEST(u.index((u[i] + u[i + 1]) / 2) == i);
  }
  BOOST_CHECK_THROW(u.index(u[0] - 1), std::runtime_error);
  BOOST_CHECK_THROW(u.index(u[u.size() - 1] + 1), std::runtime_error);
}

template <typename X, typename Y>
struct LinSplineFixture {
  std::vector<X> u {1, 2, 3, 4};
  Splider::SplineIntervals domain = Splider::SplineIntervals(u);
  Splider::SplineBuilder builder = Splider::SplineBuilder(u);
  std::vector<X> x {1.1, 2.5, 3.9};
  std::vector<Y> v {10, 20, 30, 40};
  std::vector<Y> y {11, 25, 39};
};

using RealLinSplineFixture = LinSplineFixture<double, double>;
using ComplexLinSplineFixture = LinSplineFixture<double, std::complex<double>>;

BOOST_FIXTURE_TEST_CASE(real_spline_test, RealLinSplineFixture) {
  const Splider::Spline<double> spline(domain, v);
  std::vector<double> out;
  for (const auto& e : x) {
    out.push_back(spline(e));
  }
  BOOST_TEST(out.size() == x.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    BOOST_TEST(out[i] > v[i]);
    BOOST_TEST(out[i] < v[i + 1]);
  }
}

BOOST_AUTO_TEST_CASE(real_interpolant_test) {
  const Splider::SplineBuilder b {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  const std::vector<double> v {10, 20, 30, 40};
  const auto spline = b.interpolant(v);
  std::vector<double> y;
  for (const auto& e : x) {
    y.push_back(spline(e));
  }
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i] > v[i]);
    BOOST_TEST(y[i] < v[i + 1]);
  }
}

BOOST_AUTO_TEST_CASE(real_resampler_test) {
  const Splider::SplineBuilder b {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  const std::vector<double> v {10, 20, 30, 40};
  const auto resample = b.resampler(x);
  const auto y = resample(v);
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i] > v[i]);
    BOOST_TEST(y[i] < v[i + 1]);
  }
}

BOOST_AUTO_TEST_CASE(complex_interpolant_test) {
  const Splider::SplineBuilder b {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  const std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  const auto resample = b.resampler(x);
  const auto y = resample(v);
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i].real() > v[i].real());
    BOOST_TEST(y[i].real() < v[i + 1].real());
    BOOST_TEST(y[i].imag() < v[i].imag());
    BOOST_TEST(y[i].imag() > v[i + 1].imag());
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
