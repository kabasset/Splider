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

BOOST_AUTO_TEST_CASE(scalar_real_test) {
  const Splider::SplineIntervals u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  const std::vector<double> v {10, 20, 30, 40};
  const Splider::Spline<double> spline(u, v);
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

BOOST_AUTO_TEST_CASE(vector_real_test) {
  const Splider::SplineIntervals u {1, 2, 3, 4};
  const auto x = makeSplineArgs(u, {1.1, 2.5, 3.9});
  const std::vector<double> v {10, 20, 30, 40};
  const Splider::Spline<double> spline(u, v);
  const auto y = spline(x);
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    BOOST_TEST(y[i] > v[i]);
    BOOST_TEST(y[i] < v[i + 1]);
  }
}

BOOST_AUTO_TEST_CASE(scalar_complex_test) {
  const Splider::SplineIntervals u {1, 2, 3, 4};
  const std::vector<double> x {1.1, 2.5, 3.9};
  const std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  const Splider::Spline<std::complex<double>> spline(u, v);
  std::vector<std::complex<double>> y;
  for (const auto& e : x) {
    y.push_back(spline(e));
  }
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
