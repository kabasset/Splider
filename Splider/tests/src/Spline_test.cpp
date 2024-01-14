/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Splider/Spline.h"

#include <boost/test/unit_test.hpp>
#include <complex>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Spline_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(index_test) {
  const Splider::Partition u {1, 2, 3, 4};
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

struct RealLinSplineFixture {
  std::vector<double> u {1, 2, 3, 4};
  Splider::Partition domain = Splider::Partition(u);
  std::vector<double> x {1.1, 2.5, 3.9};
  std::vector<double> v {10, 20, 30, 40};
  std::vector<double> y {11, 25, 39};
};

struct ComplexLinSplineFixture {
  std::vector<double> u {1, 2, 3, 4};
  Splider::Partition domain = Splider::Partition(u);
  std::vector<double> x {1.1, 2.5, 3.9};
  std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  std::vector<std::complex<double>> y {{11, -1.1}, {25, -2.5}, {39, -3.9}};
};

BOOST_FIXTURE_TEST_CASE(real_spline_test, RealLinSplineFixture) {
  Splider::Spline<double> spline(domain, v);
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

BOOST_FIXTURE_TEST_CASE(real_lazy_spline_test, RealLinSplineFixture) {
  Splider::Spline<double, Splider::SplineCache::Lazy> spline(domain, v);
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

BOOST_FIXTURE_TEST_CASE(real_interpolant_test, RealLinSplineFixture) {
  Splider::Spline<double> spline(domain, v);
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

BOOST_FIXTURE_TEST_CASE(real_resampler_test, RealLinSplineFixture) {
  Splider::Cospline resampler(domain, x);
  const auto out = resampler(v);
  BOOST_TEST(out.size() == x.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    BOOST_TEST(out[i] > v[i]);
    BOOST_TEST(out[i] < v[i + 1]);
  }
}

BOOST_FIXTURE_TEST_CASE(complex_interpolant_test, ComplexLinSplineFixture) {
  Splider::Cospline resampler(domain, x);
  const auto out = resampler(v);
  BOOST_TEST(out.size() == x.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    BOOST_TEST(out[i].real() > v[i].real());
    BOOST_TEST(out[i].real() < v[i + 1].real());
    BOOST_TEST(out[i].imag() < v[i].imag());
    BOOST_TEST(out[i].imag() > v[i + 1].imag());
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
