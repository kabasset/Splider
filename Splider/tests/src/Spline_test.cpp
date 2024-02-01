/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Linx/Data/Sequence.h"
#include "Splider/Spline.h"

#include <boost/test/unit_test.hpp>
#include <complex>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Spline_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(index_test)
{
  const Splider::Partition<> u {1, 2, 3, 4};
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

template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u, const V& v, const X& x)
{
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, u.size());
  std::vector<double> y;
  gsl_spline_init(spline, u.data(), v.data(), u.size());
  for (const auto& e : x) {
    y.push_back(gsl_spline_eval(spline, e, acc));
  }
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline);
  return y;
}

struct RealLinFixture {
  std::vector<double> u {1, 2, 3, 4};
  Splider::Partition<> domain = Splider::Partition<>(u);
  std::vector<double> x {1.1, 2.5, 3.9};
  std::vector<double> v {10, 20, 30, 40};
  std::vector<double> y {11, 25, 39};
};

struct ComplexLinFixture {
  std::vector<double> u {1, 2, 3, 4};
  Splider::Partition<> domain = Splider::Partition<>(u);
  std::vector<double> x {1.1, 2.5, 3.9};
  std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  std::vector<std::complex<double>> y {{11, -1.1}, {25, -2.5}, {39, -3.9}};
};

struct RealRandomFixture {
  std::vector<double> u {1, 2, 3, 4};
  Splider::Partition<> domain = Splider::Partition<>(u);
  std::vector<double> x {1.1, 2.5, 3.9};
  Linx::Sequence<double> v = Linx::Sequence<double>(u.size()).generate(Linx::UniformNoise<double>(0, 1));
};

BOOST_FIXTURE_TEST_CASE(real_lin_spline_test, RealLinFixture)
{
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
  auto expected = resample_with_gsl(u, v, x);
  BOOST_TEST(out == expected, boost::test_tools::tolerance(1.e-6) << boost::test_tools::per_element());
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
