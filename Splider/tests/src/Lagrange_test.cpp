/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Linx/Data/Sequence.h"
#include "Splider/Lagrange.h"

#include <boost/test/unit_test.hpp>
#include <complex>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Lagrange_test)

//-----------------------------------------------------------------------------

using Spline = Splider::Lagrange;

struct RealLinFixture {
  std::vector<double> u {1, 2, 3, 4};
  std::vector<double> x {1.1, 2.5, 3.9};
  std::vector<double> v {10, 20, 30, 40};
  std::vector<double> y {11, 25, 39};
};

struct ComplexLinFixture {
  std::vector<double> u {1, 2, 3, 4};
  std::vector<double> x {1.1, 2.5, 3.9};
  std::vector<std::complex<double>> v {{10, -1}, {20, -2}, {30, -3}, {40, -4}};
  std::vector<std::complex<double>> y {{11, -1.1}, {25, -2.5}, {39, -3.9}};
};

struct RealRandomFixture {
  std::vector<double> u {1, 2, 3, 4};
  std::vector<double> x {1.1, 2.5, 3.9};
  Linx::Sequence<double> v = Linx::Sequence<double>(u.size()).generate(Linx::UniformNoise<double>(0, 1));
};

BOOST_FIXTURE_TEST_CASE(real_lin_spline_test, RealLinFixture)
{
  const auto build = Spline::builder(u);
  auto spline = build.spline(v);
  auto out = spline(x);
  BOOST_TEST(out.size() == x.size());
  BOOST_TEST(out == y, boost::test_tools::tolerance(1.e-6) << boost::test_tools::per_element());
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
