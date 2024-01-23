/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "SpliderRun/Demo.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(SpliderDemo_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(default_interpolant_test) {

  //! [Default interpolant]

  Splider::Partition u {1, 2, 3, 4};
  std::vector<double> v {10, 20, 30, 40};
  Splider::Spline<double> spline(u, v);

  std::vector<double> x {1.1, 2.5, 3.9};
  auto y = spline(x);

  //! [Default interpolant]
}

BOOST_AUTO_TEST_CASE(default_resampler_test) {

  //! [Default resampler]

  Splider::Partition u {1, 2, 3, 4};
  std::vector<double> x {1.1, 2.5, 3.9};
  Splider::Cospline resampler(u, x);

  std::vector<double> v {10, 20, 30, 40};
  auto y = resampler(v);

  //! [Default resampler]
}

BOOST_AUTO_TEST_CASE(default_bivariate_resampler_test) {

  //! [Default bivariate resampler]

  Splider::Partition u0 {1, 2, 3, 4};
  Splider::Partition u1 {1, 10, 100};
  Splider::Trajectory<2> x {{1.1, 2.}, {2.5, 10.}, {2.5, 20.}, {2.5, 50.}, {3.9, 50.}};
  Splider::BiCospline<double> resampler(u0, u1, x);

  Linx::Raster<double> v({u0.size(), u1.size()}, {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400});
  auto y = resampler(v);

  //! [Default bivariate resampler]
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
