/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Splider/BiSpline.h"

#include <boost/test/unit_test.hpp>
#include <complex>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(BiSpline_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(real_resampler_test) {
  const Splider::SplineIntervals u0 {1, 2, 3, 4};
  const Splider::SplineIntervals u1 {1, 10, 100, 1000};
  const std::vector<Linx::Vector<double, 2>> x {{1.1, 2.}, {2.5, 20.}, {3.9, 50.}};
  const Linx::Raster<double> v(
      {u0.size(), u1.size()},
      {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000});
  Splider::BiSplineResampler<double> resampler(u0, u1, x.begin(), x.end());
  const auto y = resampler(v);
  BOOST_TEST(y.size() == x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    std::cout << x[i] << " - " << y[i] << "\n";
    Linx::Position<2> p;
    for (std::size_t j = 1; j < u0.size(); ++j) {
      if (u0[j] < x[i][0]) {
        p[0] = j;
      }
    }
    for (std::size_t j = 1; j < u1.size(); ++j) {
      if (u1[j] < x[i][1]) {
        p[1] = j;
      }
    }
    BOOST_TEST(y[i] > v[p]);
    BOOST_TEST(y[i] < v[p + 1]);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
