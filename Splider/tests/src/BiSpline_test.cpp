/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Splider/BiSpline.h"

#include <boost/test/unit_test.hpp>
#include <complex>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(BiSpline_test)

//-----------------------------------------------------------------------------

template <typename U, typename V, typename X>
std::vector<double> resampleWithGsl(const U& u0, const U& u1, const V& v, const X& x) {
  gsl_interp_accel* xacc = gsl_interp_accel_alloc();
  gsl_interp_accel* yacc = gsl_interp_accel_alloc();
  gsl_spline2d* spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, u0.size(), u1.size());
  gsl_spline2d_init(spline, u0.data(), u1.data(), v.data(), u0.size(), u1.size());
  std::vector<double> y;
  for (const auto& e : x) {
    y.push_back(gsl_spline2d_eval(spline, e[0], e[1], xacc, yacc));
  }
  return y;
}

template <Splider::SplineCache Cache>
struct RealLinExpSplineFixture {
  std::vector<double> u0 {1, 2, 3, 4};
  std::vector<double> u1 {1, 10, 100, 1000};
  std::vector<Linx::Vector<double, 2>> x {{1.1, 2.}, {2.5, 10.}, {2.5, 20.}, {2.5, 50.}, {3.9, 50.}};
  Linx::Raster<double> v = Linx::Raster<double>(
      {u0.size(), u1.size()},
      {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400, 1000, 2000, 3000, 4000});
  Splider::SplineIntervals domain0 = Splider::SplineIntervals(u0);
  Splider::SplineIntervals domain1 = Splider::SplineIntervals(u1);
  Splider::BiSplineResampler<double, Cache> resampler = Splider::BiSplineResampler<double, Cache>(domain0, domain1, x);
};

using EarlyRealLinExpSplineFixture = RealLinExpSplineFixture<Splider::SplineCache::Early>;
using LazyRealLinExpSplineFixture = RealLinExpSplineFixture<Splider::SplineCache::Lazy>;

BOOST_FIXTURE_TEST_CASE(real_resampler_test, EarlyRealLinExpSplineFixture) {
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

BOOST_FIXTURE_TEST_CASE(real_resampler_vs_gsl_test, EarlyRealLinExpSplineFixture) {
  const auto out = resampler(v);
  const auto gsl = resampleWithGsl(u0, u1, v, x);
  BOOST_TEST(out.size() == gsl.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    BOOST_TEST(out[i] == gsl[i]);
  }
}

BOOST_FIXTURE_TEST_CASE(real_lazy_resampler_vs_gsl_test, LazyRealLinExpSplineFixture) {
  const auto out = resampler(v);
  const auto gsl = resampleWithGsl(u0, u1, v, x);
  BOOST_TEST(out.size() == gsl.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    BOOST_TEST(out[i] == gsl[i]);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
