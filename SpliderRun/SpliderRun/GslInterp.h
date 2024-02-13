/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_GSLINTERP_H
#define _SPLIDER_GSLINTERP_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

/**
 * @brief 1D resampling with GSL.
 */
template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u, const V& v, const X& x)
{
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, u.size());
  std::vector<double> y;
  for (const auto& row : sections(v)) {
    y.clear();
    gsl_spline_init(spline, u.data(), row.data(), u.size());
    for (const auto& e : x) {
      y.push_back(gsl_spline_eval(spline, e, acc));
    }
    gsl_interp_accel_reset(acc); // FIXME needed?
  }
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline);
  return y;
}

/**
 * @brief 2D resampling with GSL.
 */
template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u0, const U& u1, const V& v, const X& x)
{
  gsl_interp_accel* xacc = gsl_interp_accel_alloc();
  gsl_interp_accel* yacc = gsl_interp_accel_alloc();
  gsl_spline2d* spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, u0.size(), u1.size());
  std::vector<double> y;
  for (const auto& plane : sections(v)) {
    y.clear();
    gsl_spline2d_init(spline, u0.data(), u1.data(), plane.data(), u0.size(), u1.size());
    for (const auto& e : x) {
      y.push_back(gsl_spline2d_eval(spline, e[0], e[1], xacc, yacc));
    }
    gsl_interp_accel_reset(xacc); // FIXME needed?
    gsl_interp_accel_reset(yacc); // FIXME needed?
  }
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  gsl_spline2d_free(spline);
  return y;
}

/**
 * @brief 2D evaluation with GSL.
 */
template <typename U, typename V, typename X>
std::vector<double> eval_with_gsl(const U& u0, const U& u1, const V& v, const X& x)
{
  gsl_interp_accel* xacc = gsl_interp_accel_alloc();
  gsl_interp_accel* yacc = gsl_interp_accel_alloc();
  gsl_spline2d* spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, u0.size(), u1.size());
  std::vector<double> y;
  gsl_spline2d_init(spline, u0.data(), u1.data(), v.data(), u0.size(), u1.size());
  for (const auto& xi : x) {
    y.push_back(gsl_spline2d_eval(spline, xi[0], xi[1], xacc, yacc));
  }
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  gsl_spline2d_free(spline);
  return y;
}

#endif
