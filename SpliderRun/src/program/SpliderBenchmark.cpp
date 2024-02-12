/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Linx/Data/Raster.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Data/Tiling.h"
#include "Linx/Run/Chronometer.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/C2.h"
#include "Splider/Cospline.h"
#include "Splider/Hermite.h"
#include "Splider/Lagrange.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <iostream>

using Duration = std::chrono::milliseconds;

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

template <typename TDuration, typename U, typename V, typename X, typename Y>
TDuration resample(const U& u, const V& v, const X& x, Y& y, const std::string& setup)
{
  Linx::Chronometer<TDuration> chrono;
  chrono.start();
  if (setup == "d") {
    using Domain = Splider::Partition<double>;
    const Domain domain(u);
    Splider::Cospline<double, Domain> cospline(domain, x);
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "l") {
    using Domain = Splider::Linspace<double>;
    const Domain domain(u[0], u[1], u.ssize());
    const Splider::Args<double> args(domain, x);
    Splider::Spline<double, Domain> spline(domain);
    for (const auto& row : sections(v)) {
      spline.assign(row);
      y = spline(args);
    }
  } else if (setup == "c2") {
    using Spline = Splider::C2;
    const auto build = Spline::builder(u);
    auto cospline = build.cospline(x);
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "c2fd") {
    using Spline = Splider::C2::FiniteDiff;
    const auto build = Spline::builder(u);
    auto cospline = build.cospline(x);
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "h") {
    using Spline = Splider::Hermite::FiniteDiff;
    const auto build = Spline::builder(u);
    auto cospline = build.cospline(x);
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "lagrange") {
    using Spline = Splider::Lagrange;
    const auto build = Spline::builder(u);
    auto cospline = build.cospline(x);
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "g") {
    y = resample_with_gsl(u, v, x);
  } else {
    throw std::runtime_error("Case not implemented");
  }
  return chrono.stop();
}

int main(int argc, const char* const argv[])
{
  Linx::ProgramOptions options("1D cospline benchmark.");
  options.named("case", "Test case: d (double), f (float), s (Spline), g (GSL)", std::string("d"));
  options.named("knots", "Number of knots", 100L);
  options.named("args", "Number of arguments", 100L);
  options.named("iters", "Numper of iterations", 1L);
  options.named("seed", "Random seed", -1L);
  options.parse(argc, argv);
  const auto setup = options.as<std::string>("case");
  const auto u_size = options.as<Linx::Index>("knots");
  const auto x_size = options.as<Linx::Index>("args");
  const auto v_iters = options.as<Linx::Index>("iters");
  const auto seed = options.as<Linx::Index>("seed");

  std::cout << "\nGenerating knots...\n" << std::endl;

  const auto u = Linx::Sequence<double>(u_size).linspace(0, Linx::pi<double>() * 4);
  auto v = Linx::Raster<double, 2>({u_size, v_iters});
  v.generate(
      [&](const auto& p) {
        auto s = std::sin(u[p[0]]);
        return s * (p[1] + 1);
      },
      v.domain());
  std::cout << "  u: " << u << std::endl;
  std::cout << "  v: " << v << std::endl;

  std::cout << "\nGenerating trajectory...\n" << std::endl;

  auto x = Linx::Sequence<double>(x_size).generate(Linx::UniformNoise<double>(u[1], u[u_size - 2], seed));
  std::vector<double> y;
  std::cout << "  x: " << x << std::endl;

  std::cout << "\nInterpolating...\n" << std::endl;

  const auto duration = resample<Duration>(u, v, x, y, setup);
  std::cout << "  y: " << Linx::Sequence<double>(y) << std::endl;

  std::cout << "  Done in " << duration.count() << "ms" << std::endl;

  std::cout << std::endl;
}
