/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Linx/Data/Raster.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Data/Tiling.h"
#include "Linx/Run/Chronometer.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/C2.h"
#include "Splider/Hermite.h"
#include "Splider/Lagrange.h"
#include "SpliderRun/GslInterp.h"

#include <iostream>

using Duration = std::chrono::milliseconds;

template <typename TSpline, typename U, typename V, typename X, typename Y>
void eval(const U& u, const V& v, const X& x, Y& y)
{
  const auto build = TSpline::Multi::builder(u, u);
  auto cospline = build.cospline(x); // FIXME use Radius = -1 for C2
  for (const auto& plane : sections(v)) {
    y = cospline(plane);
  }
}

template <typename TDuration, typename U, typename V, typename X, typename Y>
TDuration resample(const U& u, const V& v, const X& x, Y& y, const std::string& setup)
{
  Linx::Chronometer<TDuration> chrono;
  chrono.start();
  if (setup == "c2") {
    eval<Splider::C2>(u, v, x, y);
  } else if (setup == "hermite") {
    eval<Splider::Hermite::FiniteDiff>(u, v, x, y);
  } else if (setup == "lagrange") {
    eval<Splider::Lagrange>(u, v, x, y);
  } else if (setup == "gsl") {
    y = resample_with_gsl(u, u, v, x);
  } else {
    throw std::runtime_error("Case not implemented");
  }
  return chrono.stop();
}

int main(int argc, const char* const argv[])
{
  Linx::ProgramOptions options("2D cospline benchmark.");
  options.named("case", "Test case: c2, hermite, lagrange, gsl", std::string("c2"));
  options.named("knots", "Number of knots along each axis", 100L);
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
  auto v = Linx::Raster<double, 3>({u_size, u_size, v_iters});
  v.generate(
      [&](const auto& p) {
        auto s = std::sin(u[p[0]]);
        auto c = std::cos(u[p[1]]);
        return s * c * (p[2] + 1);
      },
      v.domain());
  std::cout << "  u: " << u << std::endl;
  std::cout << "  v: " << v << std::endl;

  std::cout << "\nGenerating trajectory...\n" << std::endl;

  Splider::Trajectory<2> x(x_size);
  auto si = seed;
  for (auto& xi : x) {
    if (seed != -1) {
      ++si;
    }
    xi.generate(Linx::UniformNoise<double>(u[1], u[u_size - 2], si));
  }
  std::vector<double> y;
  std::cout << "  x: " << x << std::endl;

  std::cout << "\nInterpolating...\n" << std::endl;
  const auto duration = resample<Duration>(u, v, x, y, setup);
  std::cout << "  y: " << Linx::Sequence<double>(y) << std::endl;

  // logger.debug("i\t\tx0\tx1\tf(x0, x1)\ty");
  // for (Linx::Index i = 0; i < x_size; ++i) {
  //   const auto x0 = x[i][0];
  //   const auto x1 = x[i][1];
  //   logger.debug() << i << '\t' << x0 << '\t' << x1 << '\t' << std::sin(x0) * std::cos(x1) * v_iters << '\t' << y[i];
  // }

  std::cout << "  Done in " << duration.count() << "ms" << std::endl;

  std::cout << std::endl;
}
