/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Linx/Data/Sequence.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/C2.h"
#include "Splider/CatmullRom.h"
#include "Splider/Hermite.h"
#include "Splider/Lagrange.h"
#include "SpliderRun/GslInterp.h"

#include <iostream>

template <typename TSpline, typename U, typename V, typename X>
std::vector<double> eval(const U& u, const V& v, const X& x)
{
  const auto build = TSpline::Multi::builder(u, u);
  auto cospline = build.cospline(x);
  return cospline(v);
}

int main(int argc, const char* const argv[])
{
  Linx::ProgramOptions options;
  options.positional("knots", "Number of knots", 7L);
  options.positional("args", "Number of arguments", 101L);
  options.named("seed", "Random seed", -1L);
  options.parse(argc, argv);
  const auto u_size = options.as<Linx::Index>("knots");
  const auto x_size = options.as<Linx::Index>("args");
  const auto seed = options.as<Linx::Index>("seed");
  // FIXME std::cout << options << std::endl;

  std::cout << "\nGenerating knots...\n" << std::endl;

  const auto u = Linx::Sequence<double>(u_size).linspace(0, Linx::pi<double>() * 4);
  auto v = Linx::Raster<double>({u_size, u_size});
  v.generate(
      [&](const auto& p) {
        auto s = std::sin(u[p[0]]);
        auto c = std::cos(u[p[1]]);
        return s * c;
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

  std::cout << "\nInterpolating...\n" << std::endl;

  std::vector<std::vector<double>> y;
  std::cout << "i\tx0\tx1\tgt";
  std::cout << "\tGSL";
  y.push_back(eval_with_gsl(u, u, v, x));
  std::cout << "\tC2 FD";
  y.push_back(eval<Splider::C2::FiniteDiff>(u, v, x));
  std::cout << "\tHermite FD";
  y.push_back(eval<Splider::Hermite::FiniteDiff>(u, v, x));
  std::cout << "\tCatmull-Rom";
  y.push_back(eval<Splider::Hermite::CatmullRom::Uniform>(u, v, x));
  std::cout << "\tLagrange";
  y.push_back(eval<Splider::Lagrange>(u, v, x));
  std::cout << std::endl;

  for (std::size_t i = 0; i < x.size(); ++i) {
    std::cout << i << '\t' << x[i][0] << '\t' << x[i][1];
    std::cout << '\t' << std::sin(x[i][0]) * std::cos(x[i][1]);
    for (std::size_t s = 0; s < y.size(); ++s) {
      std::cout << '\t' << y[s][i];
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;
}
