/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Linx/Data/Sequence.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/C2.h"
#include "Splider/CatmullRom.h"
#include "Splider/Hermite.h"
#include "Splider/Lagrange.h"

#include <iostream>

int main(int argc, const char* const argv[])
{
  Linx::ProgramOptions options;
  options.positional("knots", "Number of knots", 7L);
  options.positional("args", "Number of arguments", 101L);
  options.parse(argc, argv);
  const auto knot_count = options.as<Linx::Index>("knots");
  const auto arg_count = options.as<Linx::Index>("args");
  // FIXME std::cout << options << std::endl;

  std::cout << "\nGenerating knots...\n" << std::endl;

  const auto u = Linx::Sequence<double>(knot_count).linspace(0, Linx::pi<double>() * 2.5);
  const auto x = Linx::Sequence<double>(arg_count).linspace(0, max(u));
  const auto v = sin(u);
  const auto gt = sin(x);

  std::cout << "i\tu\tv" << std::endl;
  for (std::size_t i = 0; i < u.size(); ++i) {
    std::cout << i << '\t' << u[i] << '\t' << v[i] << std::endl;
  }

  std::cout << "\nInterpolating...\n" << std::endl;

  std::vector<std::vector<double>> y;
  std::cout << "i\tx\tsin(x)";
  std::cout << "\tC2 Natural";
  y.push_back(Splider::C2::eval(u, v, x));
  std::cout << "\tC2 FD";
  y.push_back(Splider::C2::FiniteDiff::eval(u, v, x));
  std::cout << "\tHermite FD";
  y.push_back(Splider::Hermite::FiniteDiff::eval(u, v, x));
  std::cout << "\tCatmull-Rom";
  y.push_back(Splider::Hermite::CatmullRom::Uniform::eval(u, v, x));
  std::cout << "\tLagrange";
  y.push_back(Splider::Lagrange::eval(u, v, x));
  std::cout << std::endl;

  for (std::size_t i = 0; i < x.size(); ++i) {
    std::cout << i << '\t' << x[i] << '\t' << gt[i];
    for (std::size_t s = 0; s < y.size(); ++s) {
      std::cout << '\t' << y[s][i];
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;
}
