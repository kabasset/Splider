/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/C2.h"
#include "Splider/CatmullRom.h"
#include "Splider/Hermite.h"
#include "Splider/Lagrange.h"

using Duration = std::chrono::milliseconds;

class SpliderBenchmark : public Elements::Program {
public:

  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override
  {
    Linx::ProgramOptions options;
    options.positional("knots", "Number of knots", 101L);
    options.positional("args", "Number of arguments", 1001L);
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override
  {
    const auto knot_count = args["knots"].as<Linx::Index>();
    const auto arg_count = args["args"].as<Linx::Index>();

    std::cout << "Generating knots..." << std::endl;
    const auto u = Linx::Sequence<double>(knot_count).linspace(0, Linx::pi<double>() * 2.5);
    const auto x = Linx::Sequence<double>(arg_count).linspace(0, max(u));
    const auto v = sin(u);
    const auto gt = sin(x);

    std::cout << "i\tu\tv" << std::endl;
    for (std::size_t i = 0; i < u.size(); ++i) {
      std::cout << i << '\t' << u[i] << '\t' << v[i] << std::endl;
    }

    std::cout << "Interpolating..." << std::endl;
    std::vector<std::vector<double>> y;
    std::cout << "i\tx\tsin(x)";
    std::cout << "\tC2";
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

    return ExitCode::OK;
  }
};

MAIN_FOR(SpliderBenchmark)
