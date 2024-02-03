/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/C2.h"
#include "Splider/CatmullRom.h"
#include "Splider/Hermite.h"
#include "Splider/Lagrange.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

static Elements::Logging logger = Elements::Logging::getLogger("SpliderSin");

using Duration = std::chrono::milliseconds;

template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u, const V& v, const X& x)
{
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, u.size());
  std::vector<double> y;
  gsl_spline_init(spline, u.data(), v.data(), u.size());
  for (const auto& e : x) {
    y.push_back(gsl_spline_eval(spline, e, acc));
  }
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline);
  return y;
}

class SpliderBenchmark : public Elements::Program {
public:

  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override
  {
    Linx::ProgramOptions options;
    options.positional("knots", "Number of knots", 101L);
    options.positional("args", "Number of arguments", 1001L);
    options.positional("case", "The test case", std::string("c2"));
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override
  {
    const auto knot_count = args["knots"].as<Linx::Index>();
    const auto arg_count = args["args"].as<Linx::Index>();
    const auto setup = args["case"].as<std::string>();

    logger.info();
    logger.info("Generating knots...");
    const auto u = Linx::Sequence<double>(knot_count).linspace(0, Linx::pi<double>() * 4);
    const auto x = Linx::Sequence<double>(arg_count).linspace(0, Linx::pi<double>() * 4);
    const auto v = sin(u);
    const auto gt = sin(x);

    logger.info("i\tu\tv");
    for (std::size_t i = 0; i < u.size(); ++i) {
      logger.info() << i << '\t' << u[i] << '\t' << v[i];
    }

    logger.info();
    logger.info("Interpolating...");
    std::vector<double> y;
    if (setup == "c2") {
      const auto build = Splider::C2::builder(u);
      auto cospline = build.cospline(x);
      y = cospline(v);
    } else if (setup == "c2fd") {
      const auto build = Splider::C2::FiniteDiff::builder(u);
      auto cospline = build.cospline(x);
      y = cospline(v);
    } else if (setup == "hermite") {
      const auto build = Splider::Hermite::FiniteDiff::builder(u);
      auto cospline = build.cospline(x);
      y = cospline(v);
    } else if (setup == "cr") {
      const auto build = Splider::CatmullRom::Uniform::builder(u);
      auto cospline = build.cospline(x);
      y = cospline(v);
    } else if (setup == "lagrange") {
      const auto build = Splider::Lagrange::builder(u);
      auto cospline = build.cospline(x);
      y = cospline(v);
    }
    const auto gsl = resample_with_gsl(u, v, x);

    logger.info("i\tx\tsin(x)\tSplider\tGSL");
    for (std::size_t i = 0; i < x.size(); ++i) {
      logger.info() << i << '\t' << x[i] << '\t' << gt[i] << '\t' << y[i] << '\t' << gsl[i];
    }

    return ExitCode::OK;
  }
};

MAIN_FOR(SpliderBenchmark)
