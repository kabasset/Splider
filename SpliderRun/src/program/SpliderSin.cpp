/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/Spline.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

static Elements::Logging logger = Elements::Logging::getLogger("SpliderSin");

using Duration = std::chrono::milliseconds;

template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u, const V& v, const X& x) {
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
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Linx::ProgramOptions options;
    options.positional("knots", "Number of knots", 101L);
    options.positional("args", "Number of arguments", 1001L);
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {
    const auto knot_count = args["knots"].as<Linx::Index>();
    const auto arg_count = args["args"].as<Linx::Index>();

    logger.info("Generating knots...");
    const auto u = Linx::Sequence<double>(knot_count).linspace(0, Linx::pi<double>() * 4);
    const auto x = Linx::Sequence<double>(arg_count).linspace(0, Linx::pi<double>() * 4);
    const auto v = sin(u);
    const auto gt = sin(x);

    logger.info("Interpolating...");
    Splider::Partition domain(u);
    Splider::Cospline cospline(domain, x);
    const auto y = cospline(v);
    const auto gsl = resample_with_gsl(u, v, x);

    logger.info("i\tx\tsin(x)\tSplider\tGSL");
    for (Linx::Index i = 0; i < x.size(); ++i) {
      logger.info() << i << '\t' << x[i] << '\t' << gt[i] << '\t' << y[i] << '\t' << gsl[i];
    }

    return ExitCode::OK;
  }
};

MAIN_FOR(SpliderBenchmark)
