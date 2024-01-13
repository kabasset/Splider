/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Run/Chronometer.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/Spline.h"

#include <algorithm>

static Elements::Logging logger = Elements::Logging::getLogger("LinxBenchmarkConvolution");

using Duration = std::chrono::milliseconds;

template <typename TDuration, typename U, typename V, typename X, typename Y>
TDuration interpolate(const U& u, const V& v, const X& x, Y& y, char setup) {
  Linx::Chronometer<TDuration> chrono;
  Splider::SplineIntervals domain(u);
  chrono.start();
  switch (setup) {
    case 'd':
      y = Splider::Spline<double>(domain, v)(x);
      break;
    case 'l':
      y = Splider::Spline<double, Splider::SplineCache::Lazy>(domain, v)(x);
      break;
    case 'r':
      y = Splider::SplineResampler(domain, x)(v);
      break;
    default:
      throw std::runtime_error("Case not implemented");
  }
  return chrono.stop();
}

class SpliderBenchmark : public Elements::Program {

public:
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Linx::ProgramOptions options;
    options.named("case", "Test case: d (default interpolant), l (lazy interpolant), r (resampler), g (GSL)", 'd');
    options.named("knots", "Number of knots", 100L);
    options.named("args", "Number of arguments", 100L);
    options.named("seed", "Random seed", -1L);
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {
    const auto setup = args["case"].as<char>();
    const auto vSize = args["knots"].as<Linx::Index>();
    const auto xSize = args["args"].as<Linx::Index>();
    const auto seed = args["seed"].as<Linx::Index>();

    logger.info("Generating knots...");
    auto u = Linx::Sequence<double>(vSize).range();
    auto v = Linx::Sequence<double>(vSize).generate(Linx::GaussianNoise<double>(0, 1, seed));
    auto x = Linx::Sequence<double>(xSize).generate(Linx::UniformNoise<double>(0, vSize - 1, seed));
    std::sort(x.begin(), x.end());
    std::vector<double> y;
    logger.info() << "  u: " << u;
    logger.info() << "  v: " << v;
    logger.info() << "  x: " << x;

    logger.info("Interpolating...");
    const auto duration = interpolate<Duration>(u, v, x, y, setup);
    logger.info() << "  y: " << Linx::Sequence<double>(y);

    logger.info() << "  Done in " << duration.count() << "ms";

    return ExitCode::OK;
  }
};

MAIN_FOR(SpliderBenchmark)