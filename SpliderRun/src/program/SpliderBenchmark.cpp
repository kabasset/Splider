/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Data/Raster.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Run/Chronometer.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/BiSpline.h"

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

static Elements::Logging logger = Elements::Logging::getLogger("LinxBenchmarkConvolution");

using Duration = std::chrono::milliseconds;

template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u0, const U& u1, const V& v, const X& x) {
  gsl_interp_accel* xacc = gsl_interp_accel_alloc();
  gsl_interp_accel* yacc = gsl_interp_accel_alloc();
  gsl_spline2d* spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, u0.size(), u1.size());
  gsl_spline2d_init(spline, u0.data(), u1.data(), v.data(), u0.size(), u1.size());
  std::vector<double> y;
  for (const auto& e : x) {
    y.push_back(gsl_spline2d_eval(spline, e[0], e[1], xacc, yacc));
  }
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  gsl_spline2d_free(spline);
  return y;
}

template <typename TDuration, typename U, typename V, typename X, typename Y>
TDuration resample(const U& u, const V& v, const X& x, Y& y, char setup) {
  Linx::Chronometer<TDuration> chrono;
  Splider::SplineIntervals domain0(u);
  Splider::SplineIntervals domain1(u);
  chrono.start();
  switch (setup) {
    case 'e':
      y = Splider::BiSplineResampler<double, Splider::SplineCache::Early>(domain0, domain1, x)(v);
      break;
    case 'l':
      y = Splider::BiSplineResampler<double, Splider::SplineCache::Lazy>(domain0, domain1, x)(v);
      break;
    case 'g':
      y = resample_with_gsl(u, u, v, x);
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
    options.named("case", "Test case: e (early), l (lazy), g (GSL)", 'e');
    options.named("knots", "Number of knots along each axis", 100L);
    options.named("args", "Number of arguments", 100L);
    options.named("seed", "Random seed", -1L);
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {
    const auto setup = args["case"].as<char>();
    const auto v_size = args["knots"].as<Linx::Index>();
    const auto x_size = args["args"].as<Linx::Index>();
    const auto seed = args["seed"].as<Linx::Index>();

    logger.info("Generating knots...");
    auto u = Linx::Sequence<double>(v_size).range();
    auto v = Linx::Raster<double>({v_size, v_size}).generate(Linx::GaussianNoise<double>(0, 1, seed));
    Linx::Sequence<Linx::Vector<double>> x(x_size);
    for (auto& xi : x) {
      xi.generate(Linx::UniformNoise<double>(0, v_size - 1, seed));
    }
    std::vector<double> y;
    logger.info() << "  u: " << u;
    logger.info() << "  v: " << v;
    logger.info() << "  x: " << x;

    logger.info("Interpolating...");
    const auto duration = resample<Duration>(u, v, x, y, setup);
    logger.info() << "  y: " << Linx::Sequence<double>(y);

    logger.info() << "  Done in " << duration.count() << "ms";

    return ExitCode::OK;
  }
};

MAIN_FOR(SpliderBenchmark)