/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Data/Raster.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Data/Tiling.h"
#include "Linx/Run/Chronometer.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/BiSpline.h"

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

static Elements::Logging logger = Elements::Logging::getLogger("SpliderBenchmark");

using Duration = std::chrono::milliseconds;

template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u0, const U& u1, const V& v, const X& x) {
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

template <typename TDuration, typename U, typename V, typename X, typename Y>
TDuration resample(const U& u, const V& v, const X& x, Y& y, char setup) {
  Linx::Chronometer<TDuration> chrono;
  chrono.start();
  if (setup == 's') {
    Splider::Partition<> domain0(u);
    Splider::Partition<> domain1(u);
    Splider::BiCospline<double> cospline(domain0, domain1, x);
    for (const auto& plane : sections(v)) {
      y = cospline(plane);
    }
  } else if (setup == 'g') {
    y = resample_with_gsl(u, u, v, x);
  } else {
    throw std::runtime_error("Case not implemented");
  }
  return chrono.stop();
}

class SpliderBenchmark : public Elements::Program {

public:
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Linx::ProgramOptions options;
    options.named("case", "Test case: s (Splider), g (GSL)", 's');
    options.named("knots", "Number of knots along each axis", 100L);
    options.named("args", "Number of arguments", 100L);
    options.named("iters", "Numper of iterations", 1L);
    options.named("seed", "Random seed", -1L);
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {
    const auto setup = args["case"].as<char>();
    const auto u_size = args["knots"].as<Linx::Index>();
    const auto x_size = args["args"].as<Linx::Index>();
    const auto v_iters = args["iters"].as<Linx::Index>();
    const auto seed = args["seed"].as<Linx::Index>();

    logger.info("Generating knots...");
    const auto u = Linx::Sequence<double>(u_size).linspace(0, Linx::pi<double>() * 4);
    auto v = Linx::Raster<double, 3>({u_size, u_size, v_iters});
    v.generate(
        [&](const auto& p) {
          auto s = std::sin(u[p[0]]);
          auto c = std::cos(u[p[1]]);
          return s * c * (p[2] + 1);
        },
        v.domain());

    logger.info("Generating trajectory...");
    Splider::Trajectory<2> x(x_size);
    auto si = seed;
    for (auto& xi : x) {
      if (seed != -1) {
        ++si;
      }
      xi.generate(Linx::UniformNoise<double>(u[1], u[u_size - 2], si));
    }
    std::vector<double> y;
    logger.info() << "  u: " << u;
    logger.info() << "  v: " << v;
    logger.info() << "  x: " << x;

    logger.info("Interpolating...");
    const auto duration = resample<Duration>(u, v, x, y, setup);
    logger.info() << "  y: " << Linx::Sequence<double>(y);

    logger.debug("i\t\tx0\tx1\tf(x0, x1)\ty");
    for (Linx::Index i = 0; i < x_size; ++i) {
      const auto x0 = x[i][0];
      const auto x1 = x[i][1];
      logger.debug() << i << '\t' << x0 << '\t' << x1 << '\t' << std::sin(x0) * std::cos(x1) * v_iters << '\t' << y[i];
    }

    logger.info() << "  Done in " << duration.count() << "ms";

    return ExitCode::OK;
  }
};

MAIN_FOR(SpliderBenchmark)