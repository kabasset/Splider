/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ElementsKernel/ProgramHeaders.h"
#include "Linx/Data/Raster.h"
#include "Linx/Data/Sequence.h"
#include "Linx/Data/Tiling.h"
#include "Linx/Run/Chronometer.h"
#include "Linx/Run/ProgramOptions.h"
#include "Splider/C2.h"
#include "Splider/Cospline.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

static Elements::Logging logger = Elements::Logging::getLogger("SpliderBenchmark");

using Duration = std::chrono::milliseconds;

template <typename U, typename V, typename X>
std::vector<double> resample_with_gsl(const U& u, const V& v, const X& x) {
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
TDuration resample(const U& u, const V& v, const X& x, Y& y, const std::string& setup) {
  Linx::Chronometer<TDuration> chrono;
  chrono.start();
  if (setup == "f") {
    using Domain = Splider::Partition<float>;
    const Domain domain(u);
    Splider::Cospline<double, Domain> cospline(domain, x); // FIXME float
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "d") {
    using Domain = Splider::Partition<double>;
    const Domain domain(u);
    Splider::Cospline<double, Domain> cospline(domain, x);
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "s") {
    using Domain = Splider::Partition<double>;
    const Domain domain(u);
    const Splider::Args<double> args(domain, x);
    Splider::Spline<double, Domain> spline(domain);
    for (const auto& row : sections(v)) {
      spline.assign(row);
      y = spline(args);
    }
  } else if (setup == "s2") {
    using Spline = Splider::Natural;
    const auto b = Spline::builder(u);
    auto spline = b.template spline<double>();
    const auto args = b.args(x);
    for (const auto& row : sections(v)) {
      spline.assign(row);
      y = spline(args);
    }
  } else if (setup == "c2") {
    using Spline = Splider::Natural;
    const auto b = Spline::builder(u);
    auto cospline = b.cospline(x);
    for (const auto& row : sections(v)) {
      y = cospline(row);
    }
  } else if (setup == "l") {
    using Domain = Splider::Linspace<double>;
    const Domain domain(u[0], u[1], u.size()); // FIXME add ssize() to Linx
    const Splider::Args<double> args(domain, x);
    Splider::Spline<double, Domain> spline(domain);
    for (const auto& row : sections(v)) {
      spline.assign(row);
      y = spline(args);
    }
  } else if (setup == "g") {
    y = resample_with_gsl(u, v, x);
  } else {
    throw std::runtime_error("Case not implemented");
  }
  return chrono.stop();
}

class SpliderBenchmark : public Elements::Program {

public:
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Linx::ProgramOptions options;
    options.named("case", "Test case: d (double), f (float), s (Spline), g (GSL)", std::string("d"));
    options.named("knots", "Number of knots", 100L);
    options.named("args", "Number of arguments", 100L);
    options.named("iters", "Numper of iterations", 1L);
    options.named("seed", "Random seed", -1L);
    return options.as_pair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {
    const auto setup = args["case"].as<std::string>();
    const auto u_size = args["knots"].as<Linx::Index>();
    const auto x_size = args["args"].as<Linx::Index>();
    const auto v_iters = args["iters"].as<Linx::Index>();
    const auto seed = args["seed"].as<Linx::Index>();

    logger.info("Generating knots...");
    const auto u = Linx::Sequence<double>(u_size).linspace(0, Linx::pi<double>() * 4);
    auto v = Linx::Raster<double, 2>({u_size, v_iters});
    v.generate(
        [&](const auto& p) {
          auto s = std::sin(u[p[0]]);
          return s * (p[1] + 1);
        },
        v.domain());

    logger.info("Generating trajectory...");
    auto x = Linx::Sequence<double>(x_size).generate(Linx::UniformNoise<double>(u[1], u[u_size - 2], seed));
    std::vector<double> y;
    logger.info() << "  u: " << u;
    logger.info() << "  v: " << v;
    logger.info() << "  x: " << x;

    logger.info("Interpolating...");
    const auto duration = resample<Duration>(u, v, x, y, setup);
    logger.info() << "  y: " << Linx::Sequence<double>(y);

    logger.debug("i\t\tx\tf(x)\ty");
    for (Linx::Index i = 0; i < x_size; ++i) {
      logger.debug() << i << '\t' << x[i] << '\t' << std::sin(x[i]) * v_iters << '\t' << y[i];
    }

    logger.info() << "  Done in " << duration.count() << "ms";

    return ExitCode::OK;
  }
};

MAIN_FOR(SpliderBenchmark)