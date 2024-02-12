/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#include "Splider/C2.h"
#include "Splider/Lagrange.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(SpliderDemo_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(default_spline_test)
{
  //! [Default spline]

  using Spline = Splider::C2;
  const auto build = Spline::builder({1, 2, 3, 4});
  auto spline = build.spline({10, 20, 30, 40});
  auto y = spline({1.1, 2.5, 3.9});

  //! [Default spline]
}

BOOST_AUTO_TEST_CASE(lagrange_spline_test)
{
  //! [Lagrange spline]

  using Spline = Splider::Lagrange;
  const auto build = Spline::builder({1, 2, 3, 4});
  auto spline = build.spline({10, 20, 30, 40});
  auto y = spline({1.1, 2.5, 3.9});

  //! [Lagrange spline]
}

BOOST_AUTO_TEST_CASE(custom_bounds_test)
{
  //! [Custom bounds]

  using Spline = Splider::C2;
  const auto build = Spline::builder<Spline::Bounds::NotAKnot>({1, 2, 3, 4});
  auto spline = build.spline({10, 20, 30, 40});
  auto y = spline({1.1, 2.5, 3.9});

  //! [Custom bounds]
}

BOOST_AUTO_TEST_CASE(eval_shortcut_test)
{
  //! [Eval shortcut]

  std::vector<double> u {1, 2, 3, 4};
  std::vector<double> v {10, 20, 30, 40};
  std::vector<double> x {1.1, 2.5, 3.9};

  using Spline = Splider::C2;
  auto y = Spline::eval<Spline::Bounds::NotAKnot>(u, v, x);

  //! [Eval shortcut]
}

BOOST_AUTO_TEST_CASE(default_cospline_test)
{
  //! [Default cospline]

  using Spline = Splider::C2;
  const auto build = Spline::builder({1, 2, 3, 4});
  auto cospline = build.cospline({1.1, 2.5, 3.9});
  auto y = cospline({10, 20, 30, 40});

  //! [Default cospline]
}

BOOST_AUTO_TEST_CASE(default_bivariate_cospline_test)
{
  //! [Default bivariate cospline]

  Linx::Sequence<double> u0 {1, 2, 3, 4};
  Linx::Sequence<double> u1 {1, 10, 100};
  Splider::Trajectory<2> x {{1.1, 2.}, {2.5, 10.}, {2.5, 20.}, {2.5, 50.}, {3.9, 50.}};
  Linx::Raster<double> v({u0.ssize(), u1.ssize()}, {1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400});

  using Spline = Splider::Lagrange;
  const auto build = Spline::Multi::builder(u0, u1);
  auto cospline = build.cospline(x);
  auto y = cospline(v);

  //! [Default bivariate cospline]
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
