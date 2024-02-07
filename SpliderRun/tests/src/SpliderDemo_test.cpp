/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

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

  using Spline = Splider::Lagrange;
  const auto build = Spline::Multi::builder({1, 2, 3, 4}, {10, 20, 30, 40});
  // auto walker = build.walker({{1.1, 11}, {2.5, 25}, {3.9, 39}});
  // FIXME walker(v)

  //! [Default bivariate cospline]
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
