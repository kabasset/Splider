/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MIXINS_BUILDER_H
#define _SPLIDER_MIXINS_BUILDER_H

#include "Splider/Builder.h"

#include <initializer_list>

namespace Splider {

/**
 * @brief Mixin for cubic spline methods.
 */
template <typename TDerived, typename TBounds>
struct BuilderMixin {
  /**
   * @brief Make a builder.
   */
  template <TBounds B = static_cast<TBounds>(0), typename TIt>
  static auto builder(TIt begin, TIt end)
  {
    using T = typename std::iterator_traits<TIt>::value_type;
    using Real = std::decay_t<T>;
    using Domain = typename TDerived::Domain<Real>;
    return Builder<Domain, TDerived, TBounds, B>(begin, end);
  }

  /**
   * @brief Make a builder.
   */
  template <TBounds B = static_cast<TBounds>(0), typename TRange>
  static auto builder(const TRange& range)
  {
    return builder<B>(std::begin(range), std::end(range));
  }

  /**
   * @brief Make a builder.
   */
  template <TBounds B = static_cast<TBounds>(0), typename T>
  static auto builder(std::initializer_list<T> list)
  {
    return builder<B>(list.begin(), list.end());
  }

  /**
   * @brief Evaluate a spline.
   */
  template <TBounds B = static_cast<TBounds>(0), typename TU, typename TV, typename TX>
  static auto eval(const TU& u, const TV& v, const TX& x)
  {
    const auto build = builder<B>(u);
    auto spline = build.spline(v);
    return spline(x);
  }
};

} // namespace Splider

#endif
