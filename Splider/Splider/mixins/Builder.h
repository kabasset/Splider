/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MIXINS_BUILDER_H
#define _SPLIDER_MIXINS_BUILDER_H

#include "Splider/Builder.h"
#include "Splider/MultiBuilder.h"

#include <initializer_list>
#include <utility> // pair

namespace Splider {

/**
 * @brief Mixin for cubic spline methods.
 */
template <typename TDerived, typename TBounds>
struct BuilderMixin {
  /**
   * @brief Multidimensional builders.
   */
  struct Multi {
    /**
     * @brief Make a ND builder.
     */
    template <TBounds B = static_cast<TBounds>(0), typename TIt0, typename... TIts>
    static auto builder(std::pair<TIt0, TIt0> u0, std::pair<TIts, TIts>... us)
    {
      using T = typename std::iterator_traits<TIt0>::value_type;
      using Real = std::decay_t<T>;
      using Domain = typename TDerived::Domain<Real>;
      static constexpr Linx::Index N = sizeof...(TIts) + 1;
      return MultiBuilder<N, Domain, TDerived, TBounds, B>(u0, us...);
    }

    /**
     * @brief Make a ND builder.
     */
    template <
        TBounds B = static_cast<TBounds>(0),
        typename TU0,
        typename std::enable_if_t<Linx::IsRange<TU0>::value>* = nullptr,
        typename... TUs>
    static auto builder(TU0&& u0, TUs&&... us)
    {
      return builder(std::make_pair(std::begin(u0), std::end(u0)), std::make_pair(std::begin(us), std::end(us))...);
    }

    /**
     * @brief Make a ND builder.
     */
    template <TBounds B = static_cast<TBounds>(0), typename TIt0, typename... TIts>
    static auto builder(std::initializer_list<TIt0> u0, std::initializer_list<TIts>... us)
    {
      return builder(std::make_pair(u0.begin(), u0.end()), std::make_pair(us.begin(), us.end())...);
    }
  };

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
