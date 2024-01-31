/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MIXINS_CUBIC_H
#define _SPLIDER_MIXINS_CUBIC_H

#include "Splider/Builder.h"
#include "Splider/Linspace.h"
#include "Splider/Partition.h"

#include <initializer_list>

namespace Splider {

/**
 * @brief Mixin for cubic spline methods.
 */
template <typename TDerived>
class CubicMixin {
public:
  /**
   * @brief The underlying method.
   */
  using Method = TDerived;

  /**
   * @brief Make a builder.
   */
  template <typename TIt>
  static auto builder(TIt begin, TIt end) {
    using T = typename std::iterator_traits<TIt>::value_type;
    using Real = std::decay_t<T>;
    return Builder<Partition<Real>, Method>(begin, end);
  }

  /**
   * @brief Make a builder.
   */
  template <typename TRange>
  static auto builder(const TRange& range) {
    return builder(std::begin(range), std::end(range));
  }

  /**
   * @brief Make a builder.
   */
  template <typename T>
  static auto builder(std::initializer_list<T> list) {
    return builder(list.begin(), list.end());
  }
};

} // namespace Splider

#endif
