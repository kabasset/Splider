/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MIXINS_HERMITE_H
#define _SPLIDER_MIXINS_HERMITE_H

#include "Splider/mixins/Builder.h"

#include <initializer_list>

namespace Splider {

/**
 * @brief Mixin for Hermite splines.
 * 
 * Hermite splines are C1.
 */
template <typename TDerived>
class HermiteMixin : public BuilderMixin<TDerived> {
protected:

  /**
   * @brief An argument.
   */
  template <typename T>
  struct Arg {
    /**
     * @brief The knot value type.
     */
    using Real = T;

    Linx::Index i; ///< The subinterval index
    Real cv0; ///< The `v[i]` coefficient
    Real cv1; ///< The `v[i + 1]` coefficient
    Real cd0; ///< The `d[i]` coefficient
    Real cd1; ///< The `d[i + 1]` coefficient
  };

  /**
   * @brief The parameters.
   */
  template <typename T>
  struct Params {
    /**
     * @brief The knot value type.
     */
    using Value = T;

    /**
     * @brief Evaluate the spline for a given argument.
     */
    template <typename TReal>
    Value operator()(const Arg<TReal>& arg)
    {
      const auto i = arg.index;
      return v[i] * arg.cv0 + v[i + 1] * arg.cv1 + d[i] * arg.cd0 + d[i + 1] * arg.cd1;
    }

    std::vector<Value> v; ///< The knot values
    std::vector<Value> d; ///< The knot derivatives
  };
};

class Hermite {
public:

  /**
   * @brief Finite difference Hermite spline.
   * 
   * Tangents are computed as finite differences from the neighboring knots.
   */
  class FiniteDiff : public HermiteMixin<FiniteDiff> {};
};

/**
 * @brief Catmull-Rom spline.
 */
class CatmullRom {
public:

  /**
   * @brief Catmull-Rom spline with uniform parametrization.
   */
  class Uniform : public HermiteMixin<Uniform> {};
};

/**
 * @brief Akima spline.
 */
class Akima : public HermiteMixin<Akima> {
  /**
 * @brief Modified Akima spline.
 * 
 * TODO
 */
  class Modified : public HermiteMixin<Modified> {};
};

/**
 * @brief Monotone cubic spline.
 */
class Monotone : public HermiteMixin<Monotone> {};

} // namespace Splider

#endif
