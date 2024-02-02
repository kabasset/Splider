/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MIXINS_HERMITE_H
#define _SPLIDER_MIXINS_HERMITE_H

#include "Splider/mixins/Builder.h"

#include <initializer_list>

namespace Splider {

/**
 * @brief The Hermite spline knot domain.
 */
template <typename TReal>
using HermiteDomain = Partition<TReal>; // FIXME

/**
 * @brief A Hermite spline argument.
 */
template <typename TDomain>
class HermiteArg {
  template <typename, typename, typename>
  friend class HermiteSplineMixin;

public:

  /**
   * @brief The knots domain type.
   */
  using Domain = TDomain;

  /**
   * @brief The argument floating point type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief Constructor.
   */
  HermiteArg(const Domain& domain, Real x)
  {
    m_i = domain.index(x);
    // FIXME
  }

private:

  Linx::Index m_i; ///< The subinterval index
  Real m_cv0; ///< The `m_v[i]` coefficient
  Real m_cv1; ///< The `m_v[i + 1]` coefficient
  Real m_cd0; ///< The `m_d[i]` coefficient
  Real m_cd1; ///< The `m_d[i + 1]` coefficient
};

/**
 * @brief Mixin for Hermite splines.
 * 
 * Hermite splines are C1.
 */
template <typename TDomain, typename TValue, typename TDerived>
class HermiteSplineMixin : public BuilderMixin<TDerived> { // FIXME Upper mixin

public:

  /**
   * @brief The knots domain type.
   */
  using Domain = TDomain;

  /**
   * @brief The abscissae floating point type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief The argument type.
   */
  using Arg = HermiteArg<Domain>;

  /**
   * @brief The knot value type.
   */
  using Value = TValue;

  /**
   * @brief Null knots constructor.
   */
  explicit HermiteSplineMixin(const Domain& u) : m_domain(u), m_v(m_domain.size()), m_d(m_domain.size()), m_valid(true)
  {}

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit HermiteSplineMixin(const Domain& u, TIt begin, TIt end) :
      m_domain(u), m_v(begin, end), m_d(m_v.size()), m_valid(false)
  {}

  /**
   * @brief Evaluate the spline for a given argument.
   */
  Value operator()(const Arg& arg)
  {
    const auto i = arg.m_i;
    static_cast<TDerived&>(*this).update(i);
    return m_v[i] * arg.m_cv0 + m_v[i + 1] * arg.m_cv1 + m_d[i] * arg.m_cd0 + m_d[i + 1] * arg.m_cd1;
  }

protected:

  const Domain& m_domain; ///< The knots domain
  std::vector<Value> m_v; ///< The knot values
  std::vector<Value> m_6s; ///< The knot derivatives
  bool m_valid; ///< Validity flag // FIXME to TDerived
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
