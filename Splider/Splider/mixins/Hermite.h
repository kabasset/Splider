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
    const auto t = (x - domain[m_i]) / domain.length(m_i);
    m_cv0 = (1 + 2 * t) * (1 - t) * (1 - t);
    m_cv1 = t * t * (3 - 2 * t);
    m_cd0 = t * (1 - t) * (1 - t);
    m_cd1 = t * t * (t - 1);
    // FIXME optimize
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
class HermiteSplineMixin { // FIXME Upper mixin

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
   * @brief Get the knots domain.
   */
  inline const Domain& domain() const
  {
    return m_domain;
  }

  /**
   * @brief Assign the knot values.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end)
  {
    m_v.assign(begin, end);
    m_valid = false; // FIXME TDerived::invalidate()
  }

  /**
   * @brief Assign the knot values.
   */
  template <typename TV, typename std::enable_if_t<Linx::IsRange<TV>::value>* = nullptr>
  void assign(const TV& v)
  {
    assign(std::begin(v), std::end(v));
  }

  /**
   * @brief Assign the knot values.
   */
  template <typename TV>
  void assign(std::initializer_list<TV> v)
  {
    assign(v.begin(), v.end());
  }

  /**
   * @brief Set a knot value.
   */
  void set(Linx::Index i, Value v)
  {
    m_v[i] = v;
    m_valid = false; // FIXME TDerived::invalidate(i)
  }

  /**
   * @brief Evaluate the spline for a given argument.
   */
  inline Value operator()(Real x)
  {
    return operator()(Arg(m_domain, x));
  }

  /**
   * @brief Evaluate the spline for a given argument.
   */
  Value operator()(const Arg& arg)
  {
    const auto i = arg.m_i;
    static_cast<TDerived&>(*this).update(i);
    return m_v[i] * arg.m_cv0 + m_v[i + 1] * arg.m_cv1 + m_d[i] * arg.m_cd0 + m_d[i + 1] * arg.m_cd1;
  }

  /**
   * @brief Evaluate the spline for multiple arguments.
   */
  template <typename TIt>
  std::vector<Value> operator()(TIt begin, TIt end)
  {
    std::vector<Value> out;
    out.reserve(std::distance(begin, end));
    for (; begin != end; ++begin) {
      out.push_back(operator()(*begin));
    }
    return out;
  }

  /**
   * @brief Evaluate the spline for multiple arguments.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  std::vector<Value> operator()(const TX& x)
  {
    return operator()(std::begin(x), std::end(x));
  }

  /**
   * @brief Evaluate the spline for multiple arguments.
   */
  template <typename TX>
  std::vector<Value> operator()(std::initializer_list<TX> x)
  {
    return operator()(x.begin(), x.end());
  }

protected:

  const Domain& m_domain; ///< The knots domain
  std::vector<Value> m_v; ///< The knot values
  std::vector<Value> m_d; ///< The knot derivatives
  bool m_valid; ///< Validity flag // FIXME to TDerived
};

} // namespace Splider

#endif
