/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#ifndef _SPLIDER_MIXINS_C2_H
#define _SPLIDER_MIXINS_C2_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Splider/Partition.h" // FIXME rm
#include "Splider/mixins/Builder.h"

#include <initializer_list>

namespace Splider {

template <typename TReal>
using C2Domain = Partition<TReal>; // FIXME

/**
 * @brief A \f$C^2\f$ spline argument.
 */
template <typename TDomain>
class C2Arg {
  template <typename, typename, typename>
  friend class C2SplineMixin;

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
  explicit C2Arg(const Domain& domain, Real x)
  {
    m_i = domain.index(x);
    const auto h = domain.length(m_i);
    const auto left = x - domain[m_i];
    const auto right = h - left;
    m_cv0 = right / h;
    m_cv1 = 1. - m_cv0;
    m_c6s0 = right * (right * m_cv0 - h);
    m_c6s1 = left * (left * m_cv1 - h);
  }

  /**
   * @brief Get the subinterval index.
   */
  inline Linx::Index index() const
  {
    return m_i;
  }

private:

  Linx::Index m_i; ///< The subinterval index
  Real m_cv0; ///< The `m_v[i]` coefficient
  Real m_cv1; ///< The `m_v[i + 1]` coefficient
  Real m_c6s0; ///< The `m_6s[i]` coefficient
  Real m_c6s1; ///< The `m_6s[i + 1]` coefficient
};

/**
 * @brief \f$C^2\f$ spline parameters.
 */
template <typename TDomain, typename TValue, typename TDerived>
class C2SplineMixin { // FIXME Upper mixin

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
  using Arg = C2Arg<Domain>;

  /**
   * @brief The knot value type.
   */
  using Value = TValue;

  /**
   * @brief Null knots constructor.
   */
  explicit C2SplineMixin(const Domain& u) : m_domain(u), m_v(m_domain.size()), m_6s(m_domain.size()), m_valid(true) {}

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit C2SplineMixin(const Domain& u, TIt begin, TIt end) :
      m_domain(u), m_v(begin, end), m_6s(m_v.size()), m_valid(false)
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
    return m_v[i] * arg.m_cv0 + m_v[i + 1] * arg.m_cv1 + m_6s[i] * arg.m_c6s0 + m_6s[i + 1] * arg.m_c6s1;
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
  std::vector<Value> m_6s; ///< The knot second derivatives times 6
  bool m_valid; ///< Validity flag // FIXME to TDerived
};

} // namespace Splider

#endif
