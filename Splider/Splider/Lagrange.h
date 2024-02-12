/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#ifndef _SPLIDER_LAGRANGE_H
#define _SPLIDER_LAGRANGE_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Splider/Partition.h" // TODO rm
#include "Splider/mixins/Builder.h"

namespace Splider {

/**
 * @brief The boundary conditions.
 */
enum class LagrangeBounds {
  NotAKnot = 0 ///< Use neighboring polynomials
};

template <typename TDomain>
class LagrangeArg {
  template <typename, typename, LagrangeBounds>
  friend class LagrangeSpline;

public:

  using Domain = TDomain;

  using Real = typename Domain::Value;

  LagrangeArg(const Domain& domain, Real x)
  {
    m_i = std::clamp(domain.index(x), 1L, domain.ssize() - 3);
    const auto u0 = domain[m_i - 1];
    const auto u1 = domain[m_i];
    const auto u2 = domain[m_i + 1];
    const auto u3 = domain[m_i + 2];
    const auto l01 = (x - u0) / (u1 - u0);
    const auto l02 = (x - u0) / (u2 - u0);
    const auto l03 = (x - u0) / (u3 - u0);
    const auto l10 = (x - u1) / (u0 - u1);
    const auto l12 = (x - u1) / (u2 - u1);
    const auto l13 = (x - u1) / (u3 - u1);
    const auto l20 = (x - u2) / (u0 - u2);
    const auto l21 = (x - u2) / (u1 - u2);
    const auto l23 = (x - u2) / (u3 - u2);
    const auto l30 = (x - u3) / (u0 - u3);
    const auto l31 = (x - u3) / (u1 - u3);
    const auto l32 = (x - u3) / (u2 - u3);
    m_l[0] = l10 * l20 * l30;
    m_l[1] = l01 * l21 * l31;
    m_l[2] = l02 * l12 * l32;
    m_l[3] = l03 * l13 * l23;
    // TODO optimize, e.g. through LagrangeDomain
  }

  /**
   * @brief Get the subinterval index.
   */
  inline Linx::Index index() const
  {
    return m_i;
  }

private:

  Linx::Index m_i;
  std::array<Real, 4> m_l;
};

/**
 * @brief The spline evaluator.
 */
template <typename TDomain, typename T, LagrangeBounds B>
class LagrangeSpline {
public:

  /**
   * @brief The knots domain.
   */
  using Domain = TDomain;

  /**
   * @brief The abscissae floating point type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief The argument type.
   */
  using Arg = LagrangeArg<Domain>;

  /**
   * @brief The knot value type.
   */
  using Value = T;

  /**
   * @brief Null knots constructor.
   */
  explicit LagrangeSpline(const Domain& u) : m_domain(u), m_v(m_domain.size()) {}

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit LagrangeSpline(const Domain& u, TIt begin, TIt end) : m_domain(u), m_v(begin, end)
  {}

  /**
   * @brief Range-based constructor.
   */
  template <typename TV, typename std::enable_if_t<Linx::IsRange<TV>::value>* = nullptr>
  explicit LagrangeSpline(const Domain& u, const TV& v) : LagrangeSpline(u, std::begin(v), std::end(v))
  {}

  /**
   * @brief List-based constructor.
   */
  template <typename TV>
  explicit LagrangeSpline(const Domain& u, std::initializer_list<TV> v) : LagrangeSpline(u, v.begin(), v.end())
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
    return std::inner_product(arg.m_l.begin(), arg.m_l.end(), &m_v[i - 1], Value());
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
};

/**
 * @ingroup builders
 * @brief Piecewise cubic Lagrange polynomials (\f$C^0\f$).
 * 
 * This spline is built by fitting a cubic Lagrange polynomial over a sliding 4-knot window.
 * It is only guaranteed to be (\f$C^0\f$).
 * 
 * In the first and last subintervals, cubic polynomials cannot be fitted,
 * and the next and previous ones are used, respectively.
 */
struct Lagrange : BuilderMixin<Lagrange, LagrangeBounds> {
  /**
   * @brief The knots domain type.
   */
  template <typename TReal>
  using Domain = Partition<TReal>; // TODO LagrangeDomain

  /**
   * @brief The argument type.
   */
  template <typename TDomain>
  using Arg = LagrangeArg<TDomain>;

  /**
   * @brief The spline evaluator.
   */
  template <typename TDomain, typename T, LagrangeBounds B>
  using Spline = LagrangeSpline<TDomain, T, B>;
};

} // namespace Splider

#endif
