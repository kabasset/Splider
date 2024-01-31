/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MIXINS_C2_H
#define _SPLIDER_MIXINS_C2_H

#include "Splider/mixins/Cubic.h"

#include <initializer_list>

namespace Splider {

/**
 * @brief Mixin for C2 splines.
 */
template <typename TDerived>
class C2Mixin : public CubicMixin<TDerived> {

public:
  /**
   * @brief An argument.
   */
  template <typename TDomain>
  class Arg {
    template <typename, typename>
    friend class TDerived::Spline;

  public:
    /**
     * @brief The domain type.
     */
    using Domain = TDomain;

    /**
     * @brief The argument floating point type.
     */
    using Real = typename Domain::Value;

    /**
     * @brief Constructor.
     */
    explicit Arg(const Domain& domain, Real x) {
      m_i = domain.index(x);
      const auto h = domain.length(m_i);
      const auto left = x - domain[m_i];
      const auto right = h - left;
      m_cv0 = right / h;
      m_cv1 = 1 - m_cv0;
      m_cs0 = right / 6. * (right * m_cv0 - h);
      m_cs1 = left / 6. * (left * m_cv1 - h);
      // FIXME rm 6 ?
    }

  private:
    Linx::Index m_i; ///< The subinterval index
    Real m_cv0; ///< The `v[i]` coefficient
    Real m_cv1; ///< The `v[i + 1]` coefficient
    Real m_cs0; ///< The `s[i]` coefficient
    Real m_cs1; ///< The `s[i + 1]` coefficient
  };

  /**
   * @brief The parameters.
   */
  template <typename TDomain, typename T>
  class Spline { // FIXME Mixin for operator()?
    friend class TDerived::Method;

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
    using Arg = typename TDerived::Arg<Domain>;

    /**
     * @brief The knot value type.
     */
    using Value = T;

    /**
     * @brief Null knots constructor.
     */
    explicit Spline(const Domain& u) : m_domain(u), m_v(m_domain.size()), m_s(m_domain.size()), m_valid(true) {}

    /**
     * @brief Iterator-based constructor.
     */
    template <typename TIt>
    explicit Spline(const Domain& u, TIt begin, TIt end) :
        m_domain(u), m_v(begin, end), m_s(m_v.size()), m_valid(false) {}

    /**
   * @brief Range-based constructor.
   */
    template <typename TV>
    explicit Spline(const Domain& u, const TV& v) : Spline(u, std::begin(v), std::end(v)) {}

    /**
   * @brief List-based constructor.
   */
    template <typename TV>
    explicit Spline(const Domain& u, std::initializer_list<TV> v) : Spline(u, v.begin(), v.end()) {}

    /**
     * @brief Get the knots domain.
     */
    inline const Domain& domain() const {
      return m_domain;
    }

    /**
     * @brief Assign the knot values.
     */
    template <typename TIt>
    void assign(TIt begin, TIt end) {
      m_v.assign(begin, end);
      m_valid = false;
    }

    /**
     * @brief Assign the knot values.
     */
    template <typename TV>
    void assign(const TV& v) {
      assign(std::begin(v), std::end(v));
    }

    /**
     * @brief Assign the knot values.
     */
    template <typename TV>
    void assign(std::initializer_list<TV> v) {
      assign(v.begin(), v.end());
    }

    /**
     * @brief Evaluate the spline for a given argument.
     */
    inline Value operator()(Real x) {
      return operator()(Arg(m_domain, x));
    }

    /**
     * @brief Evaluate the spline for a given argument.
     */
    Value operator()(const Arg& arg) {
      const auto i = arg.m_i;
      TDerived::update(*this, i);
      return m_v[i] * arg.m_cv0 + m_v[i + 1] * arg.m_cv1 + m_s[i] * arg.m_cs0 + m_s[i + 1] * arg.m_cs1;
    }

    /**
     * @brief Evaluate the spline for multiple arguments.
     */
    template <typename TIt>
    std::vector<Value> operator()(TIt begin, TIt end) {
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
    template <typename TX>
    std::vector<Value> operator()(const TX& x) {
      return operator()(std::begin(x), std::end(x));
    }

    /**
     * @brief Evaluate the spline for multiple arguments.
     */
    template <typename TX>
    std::vector<Value> operator()(std::initializer_list<TX> x) {
      return operator()(x.begin(), x.end());
    }

  private:
    const Domain& m_domain; ///< The knots domain
    std::vector<Value> m_v; ///< The knot values
    std::vector<Value> m_s; ///< The knot second derivatives
    bool m_valid; ///< Validity flag // FIXME as std::set
  };
};

} // namespace Splider

#endif
