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
    }

  private:
    Linx::Index m_i; ///< The subinterval index
    Real m_cs1; ///< The `d[i + 1]` coefficient
    Real m_cs0; ///< The `d[i]` coefficient
    Real m_cv1; ///< The `v[i + 1]` coefficient
    Real m_cv0; ///< The `v[i]` coefficient
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
     * @brief Evaluate the spline for a given argument.
     */
    Value operator()(Real x) {
      return operator()(Arg(m_domain, x));
    }

    /**
     * @brief Evaluate the spline for a given argument.
     */
    Value operator()(const Arg<Domain>& arg) {
      const auto i = arg.m_i;
      TDerived::update(*this, i);
      return m_s[i + 1] * arg.m_cs1 + m_s[i] * arg.m_cs0 + m_v[i + 1] * arg.m_cv1 + m_v[i] * arg.m_cv0;
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

/**
 * @brief Mixin for Hermite splines.
 * 
 * Hermite splines are C1.
 */
template <typename TDerived>
class HermiteMixin : public CubicMixin<TDerived> {

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
    Real cd1; ///< The `d[i + 1]` coefficient
    Real cd0; ///< The `d[i]` coefficient
    Real cv1; ///< The `v[i + 1]` coefficient
    Real cv0; ///< The `v[i]` coefficient
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
    Value operator()(const Arg<TReal>& arg) {
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

/**
 * @brief Lagrange cubic spline.
 * 
 * This spline is built by fitting a Lagrange cubic polynomial over a sliding 4-knot window.
 * It is only guaranteed to be C0.
 */
class Lagrange;

} // namespace Splider

#endif
