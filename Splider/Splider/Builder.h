/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_BUILDER_H
#define _SPLIDER_BUILDER_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Splider/Co.h"

#include <initializer_list>
#include <iterator>

/**
 * @brief Spline builder.
 */
namespace Splider {

/**
 * @brief Spline builder.
 */
template <typename TDomain, typename TMethod>
class Builder {
public:

  /**
   * @brief The knot domain type.
   */
  using Domain = TDomain;

  /**
   * @brief The abscissae floating point type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief The spline type.
   */
  using Method = TMethod;

  /**
   * @brief The spline argument type.
   */
  using Arg = typename Method::Arg<Domain>;

  /**
   * @brief Constructor.
   * @param params The domain constructor parameters
   */
  template <typename... Ts>
  Builder(Ts... params) : m_domain(LINX_FORWARD(params)...)
  {}

  /**
   * @brief Get the knots domain.
   */
  inline const Domain& domain() const
  {
    return m_domain;
  }

  /**
   * @brief Create an argument at given abscissa.
   */
  inline Arg arg(Real x) const
  {
    return Method::Arg(m_domain, x);
  }

  /**
   * @brief Create multiple arguments at given abscissae.
   */
  template <typename TIt>
  std::vector<Arg> args(TIt begin, TIt end) const
  {
    std::vector<Arg> out;
    out.reserve(std::distance(begin, end));
    for (; begin != end; ++begin) {
      out.emplace_back(m_domain, *begin);
    }
    return out;
  }

  /**
   * @brief Create multiple arguments at given abscissae.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  std::vector<Arg> args(const TX& x) const
  {
    return args(std::begin(x), std::end(x));
  }

  /**
   * @brief Create multiple arguments at given abscissae.
   */
  template <typename TX>
  std::vector<Arg> args(std::initializer_list<TX> x) const
  {
    return args(x.begin(), x.end());
  }

  /**
   * @brief Create a spline with null knots.
   */
  template <typename TV>
  auto spline() const
  {
    return typename Method::Spline<Domain, TV>(m_domain);
  }

  /**
   * @brief Create a spline with given knot values.
   */
  template <typename TIt>
  auto spline(TIt begin, TIt end) const
  {
    using T = typename std::iterator_traits<TIt>::value_type;
    using Value = std::decay_t<T>;
    return typename Method::Spline<Domain, Value>(m_domain, begin, end);
  }

  /**
   * @brief Create a spline with given knot values.
   */
  template <typename TV, typename std::enable_if_t<Linx::IsRange<TV>::value>* = nullptr>
  auto spline(const TV& v) const
  {
    return spline(std::begin(v), std::end(v));
  }

  /**
   * @brief Create a spline with given knot values.
   */
  template <typename TV>
  auto spline(std::initializer_list<TV> v) const
  {
    return spline(v.begin(), v.end());
  }

  /**
   * @brief Create a cospline with given arguments.
   */
  template <typename TIt>
  auto cospline(TIt begin, TIt end) const
  {
    using T = typename std::iterator_traits<TIt>::value_type;
    using Value = std::decay_t<T>;
    return Co<typename Method::Spline<Domain, Value>>(m_domain, begin, end);
  }

  /**
   * @brief Create a cospline with given arguments.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  auto cospline(const TX& x) const
  {
    return cospline(std::begin(x), std::end(x));
  }

  /**
   * @brief Create a cospline with given arguments.
   */
  template <typename TX>
  auto cospline(std::initializer_list<TX> x) const
  {
    return cospline(x.begin(), x.end());
  }

private:

  Domain m_domain; ///< The knot domain.
};

} // namespace Splider

#endif
