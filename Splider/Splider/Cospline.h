/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_COSPLINE_H
#define _SPLIDER_COSPLINE_H

#include "Splider/Spline.h"

#include <algorithm>
#include <vector>

namespace Splider {

/**
 * @brief Natural cubic spline resampler.
 * 
 * A resampler is parametrized with the list of knot and resampling abscissae,
 * and is evaluated on a vector of knot values.
 */
template <typename TDomain>
class Cospline {

public:
  /**
   * @brief The knot domain type.
   */
  using Domain = TDomain;

  /**
   * @brief The real number type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief The argument type.
   */
  using Arg = SplineArg<Real>;

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Cospline(const Domain& domain, TIt begin, TIt end) : m_domain(domain), m_args(m_domain, begin, end) {}

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Cospline(const Domain& u, const TRange& x) : Cospline(u, x.begin(), x.end()) {}

  /**
   * @brief List-based constructor.
   */
  template <typename T>
  explicit Cospline(const Domain& u, std::initializer_list<T> x) : Cospline(u, x.begin(), x.end()) {}

  /**
   * @brief Assign arguments from an iterator.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end) {
    m_args = Args<Real>(m_domain, begin, end);
  }

  /**
   * @brief Assign arguments from a range.
   */
  template <typename TRange>
  void assign(const TRange& x) {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Assign arguments from a list.
   */
  void assign(const std::initializer_list<Real>& x) {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Resample a spline defined by an iterator over knot values.
   */
  template <typename TIt>
  auto operator()(TIt begin, TIt end) const {
    using T = typename std::iterator_traits<TIt>::value_type;
    using Value = typename std::decay_t<T>;
    return Spline<Value, Domain>(m_domain, begin, end)(m_args);
  }

  /**
   * @brief Resample a spline defined by a range of knot values.
   */
  template <typename TRange>
  auto operator()(const TRange& v) const {
    return operator()(v.begin(), v.end());
  }

  /**
   * @brief Resample a spline defined by a list of knot values.
   */
  template <typename T>
  auto operator()(std::initializer_list<T> v) const {
    return operator()(v.begin(), v.end());
  }

private:
  const Domain& m_domain; ///< The knot abscissae
  Args<Real> m_args; ///< The resampling abscissae
};

} // namespace Splider

#endif
