/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#ifndef _SPLIDER_COSPLINE_H
#define _SPLIDER_COSPLINE_H

#include "Splider/Spline.h"

#include <algorithm>
#include <vector>

namespace Splider {

/// @cond
template <typename>
class Partition;
/// @endcond

/**
 * @brief Natural cubic spline resampler.
 * 
 * A resampler is parametrized with the list of knot and resampling abscissae,
 * and is evaluated on a vector of knot values.
 */
template <typename T, typename TDomain = Partition<double>>
class Cospline {
public:

  /**
   * @brief The knot value type.
   */
  using Value = T;

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
  explicit Cospline(const Domain& domain, TIt begin, TIt end) : m_spline(domain), m_args(domain, begin, end)
  {}

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Cospline(const Domain& u, const TRange& x) : Cospline(u, x.begin(), x.end())
  {}

  /**
   * @brief List-based constructor.
   */
  template <typename TX>
  explicit Cospline(const Domain& u, std::initializer_list<TX> x) : Cospline(u, x.begin(), x.end())
  {}

  /**
   * @brief Get the knots abscissae.
   */
  const Domain& domain() const
  {
    return m_spline.domain();
  }

  /**
   * @brief Assign arguments from an iterator.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end)
  {
    m_args = Args<Real>(m_spline.domain(), begin, end);
  }

  /**
   * @brief Assign arguments from a range.
   */
  template <typename TRange>
  void assign(const TRange& x)
  {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Assign arguments from a list.
   */
  template <typename TX>
  void assign(const std::initializer_list<TX>& x)
  {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Resample a spline defined by an iterator over knot values.
   */
  template <typename TIt>
  std::vector<Value> operator()(TIt begin, TIt end)
  {
    m_spline.assign(begin, end);
    return m_spline(m_args);
  }

  /**
   * @brief Resample a spline defined by a range of knot values.
   */
  template <typename TRange>
  std::vector<Value> operator()(const TRange& v)
  {
    return operator()(v.begin(), v.end());
  }

  /**
   * @brief Resample a spline defined by a list of knot values.
   */
  template <typename TV>
  std::vector<Value> operator()(std::initializer_list<TV> v)
  {
    return operator()(v.begin(), v.end());
  }

private:

  Spline<Value, Domain> m_spline; ///< The cached `Spline`
  Args<Real> m_args; ///< The resampling abscissae
};

} // namespace Splider

#endif
