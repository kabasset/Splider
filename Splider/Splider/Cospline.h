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
class Cospline {

public:
  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Cospline(const Partition& domain, TIt begin, TIt end) : m_domain(domain), m_args(end - begin) {
    std::transform(begin, end, m_args.begin(), [&](const auto& e) {
      return SplineArg(m_domain, e);
    });
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Cospline(const Partition& u, const TRange& x) : Cospline(u, x.begin(), x.end()) {}

  /**
   * @brief List-based constructor.
   */
  template <typename T>
  explicit Cospline(const Partition& u, std::initializer_list<T> x) : Cospline(u, x.begin(), x.end()) {}

  /**
   * @brief Assign arguments from an iterator.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end) {
    m_args.resize(std::distance(begin, end));
    std::transform(begin, end, m_args.begin(), [&](const auto& x) {
      return SplineArg(m_domain, x);
    });
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
  void assign(const std::initializer_list<double>& x) {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Resample a spline defined by an iterator over knot values.
   */
  template <typename TIt>
  auto operator()(TIt begin, TIt end) const {
    using T = typename std::iterator_traits<TIt>::value_type;
    return Spline<std::decay_t<T>>(m_domain, begin, end)(m_args);
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
  const Partition& m_domain; ///< The knot abscissae
  std::vector<SplineArg> m_args; ///< The resampling abscissae
};

} // namespace Splider

#endif
