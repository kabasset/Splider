/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_CO_H
#define _SPLIDER_CO_H

#include "Linx/Base/SeqUtils.h" // IsRange

#include <initializer_list>
#include <vector>

namespace Splider {

/// @cond
template <typename>
class Partition;
/// @endcond

/**
 * @brief Cospline.
 */
template <typename TSpline>
class Co {

public:
  /**
   * @brief The spline type.
   */
  using Method = TSpline;

  /**
   * @brief The knot domain type.
   */
  using Domain = typename Method::Domain;

  /**
   * @brief The abscissae floating point type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief The argument type.
   */
  using Arg = typename Method::Arg;

  /**
   * @brief The knot value type.
   */
  using Value = typename Method::Value;

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Co(const Domain& domain, TIt begin, TIt end) : m_spline(domain), m_args() {
    assign(begin, end);
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  explicit Co(const Domain& u, const TX& x) : Co(u, std::begin(x), std::end(x)) {}

  /**
   * @brief List-based constructor.
   */
  template <typename TX>
  explicit Co(const Domain& u, std::initializer_list<TX> x) : Co(u, x.begin(), x.end()) {}

  /**
   * @brief Get the knots abscissae.
   */
  const Domain& domain() const {
    return m_spline.domain();
  }

  /**
   * @brief Assign arguments from an iterator.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end) {
    m_args.clear();
    m_args.reserve(std::distance(begin, end));
    const auto& d = domain();
    for (; begin != end; ++begin) {
      m_args.emplace_back(d, *begin);
    }
  }

  /**
   * @brief Assign arguments from a range.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  void assign(const TX& x) {
    assign(std::begin(x), std::end(x));
  }

  /**
   * @brief Assign arguments from a list.
   */
  template <typename TX>
  void assign(const std::initializer_list<TX>& x) {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Resample a spline defined by an iterator over knot values.
   */
  template <typename TIt>
  std::vector<Value> operator()(TIt begin, TIt end) {
    m_spline.assign(begin, end);
    return m_spline(m_args);
  }

  /**
   * @brief Resample a spline defined by a range of knot values.
   */
  template <typename TV, typename std::enable_if_t<Linx::IsRange<TV>::value>* = nullptr>
  std::vector<Value> operator()(const TV& v) {
    return operator()(std::begin(v), std::end(v));
  }

  /**
   * @brief Resample a spline defined by a list of knot values.
   */
  template <typename TV>
  std::vector<Value> operator()(std::initializer_list<TV> v) {
    return operator()(v.begin(), v.end());
  }

private:
  Method m_spline; ///< The cached spline
  std::vector<Arg> m_args; ///< The resampling abscissae
};

} // namespace Splider

#endif
