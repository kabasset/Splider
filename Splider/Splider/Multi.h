/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MULTI_H
#define _SPLIDER_MULTI_H

#include "Linx/Base/SeqUtils.h" // IsRange

#include <initializer_list>
#include <vector>

namespace Splider {

/**
 * @brief Multi-dimensional spline.
 */
template <Linx::Index N, typename TSpline>
class Multi {
public:

  /**
   * @brief The domain dimension.
   */
  static constexpr Linx::Index Dimension = N;

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
  using Arg = Linx::Vector<typename Method::Arg, N>;

  /**
   * @brief The knot value type.
   */
  using Value = typename Method::Value;

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Multi(const Domain& domain, TIt begin, TIt end) : m_spline(domain), m_args()
  {
    assign(begin, end);
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  explicit Multi(const Domain& u, const TX& x) : Multi(u, std::begin(x), std::end(x))
  {}

  /**
   * @brief List-based constructor.
   */
  template <typename TX>
  explicit Multi(const Domain& u, std::initializer_list<TX> x) : Multi(u, x.begin(), x.end())
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
  void assign(const TX& x)
  {
    assign(std::begin(x), std::end(x));
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
  template <typename TV, typename std::enable_if_t<Linx::IsRange<TV>::value>* = nullptr>
  std::vector<Value> operator()(const TV& v)
  {
    return operator()(std::begin(v), std::end(v));
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

  Linx::Vector<Method, Dimension> m_splines; ///< The cached splines
  Linx::Vector<std::vector<Arg>, Dimension> m_args; ///< The resampling abscissae
};

} // namespace Splider

#endif
