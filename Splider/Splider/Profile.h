/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_PROFILE_H
#define _SPLIDER_PROFILE_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Linx/Data/Vector.h"

#include <initializer_list>
#include <vector>

namespace Splider {

/**
 * @brief Profile along a path over a multi-dimensional spline.
 */
template <typename TSpline, Linx::Index N = 2>
class Profile {
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
   * @brief The knot multidimensional domain type.
   */
  using Domain = Linx::Vector<typename Method::Domain, Dimension>;

  /**
   * @brief The abscissae floating point type.
   */
  using Real = typename Method::Domain::Value;

  /**
   * @brief The multidimensional argument type.
   */
  using Arg = Linx::Vector<typename Method::Arg, Dimension>;

  /**
   * @brief The path type.
   */
  using Path = std::vector<Arg>;

  /**
   * @brief The knot value type.
   */
  using Value = typename Method::Value;

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Profile(const Domain& domain, TIt begin, TIt end) :
      m_splines(domain.begin(), domain.end()), m_args(shape(domain))
  {
    assign(begin, end);
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  explicit Profile(const Domain& u, const TX& x) : Profile(u, std::begin(x), std::end(x))
  {}

  /**
   * @brief List-based constructor.
   */
  template <typename TX>
  explicit Profile(const Domain& u, std::initializer_list<TX> x) : Profile(u, x.begin(), x.end())
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
      m_args.emplace_back(d, *begin); // FIXME
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
