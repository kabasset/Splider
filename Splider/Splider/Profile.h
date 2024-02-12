/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#ifndef _SPLIDER_PROFILE_H
#define _SPLIDER_PROFILE_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Linx/Data/Vector.h"

#include <initializer_list>
#include <vector>

namespace Splider {

/**
 * @brief Profile along a path (sparse positions) over a multi-dimensional spline.
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
      m_splines(domain.begin(), domain.end()), m_x(std::distance(begin, end)), m_mask(shape(domain))
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
   * @brief Assign arguments from an iterator.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end)
  {
    using T = std::decay_t<typename std::iterator_traits<TIt>::value_type>;

    // Assign the arguments
    std::vector<T> x(std::distance(begin, end));
    for (Linx::Index i = 0; i < m_args.ssize(); ++i) {
      std::transform(begin, end, x.begin(), [=](const auto& xj) {
        return xj[i];
      });
      m_args[i].assign(x);
    }

    // Assign the mask
    m_mask.clear();
    Linx::Position<Dimension> p(m_args.ssize());
    for (std::size_t j = 0; j < x.size(); ++j) {
      for (Linx::Index i = 0; i < m_args.ssize(); ++i) {
        p[i] = m_args[i][j].index();
      }
      Linx::Box<Dimension> neighborhood(p - 1, p + 2); // FIXME use Method::Radius
      neighborhood &= m_mask.domain();
      m_mask(neighborhood).fill(true);
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
   * @brief Resample a spline defined by knot values data.
   */
  template <typename T>
  std::vector<Value> operator()(const T* data)
  {
    Linx::PtrRaster<Value, const T> v(m_shape, data);
    std::array<Value, 4> v; // FIXME use Method::Radius
    for (const auto& e : m_args) {
      for (Linx::Index j = -1; j < 3; ++j) {
        // FIXME
      }
    }
    m_cospline.assign(begin, end);
    return 0; // FIXME
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

  /**
   * @brief Compute the shape of a given multidimensional domain.
   */
  static Position<Dimension> shape(const Domain& u)
  {
    Position<Dimension> domain(u.size());
    domain.generate(
        [](const auto& ui) {
          return ui.ssize();
        },
        u);
    return domain;
  }

  Linx::Vector<Method, Dimension> m_splines; ///< The cached splines
  Linx::Vector<std::vector<typename Method::Arg>, Dimension> m_x; ///< The cached arguments
  // FIXME as Raster<Arg, 2>?
  Linx::Raster<bool, Dimension> m_mask; ///< The neighboring knot abscissae
};

} // namespace Splider

#endif
