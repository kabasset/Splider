/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MULTIBUILDER_H
#define _SPLIDER_MULTIBUILDER_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Linx/Data/Vector.h"
#include "Splider/BiSpline.h" // FIXME replace with Profile

#include <initializer_list>
#include <iterator>

/**
 * @brief Spline builder.
 */
namespace Splider {

/**
 * @brief Spline builder.
 */
template <Linx::Index N, typename TDomain, typename TMethod, typename TBounds, TBounds B>
class MultiBuilder {
public:

  /**
   * @brief The builder dimension.
   */
  static constexpr Linx::Index Dimension = N;

  /**
   * @brief The knot 1D domain type.
   */
  using Domain = TDomain;

  /**
   * @brief The 1D abscissae floating point type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief The 1D spline type.
   */
  using Method = TMethod;

  /**
   * @brief The spline ND argument type.
   */
  using Arg = Linx::Vector<typename Method::Arg<Domain>, Dimension>;

  /**
   * @brief Constructor.
   * @param us The domains
   */
  template <typename... TIts>
  MultiBuilder(std::pair<TIts, TIts>... us) : m_domains {Domain(us.first, us.second)...}
  {}

  /**
   * @brief Get the knots domain.
   */
  constexpr inline const Domain& domain(Linx::Index axis) const
  {
    return m_domains[axis];
  }

  /**
   * @brief Create a walker with given arguments.
   */
  template <typename TIt>
  auto walker(TIt begin, TIt end) const
  {
    static_assert(Dimension == 2, "Case N != 2 not yet implemented.");
    using T = typename std::iterator_traits<TIt>::value_type;
    using Value = typename std::decay_t<T>::value_type;
    using Spline = typename Method::Spline<Domain, Value, B>;
    return BiCospline<Spline>(m_domains[0], m_domains[1], begin, end);
  }

  /**
   * @brief Create a walker with given arguments.
   */
  template <typename TX, typename std::enable_if_t<Linx::IsRange<TX>::value>* = nullptr>
  auto walker(const TX& x) const
  {
    return walker(std::begin(x), std::end(x));
  }

  /**
   * @brief Create a walker with given arguments.
   */
  template <typename TX>
  auto walker(std::initializer_list<TX> x) const
  {
    return walker(x.begin(), x.end());
  }

private:

  std::array<Domain, Dimension> m_domains; ///< The knot domains.
  // FIXME Cannot use Linx::Vector because Domain has no operator==
};

} // namespace Splider

#endif
