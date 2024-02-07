/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_BISPLINE_H
#define _SPLIDER_BISPLINE_H

#include "Linx/Data/Mask.h"
#include "Linx/Data/Raster.h"
#include "Linx/Data/Sequence.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace Splider {

/**
 * @brief Alias for a sequence of 2D position.
 * @see `BiCospline`
 */
template <Linx::Index N, typename T = double>
using Trajectory = Linx::Sequence<Linx::Vector<T, N>>;

/**
 * @brief Bivariate natural cubic spline resampler.
 * 
 * The transform maps an input 2D rectilinear grid `u` to an output trajectory, i.e. collection of coordinates `x`.
 * The resampler is called on a 2D map of knot values `v` to output the interpolated values `y` at abscissae `x`.
 * 
 * The grid is specified as a pair of `Partition` parameters.
 * 
 * The trajectory is made of objects on which `operator[]()` is called to get the first and second components.
 * Typical classes which fulfill this requirement are raw arrays, `std::array<double, 2>` or `Linx::Vector<double, 2>`,
 * but `std::pair<double, double>` is not compatible.
 * Alias `Trajectory<2>` is defined for conciseness as a `Linx::Sequence` (see [Linx](https://github.com/kabasset/Linx) documentation).
 * 
 * The values are given as a 2D `Linx::Raster`.
 * 
 * For optimization purpose, as opposed to the 1D `Cospline`, this class is templated by the value type.
 * 
 * Similarly to `Spline`, the resampler can rely on various caching strategies:
 * see `Caching` documentation for selecting the most appropriate one.
 */
template <typename TSpline>
class BiCospline {
public:

  /**
   * @brief The dimension.
   */
  static constexpr Linx::Index Dimension = 2;

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
  BiCospline(const Domain& domain0, const Domain& domain1, TIt begin, TIt end) :
      m_splines0(domain1.size(), Method(domain0)), m_spline1(domain1), m_x(),
      m_mask(Linx::Position<Dimension>::zero(), {domain0.ssize() - 1, domain1.ssize() - 1}, false)
  {
    for (; begin != end; ++begin) {
      std::array<Arg, Dimension> xi {Arg(domain0, (*begin)[0]), Arg(domain1, (*begin)[1])};
      const auto i0 = xi[0].index();
      const auto i1 = xi[1].index();
      const auto min0 = std::max(i0 - 1, 0L);
      const auto max0 = std::min(i0 + 2, domain0.ssize() - 1);
      const auto min1 = std::max(i1 - 1, 0L);
      const auto max1 = std::min(i1 + 2, domain1.ssize() - 1);
      for (auto j1 = min1; j1 <= max1; ++j1) {
        for (auto j0 = min0; j0 <= max0; ++j0) {
          m_mask[{j0, j1}] = true;
        }
      }
      m_x.push_back(std::move(xi));
    }
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  BiCospline(const Domain& domain0, const Domain& domain1, const TRange& x) :
      BiCospline(domain0, domain1, x.begin(), x.end())
  {}

  /**
   * @brief List-based constructor.
   */
  BiCospline(const Domain& domain0, const Domain& domain1, std::initializer_list<Value> x) :
      BiCospline(domain0, domain1, x.begin(), x.end())
  {}

  /**
   * @brief Resample an input raster of knot values.
   */
  template <typename TRaster>
  std::vector<Value> operator()(const TRaster& v)
  {
    for (const auto& p : m_mask) {
      m_splines0[p[1]].set(p[0], v[p]);
    }
    std::vector<Value> y;
    y.reserve(m_x.size());
    const auto max1 = m_spline1.domain().ssize() - 1;
    for (const auto& x : m_x) {
      const auto i1 = x[1].index();
      const auto min = std::max(i1 - 1, 0L);
      const auto max = std::min(i1 + 2, max1);
      for (auto i = min; i <= max; ++i) {
        m_spline1.set(i, m_splines0[i](x[0]));
      }
      y.push_back(m_spline1(x[1]));
    }
    return y;
  }

private:

  std::vector<Method> m_splines0; ///< Splines along axis 0
  Method m_spline1; ///< Spline along axis 1
  std::vector<std::array<Arg, Dimension>> m_x; ///< The arguments
  Linx::Mask<Dimension> m_mask; ///< The neighboring knot abscissae
};

} // namespace Splider

#endif
