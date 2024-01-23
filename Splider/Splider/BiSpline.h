/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_BISPLINE_H
#define _SPLIDER_BISPLINE_H

#include "Linx/Data/Mask.h"
#include "Linx/Data/Raster.h"
#include "Linx/Data/Sequence.h"
#include "Splider/Spline.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace Splider {

/**
 * @brief Alias for a sequence of 2D position.
 * @see `BiCospline`
 */
template <Linx::Index N>
using Trajectory = Linx::Sequence<Linx::Vector<double, N>>; // FIXME Raster?

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
template <typename T> // FIXME possible to rm T?
class BiCospline {

public:
  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  BiCospline(const Partition& domain0, const Partition& domain1, TIt begin, TIt end) :
      m_domain0(domain0), m_domain1(domain1), // FIXME useful?
      m_splines0(m_domain1.size(), Spline<T>(m_domain0)), m_spline1(m_domain1), m_x(std::distance(begin, end)),
      m_mask(Linx::Position<2>::zero(), {m_domain0.size() - 1, m_domain1.size() - 1}, false) {
    for (auto it = m_x.begin(); begin != end; ++begin, ++it) {
      *it = {SplineArg(m_domain0, (*begin)[0]), SplineArg(m_domain1, (*begin)[1])};
      const auto i0 = (*it)[0].m_index;
      const auto i1 = (*it)[1].m_index;
      const auto min0 = std::max(i0 - 1, 0L);
      const auto max0 = std::min(i0 + 2, static_cast<Linx::Index>(m_domain0.size()) - 1);
      const auto min1 = std::max(i1 - 1, 0L);
      const auto max1 = std::min(i1 + 2, static_cast<Linx::Index>(m_domain1.size()) - 1);
      for (auto j1 = min1; j1 <= max1; ++j1) {
        for (auto j0 = min0; j0 <= max0; ++j0) {
          m_mask[{j0, j1}] = true;
        }
      }
    }
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  BiCospline(const Partition& domain0, const Partition& domain1, const TRange& x) :
      BiCospline(domain0, domain1, x.begin(), x.end()) {}

  /**
   * @brief List-based constructor.
   */
  BiCospline(const Partition& domain0, const Partition& domain1, std::initializer_list<T> x) :
      BiCospline(domain0, domain1, x.begin(), x.end()) {}

  /**
   * @brief Resample an input raster of knot values.
   */
  template <typename TRaster>
  auto operator()(const TRaster& v) {
    for (const auto& p : m_mask) {
      m_splines0[p[1]].v(p[0], v[p]);
    }
    std::vector<std::decay_t<typename TRaster::value_type>> y;
    y.reserve(m_x.size());
    for (const auto& x : m_x) {
      const auto i1 = x[1].m_index;
      const auto min = std::max(i1 - 1, 0L);
      const auto max = std::min(i1 + 2, static_cast<Linx::Index>(m_domain1.size()) - 1);
      for (auto i = min; i <= max; ++i) {
        m_spline1.v(i, m_splines0[i](x[0]));
      }
      y.push_back(m_spline1(x[1]));
    }
    return y;
  }

private:
  const Partition& m_domain0; ///< The intervals along axis 0
  const Partition& m_domain1; ///< The intervals along axis 1
  std::vector<Spline<T>> m_splines0; ///< Splines along axis 0
  Spline<T> m_spline1; ///< Spline along axis 1
  std::vector<std::array<SplineArg, 2>> m_x; ///< The arguments
  Linx::Mask<2> m_mask; ///< The neighboring knot abscissae
};

} // namespace Splider

#endif
