/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_BISPLINE_H
#define _SPLIDER_BISPLINE_H

#include "LinxCore/Mask.h"
#include "LinxCore/Raster.h"
#include "Splider/Spline.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace Splider {

/**
 * @brief Bivariate spline resampler.
 */
template <typename T>
class BiSplineResampler {

public:
  /**
   * @brief The dimension.
   */
  static constexpr Linx::Index Dimension = 2;

  template <typename TIt>
  BiSplineResampler(const SplineIntervals& domain0, const SplineIntervals& domain1, TIt begin, TIt end) :
      m_domain0(domain0), m_domain1(domain1), // FIXME useful?
      m_splines0(m_domain1.size(), Spline<T, SplineCache::Lazy>(m_domain0)), m_spline1(m_domain1),
      m_x(std::distance(begin, end)),
      m_mask(Linx::Position<Dimension>::zero(), {m_domain0.size() - 1, m_domain1.size() - 1}, false) {
    auto it = m_x.begin();
    for (; begin < end; ++begin, ++it) {
      *it = {SplineArg(m_domain0, (*begin)[0]), SplineArg(m_domain1, (*begin)[1])};
      const auto i0 = (*it)[0].m_index;
      const auto i1 = (*it)[1].m_index;
      for (auto j1 = i1 - 1; j1 <= i1 + 2; ++j1) {
        for (auto j0 = i0 - 1; j0 <= i0 + 2; ++j0) {
          // FIXME bounds
          m_mask[{j0, j1}] = true;
        }
      }
    }
  }

  template <typename TKnots>
  std::vector<typename TKnots::value_type> operator()(const TKnots& v) {
    for (const auto& p : m_mask) {
      m_splines0[p[1]].v(p[0], v[p]);
    }
    std::vector<typename TKnots::value_type> y;
    y.reserve(m_x.size());
    for (const auto& x : m_x) {
      const auto i1 = x[1].m_index;
      const auto min = std::max(i1 - 1, std::size_t(0));
      const auto max = std::min(i1 + 2, m_domain1.size() - 1);
      for (auto i = min; i < max; ++i) {
        m_spline1.v(i, m_splines0[i](x[0]));
      }
      y.push_back(m_spline1(x[1]));
    }
    return y;
  }

private:
  const SplineIntervals& m_domain0; ///< The intervals along axis 0
  const SplineIntervals& m_domain1; ///< The intervals along axis 1
  std::vector<Spline<T, SplineCache::Lazy>> m_splines0; ///< Splines along axis 0
  Spline<T, SplineCache::Lazy> m_spline1; ///< Spline along axis 1
  std::vector<std::array<SplineArg, Dimension>> m_x; ///< The arguments
  Linx::Mask<Dimension> m_mask; ///< The neighboring knot abscissae
};

} // namespace Splider

#endif // _SPLIDER_SPLINE_H
