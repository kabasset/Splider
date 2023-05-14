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
template <typename T> // FIXME rm?
class BiSplineResampler {

public:
  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  BiSplineResampler(const SplineIntervals& domain0, const SplineIntervals& domain1, TIt begin, TIt end) :
      m_domain0(domain0), m_domain1(domain1), // FIXME useful?
      m_splines0(m_domain1.size(), Spline<T, SplineCache::Lazy>(m_domain0)), m_spline1(m_domain1),
      m_x(std::distance(begin, end)),
      m_mask(Linx::Position<2>::zero(), {m_domain0.size() - 1, m_domain1.size() - 1}, false) {
    auto it = m_x.begin();
    for (; begin < end; ++begin, ++it) {
      *it = {SplineArg(m_domain0, (*begin)[0]), SplineArg(m_domain1, (*begin)[1])};
      const auto i0 = (*it)[0].m_index;
      const auto i1 = (*it)[1].m_index;
      const auto min0 = i0 > 0 ? i0 - 1 : 0;
      const auto max0 = std::min(i0 + 2, m_domain0.size() - 1);
      const auto min1 = i1 > 0 ? i1 - 1 : 0;
      const auto max1 = std::min(i1 + 2, m_domain1.size() - 1);
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
  BiSplineResampler(const SplineIntervals& domain0, const SplineIntervals& domain1, const TRange& x) :
      BiSplineResampler(domain0, domain1, x.begin(), x.end()) {}

  /**
   * @brief List-based constructor.
   */
  BiSplineResampler(const SplineIntervals& domain0, const SplineIntervals& domain1, std::initializer_list<T> x) :
      BiSplineResampler(domain0, domain1, x.begin(), x.end()) {}

  /**
   * @brief Resample an input raster of knot values.
   */
  template <typename TRaster>
  std::vector<typename TRaster::value_type> operator()(const TRaster& v) {
    for (const auto& p : m_mask) {
      m_splines0[p[1]].v(p[0], v[p]);
    }
    std::vector<typename TRaster::value_type> y;
    y.reserve(m_x.size());
    for (const auto& x : m_x) {
      const auto i1 = x[1].m_index;
      const auto min = i1 > 0 ? i1 - 1 : 0;
      const auto max = std::min(i1 + 2, m_domain1.size() - 1);
      for (auto i = min; i <= max; ++i) {
        printf("m_spline1.v(%li, m_splines0[%li](%li) = %f)\n", i, x[0].m_index, i, m_splines0[i](x[0]));
        m_spline1.v(i, m_splines0[i](x[0]));
      }
      printf("m_spline1(%li) = %f\n", x[1].m_index, m_spline1(x[1]));
      y.push_back(m_spline1(x[1]));
    }
    return y;
  }

private:
  const SplineIntervals& m_domain0; ///< The intervals along axis 0
  const SplineIntervals& m_domain1; ///< The intervals along axis 1
  std::vector<Spline<T, SplineCache::Lazy>> m_splines0; ///< Splines along axis 0
  Spline<T, SplineCache::Lazy> m_spline1; ///< Spline along axis 1
  std::vector<std::array<SplineArg, 2>> m_x; ///< The arguments
  Linx::Mask<2> m_mask; ///< The neighboring knot abscissae
};

} // namespace Splider

#endif
