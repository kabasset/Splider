/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MULTISPLINE_H
#define _SPLIDER_MULTISPLINE_H

#include "LinxCore/Raster.h"

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

  template <typename TDomain> // Iterable over interval refs
  BiSplineResampler(const TDomain& domain) :
      m_domain(domain.begin(), domain.end()), m_splines0(m_domain[1].size()), // FIXME Raster<Spline, Dimension-1>
      m_spline1() {}

  template <typename TPatch>
  void set(const TPatch& knots) {
    for (const auto& p : knots.domain()) {
      for (Linx::Index i = 0; i < Dimension; ++i) {
        auto q = p;
        ++q[i];
        m_splines0[q[1]].set(q[0], knots.parent()[q]);
      }
    }
  }

  T operator()(const Linx::Vector<SplineArg, Dimension>& x) {
    const auto i0 = m_domain0.index(x[0]);
    const auto i1 = m_domain1.index(x[1]);
    m_spline1.set(i1 - 1, m_splines0[i1 - 1](x[0]));
    m_spline1.set(i1 + 0, m_splines0[i1 + 0](x[0]));
    m_spline1.set(i1 + 1, m_splines0[i1 + 1](x[0]));
    m_spline1.set(i1 + 2, m_splines0[i1 + 2](x[0]));
    // FIXME out of bounds
    return m_spline1(x[1]);
  }

private:
  Linx::Vector<const SplineIntervals&, Dimension> m_domain; ///< Intervals
  std::vector<Spline<T, SplineCache::Lazy>> m_splines0; ///< Splines along axis 0
  Spline<T, SplineCache::Lazy> m_spline1; ///< Spline along axis 1
};

} // namespace Splider

#endif // _SPLIDER_SPLINE_H
