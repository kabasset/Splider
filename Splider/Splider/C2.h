/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_C2_H
#define _SPLIDER_C2_H

#include "Splider/mixins/Cubic.h"

namespace Splider {

/**
 * @brief Natural cubic C2 spline.
 * 
 * This is the only C2 cubic spline, which comes at the cost of more computing.
 */
class Natural : public C2Mixin<Natural> {

public:
  /**
   * @brief Solve the tridiagonal system using Thomas algorithm.
   */
  template <typename TSpline>
  static void update(TSpline& params, Linx::Index) {
    if (params.m_valid) {
      return;
    }

    using Real = typename TSpline::Real;
    using Value = typename TSpline::Value;
    const Linx::Index n = params.m_s.size();
    std::vector<Real> b(n);
    std::vector<Value> d(n);

    // Initialize i = 1 for merging initialization and forward pass
    auto h0 = params.m_domain.length(0);
    auto h1 = params.m_domain.length(1);
    b[1] = 2. * (h0 + h1);
    d[1] = 6. * ((params.m_v[2] - params.m_v[1]) / h1 - (params.m_v[1] - params.m_v[0]) / h0);

    for (Linx::Index i = 2; i < n - 1; ++i, h0 = h1) {
      // Initialize
      h1 = params.m_domain.length(i);
      b[i] = 2. * (h0 + h1);
      d[i] = 6. * ((params.m_v[i + 1] - params.m_v[i]) / h1 - (params.m_v[i] - params.m_v[i - 1]) / h0);

      // Forward
      const auto w = h1 / b[i - 1];
      b[i] -= w * h1;
      d[i] -= w * d[i - 1];
    }

    // Backward
    params.m_s[n - 2] = d[n - 2] / b[n - 2];
    for (auto i = n - 3; i > 0; --i) {
      params.m_s[i] = (d[i] - params.m_domain.length(i) * params.m_s[i + 1]) / b[i];
    }

    // Natutal spline // FIXME useful?
    params.m_s[0] = 0;
    params.m_s[n - 1] = 0;

    params.m_valid = true;
  }
};

namespace FiniteDiff {
class Natural : public C2Mixin<FiniteDiff::Natural> {};
} // namespace FiniteDiff

} // namespace Splider

#endif
