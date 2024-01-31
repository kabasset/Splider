/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_NATURAL_H
#define _SPLIDER_NATURAL_H

#include "Splider/mixins/C2.h"

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
    const auto& domain = params.m_domain;
    const auto& v = params.m_v;
    auto& s = params.m_s;
    const Linx::Index n = s.size();
    std::vector<Real> b(n);
    std::vector<Value> d(n);

    // Initialize i = 1 for merging initialization and forward pass
    auto h0 = domain.length(0);
    auto h1 = domain.length(1);
    auto dv0 = (v[1] - v[0]) / h0;
    auto dv1 = (v[2] - v[1]) / h1;
    b[1] = 2. * (h0 + h1);
    d[1] = 6. * (dv1 - dv0);

    // Initialization and forward pass
    for (Linx::Index i = 2; i < n - 1; ++i) {
      h0 = h1;
      h1 = domain.length(i);
      dv0 = dv1;
      dv1 = (v[i + 1] - v[i]) / h1;
      const auto w = h1 / b[i - 1];
      b[i] = 2. * (h0 + h1) - w * h1;
      d[i] = 6. * (dv1 - dv0) - w * d[i - 1];
    }

    // Backward pass
    s[n - 2] = d[n - 2] / b[n - 2];
    for (auto i = n - 3; i > 0; --i) {
      s[i] = (d[i] - domain.length(i) * s[i + 1]) / b[i];
    }

    // Natutal spline // FIXME useful?
    s[0] = 0;
    s[n - 1] = 0;

    params.m_valid = true;
  }
};

namespace FiniteDiff {
class Natural : public C2Mixin<FiniteDiff::Natural> {};
} // namespace FiniteDiff

} // namespace Splider

#endif
