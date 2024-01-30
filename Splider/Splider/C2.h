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
class Natural : public C2Mixin<Natural> {};

template <typename TDomain, typename T>
void update(Natural::Spline<TDomain, T>& params, Linx::Index) {
  using Real = typename TDomain::Value;
  using Value = T;

  const Linx::Index n = params.m_s.size();
  std::vector<Real> b(n);
  std::vector<Value> d(n);

  for (Linx::Index i = 1; i < n - 1; ++i) {
    const auto h0 = params.m_domain.length(i - 1);
    const auto h1 = params.m_domain.length(i);
    b[i] = 2. * (h0 + h1);
    d[i] = 6. * ((params.m_v[i + 1] - params.m_v[i]) / h1 - (params.m_v[i] - params.m_v[i - 1]) / h0);
  }

  // Forward
  for (Linx::Index i = 2; i < n - 1; ++i) {
    const auto h = params.m_domain.length(i);
    const auto w = h / b[i - 1];
    b[i] -= w * h;
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

namespace FiniteDiff {
class Natural : public C2Mixin<FiniteDiff::Natural> {};
} // namespace FiniteDiff

} // namespace Splider

#endif
