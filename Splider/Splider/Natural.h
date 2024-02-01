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
class Natural : public BuilderMixin<Natural> {

public:
  /**
   * @brief The argument type.
   */
  template <typename TDomain>
  using Arg = C2Arg<TDomain>;

  /**
   * @brief The spline evaluator.
   */
  template <typename TDomain, typename T>
  class Spline : public C2SplineMixin<TDomain, T, Spline<TDomain, T>> {
    using Mixin = C2SplineMixin<TDomain, T, Spline>;

  public:
    template <typename... TParams>
    Spline(TParams&&... params) : Mixin(LINX_FORWARD(params)...) {}

    /**
     * @brief Solve the tridiagonal system using Thomas algorithm.
     */
    void update(Linx::Index) {
      if (Mixin::m_valid) {
        return;
      }

      const Linx::Index n = this->m_s.size();
      std::vector<typename Mixin::Real> diag(n);
      std::vector<typename Mixin::Value> rhs(n);

      // Initialize i = 1 for merging initialization and forward pass
      auto h0 = this->m_domain.length(0);
      auto h1 = this->m_domain.length(1);
      auto dv0 = (this->m_v[1] - this->m_v[0]) / h0;
      auto dv1 = (this->m_v[2] - this->m_v[1]) / h1;
      diag[1] = 2. * (h0 + h1);
      rhs[1] = 6. * (dv1 - dv0);

      // Initialization and forward pass
      for (Linx::Index i = 2; i < n - 1; ++i) {
        h0 = h1;
        h1 = this->m_domain.length(i);
        dv0 = dv1;
        dv1 = (this->m_v[i + 1] - this->m_v[i]) / h1;
        const auto w = h1 / diag[i - 1];
        diag[i] = 2. * (h0 + h1) - w * h1;
        rhs[i] = 6. * (dv1 - dv0) - w * rhs[i - 1];
      }

      // Backward pass
      this->m_s[n - 2] = rhs[n - 2] / diag[n - 2];
      for (auto i = n - 3; i > 0; --i) {
        this->m_s[i] = (rhs[i] - this->m_domain.length(i) * this->m_s[i + 1]) / diag[i];
      }

      // Natutal spline // FIXME useful?
      this->m_s[0] = 0;
      this->m_s[n - 1] = 0;

      this->m_valid = true;
    }
  };
};

namespace FiniteDiff {
class Natural : public BuilderMixin<FiniteDiff::Natural> {};
} // namespace FiniteDiff

} // namespace Splider

#endif
