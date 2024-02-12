/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#ifndef _SPLIDER_HERMITE_H
#define _SPLIDER_HERMITE_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Splider/mixins/Hermite.h"

namespace Splider {

/**
 * @brief Cubic Hermite splines.
 */
struct Hermite {
  struct FiniteDiff;
  struct Akima;
  struct CatmullRom;
};

/**
 * @brief The finite difference Hermite splines boundary conditions.
 */
enum class FiniteDiffHermiteBounds {
  OneSided = 0, ///< One-sided finite difference
};

/**
 * @brief The Hermite spline evaluator.
 */
template <typename TDomain, typename TValue, FiniteDiffHermiteBounds B>
class FiniteDiffHermiteSpline :
    public HermiteSplineMixin<TDomain, TValue, FiniteDiffHermiteSpline<TDomain, TValue, B>> {
  using Mixin = HermiteSplineMixin<TDomain, TValue, FiniteDiffHermiteSpline>;

public:

  /**
   * @brief Constructor.
   */
  template <typename... TParams>
  FiniteDiffHermiteSpline(TParams&&... params) : Mixin(LINX_FORWARD(params)...)
  {}

  /**
   * @brief Update the derivatives.
   */
  void update(Linx::Index) // FIXME use index
  {
    if (Mixin::m_valid) {
      return;
    }

    const Linx::Index n = this->m_d.size();

    this->m_d[0] = (this->m_v[1] - this->m_v[0]) / this->m_domain.length(0);

    for (Linx::Index i = 1; i < n - 1; ++i) {
      auto h0 = this->m_domain.length(i - 1);
      auto h1 = this->m_domain.length(i);
      auto d0 = (this->m_v[i] - this->m_v[i - 1]) / h0;
      auto d1 = (this->m_v[i + 1] - this->m_v[i]) / h1;
      this->m_d[i] = (d1 + d0) * 0.5;
      // FIXME optimize
    }

    this->m_d[n - 1] = (this->m_v[n - 1] - this->m_v[n - 1]) / this->m_domain.length(n - 2);

    this->m_valid = true;
  }
};

/**
 * @ingroup builders
 * @brief Cubic Hermite spline with finite difference approximation of the derivatives.
 */
struct Hermite::FiniteDiff : BuilderMixin<Hermite::FiniteDiff, FiniteDiffHermiteBounds> {
  /**
   * @brief The boundary conditions.
   */
  using Bounds = FiniteDiffHermiteBounds;

  /**
   * @brief The knots domain type.
   */
  template <typename TReal>
  using Domain = Partition<TReal>; // FIXME HermiteDomain

  /**
   * @brief The argument type.
   */
  template <typename TDomain>
  using Arg = HermiteArg<TDomain>;

  /**
   * @brief The spline evaluator.
   */
  template <typename TDomain, typename TValue, FiniteDiffHermiteBounds B>
  using Spline = FiniteDiffHermiteSpline<TDomain, TValue, B>;
};

} // namespace Splider

#endif
