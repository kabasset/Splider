/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_C2_H
#define _SPLIDER_C2_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Splider/Partition.h" // FIXME rm
#include "Splider/mixins/C2.h"

namespace Splider {

/**
 * @brief The \f$C^2\f$ splines boundary conditions.
 */
enum class C2Bounds {
  Natural = 0, ///< Null second derivatives at bounds
  NotAKnot ///< Neighboring subinterval fitting
};

/**
 * @brief The \f$C^2\f$ spline evaluator.
 */
template <typename TDomain, typename TValue, C2Bounds B>
class C2Spline : public C2SplineMixin<TDomain, TValue, C2Spline<TDomain, TValue, B>> {
  using Mixin = C2SplineMixin<TDomain, TValue, C2Spline>;

public:

  /**
   * @brief Constructor.
   */
  template <typename... TParams>
  C2Spline(TParams&&... params) : Mixin(LINX_FORWARD(params)...)
  {}

  /**
   * @brief Solve the tridiagonal system using Thomas algorithm.
   */
  void update(Linx::Index)
  {
    if (Mixin::m_valid) {
      return;
    }

    const Linx::Index n = this->m_6s.size();
    std::vector<typename Mixin::Real> diag(n);
    std::vector<typename Mixin::Value> rhs(n);

    // Initialize i = 1 for merging initialization and forward pass
    auto h0 = this->m_domain.length(0);
    auto h1 = this->m_domain.length(1);
    auto dv0 = (this->m_v[1] - this->m_v[0]) / h0;
    auto dv1 = (this->m_v[2] - this->m_v[1]) / h1;
    diag[1] = 2. * (h0 + h1);
    rhs[1] = dv1 - dv0;

    // Initialization and forward pass
    for (Linx::Index i = 2; i < n - 1; ++i) {
      h0 = h1;
      h1 = this->m_domain.length(i);
      dv0 = dv1;
      dv1 = (this->m_v[i + 1] - this->m_v[i]) / h1;
      const auto w = h1 / diag[i - 1];
      diag[i] = 2. * (h0 + h1) - w * h1;
      rhs[i] = dv1 - dv0 - w * rhs[i - 1];
    }

    this->m_6s[n - 1] = 0;

    // Backward pass
    this->m_6s[n - 2] = rhs[n - 2] / diag[n - 2];
    for (auto i = n - 3; i > 0; --i) {
      this->m_6s[i] = (rhs[i] - this->m_domain.length(i) * this->m_6s[i + 1]) / diag[i];
    }

    this->m_6s[0] = 0;

    this->m_valid = true;
  }
};

/**
 * @ingroup builders
 * @brief \f$C^2\f$ cubic spline.
 */
struct C2 : BuilderMixin<C2, C2Bounds> {
  /**
   * @brief The boundary conditions.
   */
  using Bounds = C2Bounds;

  /**
   * @brief The knots domain type.
   */
  template <typename TReal>
  using Domain = Partition<TReal>; // FIXME C2Domain

  /**
   * @brief The argument type.
   */
  template <typename TDomain>
  using Arg = C2Arg<TDomain>;

  /**
   * @brief The spline evaluator.
   */
  template <typename TDomain, typename TValue, C2Bounds B>
  using Spline = C2Spline<TDomain, TValue, B>;

  struct FiniteDiff;
};

/**
 * @brief The \f$C^2\f$ spline evaluator.
 */
template <typename TDomain, typename TValue, C2Bounds B>
class FiniteDiffC2Spline : public C2SplineMixin<TDomain, TValue, FiniteDiffC2Spline<TDomain, TValue, B>> {
  using Mixin = C2SplineMixin<TDomain, TValue, FiniteDiffC2Spline>;

public:

  /**
   * @brief Constructor.
   */
  template <typename... TParams>
  FiniteDiffC2Spline(TParams&&... params) : Mixin(LINX_FORWARD(params)...)
  {}

  /**
   * @brief Update the second derivatives.
   */
  void update(Linx::Index) // FIXME use index
  {
    if (Mixin::m_valid) {
      return;
    }

    const Linx::Index n = this->m_6s.size();

    this->m_6s[0] = 0;

    for (Linx::Index i = 1; i < n - 1; ++i) {
      auto h0 = this->m_domain.length(i - 1);
      auto h1 = this->m_domain.length(i);
      auto d0 = (this->m_v[i] - this->m_v[i - 1]) / h0;
      auto d1 = (this->m_v[i + 1] - this->m_v[i]) / h1;
      this->m_6s[i] = (d1 - d0) * 2. / (h1 + h0);
      // FIXME optimize
    }

    this->m_6s[n - 1] = 0;

    this->m_valid = true;
  }
};

/**
 * @ingroup builders
 * @brief \f$C^2\f$ cubic spline with finite difference approximation of the second derivatives.
 * 
 * This is an approximation of `C2` which enables local evaluation of the coefficients,
 * in lieu of the global tridiagonal system solving of the latter.
 */
struct C2::FiniteDiff : BuilderMixin<C2::FiniteDiff, C2Bounds> {
  /**
   * @brief The boundary conditions.
   */
  using Bounds = C2Bounds;

  /**
   * @brief The knots domain type.
   */
  template <typename TReal>
  using Domain = Partition<TReal>; // FIXME C2Domain

  /**
   * @brief The argument type.
   */
  template <typename TDomain>
  using Arg = C2Arg<TDomain>;

  /**
   * @brief The spline evaluator.
   */
  template <typename TDomain, typename TValue, C2Bounds B>
  using Spline = FiniteDiffC2Spline<TDomain, TValue, B>; // FIXME
};

} // namespace Splider

#endif
