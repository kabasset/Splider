/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_CATMULLROM_H
#define _SPLIDER_CATMULLROM_H

#include "Linx/Base/SeqUtils.h" // IsRange
#include "Splider/Hermite.h"

namespace Splider {

/**
 * @brief Cubic Catmull-Rom splines.
 */
struct Hermite::CatmullRom {
  struct Uniform;
};

/**
 * @brief The Catmull-Rom splines boundary conditions.
 */
enum class CatmullRomBounds {
  OneSided = 0, ///< One-sided finite difference
};

/**
 * @brief The Catmull-Rom spline evaluator.
 */
template <typename TDomain, typename TValue, CatmullRomBounds B>
class UniformCatmullRomSpline :
    public HermiteSplineMixin<TDomain, TValue, UniformCatmullRomSpline<TDomain, TValue, B>> {
  using Mixin = HermiteSplineMixin<TDomain, TValue, UniformCatmullRomSpline>;

public:

  /**
   * @brief Constructor.
   */
  template <typename... TParams>
  UniformCatmullRomSpline(TParams&&... params) : Mixin(LINX_FORWARD(params)...)
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
      this->m_d[i] = (this->m_v[i + 1] - this->m_v[i - 1]) / (h0 + h1);
      // FIXME optimize
    }

    this->m_d[n - 1] = (this->m_v[n - 1] - this->m_v[n - 1]) / this->m_domain.length(n - 2);

    this->m_valid = true;
  }
};

/**
 * @ingroup builders
 * @brief Catmull-Rom spline with uniform parametrization.
 */
struct Hermite::CatmullRom::Uniform : BuilderMixin<Hermite::CatmullRom::Uniform, CatmullRomBounds> {
  /**
   * @brief The boundary conditions.
   */
  using Bounds = CatmullRomBounds;

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
  template <typename TDomain, typename TValue, CatmullRomBounds B>
  using Spline = UniformCatmullRomSpline<TDomain, TValue, B>;
};

} // namespace Splider

#endif
