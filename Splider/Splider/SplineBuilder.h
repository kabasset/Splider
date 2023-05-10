/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_SPLINEBUILDER_H
#define _SPLIDER_SPLINEBUILDER_H

#include "Splider/Spline.h"

namespace Splider {

class SplineBuilder {
public:
  template <typename... Ts>
  SplineBuilder(Ts&&... u) : m_domain(std::forward<Ts>(u)...) {}

  template <typename T, typename... Ts>
  Spline<T> interpolant(Ts&&... v) {
    return Spline<T>(std::forward<Ts>(v)...);
  }

  template <typename... Ts>
  SplineResampler resampler(Ts&&... x) {
    return SplineResampler(std::forward<Ts>(x)...);
  }

private:
  SplineIntervals m_domain;
};

} // namespace Splider

#endif // _SPLIDER_SPLINEBUILDER_H
