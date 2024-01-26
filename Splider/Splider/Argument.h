/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_ARGUMENT_H
#define _SPLIDER_ARGUMENT_H

#include "Linx/Data/Vector.h" // Index
#include "Splider/Mode.h"
#include "Splider/Partition.h"

namespace Splider {

/**
 * @brief A spline argument.
 * 
 * A spline argument is an absissa value for which a spline should be evaluated.
 * It is bound to spline intervals in order to precompute a few spline coefficients.
 * 
 * Manually instantiating a spline argument is useful when it should be reused.
 * It is not necessary for single-use: in this case, simply call the spline on a mere real value.
 */
template <typename T = double>
class SplineArg {

  template <typename, typename, Mode>
  friend class Spline;

  template <typename, typename>
  friend class BiCospline;

public:
  /**
   * @brief The real number type.
   */
  using Value = T;

  /**
   * @brief Null constructor.
   * 
   * The such-constructed object is ill-formed and should only be used for compatibility, e.g. with `std::vector`.
   */
  SplineArg() = default;

  /**
   * @brief Constructor.
   */
  explicit SplineArg(const Partition<Value>& domain, Value x) {
    m_index = domain.index(x);
    const auto h = domain.length(m_index);
    const auto left = x - domain[m_index];
    const auto right = h - left;
    m_cv0 = right / h;
    m_cv1 = left / h;
    m_cs0 = right / 6. * (right * m_cv0 - h);
    m_cs1 = left / 6. * (left * m_cv1 - h);
  }

private:
  Linx::Index m_index; ///< The interval index
  Value m_cv0; ///< The v[i] coefficient
  Value m_cv1; ///< The v[i + 1] coefficient
  Value m_cs0; ///< The s[i] coefficient
  Value m_cs1; ///< The s[i + 1] coefficient
};

} // namespace Splider

#endif
