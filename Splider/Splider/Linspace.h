/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#ifndef _SPLIDER_LINSPACE_H
#define _SPLIDER_LINSPACE_H

#include "Linx/Data/Vector.h" // Index
#include "Splider/Mode.h"

#include <stdexcept>
#include <vector>

namespace Splider {

/**
 * @brief The knot abscissae.
 * @tparam T The real number type
 * 
 * This class stores both the abscissae of the knots and precomputes some spline coefficients to speed-up spline evaluation.
 */
template <typename T = double>
class Linspace {
public:

  /**
   * @brief Declare the domain as even.
   */
  static constexpr bool IsEven = true;

  /**
   * @brief The real number type
   */
  using Value = T;

  /**
   * @brief Iterator-based constructor.
   */
  explicit Linspace(Value front, Value step, Linx::Index size) : m_front(front), m_h(step), m_ssize(size) {}

  /**
   * @brief Get the first abscissa.
   */
  inline Value front() const
  {
    return m_front;
  }

  /**
   * @brief Get the last abscissa.
   */
  inline Value back() const
  {
    return operator[](ssize() - 1);
  }

  /**
   * @brief Get the number of knots.
   */
  inline std::size_t size() const
  {
    return static_cast<std::size_t>(ssize());
  }

  /**
   * @copybrief size()
   */
  inline Linx::Index ssize() const
  {
    return m_ssize;
  }

  /**
   * @brief Get the length of any subinterval.
   */
  inline Value length(Linx::Index) const
  {
    return m_h;
  }

  /**
   * @brief Get the abscissa of the i-th knot.
   */
  inline Value operator[](Linx::Index i) const
  {
    return m_front + i * m_h;
  }

  /**
   * @brief Get the index of the interval which contains a given abscissa.
   */
  Linx::Index index(Value x) const
  {
    return Linx::Index((x - m_front) / m_h);
  }

private:

  Value m_front;
  Value m_h;
  Linx::Index m_ssize;
};

} // namespace Splider

#endif
