/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_PARTITION_H
#define _SPLIDER_PARTITION_H

#include "Linx/Data/Vector.h" // Index
#include "Splider/Mode.h"

#include <stdexcept>
#include <vector>

namespace Splider {

/**
 * @brief The knot abscissae.
 * @tparam TReal The real number type
 * 
 * This class stores both the abscissae of the knots and precomputes some spline coefficients to speed-up spline evaluation.
 */
template <typename TReal = double>
class Partition {
public:

  /**
   * @brief Declare the domain as uneven, a priori.
   */
  static constexpr bool IsEven = false;

  /**
   * @brief The real number type
   */
  using Value = TReal;

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Partition(TIt begin, TIt end) : m_u(begin, end), m_h(m_u.size() - 1)
  {
    const auto size = m_h.size();
    if (size < 2) {
      throw std::runtime_error("Not enough knots (<3).");
    }
    Value h;
    for (std::size_t i = 0; i < size; ++i) {
      h = m_u[i + 1] - m_u[i];
      if (h <= 0) {
        throw std::runtime_error("Interval bounds are not strictly increasing.");
      }
      m_h[i] = h;
    }
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Partition(const TRange& u) : Partition(u.begin(), u.end())
  {}

  /**
   * @brief List-based constructor.
   */
  Partition(std::initializer_list<Value> u) : Partition(u.begin(), u.end()) {}

  /**
   * @brief Get the number of knots.
   */
  inline std::size_t size() const
  {
    return m_u.size();
  }

  /**
   * @copybrief size()
   */
  inline Linx::Index ssize() const
  {
    return static_cast<Linx::Index>(size());
  }

  /**
   * @brief Get the abscissa of the i-th knot.
   */
  inline Value operator[](Linx::Index i) const
  {
    return m_u[i];
  }

  /**
   * @brief Get the length of the i-th subinterval.
   */
  inline Value length(Linx::Index i) const
  {
    return m_h[i];
  }

  /**
   * @brief Get the index of the interval which contains a given abscissa.
   */
  Linx::Index index(Value x) const
  {
    if (x < m_u[0]) {
      throw std::runtime_error("x is too small!");
    }
    if (x > m_u[m_u.size() - 1]) {
      throw std::runtime_error("x is too large!");
    }
    auto i = ssize() - 2;
    while (x < m_u[i]) {
      --i;
    }
    return i;
  }

private:

  std::vector<Value> m_u; ///< The knot positions
  std::vector<Value> m_h; ///< The knot spacings
};

} // namespace Splider

#endif
