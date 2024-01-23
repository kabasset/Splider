/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_PARTITION_H
#define _SPLIDER_PARTITION_H

#include <stdexcept>
#include <vector>

/**
 * @brief Spline builder.
 */
namespace Splider {

/**
 * @brief The knot abscissae.
 * 
 * This class stores both the abscissae of the knots and precomputes some spline coefficients to speed-up spline evaluation.
 */
class Partition {

  template <typename>
  friend class Spline;

public:
  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Partition(TIt begin, TIt end) : m_u(std::move(begin), std::move(end)), m_h(m_u.size() - 1) {
    const auto size = m_h.size();
    if (size < 2) {
      throw std::runtime_error("Not enough knots (<3).");
    }
    double h;
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
  explicit Partition(const TRange& u) : Partition(u.begin(), u.end()) {}

  /**
   * @brief List-based constructor.
   */
  Partition(std::initializer_list<double> u) : Partition(u.begin(), u.end()) {}

  /**
   * @brief Get the number of knots.
   */
  inline std::size_t size() const {
    return m_u.size();
  }

  /**
   * @brief Get the length of the i-th interval.
   */
  inline double length(std::size_t i) const {
    return m_h[i];
  }

  /**
   * @brief Get the abscissa of the i-th knot.
   */
  inline double operator[](std::size_t i) const {
    return m_u[i];
  }

  /**
   * @brief Get the index of the interval which contains a given abscissa.
   */
  std::size_t index(double x) const {
    if (x < m_u[0]) {
      throw std::runtime_error("x is too small!");
    }
    if (x > m_u[m_u.size() - 1]) {
      throw std::runtime_error("x is too large!");
    }
    auto i = m_u.size() - 2;
    while (x < m_u[i]) {
      --i;
    }
    return i;
  }

private:
  std::vector<double> m_u; ///< The knot positions
  std::vector<double> m_h; ///< The knot spacings
};

} // namespace Splider

#endif
