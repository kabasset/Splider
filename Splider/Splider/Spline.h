/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_SPLINE_H
#define _SPLIDER_SPLINE_H

#include "LinxCore/Raster.h"
#include "LinxCore/Tiling.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

/**
 * @brief Spline builder.
 */
namespace Splider {

/**
 * @brief The knot abscissae.
 */
class SplineIntervals {

  friend class SplineArguments;

  template <typename>
  friend class Spline;

public:
  /**
   * @brief Range-based constructor.
   */
  template <typename TIt>
  SplineIntervals(TIt begin, TIt end) : m_u(std::move(begin), std::move(end)), m_h(m_u.size() - 1) {
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
   * @brief Initializer list-based constructor.
  */
  SplineIntervals(const std::initializer_list<double>& u) : SplineIntervals(u.begin(), u.end()) {}

  /**
   * @brief Get the number of knots.
   */
  std::size_t size() const {
    return m_u.size();
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

/**
 * @brief A collection of abscissae to be passed as arguments to the spline.
 */
class SplineArguments {

  template <typename T>
  friend class Spline;

public:
  /**
   * @brief Empty collection constructor.
   */
  SplineArguments(const SplineIntervals& domain) : m_domain(domain), m_x(), m_coefficients() {}

  /**
   * @brief Range-based constructor.
   */
  template <typename TIt>
  SplineArguments(const SplineIntervals& domain, TIt begin, TIt end) :
      m_domain(domain), m_x(std::move(begin), std::move(end)), m_coefficients(m_x.size()) {
    compute();
  }

  /**
   * @brief Initializer list-based constructor.
  */
  SplineArguments(const SplineIntervals& domain, const std::initializer_list<double>& x) :
      SplineArguments(domain, x.begin(), x.end()) {}

  /**
   * @brief Get the number of arguments.
   */
  std::size_t size() const {
    return m_x.size();
  }

  /**
   * @brief Assign new abscissae.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end) {
    m_x.assign(std::move(begin), std::move(end));
    m_coefficients.resize(m_x.size());
    compute();
  }

private:
  /**
   * @brief Update coefficients given m_x.
   */
  void compute() {
    std::size_t index;
    double x, h, left, right;
    for (std::size_t i = 0; i < m_x.size(); ++i) {
      x = m_x[i];
      index = m_domain.index(x);
      h = m_domain.m_h[index];
      left = x - m_domain[index];
      right = h - left;
      m_coefficients[i] = Coefficients {
          index,
          right / h, // cv0
          left / h, // cv1
          right * right * right / (6. * h) - h * right / 6., // cs0
          left * left * left / (6. * h) - h * left / 6. // cs1
      };
    }
  }

  /**
   * @brief Spline coefficients related to x.
   */
  struct Coefficients {
    std::size_t index; ///< Index of the interval of x
    double cv0; ///< Coefficient of v[index]
    double cv1; ///< Coefficient of v[index + 1]
    double cs0; ///< Coefficient of s[index]
    double cs1; ///< Coefficient of s[index + 1]
  };

  const SplineIntervals& m_domain; ///< Knot positions
  std::vector<double> m_x; ///< Arguments
  std::vector<Coefficients> m_coefficients; ///< Spline coefficients
};

/**
 * @brief Cubic spline.
 * 
 * A spline is parametrized with the list of knot abscissae and values,
 * and is evaluated on a scalar or vector argument.
 */
template <typename T>
class Spline {

public:
  /**
   * @brief Null knots constructor.
  */
  Spline(const SplineIntervals& u) : m_domain(u), m_v(u.m_u.size()), m_s(m_v.size()) {}

  /**
   * @brief Valued knots constructors.
  */
  template <typename TKnots>
  Spline(const SplineIntervals& u, const TKnots& v) : m_domain(u), m_v(v.begin(), v.end()), m_s(m_v.size()) {
    compute();
  }

  /**
   * @brief Assign knot values.
   */
  template <typename TKnots>
  void assign(const TKnots& v) {
    m_v.assign(v.begin(), v.end());
    // FIXME check size
    compute();
  }

  /**
   * @brief Evaluate the spline on a scalar value.
   */
  T operator()(double x) const {
    const auto i = m_domain.index(x);
    const double h = m_domain.m_h[i];
    const double left = x - m_domain.m_u[i];
    const double right = h - left;
    const double cv0 = right / h;
    const double cv1 = left / h;
    const double cs0 = right * right * right / (6. * h) - h * right / 6.;
    const double cs1 = left * left * left / (6. * h) - h * left / 6.;
    return m_v[i] * cv0 + m_v[i + 1] * cv1 + m_s[i] * cs0 + m_s[i + 1] * cs1;
  }

  /**
   * @brief Evaluate the spline on a vector value.
   */
  std::vector<T> operator()(const SplineArguments& x) const {
    std::vector<T> out(x.size());
    std::transform(x.m_coefficients.begin(), x.m_coefficients.end(), out.begin(), [&](const auto& c) {
      const auto i = c.index;
      return m_v[i] * c.cv0 + m_v[i + 1] * c.cv1 + m_s[i] * c.cs0 + m_s[i + 1] * c.cs1;
    });
    return out;
  }

private:
  /**
   * @brief Update coefficients given m_v.
   */
  void compute() {
    const auto size = m_v.size();
    // FIXME check size == doman.m_u.size()
    T d0 = (m_v[1] - m_v[0]) / m_domain.m_h[0]; // Because next loop starts at 1
    T d1;
    const auto* vIt = m_v.data() + 1;
    const auto* hIt = m_domain.m_h.data() + 1;
    for (auto sIt = &m_s[0] + 1; sIt != &m_s[size - 1]; ++sIt, ++vIt, ++hIt) {
      // d[i] = (m_v[i + 1] - m_v[i]) / m_h[i];
      d1 = (*(vIt + 1) - *vIt) / *hIt;
      // m_s[i] = (m_v[i + 1] - m_v[i]) / m_h[i] - (m_v[i] - m_v[i - 1]) / m_h[i - 1];
      *sIt = d1 - d0;
      d0 = d1;
    }
    // m_s[0] and m_s[size - 1] are left at 0 for natural splines
  }

  const SplineIntervals& m_domain; ///< The knots domain
  std::vector<T> m_v; ///< The knot values
  std::vector<T> m_s; ///< The knot second derivatives
};

} // namespace Splider

#endif // _SPLIDER_SPLINE_H
