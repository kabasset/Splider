/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_SPLINE_H
#define _SPLIDER_SPLINE_H

#include "Splider/Argument.h"
#include "Splider/Partition.h"

#include <stdexcept>
#include <vector>

/**
 * @brief Spline builder.
 */
namespace Splider {

/**
 * @brief Natural cubic spline interpolant.
 * @tparam T The knot value type, which can be any arithmetic type
 * @tparam M The evaluation mode
 * 
 * A spline is parametrized with the list of knot abscissae and values,
 * and is evaluated on a scalar or vector argument.
 * 
 * For repeated use of a spline over a constant set of arguments and varying values, see `Cospline`.
 */
template <typename T, Mode M = Mode::Solve | Mode::Early>
class Spline {

public:
  /**
   * @brief Null knots constructor.
   */
  explicit Spline(const Partition& u) : m_domain(u), m_v(m_domain.size()), m_s(m_domain.size()), m_valid(true) {}

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Spline(const Partition& u, TIt begin, TIt end) :
      m_domain(u), m_v(begin, end), m_s(m_v.size()), m_valid(false) {
    early_update();
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Spline(const Partition& u, const TRange& v) : Spline(u, v.begin(), v.end()) {}

  /**
   * @brief List-based constructor.
   */
  explicit Spline(const Partition& u, std::initializer_list<T> v) : Spline(u, v.begin(), v.end()) {}

  /**
   * @brief Assign knot values from an iterator.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end) {
    m_v.assign(begin, end);
    early_update();
  }

  /**
   * @brief Assign knot values from a range.
   */
  template <typename TRange>
  void assign(const TRange& v) {
    assign(v.begin(), v.end());
  }

  /**
   * @brief Assign knot values from a list.
   */
  void assign(std::initializer_list<T> v) {
    assign(v.begin(), v.end());
  }

  /**
   * @brief Check whether the evaluation mode matches a given mode.
   */
  static constexpr bool mode_matches(Mode mode) {
    if constexpr (M == Mode::Manual && mode == Mode::Manual) {
      return true;
    } else {
      return static_cast<bool>(M & mode);
    }
  }

  /**
   * @brief Check whether the coefficients have been evaluated already.
   */
  bool is_valid(Linx::Index i) const {
    // FIXME check around i
    return m_valid;
  }

  /**
   * @brief Get the i-th knot value.
  */
  inline const T& v(std::size_t i) const {
    return m_v[i];
  }

  /**
   * @brief Set the i-th knot value.
   */
  inline void v(std::size_t i, const T& value) {
    m_v[i] = value;
    m_valid = false;
    early_update();
  }

  /**
   * @brief Get the i-th second derivative.
   */
  inline const T& dv2(std::size_t i) {
    lazy_update(i);
    return m_s[i];
  }

  /**
   * @brief Evaluate the spline.
   */
  inline T operator()(double x) {
    return operator()(SplineArg(m_domain, x));
  }

  /**
   * @brief Evaluate the spline.
   */
  inline T operator()(const SplineArg& x) {
    const auto i = x.m_index;
    lazy_update(i);
    return m_v[i] * x.m_cv0 + m_v[i + 1] * x.m_cv1 + m_s[i] * x.m_cs0 + m_s[i + 1] * x.m_cs1;
  }

  /**
   * @brief Evaluate the spline over an iterator.
   * @return The vector of interpolated values
   * 
   * The iterator can either point to `double`s or `SplineArg`s.
   */
  template <typename TIt>
  std::vector<T> operator()(TIt begin, TIt end) {
    std::vector<T> out;
    out.reserve(std::distance(begin, end));
    for (; begin != end; ++begin) {
      out.push_back(operator()(*begin));
    }
    return out;
  }

  /**
   * @brief Evaluate the spline over a range.
   * @return The vector of interpolated values
   * 
   * The range can either contain `double`s or `SplineArg`s.
   */
  template <typename TRange>
  std::vector<T> operator()(const TRange& x) {
    return operator()(x.begin(), x.end());
  }

  /**
   * @brief Evaluate the spline over a list.
   * @return The vector of interpolated values
   * 
   * The list can either contain `double`s or `SplineArg`s.
   */
  std::vector<T> operator()(std::initializer_list<T> x) {
    return operator()(x.begin(), x.end());
  }

  /**
   * @brief Solve the tridiagonal system.
   */
  void solve() {
    Linx::Index n = m_s.size();
    std::vector<double> b(n);
    const auto& c = m_domain.m_h;
    std::vector<T> d(n);

    for (Linx::Index i = 1; i < n - 1; ++i) {
      b[i] = 2. * (c[i - 1] + c[i]);
      d[i] = 6. * ((m_v[i + 1] - m_v[i]) / c[i] - (m_v[i] - m_v[i - 1]) / c[i - 1]);
    }

    // Forward
    for (Linx::Index i = 2; i < n - 1; ++i) {
      auto w = c[i - 1] / b[i - 1];
      b[i] -= w * c[i - 1];
      d[i] -= w * d[i - 1];
    }

    // Backward
    m_s[n - 2] = d[n - 2] / b[n - 2];
    for (auto i = n - 3; i > 0; --i) {
      m_s[i] = (d[i] - c[i] * m_s[i + 1]) / b[i];
    }

    // Natutal spline // FIXME useful?
    m_s[0] = 0;
    m_s[n - 1] = 0;

    m_valid = true;
  }

  /**
   * @brief Approximate the second derivatives with finite differences.
   */
  void approximate() { // FIXME take i as input
    for (Linx::Index i = 1; i < m_s.size() - 1; ++i) {
      auto d = (m_v[i + 1] - m_v[i]) / m_domain.m_h[i] - (m_v[i] - m_v[i - 1]) / m_domain.m_h[i - 1];
      m_s[i] = d * 2. / (m_domain.m_h[i] + m_domain.m_h[i - 1]);
    }
    m_valid = true;
  }

private:
  /**
   * @brief Update if early.
   */
  inline void early_update() {
    if constexpr (mode_matches(Mode::Early)) {
      update();
    }
  }

  /**
   * @brief Update if lazy.
   */
  inline void lazy_update(Linx::Index i) {
    if constexpr (mode_matches(Mode::Lazy)) {
      // FIXME update around i only
      update();
    }
  }

  /**
   * @brief Evaluate the internal coefficients.
   */
  inline void update() {
    if constexpr (mode_matches(Mode::Manual)) {
      return;
    } else {
      if (m_valid) {
        return;
      }
      if constexpr (mode_matches(Mode::Solve)) {
        solve();
      } else if constexpr (mode_matches(Mode::Approximate)) {
        approximate();
      }
    }
  }

private:
  const Partition& m_domain; ///< The knots domain
  std::vector<T> m_v; ///< The knot values
  std::vector<T> m_s; ///< The knot second derivatives
  bool m_valid; ///< Validity flags
  // FIXME local validity
};

} // namespace Splider

#endif
