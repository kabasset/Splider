/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_SPLINE_H
#define _SPLIDER_SPLINE_H

#include "Splider/Argument.h"
#include "Splider/Linspace.h"
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
 * @tparam TDomain The knot domain type
 * @tparam M The evaluation mode
 * 
 * A spline is parametrized with the list of knot abscissae and values,
 * and is evaluated on a scalar or vector argument.
 * 
 * For repeated use of a spline over a constant set of arguments and varying values, see `Cospline`.
 */
template <typename T, typename TDomain = Partition<double>, Mode M = Mode::Solve | Mode::Early>
class Spline {

public:
  /**
   * @brief The knot value type.
   */
  using Value = T;

  /**
   * @brief The knot domain type.
   */
  using Domain = TDomain;

  /**
   * @brief The real number type.
   */
  using Real = typename Domain::Value;

  /**
   * @brief The argument type.
   */
  using Arg = SplineArg<Real>;

  /**
   * @brief Null knots constructor.
   */
  explicit Spline(const Domain& u) : m_domain(u), m_v(m_domain.size()), m_s(m_domain.size()), m_valid(true) {}

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Spline(const Domain& u, TIt begin, TIt end) : m_domain(u), m_v(begin, end), m_s(m_v.size()), m_valid(false) {
    early_update();
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Spline(const Domain& u, const TRange& v) : Spline(u, v.begin(), v.end()) {}

  /**
   * @brief List-based constructor.
   */
  template <typename TV>
  explicit Spline(const Domain& u, std::initializer_list<TV> v) : Spline(u, v.begin(), v.end()) {}

  /**
   * @brief Get the knots abscissae.
   */
  const Domain& domain() const {
    return m_domain;
  }

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
  template <typename TV>
  void assign(std::initializer_list<TV> v) {
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
  inline const Value& v(Linx::Index i) const {
    return m_v[i];
  }

  /**
   * @brief Set the i-th knot value.
   */
  inline void v(Linx::Index i, const Value& value) {
    m_v[i] = value;
    m_valid = false;
    early_update();
  }

  /**
   * @brief Get the i-th second derivative.
   */
  inline const Value& dv2(Linx::Index i) {
    lazy_update(i);
    return m_s[i];
  }

  /**
   * @brief Evaluate the spline.
   */
  inline Value operator()(Real x) {
    return operator()(Arg(m_domain, x));
  }

  /**
   * @brief Evaluate the spline.
   */
  inline Value operator()(const Arg& x) {
    const auto i = x.m_index;
    lazy_update(i);
    return m_v[i] * x.m_cv0 + m_v[i + 1] * x.m_cv1 + m_s[i] * x.m_cs0 + m_s[i + 1] * x.m_cs1;
  }

  std::vector<Value> operator()(const Args<Real>& x) {
    lazy_update(0); // FIXME i
    std::vector<Value> out;
    out.reserve(x.size());
    for (const auto& arg : x.m_args) {
      const auto i = arg.m_index;
      out.push_back(m_v[i] * arg.m_cv0 + m_v[i + 1] * arg.m_cv1 + m_s[i] * arg.m_cs0 + m_s[i + 1] * arg.m_cs1);
    }
    return out;
  }

  /**
   * @brief Evaluate the spline over an iterator.
   * @return The vector of interpolated values
   * 
   * The iterator can either point to `Real`s or `Arg`s.
   */
  template <typename TIt>
  inline std::vector<Value> operator()(TIt begin, TIt end) {
    return operator()(Args<Real>(m_domain, begin, end));
  }

  /**
   * @brief Evaluate the spline over a range.
   * @return The vector of interpolated values
   * 
   * The range can either contain `Real`s or `Arg`s.
   */
  template <typename TRange>
  std::vector<Value> operator()(const TRange& x) {
    return operator()(x.begin(), x.end());
  }

  /**
   * @brief Evaluate the spline over a list.
   * @return The vector of interpolated values
   * 
   * The list can either contain `Real`s or `Arg`s.
   */
  template <typename TX>
  std::vector<Value> operator()(std::initializer_list<TX> x) {
    return operator()(x.begin(), x.end());
  }

  /**
   * @brief Solve the tridiagonal system.
   */
  void solve() {
    if constexpr (Domain::IsEven) {
      solve_even();
    } else {
      solve_uneven();
    }
  }

  /**
   * @brief Approximate the second derivatives with finite differences.
   */
  void approximate() { // FIXME take i as input
    for (Linx::Index i = 1; i < m_s.size() - 1; ++i) {
      auto d = (m_v[i + 1] - m_v[i]) * m_domain.m_g[i] - (m_v[i] - m_v[i - 1]) * m_domain.m_g[i - 1];
      m_s[i] = d * 2. / (m_domain.m_h[i] + m_domain.m_h[i - 1]);
    }
    m_valid = true;
  }

private:
  void solve_even() {
    const Linx::Index n = m_s.size();
    const auto h = m_domain.length(0);
    const auto g = 1. / h;
    std::vector<Real> b(n, 4. * h);
    std::vector<Value> d(n);

    for (Linx::Index i = 1; i < n - 1; ++i) {
      d[i] = 6. * (m_v[i + 1] - 2 * m_v[i] + m_v[i - 1]) * g;
    }

    // Forward
    for (Linx::Index i = 2; i < n - 1; ++i) {
      const auto w = h / b[i - 1];
      b[i] -= w * h;
      d[i] -= w * d[i - 1];
    }

    // Backward
    m_s[n - 2] = d[n - 2] / b[n - 2];
    for (auto i = n - 3; i > 0; --i) {
      m_s[i] = (d[i] - h * m_s[i + 1]) / b[i];
    }

    // Natutal spline // FIXME useful?
    m_s[0] = 0;
    m_s[n - 1] = 0;

    m_valid = true;
  }

  void solve_uneven() {
    const Linx::Index n = m_s.size();
    std::vector<Real> b(n);
    std::vector<Value> d(n);

    for (Linx::Index i = 1; i < n - 1; ++i) {
      const auto h0 = m_domain.length(i - 1);
      const auto h1 = m_domain.length(i);
      b[i] = 2. * (h0 + h1);
      d[i] = 6. * ((m_v[i + 1] - m_v[i]) / h1 - (m_v[i] - m_v[i - 1]) / h0);
    }

    // Forward
    for (Linx::Index i = 2; i < n - 1; ++i) {
      const auto h = m_domain.length(i);
      const auto w = h / b[i - 1];
      b[i] -= w * h;
      d[i] -= w * d[i - 1];
    }

    // Backward
    m_s[n - 2] = d[n - 2] / b[n - 2];
    for (auto i = n - 3; i > 0; --i) {
      m_s[i] = (d[i] - m_domain.length(i) * m_s[i + 1]) / b[i];
    }

    // Natutal spline // FIXME useful?
    m_s[0] = 0;
    m_s[n - 1] = 0;

    m_valid = true;
  }

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
  const Domain& m_domain; ///< The knots domain
  std::vector<Value> m_v; ///< The knot values
  std::vector<Value> m_s; ///< The knot second derivatives
  bool m_valid; ///< Validity flags
  // FIXME local validity
};

} // namespace Splider

#endif
