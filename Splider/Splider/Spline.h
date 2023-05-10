/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_SPLINE_H
#define _SPLIDER_SPLINE_H
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

  template <typename T>
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

/**
 * @brief A spline argument.
 */
class SplineArg {

  template <typename T>
  friend class Spline;

public:
  SplineArg(const SplineIntervals& domain, double x) {
    m_index = domain.index(x);
    const auto h = domain.length(m_index);
    const auto left = x - domain[m_index];
    const auto right = h - left;
    m_cv0 = right / h;
    m_cv1 = left / h;
    m_cs0 = right * right * right / (6. * h) - h * right / 6.;
    m_cs1 = left * left * left / (6. * h) - h * left / 6.;
  }

private:
  std::size_t m_index;
  double m_cv0;
  double m_cv1;
  double m_cs0;
  double m_cs1;
};

template <typename TArgs>
std::vector<SplineArg> makeSplineArgs(const SplineIntervals& domain, const TArgs& args) {
  std::vector<SplineArg> out;
  out.reserve(args.size());
  for (const auto& a : args) {
    out.push_back({domain, a});
  }
  return out;
}

template <typename T>
std::vector<SplineArg> makeSplineArgs(const SplineIntervals& domain, const std::initializer_list<T>& args) {
  return makeSplineArgs<std::initializer_list<T>>(domain, args);
}

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
    return operator()(SplineArg(m_domain, x));
  }

  /**
   * @brief Evaluate the spline on a scalar value.
   */
  T operator()(const SplineArg& x) const {
    const auto i = x.m_index;
    return m_v[i] * x.m_cv0 + m_v[i + 1] * x.m_cv1 + m_s[i] * x.m_cs0 + m_s[i + 1] * x.m_cs1;
  }

  /**
   * @brief Evaluate the spline on a vector value.
   */
  template <typename TArgs>
  std::vector<T> operator()(const TArgs& x) const {
    std::vector<T> out(x.size());
    std::transform(x.begin(), x.end(), out.begin(), [&](const auto& e) {
      const auto i = e.m_index;
      return m_v[i] * e.m_cv0 + m_v[i + 1] * e.m_cv1 + m_s[i] * e.m_cs0 + m_s[i + 1] * e.m_cs1;
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

class SplineResampler {

public:
  SplineResampler(const SplineIntervals& domain) : m_domain(domain), m_args() {}

  template <typename TArgs>
  void assign(const TArgs& x) {
    for (const auto& e : x) {
      m_args.emplace_back(m_domain, e);
    }
  }

  template <typename TKnots>
  std::vector<typename TKnots::value_type> operator()(const TKnots& v) {
    Spline<typename TKnots::value_type> spline(m_domain, v);
    return spline(m_args);
  }

private:
  const SplineIntervals& m_domain;
  std::vector<SplineArg> m_args;
};

} // namespace Splider

#endif // _SPLIDER_SPLINE_H
