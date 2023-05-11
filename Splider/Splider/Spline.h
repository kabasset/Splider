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
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit SplineIntervals(TIt begin, TIt end) : m_u(std::move(begin), std::move(end)), m_h(m_u.size() - 1) {
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
  explicit SplineIntervals(const TRange& u) : SplineIntervals(u.begin(), u.end()) {}

  /**
   * @brief List-based constructor.
   */
  SplineIntervals(std::initializer_list<double> u) : SplineIntervals(u.begin(), u.end()) {}

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
  /**
   * @brief Null constructor.
   */
  SplineArg() = default;

  /**
   * @brief Constructor.
   */
  explicit SplineArg(const SplineIntervals& domain, double x) {
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
  std::size_t m_index; ///< The interval index
  double m_cv0; ///< The v[i] coefficient
  double m_cv1; ///< The v[i + 1] coefficient
  double m_cs0; ///< The s[i] coefficient
  double m_cs1; ///< The s[i + 1] coefficient
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
  explicit Spline(const SplineIntervals& u) : m_domain(u), m_v(u.m_u.size()), m_s(m_v.size()) {}

  template <typename TIt>
  explicit Spline(const SplineIntervals& u, TIt begin, TIt end) :
      m_domain(u), m_v(std::move(begin), std::move(end)), m_s(m_v.size()) {
    compute();
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Spline(const SplineIntervals& u, const TRange& v) : Spline(u, v.begin(), v.end()) {}

  /**
   * @brief List-based constructor.
   */
  explicit Spline(const SplineIntervals& u, std::initializer_list<T> v) : Spline(u, v.begin(), v.end()) {}

  /**
   * @brief Assign knot values.
   */
  template <typename TKnots>
  void assign(const TKnots& v) { // FIXME it-based, list-based
    m_v.assign(v.begin(), v.end());
    // FIXME check size
    compute();
  }

  /**
   * @brief Evaluate the spline.
   */
  T operator()(double x) const {
    return operator()(SplineArg(m_domain, x));
  }

  /**
   * @brief Evaluate the spline with caching.
   */
  T operator()(const SplineArg& x) const {
    const auto i = x.m_index;
    return m_v[i] * x.m_cv0 + m_v[i + 1] * x.m_cv1 + m_s[i] * x.m_cs0 + m_s[i + 1] * x.m_cs1;
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
  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit SplineResampler(const SplineIntervals& domain, TIt begin, TIt end) : m_domain(domain), m_args(end - begin) {
    std::transform(std::move(begin), std::move(end), m_args.begin(), [&](const auto& e) {
      return SplineArg(m_domain, e);
    });
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit SplineResampler(const SplineIntervals& u, const TRange& x) : SplineResampler(u, x.begin(), x.end()) {}

  /**
   * @brief List-based constructor.
   */
  template <typename T>
  explicit SplineResampler(const SplineIntervals& u, std::initializer_list<T> x) :
      SplineResampler(u, x.begin(), x.end()) {}

  template <typename TRange>
  void assign(const TRange& x) { // FIXME it-based, list-based
    for (const auto& e : x) {
      m_args.emplace_back(m_domain, e);
    }
  }

  template <typename TIt>
  std::vector<typename std::iterator_traits<TIt>::value_type> operator()(TIt begin, TIt end) const {
    using T = typename std::iterator_traits<TIt>::value_type;
    Spline<T> spline(m_domain, std::move(begin), std::move(end));
    std::vector<T> out(m_args.size());
    std::transform(m_args.begin(), m_args.end(), out.begin(), [&](const auto& e) {
      return spline(e);
    });
    return out;
  }

  template <typename TRange>
  std::vector<typename TRange::value_type> operator()(const TRange& v) const {
    return operator()(v.begin(), v.end());
  }

  template <typename T>
  std::vector<T> operator()(std::initializer_list<T> v) const {
    return operator()(v.begin(), v.end());
  }

private:
  const SplineIntervals& m_domain;
  std::vector<SplineArg> m_args;
};

class SplineBuilder {
public:
  template <typename TIt>
  explicit SplineBuilder(TIt begin, TIt end) : m_domain(std::move(begin), std::move(end)) {}

  template <typename TRange>
  explicit SplineBuilder(const TRange& u) : SplineBuilder(u.begin(), u.end()) {}

  template <typename T>
  SplineBuilder(std::initializer_list<T> u) : SplineBuilder(u.begin(), u.end()) {}

  template <typename TIt>
  Spline<typename std::iterator_traits<TIt>::value_type> interpolant(TIt begin, TIt end) const {
    return Spline<typename std::iterator_traits<TIt>::value_type>(m_domain, std::move(begin), std::move(end));
  }

  template <typename TRange>
  Spline<typename TRange::value_type> interpolant(const TRange& v) const {
    return interpolant(v.begin(), v.end());
  }

  template <typename T>
  Spline<T> interpolant(std::initializer_list<T> v) const {
    return interpolant(v.begin(), v.end());
  }

  template <typename TIt>
  SplineResampler resampler(TIt begin, TIt end) const {
    return SplineResampler(m_domain, std::move(begin), std::move(end));
  }

  template <typename TRange>
  SplineResampler resampler(const TRange& x) const {
    return resampler(x.begin(), x.end());
  }

  template <typename T>
  SplineResampler resampler(std::initializer_list<T> x) const {
    return resampler(x.begin(), x.end());
  }

private:
  SplineIntervals m_domain;
};

} // namespace Splider

#endif // _SPLIDER_SPLINE_H
