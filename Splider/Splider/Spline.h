/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_SPLINE_H
#define _SPLIDER_SPLINE_H

#include <algorithm>
#include <set>
#include <stdexcept>
#include <vector>

/**
 * @brief Spline builder.
 */
namespace Splider {

/**
 * @brief Caching strategy.
 * 
 * Splines rely on knot-wise coefficients which need to be evaluated from the knot positions and values.
 * This evaluation is expensive wrt. other operations.
 * The coefficients can therefore be cached to speed-up computation.
 * Several strategies are available:
 * 
 * - When the spline is to be evaluated over many arguments, i.e. most knots are effectively used, it is best to rely on early evaluation.
 * - When knots are used sparsely, i.e. arguments lie in few of the known intervals only, lazy evaluation is faster.
 * 
 * In both of these cases, spline evaluation is guaranteed to be well-defined.
 * However, for more advanced usage, cache can be user-managed.
 * In this case, the user is responsible for updating the spline state before evaluating it... Use at your own risks.
 * 
 * In general, early and lazy options are preferrable, and should be compared through representative tests to select the fastest way.
 */
enum class SplineCache {
  Early, ///< Cache at assignment
  Lazy, ///< Cache at usage
  User ///< Do not cache automatically
};

/**
 * @brief The knot abscissae.
 * 
 * This class stores both the abscissae of the knots and precomputes some spline coefficients to speed-up spline evaluation.
 */
class Partition {

  template <typename, SplineCache>
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
 * 
 * A spline argument is an absissa value for which a spline should be evaluated.
 * It is bound to spline intervals in order to precompute a few spline coefficients.
 * 
 * Manually instantiating a spline argument is useful when it should be reused.
 * It is not necessary for single-use: in this case, simply call the spline on a mere real value.
 */
class SplineArg {

  template <typename, SplineCache>
  friend class Spline;

  template <typename, SplineCache>
  friend class BiCospline;

public:
  /**
   * @brief Null constructor.
   * 
   * The such-constructed object is ill-formed and should only be used for compatibility, e.g. with `std::vector`.
   */
  SplineArg() = default;

  /**
   * @brief Constructor.
   */
  explicit SplineArg(const Partition& domain, double x) {
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
 * @brief Natural cubic spline interpolant.
 * 
 * A spline is parametrized with the list of knot abscissae and values,
 * and is evaluated on a scalar or vector argument.
 * 
 * For repeated use of a spline over a constant set of arguments and varying values, see `Cospline`.
 */
template <typename T, SplineCache Cache = SplineCache::Early>
class Spline {

public:
  /**
   * @brief Null knots constructor.
   */
  explicit Spline(const Partition& u) : m_domain(u), m_v(m_domain.size()), m_s(m_domain.size()), m_cache() {
    if constexpr (Cache != SplineCache::Early) {
      m_cache.insert({0, m_s.size() - 1}); // Natural cubic splines
    }
  }

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Spline(const Partition& u, TIt begin, TIt end) :
      m_domain(u), m_v(std::move(begin), std::move(end)), m_s(m_v.size()), m_cache() {
    if constexpr (Cache == SplineCache::Early) {
      update();
    } else {
      m_cache.insert({0, m_s.size() - 1}); // Natural cubic splines
    }
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
    // FIXME check size
    if constexpr (Cache == SplineCache::Early) {
      update();
    } else {
      m_cache = {0, m_s.size() - 1};
    }
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
    const auto min = i > 1 ? i - 1 : 1;
    const auto max = std::min(i + 1, m_v.size() - 2);
    for (auto j = min; j <= max; ++j) { // Natural cubic spline
      if constexpr (Cache == SplineCache::Early) {
        update(j);
      } else {
        m_cache.erase(j);
      }
    }
  }

  /**
   * @brief Get the i-th second derivative.
   */
  inline const T& dv2(std::size_t i) {
    if constexpr (Cache == SplineCache::Lazy) {
      if (not valid(i)) {
        update(i);
      }
    }
    return m_s[i];
  }

  /**
   * @brief Check whether the internal coefficients of the i-th knot are valid in cache.
   * 
   * This is always true for early caching and is the method is mostly useful for user-triggered caching.
   */
  bool valid(std::size_t i) const {
    if constexpr (Cache == SplineCache::Early) {
      return true;
    } else {
      return m_cache.find(i) != m_cache.end();
    }
  }

  /**
   * @brief Update all the internal coefficients.
   */
  void update() {
    const auto size = m_v.size();
    // FIXME check size == doman.m_u.size()
    T d0 = (m_v[1] - m_v[0]) / m_domain.m_h[0]; // Because next loop starts at 1
    T d1;
    const auto* vIt = m_v.data() + 1;
    const auto* hIt = m_domain.m_h.data() + 1;
    for (auto sIt = &m_s[1]; sIt != &m_s[size - 1]; ++sIt, ++vIt, ++hIt) {
      // d[i] = (m_v[i + 1] - m_v[i]) / m_h[i];
      d1 = (*(vIt + 1) - *vIt) / *hIt;
      // m_s[i] = (m_v[i + 1] - m_v[i]) / m_h[i] - (m_v[i] - m_v[i - 1]) / m_h[i - 1];
      *sIt = d1 - d0;
      d0 = d1;
    }
    // m_s[0] and m_s[size - 1] are left at 0 for natural splines
    if constexpr (Cache != SplineCache::Early) {
      const auto max = m_s.size() - 2;
      for (std::size_t i = 1; i <= max; ++i) {
        m_cache.insert(i);
      }
    }
  }

  /**
   * @brief Update the internal coefficients associated to the i-th knot.
  */
  inline void update(std::size_t i) {
    if constexpr (Cache != SplineCache::Early) {
      m_s[i] = (m_v[i + 1] - m_v[i]) / m_domain.m_h[i] - (m_v[i] - m_v[i - 1]) / m_domain.m_h[i - 1];
      m_cache.insert(i);
    }
  }

  /**
   * @brief Evaluate the spline.
   */
  T operator()(double x) {
    return operator()(SplineArg(m_domain, x));
  }

  /**
   * @brief Evaluate the spline.
   */
  T operator()(const SplineArg& x) {
    const auto i = x.m_index;
    return m_v[i] * x.m_cv0 + m_v[i + 1] * x.m_cv1 + dv2(i) * x.m_cs0 + dv2(i + 1) * x.m_cs1;
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

private:
  const Partition& m_domain; ///< The knots domain
  std::vector<T> m_v; ///< The knot values
  std::vector<T> m_s; ///< The knot second derivatives
  std::set<std::size_t> m_cache; ///< Cached indices of m_s
};

/**
 * @brief Natural cubic spline resampler.
 * 
 * A resampler is parametrized with the list of knot and resampling abscissae,
 * and is evaluated on a vector of knot values.
 */
class Cospline {

public:
  /**
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Cospline(const Partition& domain, TIt begin, TIt end) : m_domain(domain), m_args(end - begin) {
    std::transform(std::move(begin), std::move(end), m_args.begin(), [&](const auto& e) {
      return SplineArg(m_domain, e);
    });
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Cospline(const Partition& u, const TRange& x) : Cospline(u, x.begin(), x.end()) {}

  /**
   * @brief List-based constructor.
   */
  template <typename T>
  explicit Cospline(const Partition& u, std::initializer_list<T> x) : Cospline(u, x.begin(), x.end()) {}

  /**
   * @brief Assign arguments from an iterator.
   */
  template <typename TIt>
  void assign(TIt begin, TIt end) {
    m_args.resize(std::distance(begin, end));
    std::transform(std::move(begin), std::move(end), m_args.begin(), [&](const auto& x) {
      return SplineArg(m_domain, x);
    });
  }

  /**
   * @brief Assign arguments from a range.
   */
  template <typename TRange>
  void assign(const TRange& x) {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Assign arguments from a list.
   */
  void assign(const std::initializer_list<double>& x) {
    assign(x.begin(), x.end());
  }

  /**
   * @brief Resample a spline defined by an iterator over knot values.
   */
  template <typename TIt>
  std::vector<typename std::iterator_traits<TIt>::value_type> operator()(TIt begin, TIt end) const {
    using T = typename std::iterator_traits<TIt>::value_type;
    return Spline<T>(m_domain, std::move(begin), std::move(end))(m_args);
  }

  /**
   * @brief Resample a spline defined by a range of knot values.
   */
  template <typename TRange>
  std::vector<typename TRange::value_type> operator()(const TRange& v) const {
    return operator()(v.begin(), v.end());
  }

  /**
   * @brief Resample a spline defined by a list of knot values.
   */
  template <typename T>
  std::vector<T> operator()(std::initializer_list<T> v) const {
    return operator()(v.begin(), v.end());
  }

private:
  const Partition& m_domain; ///< The knot abscissae
  std::vector<SplineArg> m_args; ///< The resampling abscissae
};

} // namespace Splider

#endif
