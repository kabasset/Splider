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
    const auto h = domain.m_h[m_index];
    const auto g = domain.m_g[m_index];
    const auto left = x - domain[m_index];
    const auto right = h - left;
    m_cv0 = right * g;
    m_cv1 = left * g;
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

/**
 * @brief Spline arguments.
*/
template <typename T = double>
class Args {

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
   * @brief Iterator-based constructor.
   */
  template <typename TIt>
  explicit Args(const Partition<Value>& domain, TIt begin, TIt end) :
      m_indices(std::distance(begin, end)), m_cv0s(m_indices.size()), m_cv1s(m_indices.size()),
      m_cs0s(m_indices.size()), m_cs1s(m_indices.size()) {
    auto iit = m_indices.begin();
    auto v0it = m_cv0s.begin();
    auto v1it = m_cv1s.begin();
    auto s0it = m_cs0s.begin();
    auto s1it = m_cs1s.begin();

    for (; begin != end; ++begin, ++iit, ++v0it, ++v1it, ++s0it, ++s1it) {
      const auto i = domain.index(*begin);
      *iit = i;
      const auto h = domain.m_h[i];
      const auto g = domain.m_g[i];
      const auto left = *begin - domain[i];
      const auto right = h - left;
      *v0it = right * g;
      *v1it = left * g;
      *s0it = right / 6. * (right * *v0it - h);
      *s1it = left / 6. * (left * *v1it - h);
    }
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TRange>
  explicit Args(const Partition<Value>& domain, const TRange& u) : Args(domain, u.begin(), u.end()) {}

  /**
   * @brief List-based constructor.
   */
  Args(const Partition<Value>& domain, std::initializer_list<Value> u) : Args(domain, u.begin(), u.end()) {}

  inline std::size_t size() const {
    return m_indices.size();
  }

private:
  std::vector<Linx::Index> m_indices; ///< The interval indices
  std::vector<Value> m_cv0s; ///< The v[i] coefficients
  std::vector<Value> m_cv1s; ///< The v[i + 1] coefficients
  std::vector<Value> m_cs0s; ///< The s[i] coefficients
  std::vector<Value> m_cs1s; ///< The s[i + 1] coefficients
};

} // namespace Splider

#endif
