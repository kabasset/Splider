/// @copyright 2023-2024, Antoine Basset (CNES)
// This file is part of Splider <github.com/kabasset/Splider>
// SPDX-License-Identifier: Apache-2.0

#ifndef _SPLIDER_ARGUMENT_H
#define _SPLIDER_ARGUMENT_H

#include "Linx/Data/Vector.h" // Index
#include "Splider/Mode.h"

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

  template <typename>
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
  template <typename TDomain>
  explicit SplineArg(const TDomain& domain, Value x)
  {
    m_index = domain.index(x);
    const auto h = domain.length(m_index);
    const auto left = x - domain[m_index];
    const auto right = h - left;
    m_cv0 = right / h;
    m_cv1 = 1. - m_cv0;
    m_c6s0 = right * (right * m_cv0 - h);
    m_c6s1 = left * (left * m_cv1 - h);
  }

  /**
   * @brief Subinterval-based constructor.
   * @param index The subinterval index
   * @param index The subinterval length
   * @param left The distance between the subinterval min and the argument
   */
  explicit SplineArg(Linx::Index index, Value length, Value left)
  {
    m_index = index;
    const auto right = length - left;
    m_cv0 = right / length;
    m_cv1 = left / length;
    m_c6s0 = right * (right * m_cv0 - length);
    m_c6s1 = left * (left * m_cv1 - length);
  }

private:

  Linx::Index m_index; ///< The interval index
  Value m_cv0; ///< The v[i] coefficient
  Value m_cv1; ///< The v[i + 1] coefficient
  Value m_c6s0; ///< The s[i] coefficient
  Value m_c6s1; ///< The s[i + 1] coefficient
};

/**
 * @brief Spline arguments.
*/
template <typename T = double>
class Args {
  template <typename, typename, Mode>
  friend class Spline;

  template <typename>
  friend class BiCospline; // TODO rm

public:

  /**
   * @brief The real number type.
   */
  using Value = T;

  /**
   * @brief Iterator-based constructor.
   */
  template <typename TDomain, typename TIt>
  explicit Args(const TDomain& domain, TIt begin, TIt end) : m_args()
  {
    m_args.reserve(std::distance(begin, end));
    for (; begin != end; ++begin) {
      m_args.emplace_back(domain, *begin);
    }
  }

  /**
   * @brief Range-based constructor.
   */
  template <typename TDomain, typename TRange>
  explicit Args(const TDomain& domain, const TRange& u) : Args(domain, u.begin(), u.end())
  {}

  /**
   * @brief List-based constructor.
   */
  template <typename TDomain>
  Args(const TDomain& domain, std::initializer_list<Value> u) : Args(domain, u.begin(), u.end())
  {}

  /**
   * @brief Get the number of arguments.
   */
  std::size_t size() const
  {
    return m_args.size();
  }

  /**
   * @copybrief size()
   */
  Linx::Index ssize() const
  {
    return static_cast<Linx::Index>(size());
  }

private:

  std::vector<SplineArg<Value>> m_args; ///< The arguments
};

} // namespace Splider

#endif
