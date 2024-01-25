/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef _SPLIDER_MODE_H
#define _SPLIDER_MODE_H

namespace Splider {

/**
 * @brief The spline coefficients evaluation mode.
 */
enum class Mode : char {
  Early = 0b0001, ///< Early evaluation
  Lazy = 0b0010, ///< Lazy evaluation
  Solve = 0b0100, ///< Global, exact solving of the tridiagonal system
  Approximate = 0b1000, ///< Local, finite difference approximation
  Manual = 0 ///< Manual calling of an evaluation function
};

constexpr Mode operator|(Mode lhs, Mode rhs) {
  // FIXME forbid Manual
  return static_cast<Mode>(static_cast<char>(lhs) | static_cast<char>(rhs));
}

constexpr Mode operator&(Mode lhs, Mode rhs) {
  return static_cast<Mode>(static_cast<char>(lhs) & static_cast<char>(rhs));
}

} // namespace Splider

#endif
