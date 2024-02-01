/// @copyright 2023, Antoine Basset
// SPDX-License-Identifier: GPL-3.0-or-later

#include "Splider/Partition.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Partition_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(index_test)
{
  const Splider::Partition u {1, 2, 3, 4};
  for (std::size_t i = 0; i < u.size() - 1; ++i) {
    BOOST_TEST(u.index(u[i]) == i);
  }
  BOOST_TEST(u.index(u[u.size() - 1]) == u.size() - 2);
  for (std::size_t i = 0; i < u.size() - 1; ++i) {
    BOOST_TEST(u.index((u[i] + u[i + 1]) / 2) == i);
  }
  BOOST_CHECK_THROW(u.index(u[0] - 1), std::runtime_error);
  BOOST_CHECK_THROW(u.index(u[u.size() - 1] + 1), std::runtime_error);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
