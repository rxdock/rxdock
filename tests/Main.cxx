#include "rxdock/Config.h"
#include "rxdock/Debug.h"
#include <gtest/gtest.h>

using namespace rxdock;

const std::string EXEVERSION = " ($Id$)";

int main(int argc, char **argv) {
  std::string strExeName(argv[0]);
  PrintStdHeader(std::cout, strExeName + EXEVERSION);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  _RBTOBJECTCOUNTER_DUMP_(cout)
  return rc;
}
