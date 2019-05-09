#include "RbtConfig.h"
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

const std::string EXEVERSION = " ($Id$)";

int main(int argc, char **argv) {
  std::string strExeName(argv[0]);
  Rbt::PrintStdHeader(std::cout, strExeName + EXEVERSION);
  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry =
      CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest(registry.makeTest());
  runner.run();
  _RBTOBJECTCOUNTER_DUMP_(std::cout)
  return 0;
}
