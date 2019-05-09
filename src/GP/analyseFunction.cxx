/***********************************************************************
 * The rDock program was developed from 1998 - 2006 by the software team
 * at RiboTargets (subsequently Vernalis (R&D) Ltd).
 * In 2006, the software was licensed to the University of York for
 * maintenance and distribution.
 * In 2012, Vernalis and the University of York agreed to release the
 * program as Open Source software.
 * This version is licensed under GNU-LGPL version 3.0 with support from
 * the University of Barcelona.
 * http://rdock.sourceforge.net/
 ***********************************************************************/

#include "../inc/RbtGPGenome.h"
#include "../inc/RbtGPPopulation.h"
#include "../inc/RbtParser.h"
#include <cstdio>
#include <fstream>
#include <sstream>

int main(int argc, char *argv[]) {
  /*std::cout << "Genome file: \n";
  std::string gfile;
  std::cin >> gfile; */
  std::string gfile(argv[1]);
  std::ifstream gstr(gfile.c_str(), std::ios::in);
  RbtGPGenome gen(gstr);
  /*std::cout << "Names file: \n";
  std::string nfile;
  std::cin >> nfile;*/
  std::string nfile(argv[2]);
  std::ifstream nstr(nfile.c_str(), std::ios::in);
  RbtParser p(gen.GetNIP(), gen.GetNIF(), gen.GetNN(), gen.GetNO());
  RbtInt i = 1;
  RbtInt pos = nstr.tellg();
  while (i != EOF) {
    gen.GetChrom()->Clear();
    nstr.close();
    nstr.open(nfile.c_str(), std::ios::in);
    std::cout << p.PrintParse(nstr, gen.GetChrom(), true, false) << std::endl;
    gen.GetChrom()->Clear();
    nstr.close();
    nstr.open(nfile.c_str(), std::ios::in);
    std::cout << p.PrintParse(nstr, gen.GetChrom(), false, false);
    i = EOF; // cin.get();
  }
  return 0;
}
