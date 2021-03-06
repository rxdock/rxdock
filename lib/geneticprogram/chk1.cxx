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

#include "rxdock/Debug.h"
#include "rxdock/Parser.h"
#include "rxdock/geneticprogram/GPFFCHK1.h"
#include "rxdock/geneticprogram/GPFitnessFunction.h"
#include "rxdock/geneticprogram/GPGenome.h"
#include "rxdock/geneticprogram/GPPopulation.h"
#include <cstdio>
#include <fstream>

int main(int argc, char *argv[]) {
  try {
    Rand &theRand = GetRand(); // ref to random number generator
    if (argc > 1)
      theRand.Seed(std::atoi(argv[2]));
    else
      theRand.SeedFromRandomDevice();
    std::cout << "Seed: " << theRand.GetSeed() << std::endl;
    GPGenome::SetStructure(56, 3, 1, 1, 7, 1, 200, 100);
    ifstream desc("descnames", ios::in);
    ContextPtr contextp(new CellContext(desc));
    GPFitnessFunctionPtr ff = new GPFFCHK1(contextp);
    std::string strTrainingFile = argv[1];
    strTrainingFile += "_training.dat";
    fstream trainingFile(strTrainingFile.c_str(), ios::in);
    if (!trainingFile)
      throw Error(_WHERE_, "can't open " + strTrainingFile);
    std::string strTestFile = argv[1];
    strTestFile += "_testing.dat";
    fstream testFile(strTestFile.c_str(), ios::in);
    if (!testFile)
      throw Error(_WHERE_, "can't open " + strTestFile);
    ReturnTypeArray ittrain, sfttrain, ittest, sfttest;
    ff->ReadTables(trainingFile, ittrain, sfttrain);
    //    if (!testFile)
    //      throw Error(_WHERE_, "can't open" + strTestFile);
    ff->ReadTables(testFile, ittest, sfttest);

    double hitlimit = 0.0;
    GPPopulation p(5, 5, ff, ittrain, sfttrain);
    p.Initialise(hitlimit, false);
    double b = 0.0, currentBest = 0.0;
    int i = 0, count;
    while ((b < 0.9) && (i < 50000)) {
      p.EPstep("", 1.0, 0.08, 0.0, 0.0, hitlimit, false);
      if ((i % 500) == 0) {
        b = ff->CalculateFitness(p.Best(), ittrain, sfttrain, false);
        if (b > currentBest) {
          currentBest = b;
          count = 0;
        } else
          count++;
        std::cout << b << "\t"
                  << ff->CalculateFitness(p.Best(), ittrain, sfttrain, hitlimit,
                                          false)
                  << std::endl;
        if (count >= 7)
          break;
      }
      i++;
    }
    std::cout << i << "\t" << *(p.Best()) << std::endl;
    std::cout << "best with training: "
              << ff->CalculateFitness(p.Best(), ittrain, sfttrain, false)
              << "\n"
              << ff->CalculateFitness(p.Best(), ittrain, sfttrain, 0.0, false)
              << std::endl;
    std::cout << "best with test:     "
              << ff->CalculateFitness(p.Best(), ittest, sfttest, false) << "\n"
              << ff->CalculateFitness(p.Best(), ittest, sfttest, 0.0, false)
              << std::endl;

  } catch (Error &e) {
    std::cout << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)
  return 0;
}
