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

#include "RbtDebug.h"
#include "RbtGPFFHSP90.h"
#include "RbtGPFitnessFunction.h"
#include "RbtGPGenome.h"
#include "RbtGPPopulation.h"
#include "RbtParser.h"
#include <cstdio>
#include <fstream>

int main(int argc, char *argv[]) {
  try {
    RbtRand &theRand = GetRbtRand(); // ref to random number generator
    if (argc > 1)
      theRand.Seed(std::atoi(argv[1]));
    else
      theRand.SeedFromRandomDevice();
    std::cout << "Seed: " << theRand.GetSeed() << std::endl;
    std::cout << "Output taken from: \n";
    std::string strInputFile;
    std::cin >> strInputFile;
    std::ifstream inputFile(strInputFile.c_str(), std::ios::in);
    if (!inputFile)
      throw RbtError(_WHERE_, "can't open " + strInputFile);
    RbtGPGenome::SetStructure(56, 3, 1, 1, 7, 1, 200, 100);
    std::ifstream desc("descnames", std::ios::in);
    RbtContextPtr contextp(new RbtCellContext(desc));
    RbtGPFitnessFunctionPtr ff = new RbtGPFFHSP90(contextp);
    std::cout << "Number of input sets: \n";
    int nInputSets;
    std::cin >> nInputSets;
    int nTestSets = nInputSets / 3;
    std::string strTrainingFile = "/tmp/beatriz/" + strInputFile + ".training";
    std::fstream trainingFile(strTrainingFile.c_str(),
                              std::ios::in | std::ios::out);
    std::string strTestFile = "/tmp/beatriz/" + strInputFile + ".test";
    std::fstream testFile(strTestFile.c_str(), std::ios::in | std::ios::out);
    std::set<unsigned int> randomList;
    std::cout << nTestSets << std::endl;
    for (int i = 0; i < nTestSets; i++) {
      bool insertb;
      do {
        unsigned int ran = theRand.GetRandomInt(nInputSets);
        insertb = (randomList.insert(ran)).second;
      } while (!insertb);
    }
    std::set<unsigned int>::iterator iter = randomList.end();
    for (iter = randomList.begin(); iter != randomList.end(); iter++)
      std::cout << *iter << "\t";
    std::cout << std::endl;
    std::string line;
    getline(inputFile, line);
    testFile << line << std::endl;
    trainingFile << line << std::endl;
    for (int i = 0; i < nInputSets; i++) {
      getline(inputFile, line);
      if (randomList.find(i) != randomList.end())
        // found it
        testFile << line << std::endl;
      else
        trainingFile << line << std::endl;
    }

    trainingFile.seekp(0, std::ios::beg);
    testFile.seekp(0, std::ios::beg);

    RbtReturnTypeArray ittrain, sfttrain, ittest, sfttest;
    ff->ReadTables(trainingFile, ittrain, sfttrain);
    //    if (!testFile)
    //      throw RbtError(_WHERE_, "can't open" + strTestFile);
    ff->ReadTables(testFile, ittest, sfttest);

    double hitlimit = 0.0;
    RbtGPPopulation p(5, 5, ff, ittrain, sfttrain);
    p.Initialise(hitlimit, false);
    double b = 0.0;
    int i = 0;
    while ((b < 0.9) && (i < 10000)) {
      p.EPstep("", 1.0, 0.08, 0.0, 0.0, hitlimit, false);
      b = ff->CalculateFitness(p.Best(), ittrain, sfttrain, false);
      if ((i % 500) == 0) {
        //        std::cout << ff->CalculateFitness(p.Best(),ittest,sfttest,
        //                                    hitlimit,false) << std::endl;
        std::cout << b << "\t"
                  << ff->CalculateFitness(p.Best(), ittrain, sfttrain, false)
                  << std::endl;
        if ((i == 5000) && (b < 0))
          break;
      }
      i++;
    }
    std::cout << i << "\t" << *(p.Best()) << std::endl;
    std::cout << "best with training: "
              << ff->CalculateFitness(p.Best(), ittrain, sfttrain, false)
              << std::endl;
    std::cout << "best with test:     "
              << ff->CalculateFitness(p.Best(), ittest, sfttest, false)
              << std::endl;
    inputFile.seekg(0, std::ios::beg);
    ff->ReadTables(inputFile, ittest, sfttest);
    std::cout << "best with all:      "
              << ff->CalculateFitness(p.Best(), ittest, sfttest, false)
              << std::endl;

  } catch (RbtError &e) {
    std::cout << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)
  return 0;
}
