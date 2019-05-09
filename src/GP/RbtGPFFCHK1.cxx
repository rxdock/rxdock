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

#include "RbtGPFFCHK1.h"
#include "RbtCellTokenIter.h"
#include "RbtDebug.h"
#include "RbtFilterExpression.h"
#include "RbtFilterExpressionVisitor.h"
#include "RbtGPGenome.h"
#include "RbtGPParser.h"
#include "RbtParser.h"
#include "RbtRand.h"
#include "RbtTokenIter.h"
#include <fstream>
#include <sstream>

std::string RbtGPFFCHK1::_CT("RbtGPFFCHK1");

void RbtGPFFCHK1::ReadTables(std::istream &in, RbtReturnTypeArray &it,
                             RbtReturnTypeArray &sft) {
  RbtReturnType value;
  int nip, nsfi;
  int i = 0, j, recordn;
  in >> nip;
  in.get();
  int nctes = 15;
  RbtGPGenome::SetNIP(nip + nctes);
  in >> nsfi;
  RbtGPGenome::SetNSFI(nsfi);
  in >> recordn; // read record number, don't store it
  inputTable.clear();
  SFTable.clear();
  char name[80];
  CreateRandomCtes(nctes);
  while (!in.eof()) {
    // read name, don't store it for now
    in.get();                  // read ','
    in.getline(name, 80, ','); // read name, don't store it
    in >> value;
    RbtReturnTypeList ivalues;
    ivalues.push_back(new RbtReturnType(value));
    for (j = 1; j < nip; j++) {
      in.get();
      in >> value;
      ivalues.push_back(new RbtReturnType(value));
    }
    for (j = 0; j < nctes; j++)
      ivalues.push_back(new RbtReturnType(ctes[j]));
    RbtReturnTypeList sfvalues;
    for (j = 0; j < nsfi; j++) {
      in.get();
      in >> value;
      sfvalues.push_back(new RbtReturnType(value));
    }
    inputTable.push_back(RbtReturnTypeList(ivalues));
    SFTable.push_back(RbtReturnTypeList(sfvalues));
    i++;
    in >> recordn;
  }
  std::cout << "Read: " << inputTable[0][0] << std::endl;
  it = inputTable;
  sft = SFTable;
}

double RbtGPFFCHK1::CalculateFitness(RbtGPGenomePtr g, RbtReturnTypeArray &it,
                                     RbtReturnTypeArray &sft, bool function) {
  double tp = 0.0; // True Positives
  double fp = 0.0; // False Positives
  double fn = 0.0; // True Negatives
  double tn = 0.0; // False Negatives
  RbtGPChromosomePtr c = g->GetChrom();
  RbtParser parse;
  RbtCellTokenIterPtr tokenIter(new RbtCellTokenIter(c, contextp));
  RbtFilterExpressionPtr fe(parse.Parse(tokenIter, contextp));
  double actualOutput, predictedOutput;
  double hitlimit = 0.0;
  for (int i = 0; i < it.size(); i++) {
    RbtReturnTypeList inputs(it[i]);
    RbtReturnType pred;
    for (int j = 0; j < inputs.size(); j++)
      contextp->Assign(j, *(inputs[j]));
    RbtReturnTypeList SFValues = sft[i];
    actualOutput = *(SFValues[SFValues.size() - 1]);
    EvaluateVisitor visitor(contextp);
    fe->Accept(visitor);
    predictedOutput = RbtReturnType(fe->GetValue());
    if (actualOutput < hitlimit) // is it a hit?
    {                            // hit
      if (predictedOutput < hitlimit)
        tp++; // hit predicted: True Positive
      else
        fn++; // miss predicted: False Negative
    } else {  // miss
      if (predictedOutput >= hitlimit)
        tn++; // miss predicted: True Negative
      else
        fp++; // hit predicted: False Positive
    }
    //        std::cout << tp
  }
  /*
  objective = 1.0 * (tp - fn) - fp;
objective = ( tp / (tp +fn) ) * ( tn / (tn + fp) );
// true positive rate * true negative rate
objective = ( tp + tn) / (tp + tn + fp + fn);
// true positive rate * true negative rate */
  objective = (tp / (tp + fp)) * (tp / (tp + fn));
  // precision * true positive rate
  fitness = objective;
  g->SetFitness(fitness);
  return fitness;
}

double RbtGPFFCHK1::CalculateFitness(RbtGPGenomePtr g, RbtReturnTypeArray &it,
                                     RbtReturnTypeArray &sft, double hitlimit,
                                     bool function) {
  double tp = 0.0; // True Positives
  double fp = 0.0; // False Positives
  double fn = 0.0; // True Negatives
  double tn = 0.0; // False Negatives
  RbtGPChromosomePtr c = g->GetChrom();
  RbtParser parse;
  RbtCellTokenIterPtr tokenIter(new RbtCellTokenIter(c, contextp));
  RbtFilterExpressionPtr fe(parse.Parse(tokenIter, contextp));
  double actualOutput, predictedOutput;
  for (int i = 0; i < it.size(); i++) {
    RbtReturnTypeList inputs(it[i]);
    RbtReturnType pred;
    for (int j = 0; j < inputs.size(); j++)
      contextp->Assign(j, *(inputs[j]));
    RbtReturnTypeList SFValues = sft[i];
    actualOutput = *(SFValues[SFValues.size() - 1]);
    EvaluateVisitor visitor(contextp);
    fe->Accept(visitor);
    predictedOutput = RbtReturnType(fe->GetValue());
    if (actualOutput < hitlimit) // is it a hit?
    {                            // hit
      if (predictedOutput < hitlimit)
        tp++; // hit predicted: True Positive
      else
        fn++; // miss predicted: False Negative
    } else {  // miss
      if (predictedOutput >= hitlimit)
        tn++; // miss predicted: True Negative
      else
        fp++; // hit predicted: False Positive
    }
  }
  objective = (tp / (tp + fp)) * (tp / (tp + fn));
  // precision * true positive rate
  std::cout << fn << "\t" << tn << "\t" << fp << "\t" << tp << std::endl;
  fitness = objective;
  // g->SetFitness(fitness);
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  return fitness;
}

void RbtGPFFCHK1::CreateRandomCtes(int nctes) {
  if (ctes.size() == 0) // it hasn't already been initialized
  {
    int a, b;
    double c;
    ctes.push_back(0.0);
    ctes.push_back(1.0);
    std::cout << "c0 \t0.0" << std::endl;
    std::cout << "c1 \t1.0" << std::endl;
    for (int i = 0; i < (nctes - 2); i++) {
      a = m_rand.GetRandomInt(200) - 100;
      b = m_rand.GetRandomInt(10) - 5;
      c = (a / 10.0) * pow(10, b);
      std::cout << "c" << i + 2 << " \t" << c << std::endl;
      ctes.push_back(c);
    }
  }
}
