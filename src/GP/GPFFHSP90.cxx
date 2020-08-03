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

#include "GPFFHSP90.h"
#include "CellTokenIter.h"
#include "Debug.h"
#include "FilterExpression.h"
#include "FilterExpressionVisitor.h"
#include "GPGenome.h"
#include "GPParser.h"
#include "Parser.h"
#include "Rand.h"
#include "TokenIter.h"

#include <loguru.hpp>

#include <cassert>
#include <fstream>
#include <sstream>

using namespace rxdock;
using namespace rxdock::geneticprogram;

std::string GPFFHSP90::_CT("GPFFHSP90");
int nInversions(int idx, ReturnTypeArray &list);

void GPFFHSP90::ReadTables(std::istream &in, ReturnTypeArray &it,
                           ReturnTypeArray &sft) {
  ReturnType value;
  int nip, nsfi;
  int i = 0, j;
  in >> nip;
  in.get();
  int nctes = 15;
  GPGenome::SetNIP(nip + nctes);
  in >> nsfi;
  GPGenome::SetNSFI(nsfi);
  inputTable.clear();
  SFTable.clear();
  char name[80];
  // read mol, don't store it for now
  in.getline(name, 80, ',');
  CreateRandomCtes(nctes);
  while (!in.eof()) {
    ReturnTypeList sfvalues;
    for (j = 0; j < nsfi; j++) {
      in >> value;
      in.get(); // read ','
      sfvalues.push_back(new ReturnType(value));
    }
    // read Ki, don't store it
    in.getline(name, 80, ',');
    // read Name, don't store it
    in.getline(name, 80, ',');
    // read RMSD, don't store it
    in.getline(name, 80, ',');
    in >> value;
    ReturnTypeList ivalues;
    ivalues.push_back(new ReturnType(value));
    for (j = 1; j < nip; j++) {
      in.get(); // read ','
      in >> value;
      ivalues.push_back(new ReturnType(value));
    }
    for (j = 0; j < nctes; j++) {
      ivalues.push_back(new ReturnType(ctes[j]));
    }
    inputTable.push_back(ReturnTypeList(ivalues));
    SFTable.push_back(ReturnTypeList(sfvalues));
    i++;
    in.getline(name, 80, ',');
  }
  it = inputTable;
  sft = SFTable;
}

double GPFFHSP90::CalculateFitness(
    GPGenomePtr g, ReturnTypeArray &it, ReturnTypeArray &sft,
    //                                         Double limit,
    bool function) {
  ReturnTypeArray o;
  //    o.push_back(new ReturnType(1.1));
  GPChromosomePtr c(g->GetChrom());
  Parser p2;
  CellTokenIterPtr ti(new CellTokenIter(c, contextp));
  FilterExpressionPtr fe(p2.Parse(ti, contextp));
  // double meanReal = 0.0, meanPred = 0.0;
  for (unsigned int i = 0; i < it.size(); i++) {
    ReturnTypeList inputs(it[i]);
    ReturnTypeList pred;
    pred.clear();
    for (unsigned int j = 0; j < inputs.size(); j++)
      contextp->Assign(j, *(inputs[j]));
    ReturnTypeList SFValues = sft[i];
    EvaluateVisitor visitor(contextp);
    fe->Accept(visitor);
    pred.push_back(new ReturnType(fe->GetValue()));
    o.push_back(pred);
    LOG_F(1, "{} {}", *(SFValues[0]), *(o[i][0]));
  }
  int ninv = nInversions(0, o);
  LOG_F(1, "n inversions {}", ninv);
  /*    Double sumDiff = 0.0,
                sumSqrReal = 0.0, sumSqrPred = 0.0;

      meanReal = meanReal / it.size();
      meanPred = meanPred / it.size();
  // LOG_F(1, "{} {} {}", it.size(), meanPred, meanReal);
  for (Int i = 0 ; i < it.size() ; i++)
      {
        Double t1 = *sft[i][0]; // - meanReal;
        Double t2 = *o[i][0]; // - meanPred;
        sumDiff += std::abs(t1 - t2);
  // LOG_F(1, "{} {} {}", sumDiff, sumSqrReal, sumSqrPred);
      }
      objective = sumDiff / it.size();

  // LOG_F(1, "{}", objective); */

  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  // objective value increasing function
  fitness = -ninv;
  g->SetFitness(fitness);
  return fitness;
}

double GPFFHSP90::CalculateFitness(GPGenomePtr g, ReturnTypeArray &it,
                                   ReturnTypeArray &sft, double limit,
                                   bool function) {
  ReturnTypeArray o;
  GPChromosomePtr c(g->GetChrom());
  Parser p2;
  CellTokenIterPtr ti(new CellTokenIter(c, contextp));
  FilterExpressionPtr fe(p2.Parse(ti, contextp));
  double meanReal = 0.0, meanPred = 0.0;
  for (unsigned int i = 0; i < it.size(); i++) {
    ReturnTypeList inputs(it[i]);
    ReturnTypeList pred;
    pred.clear();
    for (unsigned int j = 0; j < inputs.size(); j++)
      contextp->Assign(j, *(inputs[j]));
    ReturnTypeList SFValues = sft[i];
    EvaluateVisitor visitor(contextp);
    fe->Accept(visitor);
    pred.push_back(new ReturnType(fe->GetValue()));
    o.push_back(pred);
    meanPred += *(o[i][0]);
    meanReal += *(SFValues[0]);
    LOG_F(1, "{} {}", *(SFValues[0]), *(o[i][0]));
  }
  double sumDiff = 0.0, sumSqrReal = 0.0, sumSqrPred = 0.0;

  meanReal = meanReal / it.size();
  meanPred = meanPred / it.size();
  for (unsigned int i = 0; i < it.size(); i++) {
    double t1 = *sft[i][0] - meanReal;
    double t2 = *o[i][0] - meanPred;
    sumDiff += t1 * t2;
    sumSqrReal += t1 * t1;
    sumSqrPred += t2 * t2;
  }
  objective = (sumDiff * sumDiff) / (sumSqrReal * sumSqrPred);
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  // objective value increasing function
  fitness = objective;
  g->SetFitness(fitness);
  return fitness;
}

void GPFFHSP90::CreateRandomCtes(int nctes) {
  if (ctes.size() == 0) // it hasn't already been initialized
  {
    int a, b;
    double c;
    ctes.push_back(0.0);
    ctes.push_back(1.0);
    LOG_F(1, "c0=0.0");
    LOG_F(1, "c1=1.0");
    for (int i = 0; i < (nctes - 2); i++) {
      a = m_rand.GetRandomInt(200) - 100;
      b = m_rand.GetRandomInt(10) - 5;
      c = (a / 10.0) * std::pow(10, b);
      LOG_F(1, "c{}={}", i + 2, c);
      ctes.push_back(c);
    }
  }
}

int nInversions(int idx, ReturnTypeArray &list) {
  int count = 0;
  int len = list.size();
  if (idx == (len - 1))
    return 0;
  for (int i = idx + 1; i < len; i++) {
    if (*(list[idx][0]) >= *(list[i][0]))
      count++;
  }
  return (count + nInversions(idx + 1, list));
}
