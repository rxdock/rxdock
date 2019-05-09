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

#include "RbtGPFFSpike.h"
#include "RbtCellTokenIter.h"
#include "RbtDebug.h"
#include "RbtFilterExpression.h"
#include "RbtFilterExpressionVisitor.h"
#include "RbtGPGenome.h"
#include "RbtGPParser.h"
#include "RbtParser.h"
#include "RbtTokenIter.h"
#include <cassert>
#include <fstream>
#include <sstream>

std::string RbtGPFFSpike::_CT("RbtGPFFSpike");
void RbtGPFFSpike::ReadTables(std::istream &in, RbtReturnTypeArray &it,
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
  in >> recordn;
  inputTable.clear();
  SFTable.clear();
  char name[80];
  CreateRandomCtes(nctes);
  while (!in.eof()) {
    // read name, don't store it for now
    in.get();
    in.getline(name, 80, ',');
    in >> value;
    RbtReturnTypeList ivalues;
    ivalues.push_back(new RbtReturnType(value));
    for (j = 1; j < nip; j++) {
      in.get();
      in >> value;
      ivalues.push_back(new RbtReturnType(value));
    }
    for (j = 0; j < nctes; j++) {
      ivalues.push_back(new RbtReturnType(ctes[j]));
    }
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
  //    std::cout << it[0].size() << std::endl;
  sft = SFTable;
}

double RbtGPFFSpike::CalculateFitness(RbtGPGenomePtr g, RbtReturnTypeArray &it,
                                      RbtReturnTypeArray &sft, bool function) {
  //    RbtGPParser p(RbtGPGenome::GetNIP(), RbtGPGenome::GetNIF(),
  //                RbtGPGenome::GetNN(), RbtGPGenome::GetNO());
  RbtReturnTypeList o;
  o.push_back(new RbtReturnType(1.1));
  RbtReturnType oldo;
  double tot = 0.0;
  double good = 0.0;
  double bad = 0.0;
  double hitlimit = 0.0;
  RbtGPChromosomePtr c(g->GetChrom());
  RbtParser p2;
  RbtCellTokenIterPtr ti(new RbtCellTokenIter(c, contextp));
  RbtFilterExpressionPtr fe(p2.Parse(ti, contextp));
  //    RbtCellTokenIterPtr ti(new RbtCellTokenIter(c, contextp));
  for (int i = 0; i < it.size(); i++) {
    RbtReturnTypeList inputs(it[i]);
    for (int j = 0; j < inputs.size(); j++)
      contextp->Assign(j, *(inputs[j]));
    RbtReturnTypeList SFValues = sft[i];
    double hit = *(SFValues[SFValues.size() - 1]);
    //        RbtParser p2;
    //        RbtCellTokenIterPtr ti(new RbtCellTokenIter(c, contextp));
    //        RbtFilterExpressionPtr fe(p2.Parse(ti, contextp));
    EvaluateVisitor visitor(contextp);
    fe->Accept(visitor);
    *(o[0]) = fe->GetValue();
    // std::cout << *(o[0]) << std::endl;
    for (int j = 0; j < RbtGPGenome::GetNO(); j++)
      if (function) {
        std::cout << "Error, no function possible with spike\n";
        exit(1);
      } else if (*(o[j]) < 0.0) {
        if (hit < hitlimit)
          good++; // += 1.2;
        else
          bad++;
      } else if (hit >= hitlimit)
        ; // good++;
      else
        good--;
    //        c->Clear();
  }
  // objective value always an increasing function
  objective = good * 1.5 - bad;
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  fitness = objective;
  g->SetFitness(fitness);
  return fitness;
}

double RbtGPFFSpike::CalculateFitness(RbtGPGenomePtr g, RbtReturnTypeArray &it,
                                      RbtReturnTypeArray &sft, double hitlimit,
                                      bool function) {
  //    RbtGPParser p(g->GetNIP(), g->GetNIF(), g->GetNN(), g->GetNO());
  RbtReturnTypeList o;
  o.push_back(new RbtReturnType(1.1));
  RbtReturnType oldo;
  double truehits = 0.0;
  double falsehits = 0.0;
  double truemisses = 0.0;
  double falsemisses = 0.0;
  RbtGPChromosomePtr c(g->GetChrom());
  RbtParser p2;
  RbtTokenIterPtr ti(new RbtCellTokenIter(c, contextp));
  RbtFilterExpressionPtr fe(p2.Parse(ti, contextp));
  for (int i = 0; i < it.size(); i++) {
    RbtReturnTypeList inputs = it[i];
    for (int j = 0; j < inputs.size(); j++)
      contextp->Assign(j, *(inputs[j]));
    RbtReturnTypeList SFValues = sft[i];
    double hit = *(SFValues[SFValues.size() - 1]);
    //        o = p.Parse(c, inputs);
    //        oldo = *(o[0]);
    //        RbtParser p2; //new RbtCellTokenIter(c, contextp));
    //        RbtTokenIterPtr ti(new RbtCellTokenIter(c, contextp));
    //        RbtFilterExpressionPtr fe(p2.Parse(ti, contextp));
    EvaluateVisitor visitor(contextp);
    fe->Accept(visitor);
    *(o[0]) = fe->GetValue();
    //        std::cout << "** " << *(o[0]) << std::endl;
    //        assert(oldo == *(o[0]));
    double limit = 0.0;
    for (int j = 0; j < RbtGPGenome::GetNO(); j++) {
      if (hit < hitlimit)
        if (*(o[j]) < limit)
          // true hit
          truehits++;
        else
          // false miss
          falsemisses++;
      else if (*(o[j]) < limit)
        // false hit
        falsehits++;
      else
        // true miss
        truemisses++;
    }
    //        c->Clear();
  }
  // objective value always an increasing function
  std::cout << truehits << "\t" << falsehits << "\t" << truemisses << "\t"
            << falsemisses << std::endl;
  objective = truehits / (truehits + falsehits);
  fitness = objective;
  // g->SetFitness(fitness);
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  return fitness;
}

void RbtGPFFSpike::CreateRandomCtes(int nctes) {
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
