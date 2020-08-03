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

#include "GPFFGold.h"
#include "Debug.h"
#include "GPGenome.h"
#include "GPParser.h"

#include <loguru.hpp>

#include <fstream>
#include <sstream>

using namespace rxdock;
using namespace rxdock::geneticprogram;

std::string GPFFGold::_CT("GPFFGold");

void GPFFGold::ReadTables(std::istream &in, ReturnTypeArray &,
                          ReturnTypeArray &) {
  ReturnType value;
  int i = 0, j;
  in >> value;
  inputTable.clear();
  SFTable.clear();
  std::string name;
  while (!in.eof()) {
    ReturnTypeList ivalues;
    ivalues.push_back(new ReturnType(value));
    for (j = 1; j < GPGenome::GetNIP(); j++) {
      in >> value;
      ivalues.push_back(new ReturnType(value));
    }
    ReturnTypeList sfvalues;
    for (j = 0; j < GPGenome::GetNSFI(); j++) {
      in >> value;
      sfvalues.push_back(new ReturnType(value));
    }
    // read name, don't store it for now
    in >> name;
    //        if (*sfvalues[1] < 3.0)
    //        {
    inputTable.push_back(ReturnTypeList(ivalues));
    SFTable.push_back(ReturnTypeList(sfvalues));
    //        }
    i++;
    in >> value;
  }
}

double GPFFGold::CalculateFitness(GPGenomePtr g, ReturnTypeArray &it,
                                  ReturnTypeArray &sft, bool function) {
  GPParser p(GPGenome::GetNIP(), GPGenome::GetNIF(), GPGenome::GetNN(),
             GPGenome::GetNO());
  ReturnTypeList o;
  double tot = 0.0;
  // double range = 10.0;
  double good = 0.0;
  double bad = 0.0;
  double hitlimit = 2.0;
  for (unsigned int i = 0; i < it.size(); i++) {
    ReturnTypeList inputs(it[i]);
    ReturnTypeList SFValues = sft[i];
    double scoreM = *(SFValues[0]);
    double rmsd = *(SFValues[1]);
    GPChromosomePtr c(g->GetChrom());
    o = p.Parse(c, inputs);
    LOG_F(1, "parsed: ", *(o[0]));
    for (int j = 0; j < GPGenome::GetNO(); j++)
      if (function) {
        double d = 2.5 * rmsd + scoreM - *(o[j]);
        tot += d * d; // std::abs(d); // * d;
        // if (std::abs(2.5 * rmsd + scoreM - *(o[j])) < range)
        //     tot++;
      } else if (*(o[j]) < 0.0) {
        if (rmsd < hitlimit)
          good++; // += 1.2;
        else
          bad++;
      } else if (rmsd >= hitlimit)
        ; // good++;
      else
        good--;
    c->Clear();
  }
  // objective value always an increasing function
  if (function)
    objective = -tot;
  else
    //        if (good >=250 && bad <250)
    //            objective = 200;
    //        else
    objective = good * 1 - bad;
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  fitness = objective;
  g->SetFitness(fitness);
  return fitness;
}

double GPFFGold::CalculateFitness(GPGenomePtr g, ReturnTypeArray &it,
                                  ReturnTypeArray &sft, double hitlimit,
                                  bool function) {
  GPParser p(g->GetNIP(), g->GetNIF(), g->GetNN(), g->GetNO());
  ReturnTypeList o;
  double truehits = 0.0;
  double falsehits = 0.0;
  double truemisses = 0.0;
  double falsemisses = 0.0;
  for (unsigned int i = 0; i < it.size(); i++) {
    ReturnTypeList inputs = it[i];
    ReturnTypeList SFValues = sft[i];
    double scoreM = *(SFValues[0]);
    double rmsd = *(SFValues[1]);
    LOG_F(1, "inputs:");
    for (unsigned int j = 0; j < inputs.size(); j++) {
      LOG_F(1, "{}", *(inputs[j]));
    }
    GPChromosomePtr c = g->GetChrom();
    o = p.Parse(c, inputs);
    LOG_F(1, "parsed: ", *(o[0]));
    double limit = function ? scoreM + 5 : 0.0;
    for (int j = 0; j < GPGenome::GetNO(); j++) {
      if (rmsd < hitlimit)
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
    c->Clear();
  }
  // objective value always an increasing function
  LOG_F(1, "True hits: {}, false hits: {}, true misses: {}, false misses: {}",
        truehits, falsehits, truemisses, falsemisses);
  objective = truehits / (truehits + falsehits);
  fitness = objective;
  g->SetFitness(fitness);
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  return fitness;
}
