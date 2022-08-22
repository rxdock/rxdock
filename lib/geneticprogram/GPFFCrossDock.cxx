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

#include "rxdock/geneticprogram/GPFFCrossDock.h"
#include "rxdock/Debug.h"
#include "rxdock/geneticprogram/GPGenome.h"
#include "rxdock/geneticprogram/GPParser.h"

#include <loguru.hpp>

#include <fstream>
#include <sstream>

using namespace rxdock;
using namespace rxdock::geneticprogram;

const std::string GPFFCrossDock::_CT = "GPFFCrossDock";

void GPFFCrossDock::ReadTables(std::istream &in, ReturnTypeArray &it,
                               ReturnTypeArray &sft) {
  ReturnType value;
  int nip, nsfi;
  int i = 0, j, recordn;
  in >> nip;
  in.get();
  int nctes = 15;
  GPGenome::SetNIP(nip + nctes);
  in >> nsfi;
  GPGenome::SetNSFI(nsfi);
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
    ReturnTypeList ivalues;
    ivalues.push_back(new ReturnType(value));
    for (j = 1; j < nip; j++) {
      in.get();
      in >> value;
      ivalues.push_back(new ReturnType(value));
    }
    for (j = 0; j < nctes; j++)
      ivalues.push_back(new ReturnType(ctes[j]));
    ReturnTypeList sfvalues;
    for (j = 0; j < nsfi; j++) {
      in.get();
      in >> value;
      sfvalues.push_back(new ReturnType(value));
    }
    inputTable.push_back(ReturnTypeList(ivalues));
    SFTable.push_back(ReturnTypeList(sfvalues));
    i++;
    in >> recordn;
  }
  LOG_F(1, "Read: {}", *inputTable[0][0]);
  it = inputTable;
  sft = SFTable;
}

double GPFFCrossDock::CalculateFitness(GPGenomePtr g, ReturnTypeArray &it,
                                       ReturnTypeArray &sft, bool function) {
  if (function) {
    LOG_F(ERROR, "No function possible with Cross Docking");
    std::exit(1);
  }
  GPParser p(g->GetNIP(), g->GetNIF(), g->GetNN(), g->GetNO());
  ReturnTypeList o;
  double good = 0.0;
  double bad = 0.0;
  double neutral = 0.0;
  double hitlimit = 0.0;
  int ntm, nm;
  ReturnTypeList inputs;
  ReturnTypeList SFValues;
  double hit;
  unsigned int i = 0;
  while (i < it.size()) {
    inputs = it[i];
    SFValues = sft[i];
    hit = *(SFValues[SFValues.size() - 1]);
    ntm = 0;
    nm = 0;
    while (hit > hitlimit) {
      LOG_F(1, "input:");
      for (unsigned int j = 0; j < inputs.size(); j++) {
        LOG_F(1, "{}", *(inputs[j]));
      }
      GPChromosomePtr c = g->GetChrom();
      o = p.Parse(c, inputs);
      c->Clear();
      LOG_F(1, "parsed: {}", *(o[0]));
      if (*(o[0]) >= hitlimit)
        ntm++;
      nm++;
      i++;
      inputs = it[i];
      SFValues = sft[i];
      hit = *(SFValues[SFValues.size() - 1]);
    }
    GPChromosomePtr c = g->GetChrom();
    o = p.Parse(c, inputs);
    c->Clear();
    if (*(o[0]) < hitlimit) {
      //       if (ntm >= (nm / 2.0))
      if (ntm >= ((nm * 2.0) / 3.0))
        good++;
      else
        neutral++;
    } else {
      //        if (ntm >= (nm / 2.0))
      //        if (ntm >= ((nm * 2.0) / 3.0))
      bad++;
    }
    i++;
  }
  objective = good - 0.4 * neutral - bad;
  fitness = objective;
  g->SetFitness(fitness);
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  return fitness;
}

double GPFFCrossDock::CalculateFitness(GPGenomePtr g, ReturnTypeArray &it,
                                       ReturnTypeArray &sft, double hitlimit,
                                       bool function) {
  GPParser p(g->GetNIP(), g->GetNIF(), g->GetNN(), g->GetNO());
  ReturnTypeList o;
  double good = 0.0;
  double bad = 0.0;
  double neutral = 0.0;
  int ntm, nm;
  ReturnTypeList inputs;
  ReturnTypeList SFValues;
  double hit;
  unsigned int i = 0;
  while (i < it.size()) {
    inputs = it[i];
    SFValues = sft[i];
    hit = *(SFValues[SFValues.size() - 1]);
    ntm = 0;
    nm = 0;
    while (hit > hitlimit) {
      LOG_F(1, "input:");
      for (unsigned int j = 0; j < inputs.size(); j++) {
        LOG_F(1, "{}", *(inputs[j]));
      }
      GPChromosomePtr c = g->GetChrom();
      o = p.Parse(c, inputs);
      c->Clear();
      LOG_F(1, "parsed: {}", *(o[0]));
      double limit = 0.0;
      if (*(o[0]) >= limit)
        ntm++;
      nm++;
      i++;
      inputs = it[i];
      SFValues = sft[i];
      hit = *(SFValues[SFValues.size() - 1]);
    }
    GPChromosomePtr c = g->GetChrom();
    o = p.Parse(c, inputs);
    c->Clear();
    if (*(o[0]) < hitlimit) {
      if (ntm >= ((nm * 2.0) / 3.0))
        good++;
      else
        neutral++;
    } else
      bad++;
    i++;
  }
  // objective value always an increasing function
  LOG_F(1, "good={}, neutral={}, bad={}", good, neutral, bad);
  objective = good / (good + bad);
  fitness = objective;
  // g->SetFitness(fitness);
  // For now, I am using tournament selection. So the fitness
  // function doesn't need to be scaled in any way
  return fitness;
}

void GPFFCrossDock::CreateRandomCtes(int nctes) {
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
