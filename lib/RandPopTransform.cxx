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

#include "rxdock/RandPopTransform.h"
#include "rxdock/Chrom.h"
#include "rxdock/Population.h"
#include "rxdock/WorkSpace.h"

#include <loguru.hpp>

using namespace rxdock;

const std::string RandPopTransform::_CT = "RandPopTransform";
const std::string RandPopTransform::_POP_SIZE = "POP_SIZE";
const std::string RandPopTransform::_SCALE_CHROM_LENGTH = "SCALE_CHROM_LENGTH";

RandPopTransform::RandPopTransform(const std::string &strName)
    : BaseBiMolTransform(_CT, strName) {
  AddParameter(_POP_SIZE, 50);
  AddParameter(_SCALE_CHROM_LENGTH, true);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RandPopTransform::~RandPopTransform() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

////////////////////////////////////////
// Protected methods
///////////////////
void RandPopTransform::SetupReceptor() {}

void RandPopTransform::SetupLigand() {}

void RandPopTransform::SetupSolvent() {}

void RandPopTransform::SetupTransform() {
  // Construct the overall chromosome for the system
  m_chrom = new Chrom(GetWorkSpace()->GetModels());
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in BaseTransform
// Actually apply the transform
void RandPopTransform::Execute() {
  if (m_chrom.Ptr() == nullptr) {
    return;
  }
  BaseSF *pSF = GetWorkSpace()->GetSF();
  if (pSF == nullptr) {
    return;
  }
  int popSize = GetParameter(_POP_SIZE);
  bool bScale = GetParameter(_SCALE_CHROM_LENGTH);
  if (bScale) {
    int chromLength = m_chrom->GetLength();
    popSize *= chromLength;
  }
  LOG_F(2, "RandPopTransform::Execute: popSize={}", popSize);
  PopulationPtr pop = new Population(m_chrom, popSize, pSF);
  pop->Best()->GetChrom()->SyncToModel();
  GetWorkSpace()->SetPopulation(pop);
}
