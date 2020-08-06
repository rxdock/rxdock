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

#include "rxdock/ChromElement.h"
#include "rxdock/Rand.h"

#include <loguru.hpp>

using namespace rxdock;

std::string ChromElement::_CT = "ChromElement";
double ChromElement::_THRESHOLD = 1E-4;

double &ChromElement::GetThreshold() { return _THRESHOLD; }

ChromElement::eMode ChromElement::StrToMode(const std::string &modeStr) {
  eMode retVal;
  std::string modeStrUpper = modeStr;
  // std::toupper will not work https://stackoverflow.com/a/7131881
  std::transform(modeStrUpper.begin(), modeStrUpper.end(), modeStrUpper.begin(),
                 ::toupper);
  if (modeStrUpper == "FIXED") {
    retVal = FIXED;
  } else if (modeStrUpper == "TETHERED") {
    retVal = TETHERED;
  } else if (modeStrUpper == "FREE") {
    retVal = FREE;
  } else {
    throw BadArgument(_WHERE_,
                      "Unknown chromosome sampling mode (" + modeStr + ")");
  }
  return retVal;
}

std::string ChromElement::ModeToStr(ChromElement::eMode mode) {
  std::string retVal;
  switch (mode) {
  case FIXED:
    retVal = "FIXED";
    break;
  case TETHERED:
    retVal = "TETHERED";
    break;
  case FREE:
    retVal = "FREE";
    break;
  default:
    retVal = "UNDEFINED";
    break;
  }
  return retVal;
}

ChromElement::ChromElement() : m_rand(GetRandInstance()) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

ChromElement::~ChromElement() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void ChromElement::Add(ChromElement *pChromElement) {
  throw InvalidRequest(
      _WHERE_,
      "Add(ChromElement*) invalid for non-aggregate chromosome element");
}

bool ChromElement::VectorOK(const std::vector<double> &v,
                            unsigned int i) const {
  unsigned int length = GetLength();
  // if the element is empty then any vector is valid
  return (length == 0) || ((i < v.size()) && (length <= (v.size() - i)));
}

bool ChromElement::VectorOK(const XOverList &v, unsigned int i) const {
  unsigned int length = GetXOverLength();
  // if the element is empty then any vector is valid
  return (length == 0) || ((i < v.size()) && (length <= (v.size() - i)));
}

void ChromElement::CauchyMutate(double mean, double variance) {
  // Need to convert the Cauchy random variable to a positive number
  // and use this as the relative step size for mutation
  double relStepSize = std::fabs(m_rand.GetCauchyRandom(mean, variance));
  Mutate(relStepSize);
}

double ChromElement::Compare(const ChromElement &c) const {
  double retVal(0.0);
  if (GetLength() != c.GetLength()) {
    retVal = -1.0;
  } else {
    std::vector<double> v;
    int i(0);
    c.GetVector(v);
    retVal = CompareVector(v, i);
  }
  return retVal;
}

bool ChromElement::Equals(const ChromElement &c, double threshold) const {
  double cmp = Compare(c);
  return ((cmp >= 0.0) && (cmp < threshold));
}

void ChromElement::SetVector(const std::vector<double> &v) {
  int i(0);
  SetVector(v, i);
}

void ChromElement::SetVector(const XOverList &v) {
  int i(0);
  SetVector(v, i);
}

bool rxdock::operator==(const ChromElement &c1, const ChromElement &c2) {
  return c1.Equals(c2, ChromElement::_THRESHOLD);
}

bool rxdock::operator!=(const ChromElement &c1, const ChromElement &c2) {
  return !(c1 == c2);
}

std::ostream &rxdock::operator<<(std::ostream &s, const ChromElement &c) {
  c.Print(s);
  return s;
}

void rxdock::Crossover(ChromElement *pChr1, ChromElement *pChr2,
                       ChromElement *pChr3, ChromElement *pChr4) {
  // Check all chromosomes have the same crossover length
  int length1 = pChr1->GetXOverLength();
  if ((length1 != pChr2->GetXOverLength()) ||
      (length1 != pChr3->GetXOverLength()) ||
      (length1 != pChr4->GetXOverLength())) {
    throw BadArgument(_WHERE_, "Crossover: mismatch in chromosome lengths");
  }
  // Extract the gene vectors from each parent
  XOverList v1, v2;
  pChr1->GetVector(v1);
  pChr2->GetVector(v2);
  // 2-point crossover
  // In the spirit of STL, ixbegin is the first gene to crossover, ixend is one
  // after the last gene to crossover
  Rand &rand = pChr1->GetRand();
  int ixbegin = rand.GetRandomInt(length1);
  // if ixbegin is 0, we need to avoid selecting the whole chromosome
  int ixend = (ixbegin == 0)
                  ? rand.GetRandomInt(length1 - 1) + 1
                  : rand.GetRandomInt(length1 - ixbegin) + ixbegin + 1;

  LOG_F(1, "Crossover: ixbegin = {}, ixend = {}", ixbegin, ixend);
  std::swap_ranges(v1.begin() + ixbegin, v1.begin() + ixend,
                   v2.begin() + ixbegin);
  // Now we can update the two children
  pChr3->SetVector(v1);
  pChr4->SetVector(v2);
}
