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

#include "Chrom.h"

using namespace rxdock;

std::string Chrom::_CT = "Chrom";

Chrom::Chrom() : ChromElement() { _RBTOBJECTCOUNTER_CONSTR_(_CT); }

Chrom::Chrom(const ModelList &modelList)
    : ChromElement(), m_modelList(modelList) {
  for (ModelListConstIter iter = m_modelList.begin(); iter != m_modelList.end();
       ++iter) {
    if ((*iter).Ptr()) {
      Add((*iter)->GetChrom());
    }
  }
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

Chrom::~Chrom() {
  for (ChromElementListIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    delete *iter;
  }
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

// pure virtual from ChromElement
void Chrom::Reset() {
  for (ChromElementListIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->Reset();
  }
}

void Chrom::Randomise() {
  for (ChromElementListIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->Randomise();
  }
}

void Chrom::Mutate(double relStepSize) {
  for (ChromElementListIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->Mutate(relStepSize);
  }
}

void Chrom::SyncFromModel() {
  for (ChromElementListIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->SyncFromModel();
  }
}

void Chrom::SyncToModel() {
  for (ChromElementListIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->SyncToModel();
  }
  // Force an update of all the pseudo atom coords in each model
  for (ModelListIter iter = m_modelList.begin(); iter != m_modelList.end();
       ++iter) {
    (*iter)->UpdatePseudoAtoms();
  }
}

ChromElement *Chrom::clone() const {
  ChromElement *clone = new Chrom();
  for (ChromElementListConstIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    clone->Add((*iter)->clone());
  }
  return clone;
}

void Chrom::Add(ChromElement *pChromElement) {
  if (pChromElement) {
    m_elementList.push_back(pChromElement);
  }
}

int Chrom::GetLength() const {
  int retVal(0);
  for (ChromElementListConstIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    retVal += (*iter)->GetLength();
  }
  return retVal;
}

int Chrom::GetXOverLength() const {
  int retVal(0);
  for (ChromElementListConstIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    retVal += (*iter)->GetXOverLength();
  }
  return retVal;
}

void Chrom::GetVector(std::vector<double> &v) const {
  for (ChromElementListConstIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->GetVector(v);
  }
}

void Chrom::GetVector(XOverList &v) const {
  for (ChromElementListConstIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->GetVector(v);
  }
}

void Chrom::SetVector(const std::vector<double> &v, int &i) {
  if (VectorOK(v, i)) {
    for (ChromElementListIter iter = m_elementList.begin();
         iter != m_elementList.end(); ++iter) {
      (*iter)->SetVector(v, i);
    }
  } else {
    throw BadArgument(
        _WHERE_, "Index i out of range or insufficient elements remaining");
  }
}

void Chrom::SetVector(const XOverList &v, int &i) {
  if (VectorOK(v, i)) {
    for (ChromElementListIter iter = m_elementList.begin();
         iter != m_elementList.end(); ++iter) {
      (*iter)->SetVector(v, i);
    }
  } else {
    throw BadArgument(
        _WHERE_, "Index i out of range or insufficient elements remaining");
  }
}

void Chrom::GetStepVector(std::vector<double> &v) const {
  for (ChromElementListConstIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    (*iter)->GetStepVector(v);
  }
}

// Returns the maximum difference of any of the underlying chromosome elements
// or -1 if there is a mismatch in vector sizes
double Chrom::CompareVector(const std::vector<double> &v, int &i) const {
  double retVal(0.0);
  if (VectorOK(v, i)) {
    ChromElementListConstIter iter = m_elementList.begin();
    // Exit the while loop as soon as a negative value is returned
    while ((iter != m_elementList.end()) && (retVal >= 0.0)) {
      ChromElement *element = *(iter++);
      double cmp = element->CompareVector(v, i);
      retVal = (cmp < 0) ? cmp : std::max(retVal, cmp);
    }
  } else {
    retVal = -1.0;
  }
  return retVal;
}

void Chrom::Print(std::ostream &s) const {
  s << "CHROM" << std::endl;
  int i(0);
  for (ChromElementListConstIter iter = m_elementList.begin();
       iter != m_elementList.end(); ++iter) {
    s << "ELEMENT " << (i++) << ": ";
    (*iter)->Print(s);
  }
}
