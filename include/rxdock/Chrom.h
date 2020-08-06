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

// Aggregate class for managing a vector of ChromElement objects
// Also implements the ChromElement interface itself
// to delegate requests to each element in the vector
// NOTE: Chrom destructor is responsible for deleting all managed elements
#ifndef RBTCHROM_H_
#define RBTCHROM_H_

#include "ChromElement.h"
#include "Error.h"
#include "Model.h"

namespace rxdock {

class Chrom : public ChromElement {
public:
  // Class type string
  static std::string _CT;
  // Constructor for an empty chromosome
  RBTDLL_EXPORT Chrom();
  // Constructor for a combined chromosome for each model
  // in the list
  Chrom(const ModelList &modelList);
  virtual ~Chrom();

  // pure virtual from ChromElement
  virtual void Reset();
  virtual void Randomise();
  virtual void Mutate(double relStepSize);
  virtual void SyncFromModel();
  virtual void SyncToModel();
  virtual ChromElement *clone() const;
  virtual int GetLength() const;
  virtual int GetXOverLength() const;
  virtual void GetVector(std::vector<double> &v) const;
  virtual void GetVector(XOverList &v) const;
  virtual void SetVector(const std::vector<double> &v, int &i);
  virtual void SetVector(const XOverList &v, int &i);
  virtual void GetStepVector(std::vector<double> &v) const;
  virtual double CompareVector(const std::vector<double> &v, int &i) const;
  virtual void Print(std::ostream &s) const;

  // Aggregate methods
  // Appends a new chromosome element to the vector
  // Chrom destructor is responsible for deleting the new element
  // Null operation if pChromElement is nullptr
  virtual void Add(ChromElement *pChromElement);

protected:
private:
  ChromElementList m_elementList;
  // We need to store the model list so that
  // we can call UpdatePseudoAtoms() on each model following
  // a SyncToModel
  ModelList m_modelList;
};

} // namespace rxdock

#endif /*RBTCHROM_H_*/
