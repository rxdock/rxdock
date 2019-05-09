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

#ifndef _RBTNONBONDEDHHSGRID_H_
#define _RBTNONBONDEDHHSGRID_H_

#include "RbtBaseGrid.h"
#include "RbtSATypes.h"

// typedefs are in RbtSATypes.h
// rest of the class architecture from RbtNonBondedGrid.h
// but instead RbtAtom HHS_Solvation is used as operand

class RbtNonBondedHHSGrid : public RbtBaseGrid {
public:
  static std::string _CT;
  RbtNonBondedHHSGrid(const RbtCoord &gridMin, const RbtCoord &gridStep,
                      unsigned int NX, unsigned int NY, unsigned int NZ,
                      unsigned int NPad = 0);
  RbtNonBondedHHSGrid(std::istream &istr);
  ~RbtNonBondedHHSGrid();

  RbtNonBondedHHSGrid(const RbtNonBondedHHSGrid &);
  RbtNonBondedHHSGrid(const RbtBaseGrid &);
  RbtNonBondedHHSGrid &operator=(const RbtNonBondedHHSGrid &);
  RbtNonBondedHHSGrid &operator=(const RbtBaseGrid &);

  virtual void Print(std::ostream &ostr) const;
  virtual void Write(std::ostream &ostr) const;
  virtual void Read(std::istream &istr);

  const HHS_SolvationRList &GetHHSList(unsigned int iXYZ) const;
  const HHS_SolvationRList &GetHHSList(const RbtCoord &c) const;

  void SetHHSLists(HHS_Solvation *pHHS, double radius);
  void ClearHHSLists(void);

protected:
  void OwnPrint(std::ostream &ostr) const;
  void OwnWrite(std::ostream &ostr) const;
  void OwnRead(std::istream &istr);

private:
  RbtNonBondedHHSGrid();

  void CopyGrid(const RbtNonBondedHHSGrid &);
  void CreateMap();

  HHS_SolvationListMap m_hhsMap;
  const HHS_SolvationRList m_emptyList;
};

typedef SmartPtr<RbtNonBondedHHSGrid> RbtNonBondedHHSGridPtr;
#endif // _RBTNONBONDEDHHSGRID_H_
