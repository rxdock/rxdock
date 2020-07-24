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

#include "BaseGrid.h"
#include "SATypes.h"

namespace rxdock {

// typedefs are in SATypes.h
// rest of the class architecture from NonBondedGrid.h
// but instead Atom HHS_Solvation is used as operand

class NonBondedHHSGrid : public BaseGrid {
public:
  static std::string _CT;
  NonBondedHHSGrid(const Coord &gridMin, const Coord &gridStep, unsigned int NX,
                   unsigned int NY, unsigned int NZ, unsigned int NPad = 0);
  NonBondedHHSGrid(std::istream &istr);
  ~NonBondedHHSGrid();

  NonBondedHHSGrid(const NonBondedHHSGrid &);
  NonBondedHHSGrid(const BaseGrid &);
  NonBondedHHSGrid &operator=(const NonBondedHHSGrid &);
  NonBondedHHSGrid &operator=(const BaseGrid &);

  virtual void Print(std::ostream &ostr) const;
  virtual void Write(std::ostream &ostr) const;
  virtual void Read(std::istream &istr);

  const HHS_SolvationRList &GetHHSList(unsigned int iXYZ) const;
  const HHS_SolvationRList &GetHHSList(const Coord &c) const;

  void SetHHSLists(HHS_Solvation *pHHS, double radius);
  void ClearHHSLists(void);

protected:
  void OwnPrint(std::ostream &ostr) const;
  void OwnWrite(std::ostream &ostr) const;
  void OwnRead(std::istream &istr);

private:
  NonBondedHHSGrid();

  void CopyGrid(const NonBondedHHSGrid &);
  void CreateMap();

  HHS_SolvationListMap m_hhsMap;
  const HHS_SolvationRList m_emptyList;
};

typedef SmartPtr<NonBondedHHSGrid> NonBondedHHSGridPtr;

} // namespace rxdock

#endif // _RBTNONBONDEDHHSGRID_H_
