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

#ifndef _RBTINTERACTIONTEMPLATE_H_
#define _RBTINTERACTIONTEMPLATE_H_
#include "BaseGrid.h"
#include "Config.h"

#include <loguru.hpp>

namespace rxdock {

template <class T> class InteractionNode {
  T property;

public:
  InteractionNode() { _RBTOBJECTCOUNTER_CONSTR_("InteractionNode"); }
  InteractionNode(T aProperty) {
    _RBTOBJECTCOUNTER_CONSTR_("InteractionNode");
    property = aProperty;
  }
  ~InteractionNode() { _RBTOBJECTCOUNTER_DESTR_("InteractionNode"); }

  T GetProperty() const { return property; }
  void SetProperty(T aProperty) { property = aProperty; }
};

template <typename T> // instead InteractionCenterList
class InteractionNodeList : public std::vector<InteractionNode<T> *> {
}; // std:vector ctors

template <typename T> // instead InteractionListMap
class InteractionNodeListMap : public std::vector<InteractionNodeList<T>> {
}; // std:vector ctors

template <typename T>
class InteractionNodeListMapIter : public InteractionNodeListMap<T>::iterator {
};

template <typename T>
class InteractionNodeListMapConstIter
    : public InteractionNodeListMap<T>::const_iterator {};
// Most of the stuff is copied shamelessly from InteractionGrid.[cxx,h]
template <class T> class InteractionGridTemplate : public BaseGrid {
public:
  static std::string _CT;

  InteractionGridTemplate(const Coord &gridMin, const Coord &gridStep,
                          unsigned int NX, unsigned int NY, unsigned int NZ,
                          unsigned int NPad = 0) {
    _CT = "InteractionGridTemplate";
    CreateMap();
    _RBTOBJECTCOUNTER_CONSTR_("InteractionGridTemplate");
  }

  // Constructor reading params from binary stream
  InteractionGridTemplate(std::istream &istr) {
    _CT = "InteractionGridTemplate";
    CreateMap();
    _RBTOBJECTCOUNTER_CONSTR_("InteractionGridTemplate");
  }

  ~InteractionGridTemplate() {
    ClearInteractionLists();
    _RBTOBJECTCOUNTER_DESTR_("InteractionGridTemplate");
  }

  // Copy constructor
  InteractionGridTemplate(const InteractionGridTemplate &aGrid)
      : BaseGrid(aGrid) {
    CreateMap();
    CopyGrid(aGrid);
    _RBTOBJECTCOUNTER_COPYCONSTR_("InteractionGridTemplate");
  }

  // Copy constructor taking base class argument
  InteractionGridTemplate(const BaseGrid &) {
    CreateMap();
    _RBTOBJECTCOUNTER_COPYCONSTR_("InteractionGridTemplate");
  }

  // Copy assignment
  InteractionGridTemplate &operator=(const InteractionGridTemplate &aGrid) {
    if (this != &aGrid) {
      ClearInteractionLists();
      BaseGrid::Operator = (aGrid);
      CopyGrid(aGrid);
    }
    return *this;
  }

  // Copy assignment taking base class argument
  InteractionGridTemplate &operator=(const BaseGrid &aGrid) {
    if (this != &aGrid) {
      ClearInteractionLists();
      BaseGrid::Operator = (aGrid);
    }
    return *this;
  }

  virtual void Print(std::ostream &ostr) const { // Text output
    BaseGrid::Print(ostr);
    OwnPrint(ostr);
  }

  const InteractionNodeList<T> &GetInteractionList(unsigned int iXYZ) const {
    if (isValid(iXYZ)) {
      return m_intnMap[iXYZ];
    } else {
      return m_emptyList;
    }
  }

  const InteractionNodeList<T> &GetInteractionList(const Coord &c) const {
    if (isValid(c)) {
      return m_intnMap[GetIXYZ(c)];
    } else {
      LOG_F(1, "InteractionGridTemplate::GetInteractionList: {} is off grid",
            c);
      return m_emptyList;
    }
  }

  void SetInteractionLists(T *pIntn, double radius) {
    Atom *pAtom = pIntn->GetProperty().atom;
    if (nullptr == pAtom)
      return;
    const Coord &c = pAtom->GetCoords();
    std::vector<unsigned int> sphereIndices;
    GetSphereIndices(c, radius, sphereIndices);
    for (std::vector<unsigned int>::const_iterator sphereIter =
             sphereIndices.begin();
         sphereIter != sphereIndices.end(); sphereIter++) {
      m_intnMap[*sphereIter].push_back(pIntn);
    }
  }

  void ClearInteractionLists() {
    for (InteractionNodeListMap<T>::iterator iter = m_intnMap.begin();
         iter != m_intnMap.end(); iter++) {
      (*iter).clear();
    }
  }

protected:
  void OwnPrint(std::ostream &ostr) const {
    ostr << std::endl << "Class\t" << _CT << std::endl;
    ostr << "No. of entries in the map: " << m_intnMap.size() << std::endl;
  }

private:
  InteractionGridTemplate();

  void CopyGrid(const InteractionGridTemplate<T> &aGrid) {
    m_intnMap = aGrid.m_intnMap;
  }

  void CreateMap() { m_intnMap = InteractionNodeListMap<T>(GetN()); }

  InteractionNodeListMap<T> m_intnMap;      // Used to store the interaction
                                            // center lists at each grid point
  const InteractionNodeList<T> m_emptyList; // Dummy list used by GetAtomList
};

template <typename T> // instead InteractionGridPtr
class InteractionGridTemplatePtr : public SmartPtr<InteractionGridTemplate<T>> {
};

} // namespace rxdock

#endif //  _RBTINTERACTIONTEMPLATE_H_
