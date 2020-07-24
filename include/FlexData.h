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

// Abstract base class for all model flexibility types
// Provides Accept method for FlexDataVisitors
//(subclasses must implement)
#ifndef _RBTFLEXDATA_
#define _RBTFLEXDATA_

#include "FlexDataVisitor.h"
#include "ParamHandler.h"

namespace rxdock {

class Model;
class DockingSite;

class FlexData : public ParamHandler {
public:
  virtual ~FlexData() {}
  Model *GetModel() const { return m_pModel; }
  void SetModel(Model *pModel) { m_pModel = pModel; }
  DockingSite *GetDockingSite() const { return m_pDockSite; }
  virtual void Accept(FlexDataVisitor &) = 0;

protected:
  FlexData(DockingSite *pDockSite)
      : m_pModel(nullptr), m_pDockSite(pDockSite) {}

private:
  Model *m_pModel;
  DockingSite *m_pDockSite;
};

typedef SmartPtr<FlexData> FlexDataPtr;

} // namespace rxdock

#endif //_RBTFLEXDATA_
