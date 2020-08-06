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

#include "Subject.h"

#include <loguru.hpp>

#include <algorithm>

using namespace rxdock;

// Constructor -does nothing.
Subject::Subject() { _RBTOBJECTCOUNTER_CONSTR_("Subject"); }

// Destructor - notify observers of impending destruction
Subject::~Subject() {
  LOG_F(2, "Subject destructor");
  LOG_F(1, "Notifying observers of impending destruction");
  // We need to iterate using a while loop because call to Deleted will
  // trigger a call back to Detach, reducing the size of m_observers, hence
  // invalidating a conventional iterator
  while (m_observers.size() > 0) {
    Observer *pObserver = m_observers.back();
    pObserver->Deleted(this);
  }
  _RBTOBJECTCOUNTER_DESTR_("Subject");
}

////////////////////////////////////////
// Public methods
////////////////
void Subject::Attach(Observer *pObserver) {
  LOG_F(2, "Subject::Attach");
  ObserverListIter iter =
      std::find(m_observers.begin(), m_observers.end(), pObserver);
  if (iter != m_observers.end()) {
    throw BadArgument(
        _WHERE_,
        "Subject::Attach(): pObserver is already attached to this subject");
  } else {
    m_observers.push_back(pObserver);
    LOG_F(1, "Attaching new observer; #observers={}", m_observers.size());
  }
}

void Subject::Detach(Observer *pObserver) {
  LOG_F(2, "Subject::Detach");
  ObserverListIter iter =
      std::find(m_observers.begin(), m_observers.end(), pObserver);
  if (iter == m_observers.end()) {
    throw BadArgument(
        _WHERE_, "Subject::Detach(): pObserver not attached to this subject");
  } else {
    m_observers.erase(iter);
    LOG_F(1, "Detaching observer; #observers={}", m_observers.size());
  }
}

void Subject::Notify() {
  LOG_F(2, "Subject::Notify");
  LOG_F(1, "Notifying observers of change of state");
  for (ObserverListIter iter = m_observers.begin(); iter != m_observers.end();
       iter++) {
    (*iter)->Update(this);
  }
}
