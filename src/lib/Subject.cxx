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

#include <algorithm>

using namespace rxdock;

// Constructor -does nothing.
Subject::Subject() { _RBTOBJECTCOUNTER_CONSTR_("Subject"); }

// Destructor - notify observers of impending destruction
Subject::~Subject() {
#ifdef _DEBUG
  std::cout << "Subject::~Subject: Notifying observers of impending destruction"
            << std::endl;
#endif //_DEBUG
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
  ObserverListIter iter =
      std::find(m_observers.begin(), m_observers.end(), pObserver);
  if (iter != m_observers.end()) {
    throw BadArgument(
        _WHERE_,
        "Subject::Attach(): pObserver is already attached to this subject");
  } else {
    m_observers.push_back(pObserver);
#ifdef _DEBUG
    std::cout << "Subject::Attach: Attaching new observer; #observers="
              << m_observers.size() << std::endl;
#endif //_DEBUG
  }
}

void Subject::Detach(Observer *pObserver) {
  ObserverListIter iter =
      std::find(m_observers.begin(), m_observers.end(), pObserver);
  if (iter == m_observers.end()) {
    throw BadArgument(
        _WHERE_, "Subject::Detach(): pObserver not attached to this subject");
  } else {
    m_observers.erase(iter);
#ifdef _DEBUG
    std::cout << "Subject::Detach: Detaching observer; #observers="
              << m_observers.size() << std::endl;
#endif //_DEBUG
  }
}

void Subject::Notify() {
#ifdef _DEBUG
  std::cout << "Subject::Notify: Notifying observers of change of state"
            << std::endl;
#endif //_DEBUG
  for (ObserverListIter iter = m_observers.begin(); iter != m_observers.end();
       iter++) {
    (*iter)->Update(this);
  }
}
