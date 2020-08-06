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

#ifndef _RBTNMSTATE_H_
#define _RBTNMSTATE_H_

namespace rxdock {

namespace neldermead {

/**
 * Optimization state for an optimizer
 */
template <class DataType, class ParameterType> struct State {
  ParameterType bestParameters;
  DataType bestValue;

  ParameterType currentParameters;
  DataType currentValue;

  ParameterType formerParameters;
  DataType formerValue;

  long iteration;
};

} // namespace neldermead

} // namespace rxdock

#endif /* _RBTNMSTATE_H_ */
