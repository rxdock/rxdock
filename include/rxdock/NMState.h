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

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace rxdock {

namespace neldermead {

/**
 * Optimization state for an optimizer
 */
template <class DataType, class ParameterType> struct State {
  friend void to_json(json &j, const State &state) {
    j = json{{"best-parameters", state.bestParameters},
             {"best-value", state.bestValue},
             {"current-parameters", state.currentParameters},
             {"current-value", state.currentValue},
             {"former-parameters", state.formerParameters},
             {"former-value", state.formerValue},
             {"iteration", state.iteration}};
  }

  friend void from_json(const json &j, State &state) {
    j.at("best-parameters").get_to(state.bestParameters);
    j.at("best-value").get_to(state.bestValue);
    j.at("current-parameters").get_to(state.currentParameters);
    j.at("current-value").get_to(state.currentValue);
    j.at("former-parameters").get_to(state.formerParameters);
    j.at("former-value").get_to(state.formerValue);
    j.at("iteration").get_to(state.iteration);
  }

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
