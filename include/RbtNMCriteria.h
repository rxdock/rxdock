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

#ifndef _RBTNMCRITERIA_H_
#define _RBTNMCRITERIA_H_

#include <cmath>
#include <limits>

namespace RbtNM {
class IterationCriterion {
private:
  int m_maxIterations;

public:
  IterationCriterion(long max_iterations) : m_maxIterations(max_iterations) {}

  template <class State> bool operator()(const State &state) {
    bool result = state.iteration < m_maxIterations;
#ifdef _DEBUG
    std::cout << "IterationCriteria: " << result << std::endl;
#endif // _DEBUG
    return result;
  }
};

template <class DataType> class RelativeValueCriterion {
private:
  DataType m_ftol;

public:
  RelativeValueCriterion(DataType ftol) : m_ftol(ftol) {}

  template <class State> bool operator()(const State &state) {
    bool result =
        2 * std::abs(state.currentValue - state.formerValue) /
            (std::abs(state.currentValue) + std::abs(state.formerValue) +
             std::numeric_limits<DataType>::epsilon()) >
        m_ftol;
#ifdef _DEBUG
    std::cout << "RelativeValueCriterion: " << result << std::endl;
#endif // _DEBUG
    return result;
  }
};

template <class Criteria1, class Criteria2> class AndCriteria {
private:
  Criteria1 m_criteria1;
  Criteria2 m_criteria2;

public:
  AndCriteria(const Criteria1 &criteria1, const Criteria2 &criteria2)
      : m_criteria1(criteria1), m_criteria2(criteria2) {}

  template <class State> bool operator()(const State &state) {
    return m_criteria1(state) && m_criteria2(state);
  }
};

template <class Criteria1, class Criteria2>
AndCriteria<Criteria1, Criteria2>
CreateAndCriteria(const Criteria1 &criteria1, const Criteria2 &criteria2) {
  return AndCriteria<Criteria1, Criteria2>(criteria1, criteria2);
}
} // namespace RbtNM

#endif /* _RBTNMCRITERIA_H_ */
