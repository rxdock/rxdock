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

#ifndef _RBTSIMPLEXCOSTFUNCTION_H
#define _RBTSIMPLEXCOSTFUNCTION_H

#ifdef __PGI
#define EIGEN_DONT_VECTORIZE
#endif
#include <Eigen/Core>

#include "BaseSF.h"
#include "ChromElement.h"

namespace rxdock {

// Rosenbrock cost function
class SimplexCostFunction {
public:
  SimplexCostFunction(BaseSF *pSF, ChromElementPtr chrom)
      : m_pSF(pSF), m_chrom(chrom) {}
  typedef double DataType;
  typedef Eigen::VectorXd ParameterType;
  double operator()(const ParameterType &parameters) { // const
    nCalls++;
    std::vector<double> vec(parameters.data(),
                            parameters.data() + parameters.size());
    m_chrom->SetVector(vec);
    m_chrom->SyncToModel();
    return m_pSF->Score();
  }
  unsigned long int nCalls = 0;

private:
  BaseSF *m_pSF;
  ChromElementPtr m_chrom;
};

} // namespace rxdock

#endif // _RBTSIMPLEXCOSTFUNCTION_H
