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

// Visitors will be created to
// execute the functions represented by the expressions, print them etc.

#ifndef _RBT_FILTEREXPRESSIONVISITOR_H_
#define _RBT_FILTEREXPRESSIONVISITOR_H_

#include "rxdock/Context.h"
#include "rxdock/FilterExpression.h"

namespace rxdock {

class FilterVbleExp;
class FilterAddExp;
class FilterSubExp;
class FilterMulExp;
class FilterDivExp;
class FilterAndExp;
class FilterLogExp;
class FilterExpExp;
class FilterIfExp;

class FilterExpressionVisitor {
public:
  virtual void VisitVbleExp(FilterVbleExp *) = 0;
  virtual void VisitAddExp(FilterAddExp *) = 0;
  virtual void VisitSubExp(FilterSubExp *) = 0;
  virtual void VisitMulExp(FilterMulExp *) = 0;
  virtual void VisitDivExp(FilterDivExp *) = 0;
  virtual void VisitAndExp(FilterAndExp *) = 0;
  virtual void VisitLogExp(FilterLogExp *) = 0;
  virtual void VisitExpExp(FilterExpExp *) = 0;
  virtual void VisitIfExp(FilterIfExp *) = 0;
  ///////////////////
  // Destructor
  //////////////////
  virtual ~FilterExpressionVisitor();

  ////////////////////
  // Protected methods
  ////////////////////

protected:
  FilterExpressionVisitor(); // Default constructor disabled
};

class PrintVisitor : public FilterExpressionVisitor {
public:
  PrintVisitor(ContextPtr);
  ~PrintVisitor();
  virtual void VisitVbleExp(FilterVbleExp *);
  virtual void VisitAddExp(FilterAddExp *);
  virtual void VisitSubExp(FilterSubExp *);
  virtual void VisitMulExp(FilterMulExp *);
  virtual void VisitDivExp(FilterDivExp *);
  virtual void VisitAndExp(FilterAndExp *);
  virtual void VisitLogExp(FilterLogExp *);
  virtual void VisitExpExp(FilterExpExp *);
  virtual void VisitIfExp(FilterIfExp *);

private:
  ContextPtr contextp;
  std::ostringstream strbuf;
};

class PrettyPrintVisitor : public FilterExpressionVisitor {
public:
  PrettyPrintVisitor(ContextPtr);
  ~PrettyPrintVisitor();
  virtual void VisitVbleExp(FilterVbleExp *);
  virtual void VisitAddExp(FilterAddExp *);
  virtual void VisitSubExp(FilterSubExp *);
  virtual void VisitMulExp(FilterMulExp *);
  virtual void VisitDivExp(FilterDivExp *);
  virtual void VisitAndExp(FilterAndExp *);
  virtual void VisitLogExp(FilterLogExp *);
  virtual void VisitExpExp(FilterExpExp *);
  virtual void VisitIfExp(FilterIfExp *);

private:
  ContextPtr contextp;
  std::ostringstream strbuf;
};

class EvaluateVisitor : public FilterExpressionVisitor {
public:
  EvaluateVisitor(ContextPtr);
  virtual void VisitVbleExp(FilterVbleExp *);
  virtual void VisitAddExp(FilterAddExp *);
  virtual void VisitSubExp(FilterSubExp *);
  virtual void VisitMulExp(FilterMulExp *);
  virtual void VisitDivExp(FilterDivExp *);
  virtual void VisitAndExp(FilterAndExp *);
  virtual void VisitLogExp(FilterLogExp *);
  virtual void VisitExpExp(FilterExpExp *);
  virtual void VisitIfExp(FilterIfExp *);

private:
  ContextPtr contextp;
};

} // namespace rxdock

#endif
