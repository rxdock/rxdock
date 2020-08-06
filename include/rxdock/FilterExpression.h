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

// Takes a list of tokens and parsers them to create a Filter Expression

#ifndef _RBT_FILTEREXPRESSION_H_
#define _RBT_FILTEREXPRESSION_H_

#include "rxdock/Error.h"
#include "rxdock/FilterExpressionVisitor.h"
#include "rxdock/TokenIter.h"
#include "rxdock/geneticprogram/GPTypes.h"

#include <iostream>

namespace rxdock {

class FilterExpressionVisitor;

class FilterExpression {
public:
  static std::string _CT;
  virtual void Print() = 0;
  virtual void Accept(FilterExpressionVisitor &) = 0;
  virtual int GetNOps() = 0;
  virtual SmartPtr<FilterExpression> GetOp(int) = 0;
  ///////////////////
  // Destructor
  //////////////////
  virtual ~FilterExpression();
  ReturnType GetValue() { return value; }
  void SetValue(ReturnType v) { value = v; }

  ////////////////////
  // Private methods
  ////////////////////

protected:
  FilterExpression(); // Default constructor disabled
  ReturnType value;
};

// Useful typedefs
typedef SmartPtr<FilterExpression> FilterExpressionPtr; // Smart pointer
typedef std::vector<FilterExpressionPtr> FilterExpressionList;
// Vector of smart pointers
typedef FilterExpressionList::iterator FilterExpressionListIter;
typedef FilterExpressionList::const_iterator FilterExpressionListConstIter;

class FilterVbleExp : public FilterExpression {
public:
  FilterVbleExp(const Vble &);
  virtual ~FilterVbleExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() { std::cout << "v" << vble.GetName() << " "; }
  int GetNOps() { return 0; }
  const Vble &GetVble() const { return vble; }
  FilterExpressionPtr GetOp(int i) {
    throw Error(_WHERE_, "Vble Expressions don't have operands");
  }

private:
  const Vble &vble;
};

class FilterLogExp : public FilterExpression {
public:
  FilterLogExp(FilterExpressionPtr);
  virtual ~FilterLogExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "log ";
    operand->Print();
  }
  virtual int GetNOps() { return 1; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand;
    else
      throw Error(_WHERE_, "Log Expressions only have 1 operand");
  }

private:
  FilterExpressionPtr operand;
};

class FilterExpExp : public FilterExpression {
public:
  FilterExpExp(FilterExpressionPtr);
  virtual ~FilterExpExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "exp ";
    operand->Print();
  }
  virtual int GetNOps() { return 1; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand;
    else
      throw Error(_WHERE_, "Exp Expressions only have 1 operand");
  }

private:
  FilterExpressionPtr operand;
};

class FilterAddExp : public FilterExpression {
public:
  FilterAddExp(FilterExpressionPtr, FilterExpressionPtr);
  virtual ~FilterAddExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "add ";
    operand1->Print();
    operand2->Print();
  }
  virtual int GetNOps() { return 2; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand1;
    if (i == 1)
      return operand2;
    else
      throw Error(_WHERE_, "Add Expressions only have 2 operands");
  }

private:
  FilterExpressionPtr operand1, operand2;
};

class FilterSubExp : public FilterExpression {
public:
  FilterSubExp(FilterExpressionPtr, FilterExpressionPtr);
  virtual ~FilterSubExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "sub ";
    operand1->Print();
    operand2->Print();
  }
  virtual int GetNOps() { return 2; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand1;
    if (i == 1)
      return operand2;
    else
      throw Error(_WHERE_, "Sub Expressions only have 2 operands");
  }

private:
  FilterExpressionPtr operand1, operand2;
};

class FilterMulExp : public FilterExpression {
public:
  FilterMulExp(FilterExpressionPtr, FilterExpressionPtr);
  virtual ~FilterMulExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "mul ";
    operand1->Print();
    operand2->Print();
  }
  virtual int GetNOps() { return 2; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand1;
    if (i == 1)
      return operand2;
    else
      throw Error(_WHERE_, "Mul Expressions only have 2 operands");
  }

private:
  FilterExpressionPtr operand1, operand2;
};

class FilterDivExp : public FilterExpression {
public:
  FilterDivExp(FilterExpressionPtr, FilterExpressionPtr);
  virtual ~FilterDivExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "div ";
    operand1->Print();
    operand2->Print();
  }
  virtual int GetNOps() { return 2; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand1;
    if (i == 1)
      return operand2;
    else
      throw Error(_WHERE_, "Div Expressions only have 2 operands");
  }

private:
  FilterExpressionPtr operand1, operand2;
};

class FilterAndExp : public FilterExpression {
public:
  FilterAndExp(FilterExpressionPtr, FilterExpressionPtr);
  virtual ~FilterAndExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "and ";
    operand1->Print();
    operand2->Print();
  }
  virtual int GetNOps() { return 2; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand1;
    if (i == 1)
      return operand2;
    else
      throw Error(_WHERE_, "And Expressions only have 2 operands");
  }

private:
  FilterExpressionPtr operand1, operand2;
};

class FilterIfExp : public FilterExpression {
public:
  FilterIfExp(FilterExpressionPtr, FilterExpressionPtr, FilterExpressionPtr);
  virtual ~FilterIfExp();
  void Accept(FilterExpressionVisitor &visitor);
  void Print() {
    std::cout << "if ";
    operand1->Print();
    operand2->Print();
    operand3->Print();
  }
  virtual int GetNOps() { return 3; }
  FilterExpressionPtr GetOp(int i) {
    if (i == 0)
      return operand1;
    if (i == 1)
      return operand2;
    if (i == 2)
      return operand3;
    else
      throw Error(_WHERE_, "Add Expressions only have 3 operands");
  }

private:
  FilterExpressionPtr operand1, operand2, operand3;
};

} // namespace rxdock

#endif //_FilterExpression_H_
