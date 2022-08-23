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

#include "rxdock/FilterExpression.h"

using namespace rxdock;

const std::string FilterExpression::_CT = "FilterExpression";

FilterExpression::FilterExpression() { _RBTOBJECTCOUNTER_CONSTR_(_CT); }

FilterExpression::~FilterExpression() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

FilterVbleExp::FilterVbleExp(const Vble &vn) : vble(vn) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterVbleExp::~FilterVbleExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterVbleExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitVbleExp(this);
}

FilterLogExp::FilterLogExp(FilterExpressionPtr op1) : operand(op1) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterLogExp::~FilterLogExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterLogExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitLogExp(this);
}

FilterExpExp::FilterExpExp(FilterExpressionPtr op1) : operand(op1) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterExpExp::~FilterExpExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterExpExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitExpExp(this);
}

FilterAddExp::FilterAddExp(FilterExpressionPtr op1, FilterExpressionPtr op2)
    : operand1(op1), operand2(op2) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterAddExp::~FilterAddExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterAddExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitAddExp(this);
}

FilterSubExp::FilterSubExp(FilterExpressionPtr op1, FilterExpressionPtr op2)
    : operand1(op1), operand2(op2) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterSubExp::~FilterSubExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterSubExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitSubExp(this);
}

FilterMulExp::FilterMulExp(FilterExpressionPtr op1, FilterExpressionPtr op2)
    : operand1(op1), operand2(op2) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterMulExp::~FilterMulExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterMulExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitMulExp(this);
}

FilterDivExp::FilterDivExp(FilterExpressionPtr op1, FilterExpressionPtr op2)
    : operand1(op1), operand2(op2) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterDivExp::~FilterDivExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterDivExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitDivExp(this);
}

FilterAndExp::FilterAndExp(FilterExpressionPtr op1, FilterExpressionPtr op2)
    : operand1(op1), operand2(op2) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterAndExp::~FilterAndExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterAndExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitAndExp(this);
}

FilterIfExp::FilterIfExp(FilterExpressionPtr op1, FilterExpressionPtr op2,
                         FilterExpressionPtr op3)
    : operand1(op1), operand2(op2), operand3(op3) {
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

FilterIfExp::~FilterIfExp() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void FilterIfExp::Accept(FilterExpressionVisitor &visitor) {
  visitor.VisitIfExp(this);
}

void rxdock::to_json(json &j, const FilterExpression &filterExpression) {
  j = json{{"value", filterExpression.value}};
}

void rxdock::from_json(const json &j, FilterExpression &filterExpression) {
  j.at("value").get_to(filterExpression.value);

  // subclasses missing
}