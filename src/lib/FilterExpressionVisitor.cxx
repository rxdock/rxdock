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

#include "FilterExpressionVisitor.h"

#include <loguru.hpp>

#include <cerrno>
#include <climits>
#include <cmath>

using namespace rxdock;

void PrintVisitor::VisitVbleExp(FilterVbleExp *fe) {
  std::string name = fe->GetVble().GetName();
  if (name[0] == 'c') // constant
    strbuf << fe->GetVble().GetValue() << " ";
  else
    strbuf << fe->GetVble().GetName() << " ";
}

void PrettyPrintVisitor::VisitVbleExp(FilterVbleExp *fe) {
  //	  strbuf << "Vble ";
  strbuf << fe->GetVble().GetName() << " ";
}

void EvaluateVisitor::VisitVbleExp(FilterVbleExp *fe) {
  fe->SetValue(fe->GetVble().GetValue());
}

void PrintVisitor::VisitAddExp(FilterAddExp *fe) {
  strbuf << "+ ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitAddExp(FilterAddExp *fe) {
  strbuf << "(";
  fe->GetOp(0)->Accept(*this);
  strbuf << "+ ";
  fe->GetOp(1)->Accept(*this);
  strbuf << ")";
}

void EvaluateVisitor::VisitAddExp(FilterAddExp *fe) {
  fe->GetOp(0)->Accept(*this);
  fe->GetOp(1)->Accept(*this);
  fe->SetValue(fe->GetOp(0)->GetValue() + fe->GetOp(1)->GetValue());
  // strbuf << "AddExp: " << fe->GetOp(0)->GetValue() << "\t" <<
  // fe->GetOp(1)->GetValue() << std::endl;
}

void PrintVisitor::VisitSubExp(FilterSubExp *fe) {
  strbuf << "- ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitSubExp(FilterSubExp *fe) {
  strbuf << "(";
  fe->GetOp(0)->Accept(*this);
  strbuf << "- ";
  fe->GetOp(1)->Accept(*this);
  strbuf << ")";
}

void EvaluateVisitor::VisitSubExp(FilterSubExp *fe) {
  fe->GetOp(0)->Accept(*this);
  fe->GetOp(1)->Accept(*this);
  fe->SetValue(fe->GetOp(0)->GetValue() - fe->GetOp(1)->GetValue());
  // strbuf << "sub " << fe->GetOp(0)->GetValue() <<"\t" <<
  // fe->GetOp(1)->GetValue() << std::endl;
}

void PrintVisitor::VisitMulExp(FilterMulExp *fe) {
  strbuf << "* ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitMulExp(FilterMulExp *fe) {
  strbuf << "(";
  fe->GetOp(0)->Accept(*this);
  strbuf << "* ";
  fe->GetOp(1)->Accept(*this);
  strbuf << ")";
}

void EvaluateVisitor::VisitMulExp(FilterMulExp *fe) {
  fe->GetOp(0)->Accept(*this);
  fe->GetOp(1)->Accept(*this);
  fe->SetValue(fe->GetOp(0)->GetValue() * fe->GetOp(1)->GetValue());
  // strbuf << "Mul " << fe->GetOp(0)->GetValue() <<"\t" <<
  // fe->GetOp(1)->GetValue() << std::endl;
}

void PrintVisitor::VisitDivExp(FilterDivExp *fe) {
  strbuf << "/ ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitDivExp(FilterDivExp *fe) {
  strbuf << "(";
  fe->GetOp(0)->Accept(*this);
  strbuf << "div ";
  fe->GetOp(1)->Accept(*this);
  strbuf << ")";
}

void EvaluateVisitor::VisitDivExp(FilterDivExp *fe) {
  fe->GetOp(1)->Accept(*this);
  fe->GetOp(0)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  ReturnType v1 = fe->GetOp(1)->GetValue();
  // strbuf << "div " << v0 <<"\t" << v1 << std::endl;
  if (std::fabs(v1) < 0.000001)
    fe->SetValue(v0);
  else
    fe->SetValue(v0 / v1);
}

void PrintVisitor::VisitAndExp(FilterAndExp *fe) {
  strbuf << "and ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitAndExp(FilterAndExp *fe) {
  strbuf << "( (";
  fe->GetOp(0)->Accept(*this);
  strbuf << " > 0) and (";
  fe->GetOp(1)->Accept(*this);
  strbuf << " > 0) )";
}

void EvaluateVisitor::VisitAndExp(FilterAndExp *fe) {
  fe->GetOp(1)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  if (v0 <= 0.0)
    fe->SetValue(0.0); // false
  else {
    fe->GetOp(0)->Accept(*this);
    ReturnType v1 = fe->GetOp(1)->GetValue();
    if (v1 <= 0.0)
      fe->SetValue(0.0); // false
    else
      fe->SetValue(1.0); // true
  }
}

void PrintVisitor::VisitLogExp(FilterLogExp *fe) {
  strbuf << "log ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitLogExp(FilterLogExp *fe) {
  strbuf << "(log ";
  fe->GetOp(0)->Accept(*this);
  strbuf << ")";
}

void EvaluateVisitor::VisitLogExp(FilterLogExp *fe) {
  fe->GetOp(0)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  if (std::fabs(v0) < 0.000001)
    fe->SetValue(0.0);
  else
    fe->SetValue(std::log(std::fabs(v0)));
  //  strbuf << "log " << v0 << std::endl;
}

void PrintVisitor::VisitExpExp(FilterExpExp *fe) {
  strbuf << "exp ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitExpExp(FilterExpExp *fe) {
  strbuf << "(exp ";
  fe->GetOp(0)->Accept(*this);
  strbuf << ")";
}

void EvaluateVisitor::VisitExpExp(FilterExpExp *fe) {
  fe->GetOp(0)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  errno = 0;
  if (v0 > 20)

  {
    //			float m = FLOAT_MAX; //numeric_limits<float>::max();
    // strbuf << "greater " << m << "\n";
    fe->SetValue(std::exp(20));
  } else if (v0 < -200)
    fe->SetValue(0.0);
  else
    fe->SetValue(std::exp(v0));
  //  strbuf << "exp " << v0 <<  "\t" << fe->GetValue() <<endl;
}

void PrintVisitor::VisitIfExp(FilterIfExp *fe) {
  strbuf << "if ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitIfExp(FilterIfExp *fe) {
  strbuf << "\nif ";
  fe->GetOp(0)->Accept(*this);
  strbuf << " > 0 then\n\t";
  fe->GetOp(1)->Accept(*this);
  strbuf << "\nelse\n\t";
  fe->GetOp(2)->Accept(*this);
  strbuf << "\nend\n";
}

void EvaluateVisitor::VisitIfExp(FilterIfExp *fe) {
  fe->GetOp(0)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  if (v0 > 0.0) {
    fe->GetOp(1)->Accept(*this);
    fe->SetValue(fe->GetOp(1)->GetValue());
  } else {
    fe->GetOp(2)->Accept(*this);
    fe->SetValue(fe->GetOp(2)->GetValue());
  }
}

FilterExpressionVisitor::~FilterExpressionVisitor() {}

FilterExpressionVisitor::FilterExpressionVisitor() {}

PrintVisitor::PrintVisitor(ContextPtr co) : contextp(co) {}

PrintVisitor::~PrintVisitor() { LOG_F(INFO, "{}", strbuf.str()); }

PrettyPrintVisitor::PrettyPrintVisitor(ContextPtr co) : contextp(co) {}

PrettyPrintVisitor::~PrettyPrintVisitor() { LOG_F(INFO, "{}", strbuf.str()); }

EvaluateVisitor::EvaluateVisitor(ContextPtr co) : contextp(co) {}
