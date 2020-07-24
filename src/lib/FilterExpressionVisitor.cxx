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
#include <cerrno>
#include <climits>
#include <cmath>

using namespace rxdock;

void PrintVisitor::VisitVbleExp(FilterVbleExp *fe) {
  std::string name = fe->GetVble().GetName();
  if (name[0] == 'c') // constant
    std::cout << fe->GetVble().GetValue() << " ";
  else
    std::cout << fe->GetVble().GetName() << " ";
}

void PrettyPrintVisitor::VisitVbleExp(FilterVbleExp *fe) {
  //	  std::cout << "Vble ";
  std::cout << fe->GetVble().GetName() << " ";
}

void EvaluateVisitor::VisitVbleExp(FilterVbleExp *fe) {
  fe->SetValue(fe->GetVble().GetValue());
}

void PrintVisitor::VisitAddExp(FilterAddExp *fe) {
  std::cout << "+ ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitAddExp(FilterAddExp *fe) {
  std::cout << "(";
  fe->GetOp(0)->Accept(*this);
  std::cout << "+ ";
  fe->GetOp(1)->Accept(*this);
  std::cout << ")";
}

void EvaluateVisitor::VisitAddExp(FilterAddExp *fe) {
  fe->GetOp(0)->Accept(*this);
  fe->GetOp(1)->Accept(*this);
  fe->SetValue(fe->GetOp(0)->GetValue() + fe->GetOp(1)->GetValue());
  // std::cout << "AddExp: " << fe->GetOp(0)->GetValue() << "\t" <<
  // fe->GetOp(1)->GetValue() << std::endl;
}

void PrintVisitor::VisitSubExp(FilterSubExp *fe) {
  std::cout << "- ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitSubExp(FilterSubExp *fe) {
  std::cout << "(";
  fe->GetOp(0)->Accept(*this);
  std::cout << "- ";
  fe->GetOp(1)->Accept(*this);
  std::cout << ")";
}

void EvaluateVisitor::VisitSubExp(FilterSubExp *fe) {
  fe->GetOp(0)->Accept(*this);
  fe->GetOp(1)->Accept(*this);
  fe->SetValue(fe->GetOp(0)->GetValue() - fe->GetOp(1)->GetValue());
  // std::cout << "sub " << fe->GetOp(0)->GetValue() <<"\t" <<
  // fe->GetOp(1)->GetValue() << std::endl;
}

void PrintVisitor::VisitMulExp(FilterMulExp *fe) {
  std::cout << "* ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitMulExp(FilterMulExp *fe) {
  std::cout << "(";
  fe->GetOp(0)->Accept(*this);
  std::cout << "* ";
  fe->GetOp(1)->Accept(*this);
  std::cout << ")";
}

void EvaluateVisitor::VisitMulExp(FilterMulExp *fe) {
  fe->GetOp(0)->Accept(*this);
  fe->GetOp(1)->Accept(*this);
  fe->SetValue(fe->GetOp(0)->GetValue() * fe->GetOp(1)->GetValue());
  // std::cout << "Mul " << fe->GetOp(0)->GetValue() <<"\t" <<
  // fe->GetOp(1)->GetValue() << std::endl;
}

void PrintVisitor::VisitDivExp(FilterDivExp *fe) {
  std::cout << "/ ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitDivExp(FilterDivExp *fe) {
  std::cout << "(";
  fe->GetOp(0)->Accept(*this);
  std::cout << "div ";
  fe->GetOp(1)->Accept(*this);
  std::cout << ")";
}

void EvaluateVisitor::VisitDivExp(FilterDivExp *fe) {
  fe->GetOp(1)->Accept(*this);
  fe->GetOp(0)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  ReturnType v1 = fe->GetOp(1)->GetValue();
  // std::cout << "div " << v0 <<"\t" << v1 << std::endl;
  if (std::fabs(v1) < 0.000001)
    fe->SetValue(v0);
  else
    fe->SetValue(v0 / v1);
}

void PrintVisitor::VisitAndExp(FilterAndExp *fe) {
  std::cout << "and ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitAndExp(FilterAndExp *fe) {
  std::cout << "( (";
  fe->GetOp(0)->Accept(*this);
  std::cout << " > 0) and (";
  fe->GetOp(1)->Accept(*this);
  std::cout << " > 0) )";
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
  std::cout << "log ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitLogExp(FilterLogExp *fe) {
  std::cout << "(log ";
  fe->GetOp(0)->Accept(*this);
  std::cout << ")";
}

void EvaluateVisitor::VisitLogExp(FilterLogExp *fe) {
  fe->GetOp(0)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  if (std::fabs(v0) < 0.000001)
    fe->SetValue(0.0);
  else
    fe->SetValue(std::log(std::fabs(v0)));
  //  std::cout << "log " << v0 << std::endl;
}

void PrintVisitor::VisitExpExp(FilterExpExp *fe) {
  std::cout << "exp ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitExpExp(FilterExpExp *fe) {
  std::cout << "(exp ";
  fe->GetOp(0)->Accept(*this);
  std::cout << ")";
}

void EvaluateVisitor::VisitExpExp(FilterExpExp *fe) {
  fe->GetOp(0)->Accept(*this);
  ReturnType v0 = fe->GetOp(0)->GetValue();
  errno = 0;
  if (v0 > 20)

  {
    //			float m = FLOAT_MAX; //numeric_limits<float>::max();
    // std::cout << "greater " << m << "\n";
    fe->SetValue(std::exp(20));
  } else if (v0 < -200)
    fe->SetValue(0.0);
  else
    fe->SetValue(std::exp(v0));
  //  std::cout << "exp " << v0 <<  "\t" << fe->GetValue() <<endl;
}

void PrintVisitor::VisitIfExp(FilterIfExp *fe) {
  std::cout << "if ";
  for (int i = 0; i < fe->GetNOps(); i++)
    fe->GetOp(i)->Accept(*this);
}

void PrettyPrintVisitor::VisitIfExp(FilterIfExp *fe) {
  std::cout << "\nif ";
  fe->GetOp(0)->Accept(*this);
  std::cout << " > 0 then\n\t";
  fe->GetOp(1)->Accept(*this);
  std::cout << "\nelse\n\t";
  fe->GetOp(2)->Accept(*this);
  std::cout << "\nend\n";
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

PrettyPrintVisitor::PrettyPrintVisitor(ContextPtr co) : contextp(co) {}

EvaluateVisitor::EvaluateVisitor(ContextPtr co) : contextp(co) {}
