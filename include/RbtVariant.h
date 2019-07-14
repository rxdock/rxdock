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

// Variant class along the lines of the Visual Basic variant data type.
// Handles int, double, const char*, string and std::vector<string> types
// Underlying value is stored as double and std::vector<string>

#ifndef _RBTVARIANT_H_
#define _RBTVARIANT_H_

#include <sstream>

#include "RbtContainers.h"
#include "RbtCoord.h"

// Strings representing TRUE and FALSE
const std::string _TRUE("TRUE");
const std::string _FALSE("FALSE");

class RbtVariant {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  RbtVariant() : m_d(0.0) {}
  RbtVariant(int i) { SetDouble(i); }
  RbtVariant(double d) { SetDouble(d); }
  RbtVariant(const std::string &s) { SetString(s); }
  RbtVariant(const char *c) { SetString(std::string(c)); }
  RbtVariant(const RbtStringList &sl) { SetStringList(sl); }
  RbtVariant(bool b) { SetBool(b); }
  RbtVariant(const RbtCoord &c) { SetCoord(c); }
  // Renders a vector of doubles into comma-separated strings of maxCols in
  // length
  RbtVariant(const RbtDoubleList &dl, int maxCols, int precision) {
    SetDoubleList(dl, maxCols, precision);
  }
  virtual ~RbtVariant() { m_sl.clear(); }

  ///////////////////////////////////////////////
  // Stream functions
  ///////////////////////////////////////////////

  // Insertion operator
  friend std::ostream &operator<<(std::ostream &s, const RbtVariant &v) {
    // Try and guess the most appropriate output format
    // Empty variant
    if (v.isEmpty())
      s << "Undefined variant";
    // String/double
    else if (v.Size() == 1)
      s << v.m_sl.front();
    // String list
    else {
      s << std::endl;
      for (RbtStringListConstIter iter = v.m_sl.begin(); iter != v.m_sl.end();
           iter++)
        s << *iter << std::endl;
    }
    return s;
  }

  ///////////////////////////////////////////////
  // Operator functions:
  /////////////////////

  // Set methods
  RbtVariant &operator=(int i) {
    SetDouble(i);
    return *this;
  }
  RbtVariant &operator=(double d) {
    SetDouble(d);
    return *this;
  }
  RbtVariant &operator=(const std::string &s) {
    SetString(s);
    return *this;
  }
  RbtVariant &operator=(const char *c) {
    SetString(std::string(c));
    return *this;
  }
  RbtVariant &operator=(const RbtStringList &sl) {
    SetStringList(sl);
    return *this;
  }
  RbtVariant &operator=(bool b) {
    SetBool(b);
    return *this;
  }
  RbtVariant &operator=(const RbtCoord &c) {
    SetCoord(c);
    return *this;
  }

  void operator+=(const RbtVariant &v) {
    std::copy(v.m_sl.begin(), v.m_sl.end(), std::back_inserter(m_sl));
    m_d = std::atof(String().c_str());
  }

  // Conversion operators to convert back to basic types
  operator int() const { return int(Double()); }
  operator unsigned int() const { return (unsigned int)(Double()); }
  operator double() const { return Double(); }
  operator std::string() const { return String(); }
  operator RbtStringList() const { return StringList(); }
  operator bool() const { return Bool(); }
  operator RbtCoord() const { return Coord(); }

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Get methods
  double Double() const { return m_d; }
  std::string String() const {
    return m_sl.empty() ? std::string() : m_sl.front();
  }
  RbtStringList StringList() const { return m_sl; }
  bool Bool() const { return m_d != 0.0 || String() == _TRUE; }
  RbtCoord Coord() const {
    RbtCoord c;
    if (m_sl.empty())
      return c;
    else {
      std::istringstream(m_sl.front()) >> c;
      return c;
    }
  }

  unsigned int Size() const { return m_sl.size(); }
  bool isEmpty() const { return m_sl.empty(); }

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  void SetDouble(double d) {
    m_d = d;
    m_sl.clear();
    std::ostringstream ostr;
    // Don't need "ends" with std::ostringstream apparently
    //(was introducing a non-ASCII \0 char into log files
    // ostr << d << ends;
    ostr << d;
    std::string s(ostr.str());
    // delete ostr.str();
    m_sl.push_back(s);
  }

  void SetString(const std::string &s) {
    m_sl.clear();
    if (!s.empty())
      m_sl.push_back(s);
    m_d = std::atof(String().c_str());
  }

  void SetStringList(const RbtStringList &sl) {
    m_sl = sl;
    m_d = std::atof(String().c_str());
  }

  void SetBool(bool b) {
    if (b) {
      SetString(_TRUE);
    } else {
      SetString(_FALSE);
    }
  }

  void SetCoord(const RbtCoord &c) {
    m_d = c.Length();
    m_sl.clear();
    std::ostringstream ostr;
    // Don't need "ends" with std::ostringstream apparently
    //(was introducing a non-ASCII \0 char into log files
    // ostr << c << ends;
    ostr << c;
    std::string s(ostr.str());
    // delete ostr.str();
    m_sl.push_back(s);
  }

  void SetDoubleList(const RbtDoubleList &dl, int maxColumns, int precision) {
    int nValues = dl.size();
    m_d = (nValues > 0) ? dl[0] : 0.0;
    m_sl.clear();
    std::ostringstream ostr;
    ostr.precision(precision);
    ostr.setf(std::ios_base::fixed, std::ios_base::floatfield);
    int lastIndex = nValues - 1;
    for (int i = 0; i < nValues; ++i) {
      ostr << dl[i];
      if ((ostr.tellp() >= maxColumns) || (i == lastIndex)) {
        m_sl.push_back(ostr.str());
        ostr.str("");
      } else {
        ostr << ",";
      }
    }
  }

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  double m_d;
  RbtStringList m_sl;
};

// Useful typedefs
// typedef SharedPtr<RbtVariant> RbtVariantPtr;//Smart pointer
typedef std::vector<RbtVariant> RbtVariantList; // Vector of variants
typedef RbtVariantList::iterator RbtVariantListIter;
typedef RbtVariantList::const_iterator RbtVariantListConstIter;

// Map of (key=string, value=variant)
typedef std::map<std::string, RbtVariant> RbtStringVariantMap;
typedef RbtStringVariantMap::iterator RbtStringVariantMapIter;
typedef RbtStringVariantMap::const_iterator RbtStringVariantMapConstIter;

#endif //_RBTVARIANT_H_
