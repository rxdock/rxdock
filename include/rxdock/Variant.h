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

#include "rxdock/Coord.h"

namespace rxdock {

// Strings representing TRUE and FALSE
const std::string _TRUE("TRUE");
const std::string _FALSE("FALSE");

class Variant {
public:
  ////////////////////////////////////////
  // Constructors/destructors
  Variant() : m_d(0.0) {}
  Variant(int i) { SetDouble(i); }
  Variant(double d) { SetDouble(d); }
  Variant(const std::string &s) { SetString(s); }
  Variant(const char *c) { SetString(std::string(c)); }
  Variant(const std::vector<std::string> &sl) { SetStringList(sl); }
  Variant(bool b) { SetBool(b); }
  Variant(const Coord &c) { SetCoord(c); }
  // Renders a vector of doubles into comma-separated strings of maxCols in
  // length
  Variant(const std::vector<double> &dl, int maxCols, int precision) {
    SetDoubleList(dl, maxCols, precision);
  }
  virtual ~Variant() { m_sl.clear(); }

  ///////////////////////////////////////////////
  // Stream functions
  ///////////////////////////////////////////////

  // Insertion operator
  friend std::ostream &operator<<(std::ostream &s, const Variant &v) {
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
      for (std::vector<std::string>::const_iterator iter = v.m_sl.begin();
           iter != v.m_sl.end(); iter++)
        s << *iter << std::endl;
    }
    return s;
  }

  ///////////////////////////////////////////////
  // Operator functions:
  /////////////////////

  // Set methods
  Variant &operator=(int i) {
    SetDouble(i);
    return *this;
  }
  Variant &operator=(double d) {
    SetDouble(d);
    return *this;
  }
  Variant &operator=(const std::string &s) {
    SetString(s);
    return *this;
  }
  Variant &operator=(const char *c) {
    SetString(std::string(c));
    return *this;
  }
  Variant &operator=(const std::vector<std::string> &sl) {
    SetStringList(sl);
    return *this;
  }
  Variant &operator=(bool b) {
    SetBool(b);
    return *this;
  }
  Variant &operator=(const Coord &c) {
    SetCoord(c);
    return *this;
  }

  void operator+=(const Variant &v) {
    std::copy(v.m_sl.begin(), v.m_sl.end(), std::back_inserter(m_sl));
    m_d = std::atof(GetString().c_str());
  }

  // Conversion operators to convert back to basic types
  operator int() const { return int(GetDouble()); }
  operator unsigned int() const { return (unsigned int)(GetDouble()); }
  operator double() const { return GetDouble(); }
  operator std::string() const { return GetString(); }
  operator std::vector<std::string>() const { return GetStringList(); }
  operator bool() const { return GetBool(); }
  operator Coord() const { return GetCoord(); }

  friend void to_json(json &j, const Variant &v) {
    j = json{{"double", v.m_d}, {"string-list", v.m_sl}};
  }

  friend void from_json(const json &j, Variant &v) {
    j.at("double").get_to(v.m_d);
    j.at("string-list").get_to(v.m_sl);
  }

  ////////////////////////////////////////
  // Public methods
  ////////////////

  // Get methods
  double GetDouble() const { return m_d; }
  std::string GetString() const {
    return m_sl.empty() ? std::string() : m_sl.front();
  }
  std::vector<std::string> GetStringList() const { return m_sl; }
  bool GetBool() const {
    return m_d != 0.0 || GetString() == _TRUE || GetString() == "true";
  }
  Coord GetCoord() const {
    Coord c;
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
    m_d = std::atof(GetString().c_str());
  }

  void SetStringList(const std::vector<std::string> &sl) {
    m_sl = sl;
    m_d = std::atof(GetString().c_str());
  }

  void SetBool(bool b) {
    if (b) {
      SetString(_TRUE);
    } else {
      SetString(_FALSE);
    }
  }

  void SetCoord(const Coord &c) {
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

  void SetDoubleList(const std::vector<double> &dl, int maxColumns,
                     int precision) {
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
  std::vector<std::string> m_sl;
};

// Useful typedefs
// typedef SharedPtr<Variant> VariantPtr;//Smart pointer
typedef std::vector<Variant> VariantList; // Vector of variants
typedef VariantList::iterator VariantListIter;
typedef VariantList::const_iterator VariantListConstIter;

// Map of (key=string, value=variant)
typedef std::map<std::string, Variant> StringVariantMap;
typedef StringVariantMap::iterator StringVariantMapIter;
typedef StringVariantMap::const_iterator StringVariantMapConstIter;

} // namespace rxdock

#endif //_RBTVARIANT_H_
