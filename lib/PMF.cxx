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

#include "PMF.h"

using namespace rxdock;

std::string rxdock::PMFType2Str(PMFType aType) {
  switch (aType) {
  case CF:
    return std::string("CF");
    break;
  case CP:
    return std::string("CP");
    break;
  case cF:
    return std::string("cF");
    break;
  case cP:
    return std::string("cP");
    break;
  case C3:
    return std::string("C3");
    break;
  case CW:
    return std::string("CW");
    break;
  case CO:
    return std::string("CO");
    break;
  case CN:
    return std::string("CN");
    break;
  case C0:
    return std::string("C0");
    break;
  case NC:
    return std::string("NC");
    break;
  case NP:
    return std::string("NP");
    break;
  case NA:
    return std::string("NA");
    break;
  case ND:
    return std::string("ND");
    break;
  case NR:
    return std::string("NR");
    break;
  case N0:
    return std::string("N0");
    break;
  case NS:
    return std::string("NS");
    break;
  case OC:
    return std::string("OC");
    break;
  case OA:
    return std::string("OA");
    break;
  case OE:
    return std::string("OE");
    break;
  case OR:
    return std::string("OR");
    break;
  case OS:
    return std::string("OS");
    break;
  case OD:
    return std::string("OD");
    break;
  case OW:
    return std::string("OW");
    break;
  case P:
    return std::string("P");
    break;
  case SA:
    return std::string("SA");
    break;
  case SD:
    return std::string("SD");
    break;
  case HL:
    return std::string("HL");
    break;
  case HH:
    return std::string("HH");
    break;
  case Zn:
    return std::string("Zn");
    break;
  case CL:
    return std::string("CL");
    break;
  case Mn:
    return std::string("Mn");
    break;
  case Mg:
    return std::string("Mg");
    break;
  case F:
    return std::string("F");
    break;
  case Fe:
    return std::string("Fe");
    break;
  case Br:
    return std::string("Br");
    break;
  case V:
    return std::string("V");
    break;
  case PMF_UNDEFINED:
    return std::string("PMF_UNDEFINED");
    break;
  default:
    return std::string("NO SUCH TYPE");
    break;
  }
}

PMFType rxdock::PMFStr2Type(std::string aStr) {
  if (!aStr.compare("CF"))
    return CF;
  else if (!aStr.compare("CP"))
    return CP;
  else if (!aStr.compare("cF"))
    return cF;
  else if (!aStr.compare("cP"))
    return cP;
  else if (!aStr.compare("C3"))
    return C3;
  else if (!aStr.compare("CW"))
    return CW;
  else if (!aStr.compare("CO"))
    return CO;
  else if (!aStr.compare("CN"))
    return CN;
  else if (!aStr.compare("C0"))
    return C0;
  else if (!aStr.compare("NC"))
    return NC;
  else if (!aStr.compare("NP"))
    return NP;
  else if (!aStr.compare("NA"))
    return NA;
  else if (!aStr.compare("ND"))
    return ND;
  else if (!aStr.compare("NR"))
    return NR;
  else if (!aStr.compare("N0"))
    return N0;
  else if (!aStr.compare("NS"))
    return NS;
  else if (!aStr.compare("OC"))
    return OC;
  else if (!aStr.compare("OA"))
    return OA;
  else if (!aStr.compare("OE"))
    return OE;
  else if (!aStr.compare("OR"))
    return OR;
  else if (!aStr.compare("OS"))
    return OS;
  else if (!aStr.compare("OD"))
    return OD;
  else if (!aStr.compare("OW"))
    return OW;
  else if (!aStr.compare("P_") || !aStr.compare("P"))
    return P;
  else if (!aStr.compare("SA"))
    return SA;
  else if (!aStr.compare("SD"))
    return SD;
  else if (!aStr.compare("HL"))
    return HL;
  else if (!aStr.compare("HH"))
    return HH;
  else if (!aStr.compare("Zn"))
    return Zn;
  else if (!aStr.compare("CL"))
    return CL;
  else if (!aStr.compare("Mn"))
    return Mn;
  else if (!aStr.compare("Mg"))
    return Mg;
  else if (!aStr.compare("F_") || !aStr.compare("F"))
    return F;
  else if (!aStr.compare("Fe"))
    return Fe;
  else if (!aStr.compare("Br"))
    return Br;
  else if (!aStr.compare("V_") || !aStr.compare("V"))
    return V;

  else
    return PMF_UNDEFINED;
}
