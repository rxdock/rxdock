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

#include <fstream>

#include "RbtFileError.h"
#include "RbtVdwGridSF.h"

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  // Command line arguments and default values
  std::string strInputFile;
  std::string strOutputFile("insight.grid");
  int iGrid = 999;

  // Brief help message
  if (argc == 1) {
    std::cout
        << std::endl
        << "rbconvgrid - converts RbtVdwGridSF binary grid file to InsightII "
           "ascii grid file"
        << std::endl;
    std::cout
        << std::endl
        << "Usage:\trbconvgrid -i<InputFile> [-o<OutputFile>] [-n<GridNum>]"
        << std::endl;
    std::cout
        << std::endl
        << "Options:\t-i<InputFile> - input RbtVdwGridSF binary grid filename"
        << std::endl;
    std::cout << "\t\t-o<OutputFile> - output InsightII ascii grid filename"
              << std::endl;
    std::cout
        << "\t\t-n<GridNum> - grid number to convert (default = list grids)"
        << std::endl;
    return 1;
  }

  // Check command line arguments
  std::cout << std::endl << "Command line args:" << std::endl;
  for (int iarg = 1; iarg < argc; iarg++) {
    std::cout << argv[iarg];
    std::string strArg(argv[iarg]);
    if (strArg.find("-i") == 0)
      strInputFile = strArg.substr(2);
    else if (strArg.find("-o") == 0)
      strOutputFile = strArg.substr(2);
    else if (strArg.find("-n") == 0) {
      std::string strGridNum = strArg.substr(2);
      iGrid = std::atoi(strGridNum.c_str());
    } else {
      std::cout << " ** INVALID ARGUMENT" << std::endl;
      return 1;
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;

  try {
    // Read the grid file header
#if defined(__sgi) && !defined(__GNUC__)
    std::ifstream istr(strInputFile.c_str(), std::ios_base::in);
#else
    std::ifstream istr(strInputFile.c_str(),
                       std::ios_base::in | std::ios_base::binary);
#endif
    if (istr) {
      std::cout << strInputFile << " opened OK" << std::endl;
    }
    // Read header string
    int length;
    Rbt::ReadWithThrow(istr, (char *)&length, sizeof(length));
    char *header = new char[length + 1];
    Rbt::ReadWithThrow(istr, header, length);
    // Add null character to end of string
    header[length] = '\0';
    // Compare title with
    bool match = (RbtVdwGridSF::_CT == header);
    delete[] header;
    if (!match) {
      throw RbtFileParseError(_WHERE_,
                              "Invalid title string in " + strInputFile);
    }

    // Skip the appropriate number of grids
    int nGrids;
    Rbt::ReadWithThrow(istr, (char *)&nGrids, sizeof(nGrids));
    std::cout << "File contains " << nGrids << " grids..." << std::endl;
    if ((iGrid > nGrids) || (iGrid < 1)) {
      std::cout << "Listing grids..." << std::endl;
    } else {
      std::cout << "Locating grid# " << iGrid << "..." << std::endl;
    }
    RbtRealGridPtr spGrid;
    for (int i = 1; (i <= nGrids) && (i <= iGrid); i++) {
      // Read the atom type string
      Rbt::ReadWithThrow(istr, (char *)&length, sizeof(length));
      char *szType = new char[length + 1];
      Rbt::ReadWithThrow(istr, szType, length);
      // Add null character to end of string
      szType[length] = '\0';
      std::string strType(szType);
      delete[] szType;
      RbtTriposAtomType triposType;
      RbtTriposAtomType::eType aType = triposType.Str2Type(strType);
      std::cout << "Grid# " << i << "\t"
                << "atom type=" << strType << " (type #" << aType << ")"
                << std::endl;
      spGrid = RbtRealGridPtr(new RbtRealGrid(istr));
    }
    istr.close();
    // If we are not in listing mode, write the grid
    if ((iGrid <= nGrids) && (iGrid >= 1)) {
      std::cout << "Writing grid# " << iGrid << " to " << strOutputFile << "..."
                << std::endl;
      std::ofstream ostr(strOutputFile.c_str());
      spGrid->PrintInsightGrid(ostr);
      ostr.close();
    }
  } catch (RbtError &e) {
    std::cout << e << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
