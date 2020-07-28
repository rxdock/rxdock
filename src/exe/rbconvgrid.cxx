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

#include "FileError.h"
#include "VdwGridSF.h"

using namespace rxdock;

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  // Command line arguments and default values
  std::string strInputFile;
  std::string strOutputFile("insight.grid");
  int iGrid = 999;

  // Strip off the path to the executable, leaving just the file name
  std::string strExeName(argv[0]);
  std::string::size_type i = strExeName.rfind("/");
  if (i != std::string::npos)
    strExeName.erase(0, i + 1);

  // Print a standard header
  PrintStdHeader(std::cout, strExeName);

  // Brief help message
  if (argc == 1) {
    std::cout
        << std::endl
        << "rbconvgrid - converts VdwGridSF binary grid file to InsightII "
           "ascii grid file"
        << std::endl;
    std::cout
        << std::endl
        << "Usage:\trbconvgrid -i<InputFile> [-o<OutputFile>] [-n<GridNum>]"
        << std::endl;
    std::cout << std::endl
              << "Options:\t-i<InputFile> - input VdwGridSF TOML grid filename"
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
    // Read the grid file
    std::ifstream inputFile(strInputFile);
    json vdwGrids;
    inputFile >> vdwGrids;
    inputFile.close();

    // Skip the appropriate number of grids
    std::cout << "File contains " << vdwGrids.size() << " grids..."
              << std::endl;
    if ((iGrid > vdwGrids.size()) || (iGrid < 1)) {
      std::cout << "Listing grids..." << std::endl;
    } else {
      std::cout << "Locating grid# " << iGrid << "..." << std::endl;
    }
    RealGridPtr spGrid;
    for (int i = 1; (i <= vdwGrids.size()) && (i <= iGrid); i++) {
      // Read the atom type string
      std::string strType;
      vdwGrids.at(i).at("vdw-grid").at("tripos-type").get_to(strType);
      TriposAtomType triposType;
      TriposAtomType::eType aType = triposType.Str2Type(strType);
      std::cout << "Grid# " << i << "\t"
                << "atom type=" << strType << " (type #" << aType << ")"
                << std::endl;
      spGrid = RealGridPtr(
          new RealGrid(vdwGrids.at(i).at("vdw-grid").at("real-grid")));
    }
    // If we are not in listing mode, write the grid
    if ((iGrid <= vdwGrids.size()) && (iGrid >= 1)) {
      std::cout << "Writing grid# " << iGrid << " to " << strOutputFile << "..."
                << std::endl;
      std::ofstream ostr(strOutputFile.c_str());
      spGrid->PrintInsightGrid(ostr);
      ostr.close();
    }
  } catch (Error &e) {
    std::cout << e.what() << std::endl;
  } catch (...) {
    std::cout << "Unknown exception" << std::endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(std::cout)

  return 0;
}
