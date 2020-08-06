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

#include <cstdlib>
#include <fstream>

#include "rxdock/MOEGrid.h"

#include "rxdock/ParameterFileSource.h"
#include "rxdock/RealGrid.h"

#include <loguru.hpp>

using namespace rxdock;

const std::string GRID_SIGNATURE = "SVLgrid\n";
const int VERSION = 100;

// utility SVL-related functions that are not really OO and are not strictly
// grid-related only. still the main purpose of these to export numbers into
// format that is readable by MOE. If doubt, use
// $ hexdump -dc test.grd |head -30
// to make sure the generated grid is the same as the "official" MOE version

namespace svl {
void write_int4m(std::ofstream &ofstr, int an_int) {
  char *i_string = (char *)&an_int;
  ofstr.put(i_string[3]);
  ofstr.put(i_string[2]);
  ofstr.put(i_string[1]);
  ofstr.put(i_string[0]);
}

void write_IEEE8m(std::ofstream &ofstr, double a_double) {
  char *d_string = (char *)&a_double;
  ofstr.put(d_string[7]);
  ofstr.put(d_string[6]);
  ofstr.put(d_string[5]);
  ofstr.put(d_string[4]);
  ofstr.put(d_string[3]);
  ofstr.put(d_string[2]);
  ofstr.put(d_string[1]);
  ofstr.put(d_string[0]);
}
} // namespace svl

//////////////////////////////////////////////////////////////////////////////////////
// MOEGridPoint
// singe grid point related methods
//////////////////////////////////////////////////////////////////////////////////////
//
// default ctor
MOEGridPoint::MOEGridPoint() {
  value = 1.0;
  coords.push_back(0.0); // it is a 1D single point
}

// parameterised ctor receives an n-dimensional vector and a value as
// input
MOEGridPoint::MOEGridPoint(std::vector<double> a_coord, double a_value) {
  // copy everything
  std::vector<double>::iterator crditer;
  for (crditer = a_coord.begin(); crditer != a_coord.end(); ++crditer)
    coords.push_back(*crditer);
  value = a_value;
}

// convenience method for Coord
MOEGridPoint::MOEGridPoint(Coord &a_coord, double a_value) {
  coords.push_back(ExtractXCoord(a_coord));
  coords.push_back(ExtractYCoord(a_coord));
  coords.push_back(ExtractZCoord(a_coord));
  value = a_value;
}
//////////////////////////////////////////////////////////////////////////////////////
// MOEGridShape
// SVL grid shape vector related methods
//////////////////////////////////////////////////////////////////////////////////////
// default ctor allocates a 1D grid starting at 0.0 to 1.0
// with 1.0 stepsize (2 elements at least needed for MOE)
MOEGridShape::MOEGridShape() {
  std::vector<double> grid_origin(1, 0.0);
  std::vector<double> grid_extents(1, 1.0);
  double grid_step = 1.0;
  MOEGridShape(grid_origin, grid_extents, grid_step);
}

// parametrized ctor receives grid as the origin vector, the extent
// vector and stepsize
MOEGridShape::MOEGridShape(std::vector<double> grid_origin,
                           std::vector<double> grid_extents,
                           double grid_steps) {
  // grid origin should to be a N-dimensional point with the
  // same dimensions as the extends
  if (grid_origin.size() != grid_extents.size()) {
    LOG_F(WARNING, "Grid origin and extents dimension mismatch.");
  }
  // copying ctor parameter vectors (it is not guaranteed they will
  // live forever during the code)
  std::vector<double>::iterator o_iter,
      e_iter; // for the O_rigin and the E_xtents
  for (o_iter = grid_origin.begin(); o_iter != grid_origin.end(); ++o_iter)
    origin.push_back(*o_iter);
  for (e_iter = grid_extents.begin(); e_iter != grid_extents.end(); ++e_iter)
    extents.push_back(*e_iter);
  step_size = grid_steps;

  // calculate sizes for all vertices of the grid cube
  for (o_iter = origin.begin(), e_iter = extents.begin();
       o_iter != origin.end() && e_iter != extents.end(); ++o_iter, ++e_iter) {
    i_extents.push_back((int)(((*e_iter) - (*o_iter)) / step_size) + 1);
  }
  // calculate the total number of points (initialize multiplication by 1)
  data_size = 1;
  for (std::vector<int>::const_iterator ext_iter = i_extents.begin();
       ext_iter != i_extents.end(); ++ext_iter)
    data_size *= (*ext_iter);
}

//////////////////////////////////////////////////////////////////////////////////////
// MOEGrid
// grid IO methods involving both the data and the shape vectors
//////////////////////////////////////////////////////////////////////////////////////

// default ctor
MOEGrid::MOEGrid() {
  // default stream name (should be overridden by user in command-line)
  stream_name = "moegrid.grd";
}

// parameterised ctor needs the shape and the data vector
MOEGrid::MOEGrid(MOEGridShape &a_shape, MOEGridData &a_data,
                 std::string a_strfName = "moegrid.grd") {
  stream_name = a_strfName;
  myShape = a_shape;
  myData = a_data;
}

// read a single rDock .as /.prm file to get the docking site extents
void MOEGrid::GetDockingSiteExtents(std::string &a_strPrmFile) {
  // Create a bimolecular workspace
  BiMolWorkSpacePtr spWS(new BiMolWorkSpace());
  // Set the workspace name to the root of the receptor .prm filename
  std::vector<std::string> componentList =
      ConvertDelimitedStringToList(a_strPrmFile, ".");
  std::string wsName = componentList.front();
  spWS->SetName(wsName);
  // Read docking site from file
  std::string strDockingSiteFile = spWS->GetName() + ".json";
  std::string strInputFile = GetDataFileName("data/grids", strDockingSiteFile);
  std::ifstream inputFile(strInputFile.c_str());
  json siteData;
  inputFile >> siteData;
  inputFile.close();
  DockingSitePtr dSite(new DockingSite(siteData.at("docking-site")));
  // we need only these two values
  c_min = dSite->GetMinCoord();
  c_max = dSite->GetMaxCoord();
}
// reads all the .prm files and calculates the common extents for all
// cavities to get a grid that includes them all
void MOEGrid::CalculateCommonExtents(std::vector<std::string> strPrmFiles) {
  LOG_F(1, "Receptors:");
  std::vector<std::string>::iterator strPrmFilesIter;
  for (strPrmFilesIter = strPrmFiles.begin();
       strPrmFilesIter != strPrmFiles.end(); ++strPrmFilesIter) {
    LOG_F(1, "{}", *strPrmFilesIter);
    // initialize values with the very first
    if (strPrmFilesIter == strPrmFiles.begin()) {
      GetDockingSiteExtents((*strPrmFilesIter));
      min = c_min;
      max = c_max;
    } // otherwise go through on the rest of the prms and calculate extremes
    else {
      GetDockingSiteExtents((*strPrmFilesIter));
      // change min extents to current values if required
      if (c_min.xyz(0) < min.xyz(0))
        min.xyz(0) = c_min.xyz(0);
      if (c_min.xyz(1) < min.xyz(1))
        min.xyz(1) = c_min.xyz(1);
      if (c_min.xyz(2) < min.xyz(2))
        min.xyz(2) = c_min.xyz(2);
      // change max extents
      if (c_max.xyz(0) > max.xyz(0))
        max.xyz(0) = c_max.xyz(0);
      if (c_max.xyz(1) > max.xyz(1))
        max.xyz(1) = c_max.xyz(1);
      if (c_max.xyz(2) > max.xyz(2))
        max.xyz(2) = c_max.xyz(2);
    }
  }
}

// writes the (binary format) grid
// and returns the number of gridpoints written
// that should be the same as the number of datapoints (a_data.size() above)
// mode stands for std::ios_base::trunc or ios:app
long MOEGrid::WriteGrid(std::ios_base::openmode mode) {
  // file to save (binary, overidden)
  std::ofstream ofstr;
  ofstr.open(stream_name.c_str(),
             std::ios_base::out | mode | std::ios_base::binary);
  long written = 0;                                  // data written
  long ndata = (long)(myData.size());                // data vector size
  long dimension = (long)myShape.GetOrigin().size(); // space dim
  // writing stuff
  // if we are starting a new grid, write signature
  if (std::ios_base::trunc == mode)
    ofstr.write(GRID_SIGNATURE.c_str(), GRID_SIGNATURE.size()); // signature
  svl::write_int4m(ofstr, VERSION);                             // VERSION
  svl::write_int4m(ofstr, ndata);     // size of data vector
  svl::write_int4m(ofstr, dimension); // size of shape vector/origin
  svl::write_int4m(ofstr, 0);         // the number of skipped bytes (is zero)
  std::vector<int> i_ext = myShape.GetIExtents();
  for (std::vector<int>::iterator i_ext_iter = i_ext.begin();
       i_ext_iter != i_ext.end(); ++i_ext_iter) {
    svl::write_int4m(ofstr, (*i_ext_iter)); // write shape vector extents
  }
  std::vector<double> grid_origin = myShape.GetOrigin(); // cat shape :
  std::vector<double> grid_extents =
      myShape.GetExtents(); // list all the values in the shape vector
  std::vector<double>::iterator o_iter, e_iter; // origin and extents
  for (o_iter = grid_origin.begin(), e_iter = grid_extents.begin();
       o_iter != grid_origin.end() && e_iter != grid_extents.end();
       ++o_iter, ++e_iter) {
    for (double grid_coord = (*o_iter); grid_coord <= (*e_iter);
         grid_coord += myShape.GetStepSize()) {
      svl::write_IEEE8m(ofstr, (grid_coord)); // write shape vector data
    }
  }
  // write data vector
  for (MOEGridData::iterator d_iter = myData.begin(); d_iter != myData.end();
       ++d_iter) {
    svl::write_IEEE8m(ofstr, (*d_iter).GetValue());
    ++written;
  }
  ofstr.close();
  return written;
}

long MOEGrid::ReadGrid() {
  std::ifstream ifstr(stream_name.c_str()); // file to read
  return 0;
}
