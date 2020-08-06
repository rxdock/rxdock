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

// MOE grid description at ${MOE}/html/svl/manual/gridman.htm
// here you can read:
// [...]
// Creating a grid consists of creating both a shape vector and a data vector.
// The shape vector stores gridpoint positions along each dimension. The number
// of dimensions must be between 2 and 24 with each dimension containing at
// least 2 points. Each element of the shape vector is a flat vector of strictly
// increasing grid positions along the corresponding dimension. The values of
// the field at the gridpoints are stored in row-major order, i.e., with the
// index of the last dimension changing most rapidly.
// [...]
// Generating a shape vector consists of creating coordinate vectors along each
// of the n dimensions of an n-dimensional lattice. The intervals between the
// coordinates need not be regular. To generate the shape vector, you need the
// origin of the lattice and the length of each side of the lattice, along with
// the coordinate spacings in each dimension.
// [...]
// The values of the field at the gridpoints are stored in a data vector. There
// must be one data value for every gridpoint. The length of the data vector
// must therefore equal the product of the number of gridpoints along each
// dimension.

#ifndef _RBTMOEGRID_H_
#define _RBTMOEGRID_H_

#include "rxdock/Coord.h"
#include "rxdock/Rbt.h"
#include <cstdlib>
#include <fstream>
#include <vector>

#include "rxdock/BiMolWorkSpace.h"

namespace rxdock {

// single point in a multi-dimensional grid
// the coords contains N values where N is the dimension of the grid
class MOEGridPoint {
  std::vector<double> coords;
  double value;

public:
  // default ctor with a 1D gridpoint at 0.0 value 1.0
  // see also default ctor for the shape and data vector below
  RBTDLL_EXPORT MOEGridPoint();
  // parametrized ctor when supplying data as a vector + value pair
  MOEGridPoint(std::vector<double>, double);
  // parametrized ctor getting coord values from Coord (default 3D
  // grid for convenience)
  MOEGridPoint(Coord &, double);
  // data access methods
  double GetValue() { return value; }
  const std::vector<double> &GetCoords() { return coords; }
  void SetCoords(double a_val) { coords.push_back(a_val); }
  void SetValue(double a_val) { value = a_val; }
};

// MOE grid shape vector
class MOEGridShape {
  // data that is provided usually in the SVL code
  // a shape vector is defined by
  //		i)		starting coordinates
  //		ii)		X,Y,Z extension (in absolute coordinates)
  //		iii)	grid stepsize
  // data vector defined by
  //		i)		a 1D vector containing the same number
  //				of elements than the total number of the
  //				shape vector

  std::vector<double> origin;  // the starting point
  std::vector<double> extents; // the final points
  std::vector<int> i_extents;  // extents in integers: number of
                               // points at the edges (ie 20 X 30 X 40 grid)
  double step_size;            // grid step size
  long data_size;              // total number of points

public:
  // default ctor allocates a 1D grid starting at 0.0 to 1.0
  // with 1.0 stepsize (2 elements at least needed for MOE)
  MOEGridShape();

  // parametrized ctor receives grid as the origin vector, the extent
  // vector and stepsize
  RBTDLL_EXPORT MOEGridShape(std::vector<double>, std::vector<double>, double);

  const std::vector<double> &GetOrigin() { return origin; }
  const std::vector<double> &GetExtents() { return extents; }
  const std::vector<int> &GetIExtents() { return i_extents; }
  double GetStepSize() { return step_size; }
  long GetDataSize() { return data_size; }
};

typedef std::vector<MOEGridPoint> MOEGridData;

// MOE grid data class for IO
//
// the MOE grid file format (in the present stage at 06/03/2003)
class MOEGrid {
  MOEGridShape myShape;
  MOEGridData myData;
  std::string stream_name; // filename
  Coord min, max;          // common extent min and max coords
  Coord c_min, c_max;      // current extent min and max coords
public:
  // default ctor
  RBTDLL_EXPORT MOEGrid();
  // parameterised ctor for a single grid
  // pass shape, data and outfilename
  MOEGrid(MOEGridShape &, MOEGridData &, std::string);
  // or pass parameters using these functions (when it is about
  // to save multiple grids into a single file)
  void SetShape(MOEGridShape &a_shape) { myShape = a_shape; }
  void SetData(MOEGridData &a_data) { myData = a_data; }
  void SetOutputFileName(std::string &a_filename) { stream_name = a_filename; }
  // methods to calculate grid extents
  RBTDLL_EXPORT void
  CalculateCommonExtents(std::vector<std::string> strPrmFiles); // for all
  void GetDockingSiteExtents(std::string &a_strPrmFile);        // for a single
  // get methods
  Coord GetMinExtents() { return min; }
  Coord GetMaxExtents() { return max; }

  // iofstream methods
  RBTDLL_EXPORT long WriteGrid(
      std::ios_base::openmode mode); // returns the number of points written
  long ReadGrid();                   // ditto (fields read)
};

// the relevant SVL code in gridshow.svl
//
// basic double should be the IEEE 8 byte double MSB first in Linux
// on big-endian try swab with dd and/or cpio
// MOE grid format:
//
// Note: the internal grid format (subject to change without notice!)
// Header for the grid-file only:
//  GRID_SIGNATURE          8 bytes
// Each grid:
//  version number          	int4m
//  ndata (= length data)       int4m
//  ndim (= length size)        int4m
//  skip                int4m
//  ... (reserved for future use)   byte * skip
//  size (= app length shape)   int4m * ndim
//  shape               ieee8m * (add size)
//  data                ieee8m * ndata
//
// On reading, we check that ndata, ndim, and size are 24-bit non-negatives,
// but we don't check "ndata == (mul size)"
//
// snippets from '$MOE/lib/svl/run/gridshow.svl'
// const GRID_SIGNATURE = "SVLgrid\n";
// const VERSION = 100;
// const SKIP = 0;
//
// Check file signature and return non-zero if the file of file name fn
// is a "grid file", i.e. a special file that contains exactly one grid
// and nothing else
//
// local function read_gridsign fn =
//    freadb [fn, 'char', length GRID_SIGNATURE] === GRID_SIGNATURE;

// Write grid to a file of given filename or filehandle
//
// local function writeb_fnum_grid [fn, [data, shape]]
//    local ndata = length data;
//    local ndim = length shape;
//    local size = app length shape;
//
//    fwriteb [fn, 'int4m', [VERSION, ndata, ndim, SKIP]];
//    fwriteb [fn, 'char', rep ["\0", SKIP]];
//    fwriteb [fn, 'int4m', size];
//    fwriteb [fn, 'IEEE8m', cat shape];
//    fwriteb [fn, 'IEEE8m', data];
// endfunction
//
// function grid_fwriteb [fn, grid]
//    if type fn === 'num' then writeb_fnum_grid [fn, grid]; return; endif
//    fn = fopenw fn;
//    writeb_fnum_grid [fn, grid];
//    fclose fn;
// endfunction
//
// Create/rewrite "grid-file" of given name that willcontain a given grid
//
// function grid_writefile [fn, grid]
//    fn = fopenw fn;
//    fwriteb [fn, 'char', GRID_SIGNATURE];
//    writeb_fnum_grid [fn, grid];
//    fclose fn;
// endfunction

//    fwriteb [fn, 'char', GRID_SIGNATURE];
//    fwriteb [fn, 'int4m', [VERSION, ndata, ndim, SKIP]];
//    fwriteb [fn, 'char', rep ["\0", SKIP]];
//    fwriteb [fn, 'int4m', size];
//    fwriteb [fn, 'IEEE8m', cat shape];
//    fwriteb [fn, 'IEEE8m', data];

} // namespace rxdock

#endif //_RBTMOEGRID_H_
