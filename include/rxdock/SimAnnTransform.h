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

// Simple, non-adaptive simulated annealing protocol
#ifndef _RBTSIMANNTRANSFORM_H_
#define _RBTSIMANNTRANSFORM_H_

#include "rxdock/BaseBiMolTransform.h"
#include "rxdock/ChromElement.h"
#include "rxdock/Rand.h"

namespace rxdock {

// Simple class to keep track of Monte Carlo sampling statistics
class MCStats {
public:
  MCStats();
  void Init(double score);
  void InitBlock(double score);
  void Accumulate(double score, bool bAccepted);
  double Mean() const;
  double Variance() const;
  double AccRate() const;
  double _total;
  double _total2;
  double _blockInitial;
  double _blockFinal;
  double _blockMin;
  double _blockMax;
  double _initial;
  double _final;
  double _min;
  double _max;
  int _steps;
  int _accepted;
};
typedef SmartPtr<MCStats> MCStatsPtr; // Smart pointer

class SimAnnTransform : public BaseBiMolTransform {
public:
  // Static data member for class type
  static const std::string _CT;
  // Parameter names
  static const std::string _START_T;
  static const std::string _FINAL_T;
  static const std::string _BLOCK_LENGTH;
  static const std::string _SCALE_CHROM_LENGTH;
  static const std::string _NUM_BLOCKS;
  static const std::string _STEP_SIZE;
  static const std::string _MIN_ACC_RATE;
  static const std::string _PARTITION_DIST;
  static const std::string _PARTITION_FREQ;
  static const std::string _HISTORY_FREQ;

  RBTDLL_EXPORT static const std::string &GetStartT();
  RBTDLL_EXPORT static const std::string &GetFinalT();
  RBTDLL_EXPORT static const std::string &GetBlockLength();
  RBTDLL_EXPORT static const std::string &GetNumBlocks();
  RBTDLL_EXPORT static const std::string &GetStepSize();
  RBTDLL_EXPORT static const std::string &GetPartitionDist();
  RBTDLL_EXPORT static const std::string &GetPartitionFreq();

  ////////////////////////////////////////
  // Constructors/destructors
  RBTDLL_EXPORT SimAnnTransform(const std::string &strName = "SIMANN");
  virtual ~SimAnnTransform();

  ////////////////////////////////////////
  // Public methods
  ////////////////

protected:
  ////////////////////////////////////////
  // Protected methods
  ///////////////////
  virtual void SetupReceptor(); // Called by Update when receptor is changed
  virtual void SetupLigand();   // Called by Update when ligand is changed
  virtual void
  SetupTransform(); // Called by Update when either model has changed
  void MC(double t, int blockLen, double stepSize);
  virtual void Execute();

private:
  ////////////////////////////////////////
  // Private methods
  /////////////////
  SimAnnTransform(
      const SimAnnTransform &); // Copy constructor disabled by default
  SimAnnTransform &
  operator=(const SimAnnTransform &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  MCStatsPtr m_spStats;
  Rand &m_rand; // keep a reference to the singleton random number generator
  RequestPtr m_spPartReq;  // Partitioning request
  ChromElementPtr m_chrom; // Current chromosome
  std::vector<double>
      m_minVector; // Chromosome vector corresponding to overall minimum score
  std::vector<double> m_lastGoodVector; // Saved chromosome before each MC
                                        // mutation (to allow revert)
};

// Useful typedefs
typedef SmartPtr<SimAnnTransform> SimAnnTransformPtr; // Smart pointer

} // namespace rxdock

#endif //_RBTSIMANNTRANSFORM_H_
