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

#include "RbtBaseBiMolTransform.h"
#include "RbtChromElement.h"
#include "RbtRand.h"

// Simple class to keep track of Monte Carlo sampling statistics
class RbtMCStats {
public:
  RbtMCStats();
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
typedef SmartPtr<RbtMCStats> RbtMCStatsPtr; // Smart pointer

class RbtSimAnnTransform : public RbtBaseBiMolTransform {
public:
  // Static data member for class type
  static std::string _CT;
  // Parameter names
  static std::string _START_T;
  static std::string _FINAL_T;
  static std::string _BLOCK_LENGTH;
  static std::string _SCALE_CHROM_LENGTH;
  static std::string _NUM_BLOCKS;
  static std::string _STEP_SIZE;
  static std::string _MIN_ACC_RATE;
  static std::string _PARTITION_DIST;
  static std::string _PARTITION_FREQ;
  static std::string _HISTORY_FREQ;

  RBTDLL_EXPORT static std::string &GetStartT();
  RBTDLL_EXPORT static std::string &GetFinalT();
  RBTDLL_EXPORT static std::string &GetBlockLength();
  RBTDLL_EXPORT static std::string &GetNumBlocks();
  RBTDLL_EXPORT static std::string &GetStepSize();
  RBTDLL_EXPORT static std::string &GetPartitionDist();
  RBTDLL_EXPORT static std::string &GetPartitionFreq();

  ////////////////////////////////////////
  // Constructors/destructors
  RBTDLL_EXPORT RbtSimAnnTransform(const std::string &strName = "SIMANN");
  virtual ~RbtSimAnnTransform();

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
  RbtSimAnnTransform(
      const RbtSimAnnTransform &); // Copy constructor disabled by default
  RbtSimAnnTransform &
  operator=(const RbtSimAnnTransform &); // Copy assignment disabled by default

protected:
  ////////////////////////////////////////
  // Protected data
  ////////////////

private:
  ////////////////////////////////////////
  // Private data
  //////////////
  RbtMCStatsPtr m_spStats;
  RbtRand &m_rand; // keep a reference to the singleton random number generator
  RbtRequestPtr m_spPartReq;  // Partitioning request
  RbtChromElementPtr m_chrom; // Current chromosome
  std::vector<double>
      m_minVector; // Chromosome vector corresponding to overall minimum score
  std::vector<double> m_lastGoodVector; // Saved chromosome before each MC
                                        // mutation (to allow revert)
};

// Useful typedefs
typedef SmartPtr<RbtSimAnnTransform> RbtSimAnnTransformPtr; // Smart pointer

#endif //_RBTSIMANNTRANSFORM_H_
