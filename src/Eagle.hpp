/*
   This file is part of the Eagle haplotype phasing software package
   developed by Po-Ru Loh.  Copyright (C) 2015-2018 Harvard University.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EAGLE_HPP
#define EAGLE_HPP

//#define RDTSC_TIMING

#include <vector>
#include <string>
#include <set>
#include <unordered_map>

#include "Types.hpp"
#include "GenoData.hpp"
#include "HapHedge.hpp"
#include "StaticMultimap.hpp"

namespace EAGLE {

  struct Match {
    uint n, m64jStart, m64jEnd, m64jStartCons, m64jEndCons; // [m64jStart..m64jEnd]
    double logBF, cMlenInit;
    bool operator < (const Match &m2) const {
      return logBF > m2.logBF;
    }
    Match() : n(0), m64jStart(0), m64jEnd(0), m64jStartCons(0), m64jEndCons(1<<30), logBF(0) {}
    Match(uint _n, uint _m64jStart, uint _m64jEnd, double _logBF) :
      n(_n), m64jStart(_m64jStart), m64jEnd(_m64jEnd), m64jStartCons(0), m64jEndCons(1<<30), 
      logBF(_logBF) {}
    static bool greaterEnd(const Match &m1, const Match &m2) {
      return m1.m64jEnd > m2.m64jEnd || (m1.m64jEnd == m2.m64jEnd && m1.logBF > m2.logBF);
    }
    static bool greaterLen(const Match &m1, const Match &m2) {
      return m1.cMlenInit > m2.cMlenInit ||
	(m1.cMlenInit == m2.cMlenInit && m1.logBF > m2.logBF);
    }
  };

  struct DPState {
    uint score;
    std::pair <uint, uint> from;
    DPState() : score(0), from(std::make_pair(0U, 0U)) {}
    DPState(uint _score, std::pair <uint, uint> _from) : score(_score), from(_from) {}
    bool operator < (const DPState &state2) const {
      return score < state2.score
	|| (score==state2.score && from < state2.from);
    }
  };

  class Eagle {
  public:
    mutable uint64 totTicks, extTicks, diphapTicks, lshTicks, lshCheckTicks, dpTicks, dpStaticTicks, dpSwitchTicks, dpSortTicks, dpUpdateTicks, dpUpdateCalls, blipFixTicks, blipPopTicks, blipVoteTicks, blipLshTicks;
    
  private:

    static const uint homErrCost = 1, hetErrCost = 2, switchCost = 3;
    static const uint switchScoreLutBits = 5;
    char switchScoreLut[1<<(3*switchScoreLutBits)][2];
    const uint64 N, Nref; // Nref = 0 if not in ref-mode
    const uint64 Mseg64; // number of <=64-SNP chunks
    const uint64_masks *genoBits; // [[MATRIX]]: M64 x N (is0 and is2 64-bit masks)
    const std::vector <std::vector <double> > seg64cMvecs;
    uchar *maskSnps64j; // M64x64 binary vector indicating SNPs to use
    double *cMs64j; // M64x64+1 cM coordinates
    uchar *phaseConfs, *phaseConfs2; // [[MATRIX]]: 2N x M64x64
    uint64 *haploBits; // [[MATRIX]]: M64 x 2N (is1 for hard calls)
    uint64 *haploBitsT; // [[MATRIX]]: 2N x M64 (is1 for hard calls)
    uint64 *tmpHaploBitsT; // [[MATRIX]]: 2Ntarget x M64 (temp storage for target haps in ref-mode)
    uchar *segConfs; // [[MATRIX]]: 2N x M64
    std::vector < std::vector <StaticMultimap> > hashLookups;
    std::vector < std::vector < std::vector <uint64> > > hashBits;
    const AlleleFreqs *seg64logPs;
    const std::vector <double> invLD64j; // M64x64 LD-based weights for evaluating match evidence
    const std::vector <IndivInfoX> indivs;
    const std::vector <SnpInfoX> snps; // M-vector
    std::vector <uchar> maskIndivs; // N-vector: 0 to ignore indivs (e.g., relatives)
    std::vector <bool> isFlipped64j; // in non-ref mode, SNPs are internally flipped to 0=A2=major
    const double logPerr; // genotype error rate
    
    void init(void);
    uint getGeno0123(uint64 m64j, uint64 n) const;
    void retractMatch(uint n0, Match &match, double memoLogBF[][4]) const;
    Match computeDuoLogBF(double memoLogBF[][4], double workLogBF[], uint64 n0, uint64 n1, uint64 m64cur) const;
    void trim(Match &match, const Match &ref, uint64 n0, int orientation, uint64 trimStart,
	      int inc, double workLogBF[]) const;
    std::string computePhaseString(uint64 n0, uint64 nF1, uint64 nF2,
				   const std::vector <Match> &matches,
				   const std::vector <int> &signs, uint64 start, double cMend,
				   bool cons)
      const;
    void printMatch(uint64 n0, uint64 nF1, uint64 nF2, const Match &duoMatch,
		    double memoLogBF[][4]) const;
    void findLongHalfIBD(uint64 n0, std::vector <uint> topInds[2],
			 std::vector <uint> topIndsLens[2], uint K) const;
    std::vector <uint> findMinErrDipHap(uint64 n0, uint K, bool useTargetHaps) const;
    void findLongDipHap(uint64 n0, std::vector <uint> topInds[2],
			std::vector <uint> topIndsLens[2], uint K, uint errStart) const;
    void computePhaseConfs(uint64 n0, const std::vector <Match> &matches,
			   const std::vector <int> &signs, bool cons);
    std::vector <int> trioRelPhase(uint64 n0, uint64 nF1, uint64 nF2) const;
    void checkPhase(uint64 n0, uint64 nF1, uint64 nF2, double thresh) const;
    std::vector <bool> checkPhaseConfsPhase(uint64 n0, uint64 nF1, uint64 nF2) const;
    void checkHapPhase(uint64 n0, uint64 nF1, uint64 nF2, const uint64 curHaploBitsT[], uint64 m64,
		       uint64 side, std::vector < std::vector <int> > votes=std::vector < std::vector <int> > ()) const;
    std::vector <bool> checkHapPhase1(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap,
				      uint64 m64start, uint64 m64end, int sign=1) const;
    std::vector <bool> checkHapPhase1j(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap,
				       uint64 m64jStart, uint64 m64jEnd, int sign=1) const;
    std::vector <bool> checkHapPhase1jCall(uint64 n0, uint64 nF1, uint64 nF2, uint64 callBitsT[],
					   uint64 m64jStart, uint64 m64jEnd, bool print, int sign=1) const;
    int checkHapPhase2(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap,
				      uint64 n2hap, uint64 n3hap, uint64 m64, int sign) const;
    std::vector <bool> checkHaploBits(uint64 n0, uint64 nF1, uint64 nF2, uint64 hapBits,
				      uint64 m64, int pad=0) const;
    std::pair <uint64, uint64> phaseSegHMM(uint64 n0, uint64 n1hap, uint64 n2hap, uint64 n3hap,
					   uint64 m64, uint64 &hetErrMask) const;
    std::vector <bool> checkSegPhase(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap, uint64 n2hap,
				     int sign, uint64 m64) const;
    void computeSegPhaseConfs(uint64 n0, uint64 n1hap, uint64 n2hap, int sign, uint64 m64,
			      int err);
    int numDipHapWrongBits(uint64 m64, uint64 n0, uint64 n1hap) const;
    int firstDipHapGoodBit(uint64 m64, uint64 n0, uint64 n1hap) const;
    int firstDipHapWrongBit(uint64 m64, uint64 n0, uint64 n1hap) const;
    uint computeHash(const uint64 *curHaploBitsT, const uint64 *curHashBits, uint B) const;
    uint computeHash(const uint64 *curHaploBitsT, const std::vector <uint64> &curHashBits) const;

    bool updateHelper(std::unordered_map <uint64, DPState> &dpTab, uint &dpBestScore,
		      std::pair <uint, uint> cur, std::pair <uint, uint> next, uint score) const;
    uint computeStaticScore(uint n0, uint n1hap, uint n2hap, uint64 m64) const;
    uint computeSwitchScore(uint n0, uint n1hap, uint n2hapA, uint n2hapB, uint64 m64) const;
    void updateTable(std::unordered_map <uint64, DPState> dpTable[], uint dpBestScores[],
		     uint64 m64, uint64 dist, uint n0, uint n1hapA, uint n2hapA, uint n1hapB,
		     uint n2hapB, uint score) const;
    void safeInsert(std::set <uint> &refHapSet, uint n1hap, uint n0) const;
    std::vector < std::pair <uint64, uint64> > findGoodSegs(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap) const;
    void updateFarHaps(std::vector < std::pair <uint, uint> > &farHaps, uint n1hap, uint m64jStart, uint m64jEnd) const;
    double computeLogHetP(uint64 m64j) const;

  public:
    Eagle(uint64 _N, uint64 _Mseg64, const uint64_masks *_genoBits,
	  std::vector < std::vector <double> > _seg64cMvecs, const AlleleFreqs *_seg64freqs,
	  std::vector <double> _invLD64j, const std::vector <IndivInfoX> &_indivs,
	  const std::vector <SnpInfoX> &_snps, const std::string &maskFile,
	  const std::vector <bool> &isFlipped64j, double _pErr, int runStep2);
    // constructor for ref-mode
    Eagle(uint64 _Nref, uint64 _Ntarget, uint64 _Mseg64, const uint64_masks *_genoBits,
	  std::vector < std::vector <double> > _seg64cMvecs, double _pErr);

    void reallocLRPtoPBWT(void);
    ~Eagle();
    
    void checkTrioErrorRate(uint64 n0, uint64 nF1, uint64 nF2) const;
    void randomlyPhaseTmpHaploBitsT(uint64 n0);
    std::pair <double, std::vector <double> > findLongDipMatches(uint64 n0, uint64 nF1,
								 uint64 nF2);
    double findLongHapMatches(uint64 n0, uint64 nF1, uint64 nF2, int iter);
    double runHMM(uint64 n0, uint64 nF1, uint64 nF2, int iter, uint beamWidth, uint maxHapStates);
    std::vector <bool> computeRefIsMono(const std::vector <uint> &bestHaps) const;
    float runPBWT(uint64 n0, uint64 nF1, uint64 nF2, int Kpbwt, double cMexpect, double histFactor,
		  bool runReverse, bool useTargetHaps, bool impMissing, bool isChrX);
    float runPBWT(uint64 n0, uint64 nF1, uint64 nF2, int Kpbwt, double cMexpect, double histFactor,
		  bool runReverse, bool useTargetHaps, bool impMissing, int usePS,
		  const std::vector < std::pair <int, int> > &conPS, bool isChrX);
    void imputeMissing(const HapHedge &hapHedge, uint64 n0);
    void writePhaseConfs(const std::string &tmpPhaseFile) const;
    void readPhaseConfs(const std::string &tmpPhaseFile);
    void cpPhaseConfs(uint64 n0start, uint64 n0end);
    void cpTmpHaploBitsT(uint64 n0start, uint64 n0end);
    void outputSE(const std::vector <uint> &children, const std::vector <uint> &nF1s,
		  const std::vector <uint> &nF2s, int step) const;
    void writeHapsGzSample(const std::string &prefix) const;
    void writeVcf(const std::string &tmpFile, const std::vector <bool> &isTmpPhased,
		  const std::string &outFile, int chromX, double bpStart, double bpEnd,
		  const std::string &writeMode, bool noImpMissing, bool keepMissingPloidyX,
		  int argc, char**argv) const;
    void writeVcfNonRef(const std::string &vcfFile, const std::string &outFile, int inputChrom,
			int chromX, double bpStart, double bpEnd, const std::string &writeMode,
			bool noImpMissing, bool keepMissingPloidyX, int argc, char **argv) const;
    void makeHardCalls(uint64 n0start, uint64 n0end, uint seed);
    void initRefIter(int refIter);
    void buildHashTables(int iter, int batch, int seed);
    const uint64 *getHaploBitsT(void) const;
    uint64 getNlib(int iter) const;
    uint64 getMseg64(void) const;
    const uchar *getMaskSnps64j(void) const;
    double computeHetRate(void) const;

    static int countSE(const std::vector <bool> &phaseVec);
    static int countMajorSE(const std::vector <bool> &phaseVec);

  };
}

#endif
