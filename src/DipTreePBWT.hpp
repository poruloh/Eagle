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

#ifndef DIPTREEPBWT_HPP
#define DIPTREEPBWT_HPP

#include <vector>

#include <boost/random.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_01.hpp>

#include "HapHedge.hpp"

namespace EAGLE {

  struct HapPathSplit {
    int t; // location of most recent start = tree index (note HapPath goes to tCur >= t)
    float relProbLastStop; // cumP for path ending just before t, relative to
                           // cumP = exp(cumLogP) for full HapPath ending just before tCur
    int hapPrefixInd; // index (in HapWaves::hapPrefixes[curMod][.])
                      // of the HapPrefix starting from t that ends the HapPath
    int hapPrefixTo[2]; // indices (in HapWaves::hapPrefixes[nextMod][.])
                        // of the extended HapPrefixes starting from t that have 0 (resp. 1)
                        // at split site tCur (i.e, bit 2*tCur) and no err (i.e., 0) at 2*tCur+1
    HapPathSplit(void);
    HapPathSplit(int _t); // split corresponding to new start at tCur=_t (i.e., root state)
    HapPathSplit(int _t, float relProb, int ind);
  };

  struct HapPath {
    float cumLogP; // log(cumP) for this path ending just before tCur
    int splitListLength;
    HapPathSplit *splitList; // [max size = histLength] list of prev split points, probs, prefixes
    int to[2]; // index of extension to bit = 0 (resp. 1) at split site tCur (i.e., HapPath tCur+1)
    float toCumLogP[2]; // log(cumP) for extended path
    HapPath(void);
  };

  struct HapPrefix {
    HapTreeState state; // state in tree getHapTreeMulti(split.t)
    //int to[2];
    float toHetOnlyProb[2]; // probability of extension to bit = 0 (resp. 1) at split site tCur
    HapPrefix(void);
    HapPrefix(const HapTreeState &_state);
  };

  struct RefHap {
    uint refSeq;
    short tLength;
    bool isEnd;
    short tMaskFwd, tMaskRev;
  };
  struct HapPair {
    RefHap haps[2];
  };

// history length to save for sampling ref haps (for in-sample imputation):
// needs to be a few splits longer than callLength passed to sampleRefs()
#define HAPWAVES_HIST 25

  class HapWaves {

    boost::lagged_fibonacci607 rng;
    boost::variate_generator<boost::lagged_fibonacci607&, boost::uniform_01<> > rand01;

    const HapHedgeErr &hapHedge;
    const std::vector <double> &cMcoords;
    const double cMexpect;
    const int histLength, beamWidth;
    const float pErr;
    const int maxHapPaths, maxHapPrefixes;
    int tCur, curMod, nextMod;
    int hapPathSizes[HAPWAVES_HIST];
    HapPath *hapPaths[HAPWAVES_HIST]; // [max size = 2*beamWidth each]
    int hapPrefixSizes[HAPWAVES_HIST];
    HapPrefix *hapPrefixes[HAPWAVES_HIST]; // [max size = 2*beamWidth * histLength * 2 each]

  public:
    
    HapWaves(const HapHedgeErr &_hapHedge, const std::vector <double> &_cMcoords, double cMexpect,
	     int _histLength, int _beamWidth, float _logPerr, int _tCur);
    ~HapWaves(void);

    float recombP(int tCur, int tSplit) const;

    // populate hapPrefixes[nextMod]
    // populate toCumLogP[] in hapPaths[curMod] (but don't populate hapPaths[nextMod])
    void computeAllExtensions(const std::vector <uchar> &nextPossibleBits);

    float getToCumLogProb(int ind, int nextBit) const;

    // look up/create extension of hapPaths[curMod][ind] in hapPaths[nextMod]
    // return index in hapPaths[nextMod]
    int extendPath(int ind, int nextBit);
    
    void advance(void);

    void sampleLastPrefix(int &tStart, HapTreeState &state, int t, int hapPathInd, int tBit);
  };

  struct DipTreeNode {
    int from;
    char unequalAnc, hapMat, hapPat;
    uint64 histMat, histPat;
    float boostLogP;
    float logP;
    int numErr;
    int hapPathInds[2];
    bool operator < (const DipTreeNode &dNode) const;
  };

  // constraint encoding:
  const char OPP_CONSTRAINT = -2; // require no het err, i.e., 0|1 or 1|0
  const char NO_CONSTRAINT = -1;
  // relative phase contraints are encoded as (dist=num_splits_to_ref_het<<1)|(rel_phase)
  // no hom err constraints (i.e., 0|0 at 0, 1|1 at 2) are encoded as 0 or 1 (i.e., dist=0 above)

  class DipTree {

    boost::lagged_fibonacci607 rng;
    boost::variate_generator<boost::lagged_fibonacci607&, boost::uniform_01<> > rand01;

    const HapHedgeErr &hapHedge;
    HapWaves hapWaves;
    const std::vector <uchar> &genos;
    const char *constraints;
    const int histLength, beamWidth;
    const float logPerr;
    int tCur; const int T;
    std::vector < std::vector <DipTreeNode> > nodes;
    std::vector < std::vector <float> > normProbs;

    void traceNode(int t, int i);
    void advance(void);

  public:

    DipTree(const HapHedgeErr &_hapHedge, const std::vector <uchar> &_genos,
	    const char *_constraints, const std::vector <double> &_cMcoords, double cMexpect,
	    int _histLength, int _beamWidth, float _logPerr, int _tCur);

    // compute probability of AA at hets tCallLoc1 and tCallLoc2
    float callProbAA(int tCallLoc1, int tCallLoc2, int callLength);
    // compute diploid dosage at tCallLoc
    float callDosage(int tCallLoc, int callLength);

    std::vector <HapPair> sampleRefs(int tCallLoc, int callLength, int samples,
				     const std::vector <uint> &bestHaps, bool isFwd);
  };

}

#endif
