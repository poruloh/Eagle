/*
   This file is part of the Eagle haplotype phasing software package
   developed by Po-Ru Loh.  Copyright (C) 2015-2016 Harvard University.

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

#include "HapHedge.hpp"

namespace EAGLE {

  struct HapPathSplit {
    int t; // new start
    float relProbLastStop; // relative to cumLogP in HapPath
    int hapPrefixInd;
    int hapPrefixTo[2];
    HapPathSplit(void);
    HapPathSplit(int _t);
    HapPathSplit(int _t, float relProb, int ind);
  };

  struct HapPath {
    float cumLogP;
    int splitListLength;
    HapPathSplit *splitList; // [max size = histLength]
    int to[2];
    float toCumLogP[2];
    HapPath(void);
  };

  struct HapPrefix {
    HapTreeState state;
    //int to[2];
    float toHetOnlyProb[2];
    HapPrefix(void);
    HapPrefix(const HapTreeState &_state);
  };

  class HapWaves {

    const HapHedgeErr &hapHedge;
    const std::vector <float> &cMcoords;
    const int histLength, beamWidth;
    const float pErr;
    const int maxHapPaths, maxHapPrefixes;
    int tCur, curParity, nextParity;
    int hapPathSizes[2];
    HapPath *hapPaths[2]; // [max size = 2*beamWidth each]
    int hapPrefixSizes[2];
    HapPrefix *hapPrefixes[2]; // [max size = 2*beamWidth * histLength * 2 each]

  public:
    
    HapWaves(const HapHedgeErr &_hapHedge, const std::vector <float> &_cMcoords,
	     int _histLength, int _beamWidth, float _logPerr, int _tCur);
    ~HapWaves(void);

    // populate hapPrefixes[nextParity]
    // populate toCumLogP[] in hapPaths[curParity] (but don't populate hapPaths[nextParity])
    void computeAllExtensions(const std::vector <uchar> &nextPossibleBits);

    float getToCumLogProb(int ind, int nextBit) const;

    // look up/create extension of hapPaths[curParity][ind] in hapPaths[nextParity]
    // return index in hapPaths[nextParity]
    int extendPath(int ind, int nextBit);
    
    void advance(void);

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
	    const char *_constraints, const std::vector <float> &_cMcoords,
	    int _histLength, int _beamWidth, float _logPerr, int _tCur);

    // compute probability of AA at hets tCallLoc1 and tCallLoc2
    float callProbAA(int tCallLoc1, int tCallLoc2, int callLength);
    // compute diploid dosage at tCallLoc
    float callDosage(int tCallLoc, int callLength);
  };

}

#endif
