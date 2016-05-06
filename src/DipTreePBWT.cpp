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

#include <vector>
#include <iostream>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "HapHedge.hpp"
#include "NumericUtils.hpp"
#include "Timer.hpp"
#include "DipTreePBWT.hpp"

namespace EAGLE {

  using std::vector;
  using std::cout;
  using std::endl;

  const int TO_UNKNOWN = -2, TO_NONE = -1; // TO_NONE used in HapPathSplit and HapPrefix


  // struct HapPathSplit

  HapPathSplit::HapPathSplit(void) {};
  HapPathSplit::HapPathSplit(int _t) : t(_t), relProbLastStop(1), hapPrefixInd(0) {
    hapPrefixTo[0] = hapPrefixTo[1] = TO_UNKNOWN;
  };
  HapPathSplit::HapPathSplit(int _t, float relProb, int ind)
    : t(_t), relProbLastStop(relProb), hapPrefixInd(ind) {
    hapPrefixTo[0] = hapPrefixTo[1] = TO_UNKNOWN;
  };


  // struct HapPath

  HapPath::HapPath(void) {};


  // struct HapPrefix

  HapPrefix::HapPrefix(void) {};
  HapPrefix::HapPrefix(const HapTreeState &_state) {
    state = _state;
    //to[0] = to[1] = TO_UNKNOWN;
  };

  // class HapWaves

  HapWaves::HapWaves(const HapHedgeErr &_hapHedge, const vector <float> &_cMcoords,
		     int _histLength, int _beamWidth, float _logPerr, int _tCur) :
    hapHedge(_hapHedge), cMcoords(_cMcoords), histLength(_histLength),
    beamWidth(_beamWidth), pErr(expf(_logPerr)), maxHapPaths(2*beamWidth),
    maxHapPrefixes(maxHapPaths*histLength*2+1), tCur(_tCur) {

    curParity = tCur&1; nextParity = 1-curParity;

    for (int p = 0; p < 2; p++) {
      hapPathSizes[p] = 0;
      hapPaths[p] = new HapPath[maxHapPaths];
      for (int i = 0; i < maxHapPaths; i++)
	hapPaths[p][i].splitList = new HapPathSplit[histLength];
      hapPrefixes[p] = new HapPrefix[maxHapPrefixes];
    }

    // add root of 0th HapTree as cur HapPath
    hapPaths[curParity][0].cumLogP = 0;
    hapPaths[curParity][0].splitListLength = 1;
    hapPaths[curParity][0].splitList[0] = HapPathSplit(tCur);
    hapPaths[curParity][0].to[0] = hapPaths[0][0].to[1] = TO_UNKNOWN;
    hapPathSizes[curParity] = 1;
      
    hapPrefixes[curParity][0] = HapPrefix(hapHedge.getHapTreeMulti(tCur).getRootState());
    hapPrefixSizes[curParity] = 1;

  }

  HapWaves::~HapWaves(void) {
    for (int p = 0; p < 2; p++) {
      for (int i = 0; i < maxHapPaths; i++)
	delete[] hapPaths[p][i].splitList;
      delete[] hapPaths[p];
      delete[] hapPrefixes[p];
    }
  }

  // populate hapPrefixes[nextParity]
  // populate toCumLogP[] in hapPaths[curParity] (but don't populate hapPaths[nextParity])
  void HapWaves::computeAllExtensions(const vector <uchar> &nextPossibleBits) {
    // add root of next (= new cur) HapTree as beginning of HapPrefix list
    if (tCur+1 < (int) cMcoords.size()) {
      hapPrefixes[nextParity][0] = HapPrefix(hapHedge.getHapTreeMulti(tCur+1).getRootState());
      hapPrefixSizes[nextParity] = 1;
    }
    
    float mult = hapHedge.getHapTreeMulti(tCur).getInvNhaps();

    // iterate over paths
    for (int i = 0; i < hapPathSizes[curParity]; i++) {
      float relProbStopNext[2] = {0, 0};
      // iterate over splits
      for (int j = 0; j < hapPaths[curParity][i].splitListLength; j++) {
	HapPathSplit &split = hapPaths[curParity][i].splitList[j];
	// iterate over next possible bits
	for (int b = 0; b < 2; b++) {
	  if (!((nextPossibleBits[i]>>b)&1)) continue;
	  HapPrefix &hapPrefix = hapPrefixes[curParity][split.hapPrefixInd];
	  // if extension of hap prefix hasn't been attempted, attempt to perform extension
	  if (split.hapPrefixTo[b] == TO_UNKNOWN) {
	    split.hapPrefixTo[b] = TO_NONE; // default: can't extend (overwrite if path found)
	    hapPrefix.toHetOnlyProb[b] = 0;
	    // try to extend hap prefix:
	    // fill in split.hapPrefixTo[b], hapPrefixes[curParity][split.hapPrefixInd].to*[b]
	    const HapTreeMulti &hapTree = hapHedge.getHapTreeMulti(split.t);

	    HapTreeState state = hapPrefix.state;
	    if (hapTree.next(2*tCur, state, b)) { // can extend to match at het
	      hapPrefix.toHetOnlyProb[b] += mult * state.count;
	      if (hapTree.next(2*tCur+1, state, 0)) { // no err in inter-het region
		// create and link new HapPrefix node in hapPrefixes[nextParity]; link
		split.hapPrefixTo[b] = hapPrefixSizes[nextParity]++;
		hapPrefixes[nextParity][split.hapPrefixTo[b]].state = state;
	      }
	    }
	  }
	  float recombP;
	  if (tCur+1 == (int) cMcoords.size()) recombP = 1.0f;
	  else {
	    const float cMpseudo = 2.0f, minRecombP = 0.000001f, maxRecombP = 1.0f;//pErr;
	    recombP = std::max(std::min(3 * (cMcoords[tCur+1]-cMcoords[tCur])
					/ (cMpseudo + 2*(cMcoords[tCur+1]-cMcoords[split.t])),
					maxRecombP), minRecombP);
	  }
	  relProbStopNext[b] += split.relProbLastStop * hapPrefix.toHetOnlyProb[b] * recombP;
	}
      }
      for (int b = 0; b < 2; b++) {
	if (!((nextPossibleBits[i]>>b)&1)) continue;
	float relLogP = -1000;
	if (relProbStopNext[b] != 0) relLogP = logf(relProbStopNext[b]);
	hapPaths[curParity][i].toCumLogP[b] =
	  hapPaths[curParity][i].cumLogP + relLogP;// + recombLogPs[tCur];
      }
    }
  }

  float HapWaves::getToCumLogProb(int ind, int nextBit) const {
    return hapPaths[curParity][ind].toCumLogP[nextBit];
  }

  // look up/create extension of hapPaths[curParity][ind] in hapPaths[nextParity]
  // return index in hapPaths[nextParity]
  int HapWaves::extendPath(int ind, int nextBit) {
    HapPath &curHapPath = hapPaths[curParity][ind];
    if (curHapPath.to[nextBit] == TO_UNKNOWN) {
      int nextInd = hapPathSizes[nextParity]++;
      assert(hapPathSizes[nextParity]<=maxHapPaths);
      curHapPath.to[nextBit] = nextInd;
      HapPath &nextHapPath = hapPaths[nextParity][nextInd];
      nextHapPath.cumLogP = curHapPath.toCumLogP[nextBit];
      float calibP = expf(curHapPath.cumLogP - nextHapPath.cumLogP);
      int &nSplit = nextHapPath.splitListLength; nSplit = 0;
      nextHapPath.to[0] = nextHapPath.to[1] = TO_UNKNOWN;
      for (int j = (curHapPath.splitList[0].t + histLength == tCur+1 ? 1 : 0);
	   j < curHapPath.splitListLength; j++) {
	const HapPathSplit &curSplit = curHapPath.splitList[j];
	if (curSplit.hapPrefixTo[nextBit] != TO_NONE) {
	  nextHapPath.splitList[nSplit++] = HapPathSplit(curSplit.t,
							 curSplit.relProbLastStop * calibP,
							 curSplit.hapPrefixTo[nextBit]);
	}
      }
      nextHapPath.splitList[nSplit++] = HapPathSplit(tCur+1); // restart
    }
    return curHapPath.to[nextBit];
  }
    
  void HapWaves::advance(void) {
    tCur++; curParity = tCur&1; nextParity = 1-curParity;
    hapPathSizes[nextParity] = 0;
    hapPrefixSizes[nextParity] = 0;
  }


  // struct DipTreeNode

  bool DipTreeNode::operator < (const DipTreeNode &dNode) const {
    return logP+boostLogP > dNode.logP+dNode.boostLogP;
  }


  // class DipTree

  void DipTree::traceNode(int t, int i) {
    int from = nodes[t][i].from;
    if (t>1) traceNode(t-1, from);
    cout << "(" << (int) nodes[t][i].hapMat << "," << (int) nodes[t][i].hapPat << ") ";
  }

  std::pair <uint64, uint64> truncPair(uint64 histMat, uint64 histPat, uint64 histBits) {
    uint64 mask = histBits>=64ULL ? -1ULL : (1ULL<<histBits)-1;
    uint64 x = histMat&mask, y = histPat&mask;
    return x<y ? std::make_pair(x, y) : std::make_pair(y, x);
  }

  void DipTree::advance(void) {

    bool isOppConstrained = constraints[tCur]==OPP_CONSTRAINT; // constrained to be 0|1 or 1|0
    bool isFullyConstrained = !isOppConstrained && constraints[tCur]!=NO_CONSTRAINT;
    
    // populate next possible bits: nextPossibleBits[i] corresponds to hapPaths[curParity][i]
    //                              for i = dNode.hapPathInds[0], dNode.hapPathInds[1]
    vector <uchar> nextPossibleBits(2*beamWidth);
    int checkWidth = std::min((int) nodes[tCur].size(), beamWidth);
    vector <char> reqMats(checkWidth), reqPats(checkWidth);
    const float logPthresh = 2*logPerr;//logf(0.000001f);
    for (int i = 0; i < checkWidth; i++) {
      const DipTreeNode &dNode = nodes[tCur][i];
      if (dNode.logP+dNode.boostLogP < nodes[tCur][0].logP+nodes[tCur][0].boostLogP + logPthresh) {
	checkWidth = i;
	break;
      }
      assert(dNode.hapPathInds[0] < (int) nextPossibleBits.size());
      assert(dNode.hapPathInds[1] < (int) nextPossibleBits.size());
      if (isFullyConstrained) {
	char &reqMat = reqMats[i], &reqPat = reqPats[i];
	if ((constraints[tCur]>>1) == 0) // no-hom-err constraint
	  reqMat = reqPat = constraints[tCur]&1;
	else { // rel phase constraint
	  int t = tCur, ind = i;
	  for (int d = 0; d < (constraints[tCur]>>1)-1; d++)
	    ind = nodes[t--][ind].from;
	  reqMat = nodes[t][ind].hapMat ^ (constraints[tCur]&1);
	  reqPat = nodes[t][ind].hapPat ^ (constraints[tCur]&1);
	}
	nextPossibleBits[dNode.hapPathInds[0]] |= 1<<reqMat;
	nextPossibleBits[dNode.hapPathInds[1]] |= 1<<reqPat;
      }
      else {
	nextPossibleBits[dNode.hapPathInds[0]] = 3;
	nextPossibleBits[dNode.hapPathInds[1]] = 3;
      }
    }
    // extend hap paths (part 1)
    hapWaves.computeAllExtensions(nextPossibleBits);
    
    // extend dip paths
    vector <DipTreeNode> nextNodes;
    for (int i = 0; i < checkWidth; i++) {
      const DipTreeNode &dNode = nodes[tCur][i];
      for (char hapMat = 0; hapMat < 2; hapMat++)
	for (char hapPat = 0; hapPat < 2; hapPat++) {
	  if (!dNode.unequalAnc && hapMat > hapPat) continue;
	  if (isFullyConstrained && (hapMat != reqMats[i] || hapPat != reqPats[i])) continue;
	  if (isOppConstrained && hapMat==hapPat) continue;
	  DipTreeNode nextNode;
	  nextNode.from = i;
	  nextNode.unequalAnc = dNode.unequalAnc || (hapMat != hapPat);
	  nextNode.hapMat = hapMat;
	  nextNode.hapPat = hapPat;
	  nextNode.numErr = dNode.numErr + (genos[tCur]<=2 && hapMat+hapPat != genos[tCur]);
	  nextNode.logP = hapWaves.getToCumLogProb(dNode.hapPathInds[0], hapMat) +
	    hapWaves.getToCumLogProb(dNode.hapPathInds[1], hapPat) + nextNode.numErr * logPerr;
	  nextNode.boostLogP = dNode.boostLogP;
	  if (isFullyConstrained) {
	    nextNode.histMat = dNode.histMat;
	    nextNode.histPat = dNode.histPat;
	  }
	  else {
	    nextNode.histMat = (dNode.histMat<<1ULL) | hapMat;
	    nextNode.histPat = (dNode.histPat<<1ULL) | hapPat;
	  }
	  nextNodes.push_back(nextNode);
	}
    }

    if (!isFullyConstrained) {
      // compute number of bits of history to use (histLength minus # of fully constrained sites)
      int histBits = 0;
      for (int t = tCur; t > std::max(tCur-histLength, 0); t--)
	if (constraints[t]==OPP_CONSTRAINT || constraints[t]==NO_CONSTRAINT)
	  histBits++;

      // aggregate DipTree paths that agree exactly in past histLength
      std::sort(nextNodes.begin(), nextNodes.end());
      std::map < std::pair <uint64, uint64>, int > histToInd;
      for (int i = 0; i < (int) nextNodes.size(); i++) {
	const DipTreeNode &nextNode = nextNodes[i];
	std::pair <uint64, uint64> histPair =
	  truncPair(nextNode.histMat, nextNode.histPat, histBits);
	std::map < std::pair <uint64, uint64>, int >::iterator it = histToInd.find(histPair);
	if (it == histToInd.end()) {
	  histToInd[histPair] = nodes[tCur+1].size();
	  nodes[tCur+1].push_back(nextNode);
	}
	else {
	  int j = it->second;
	  float sumLogPj = nodes[tCur+1][j].logP + nodes[tCur+1][j].boostLogP;
	  float sumLogPi = nextNode.logP + nextNode.boostLogP;
	  NumericUtils::logSumExp(sumLogPi, sumLogPj); // prob i += prob existing tCur+1 node j
	  nodes[tCur+1][j].boostLogP += sumLogPi - sumLogPj; // augment boost for existing node j
	}
      }
      //cout << " " << nodes[tCur+1].size() << "/" << nextNodes.size() << std::flush;
    }
    else
      nodes[tCur+1] = nextNodes;
    
    // extend hap paths of top beamWidth DipTree nodes (part 2)
    for (int i = 0; i < std::min((int) nodes[tCur+1].size(), beamWidth); i++) {
      DipTreeNode &nextNode = nodes[tCur+1][i];
      const DipTreeNode &dNode = nodes[tCur][nextNode.from];
      nextNode.hapPathInds[0] = hapWaves.extendPath(dNode.hapPathInds[0], nextNode.hapMat);
      nextNode.hapPathInds[1] = hapWaves.extendPath(dNode.hapPathInds[1], nextNode.hapPat);
    }

    hapWaves.advance();
    tCur++;

    float totLogP = nodes[tCur][0].logP + nodes[tCur][0].boostLogP;
    for (int i = 1; i < (int) nodes[tCur].size(); i++)
      NumericUtils::logSumExp(totLogP, nodes[tCur][i].logP + nodes[tCur][i].boostLogP);
    for (int i = 0; i < (int) nodes[tCur].size(); i++) {
      normProbs[tCur].push_back(expf(nodes[tCur][i].logP + nodes[tCur][i].boostLogP - totLogP));
      //traceNode(tCur, i); cout << normProbs[tCur].back() << endl;
    }
  }

  DipTree::DipTree(const HapHedgeErr &_hapHedge, const vector <uchar> &_genos,
		   const char *_constraints, const vector <float> &_cMcoords,
		   int _histLength, int _beamWidth, float _logPerr, int _tCur) :
    hapWaves(_hapHedge, _cMcoords, _histLength, _beamWidth, _logPerr, _tCur),
    genos(_genos), constraints(_constraints), histLength(_histLength), beamWidth(_beamWidth),
    logPerr(_logPerr), tCur(_tCur), T(_hapHedge.getNumTrees()), nodes(T+1), normProbs(T+1) {

    DipTreeNode dNode;
    dNode.from = -1; dNode.unequalAnc = 0; dNode.logP = 0; dNode.numErr = 0;
    dNode.hapPathInds[0] = dNode.hapPathInds[1] = 0;
    dNode.histMat = 0; dNode.histPat = 0; dNode.boostLogP = 0;
    nodes[tCur].push_back(dNode); // root of DipTree
  }

  // compute probability of AA at hets tCallLoc1 and tCallLoc2
  float DipTree::callProbAA(int tCallLoc1, int tCallLoc2, int callLength) {
    assert(tCallLoc1>0 && tCallLoc2<T);
    int tFront = std::min(T, tCallLoc2 + callLength);
    while (tCur < tFront)
      advance();
    float probAA = 0, probAB = 0;
    for (int i = 0; i < (int) nodes[tFront].size(); i++) {
      int t = tFront, ind = i;
      while (t != tCallLoc2+1)
	ind = nodes[t--][ind].from;
      char hapMat2 = nodes[t][ind].hapMat, hapPat2 = nodes[t][ind].hapPat; // alleles at tCallLoc2
      while (t != tCallLoc1+1)
	ind = nodes[t--][ind].from;
      char hapMat1 = nodes[t][ind].hapMat, hapPat1 = nodes[t][ind].hapPat; // alleles at tCallLoc1
      if (hapMat2 != hapPat2 && hapMat1 != hapPat1) {
	if (hapMat1 == hapMat2)
	  probAA += normProbs[tFront][i];
	else
	  probAB += normProbs[tFront][i];
      }
      else {
	probAA += normProbs[tFront][i] / 2;
	probAB += normProbs[tFront][i] / 2;
      }
    }
    if (probAA + probAB == 0) return 0.5;
    return probAA / (probAA + probAB);
  }

  // compute diploid dosage at tCallLoc
  float DipTree::callDosage(int tCallLoc, int callLength) {
    assert(tCallLoc>0 && tCallLoc<T);
    int tFront = std::min(T, tCallLoc + callLength);
    while (tCur < tFront)
      advance();
    float prob1 = 0, probTot = 0;
    for (int i = 0; i < (int) nodes[tFront].size(); i++) {
      int t = tFront, ind = i;
      while (t != tCallLoc+1)
	ind = nodes[t--][ind].from;
      char hapMat = nodes[t][ind].hapMat, hapPat = nodes[t][ind].hapPat; // alleles at tCallLoc
      prob1 += (hapMat+hapPat) * normProbs[tFront][i];
      probTot += normProbs[tFront][i];
    }
    if (probTot == 0) return 1.0;
    return prob1 / probTot;
  }

}
