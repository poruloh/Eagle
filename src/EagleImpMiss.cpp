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

#include <vector>
#include <iostream>
#include <cmath>

#include "HapHedge.hpp"
#include "MemoryUtils.hpp"
#include "NumericUtils.hpp"
#include "Eagle.hpp"


namespace EAGLE {

  using std::cout;
  using std::endl;
  using std::vector;
  using std::string;

#define HAP_BEAM_WIDTH 16
  struct ProbMaskBundle {
    float logTotProb;
    int numTop;
    float probs[HAP_BEAM_WIDTH];
    uint64 masks[HAP_BEAM_WIDTH];
  };

  struct MultMaskState {
    float mult; // 1/Nhaps * pErr^numErrs
    uint64 mask;
    HapTreeState state;
    char rmActive[2];
    inline int count(void) const {
      return state.count - rmActive[0] - rmActive[1];
    }
    inline float multCount(void) const {
      if (mult == 0) return 0;
      else return mult * count();
    }
    bool operator < (const MultMaskState &mms) const {
      return multCount() > mms.multCount();
    }
  };

  void advance(const HapTree &hapTree, const uchar rmHaps[], const uchar hap[],
	       const uchar missing[], MultMaskState states[HAP_BEAM_WIDTH], int m, float pErr) {
    MultMaskState nextStates[2*HAP_BEAM_WIDTH];
    int numNext = 0;
    // impose pruning thresh (TODO: test)
    const float minP = states[0].multCount() * pErr * pErr;
    for (int k = 0; k < HAP_BEAM_WIDTH && states[k].multCount() > minP; k++)
      for (int b = 0; b < 2; b++) {
	HapTreeState nextState = states[k].state;
	if (hapTree.next(m, nextState, b)) {
	  if (missing[m]) {
	    nextStates[numNext].mult = states[k].mult;
	    nextStates[numNext].mask = (states[k].mask<<1)|b;
	  }
	  else {
	    nextStates[numNext].mult = states[k].mult * (b==hap[m] ? 1 : pErr);
	    nextStates[numNext].mask = states[k].mask;
	  }
	  nextStates[numNext].state = nextState;
	  nextStates[numNext].rmActive[0] = states[k].rmActive[0] && ((rmHaps[m]&1)==b);
	  nextStates[numNext].rmActive[1] = states[k].rmActive[1] && ((rmHaps[m]>>1)==b);
	  numNext++;
	}
      }
    std::sort(nextStates, nextStates + numNext);
    memcpy(states, nextStates, std::min(numNext, HAP_BEAM_WIDTH) * sizeof(states[0])); // copy best
    if (numNext < HAP_BEAM_WIDTH) states[numNext].mult = 0;
  }

  string logProbToStr(float f) {
    f /= log(10);
    char buf[100];
    sprintf(buf, "%.2fe%d", pow(10, f - floor(f)), (int) floor(f));
    return buf;
  }

  // rmHaps: 2-bit phased genotypes for haplotype pair to ignore... or NULL if none (ref mode)
  // recombLogPs[T+1]: [0,T] = 0; [1..T-1] = logP for recomb in (t-0.5,t+0.5)
  void impMissing(const HapHedge &hapHedge, const uchar *rmHaps, uchar hap[],
		  const uchar missing[], const float recombLogPs[], float pErr) {
    const int maxExt = 500;
    const int M = hapHedge.getM();
    const int skip = hapHedge.getSkip();
    const int T = hapHedge.getNumTrees();
    const int dtMax = maxExt / skip;

    ProbMaskBundle *topPrefixes = (ProbMaskBundle *) calloc(T * dtMax, sizeof(topPrefixes[0]));
    // calloc => numTop initialized to 0

    vector <float> fwdLogProbs(T+1, -1000000), bwdLogProbs(T+1, -1000000);
    fwdLogProbs[0] = 0; bwdLogProbs[T] = 0;
  
    MultMaskState states[HAP_BEAM_WIDTH];
    // compute haplotype prefix beams; compute fwdLogProbs
    for (int t = 0; t < T; t++) {
      const HapTree &hapTree = hapHedge.getHapTree(t);
      states[0].mask = 0;
      states[0].state = hapTree.getRootState();
      if (rmHaps == NULL) {
	states[0].mult = hapTree.getInvNhaps();
	states[0].rmActive[0] = states[0].rmActive[1] = 0;
      }
      else {
	states[0].mult = 1 / (1/hapTree.getInvNhaps()-2); // removing 2 haps
	states[0].rmActive[0] = states[0].rmActive[1] = 1;
      }
      for (int k = 1; k < HAP_BEAM_WIDTH; k++) states[k].mult = 0;
      int m = t*skip;
      int dtMiss = 0;
      for (int dt = 0; dt<dtMax && t+dt<T; dt++) {
	int mEnd = std::min(m + skip, M);
	for (; m < mEnd; m++) {
	  if (missing[m]) dtMiss++;
	  advance(hapTree, rmHaps, hap, missing, states, m, pErr);
	}
	if (dtMiss > 64) break;
	// impose pruning threshold (TODO: test)
	float dtPrevBestProb = (dt==0 ? 0 : topPrefixes[t*dtMax + (dt-1)].probs[0]);
	float minP = dtPrevBestProb * pErr * pErr * expf(recombLogPs[t+dt]); // rel to previous
	float totProb = 0;
	ProbMaskBundle &bundle = topPrefixes[t*dtMax + dt];
	bundle.numTop = 0;
	//cout << "t = " << t << ", dt = " << dt << ":" << endl;
	for (int k = 0; k < HAP_BEAM_WIDTH && states[k].multCount() > minP; k++) {
	  float prob = states[k].multCount();
	  totProb += prob;
	  bundle.numTop++;
	  bundle.probs[k] = prob;
	  bundle.masks[k] = states[k].mask;
	  //cout << "  prob = " << prob << ", mask = " << states[k].mask << endl;
	}
	if (bundle.numTop == 0) break;
	bundle.logTotProb = logf(totProb);
	// compute fwdLogProbs
	NumericUtils::logSumExp(fwdLogProbs[t+dt+1],
				fwdLogProbs[t] + bundle.logTotProb + recombLogPs[t+dt+1]);
	//cout << endl;
      }
      //cout << endl;
    }

    // count missing sites
    int numMiss = 0;
    vector <int> tMiss;
    for (int t = 0; t < T; t++) {
      tMiss.push_back(numMiss);
      int m = t*skip;
      int mEnd = std::min(m + skip, M);
      for (; m < mEnd; m++)
	if (missing[m])
	  numMiss++;
    }
    // initialize log probs for missing sites
    float allLogProb01s[numMiss][2];
    for (int i = 0; i < numMiss; i++)
      for (int b = 0; b < 2; b++)
	allLogProb01s[i][b] = -1000000;

    // compute bwdLogProbs; compute log probs at missing sites (using saved haplotype prefix beams)
    for (int t = T-1; t >= 0; t--) {
      int m = t*skip;
      int dtMiss = 0;
      for (int dt = 0; dt<dtMax && t+dt<T; dt++) {
	// compute bwdLogProbs
	const ProbMaskBundle &bundle = topPrefixes[t*dtMax + dt];
	if (bundle.numTop == 0) break;
	NumericUtils::logSumExp(bwdLogProbs[t],
				recombLogPs[t] + bundle.logTotProb + bwdLogProbs[t+dt+1]);

	int mEnd = std::min(m + skip, M);
	for (; m < mEnd; m++)
	  if (missing[m])
	    dtMiss++;
	float prob01s[dtMiss][2]; memset(prob01s, 0, dtMiss*2*sizeof(prob01s[0][0]));
	for (int k = 0; k < bundle.numTop; k++)
	  for (int i = 0; i < dtMiss; i++)
	    prob01s[i][(bundle.masks[k]>>(dtMiss-1-i))&1] += bundle.probs[k];
	for (int i = 0; i < dtMiss; i++)
	  for (int b = 0; b < 2; b++)
	    if (prob01s[i][b] != 0)
	      NumericUtils::logSumExp(allLogProb01s[tMiss[t]+i][b],
				      fwdLogProbs[t] + logf(prob01s[i][b]) + bwdLogProbs[t+dt+1]);
      }
    }
  
    // impute missing sites
    numMiss = 0;
    for (int m = 0; m < M; m++)
      if (missing[m]) {
	hap[m] = (allLogProb01s[numMiss][1] > allLogProb01s[numMiss][0]);
	numMiss++;
      }
    /*
    for (int t = 0; t <= T; t++)
      cout << "fwdProbs[" << t << "]: " << logProbToStr(fwdLogProbs[t]) << endl;
    for (int t = 0; t <= T; t++)
      cout << "bwdProbs[" << t << "]: " << logProbToStr(bwdLogProbs[t]) << endl;
    for (int i = 0; i < numMiss; i++) {
      cout << "missing site " << i << ":";
      for (int b = 0; b < 2; b++)
	cout << " " << logProbToStr(allLogProb01s[i][b]);
      cout << endl;
    }
    */
    free(topPrefixes);
  }




  void Eagle::imputeMissing(const HapHedge &hapHedge, uint64 n0) {
    int M = hapHedge.getM(), skip = hapHedge.getSkip(), T = hapHedge.getNumTrees();
    const HapBitsT &hapBitsT = hapHedge.getHapBitsT();

    uchar *rmHaps = NULL;
    if (2*n0 < (uint64) hapBitsT.getNhaps()) { // set rmHaps
      rmHaps = ALIGNED_MALLOC_UCHARS(M * sizeof(rmHaps[0]));
      for (int m = 0; m < M; m++)
	rmHaps[m] = hapBitsT.getBit(2*n0, m) | (hapBitsT.getBit(2*n0+1, m)<<1);
    }

    vector <uchar> missing(M);
    int m = 0;
    vector <double> cMtreeStarts;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
      if (maskSnps64j[m64j]) {
	if (m % skip == 0)
	  cMtreeStarts.push_back(cMs64j[m64j]);
	missing[m++] = (genoBits[m64j/64 * N + n0].is9>>(m64j&63))&1;
      }
    cMtreeStarts.push_back(cMs64j[Mseg64*64]);

    vector <float> recombLogPs(T+1);
    const double cMswitch = 2.0;
    for (int t = 1; t < T; t++) {
      double cMdelta = (cMtreeStarts[t+1] - cMtreeStarts[t-1]) / 2;
      recombLogPs[t] = log(std::max(1 - exp(-cMdelta / cMswitch), 1e-6));
      //cout << exp(recombLogPs[t]) << " " << std::flush;
    }
    //cout << endl;

    for (uint64 nHap = 2*(n0-Nref); nHap <= 2*(n0-Nref)+1; nHap++) {
      vector <uchar> hap(M);
      m = 0;
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
	if (maskSnps64j[m64j])
	  hap[m++] = (tmpHaploBitsT[nHap*Mseg64 + m64j/64]>>(m64j&63))&1;

      impMissing(hapHedge, rmHaps, &hap[0], &missing[0], &recombLogPs[0], pow(10.0, logPerr));

      m = 0;
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
	if (maskSnps64j[m64j]) {
	  uint64 bit = 1ULL<<(m64j&63);
	  if (hap[m])
	    tmpHaploBitsT[nHap*Mseg64 + m64j/64] |= bit;
	  else
	    tmpHaploBitsT[nHap*Mseg64 + m64j/64] &= ~bit;
	  m++;
	}
    }

    if (rmHaps != NULL)
      ALIGNED_FREE(rmHaps);    
  }

}
