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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "DipTreePBWT.hpp"
#include "HapHedge.hpp"
#include "NumericUtils.hpp"
#include "Timer.hpp"
#include "Types.hpp"
#include "Eagle.hpp"

namespace EAGLE {

  using std::vector;
  using std::cout;
  using std::endl;

  struct ProbInd {
    float prob;
    int ind1, ind2;
    ProbInd(float _prob=0, int _ind1=0, int _ind2=0) : prob(_prob), ind1(_ind1), ind2(_ind2) {}
    bool operator < (const ProbInd &pi) const {
      return std::min(prob, 1-prob) < std::min(pi.prob, 1-pi.prob);
    }
  };

  inline int popcount64_012(uint64 i) {
    if (i == 0) return 0;
    else if ((i & (i-1ULL)) == 0) return 1;
    else return 2;
  }

  vector <uint> Eagle::findMinErrDipHap(uint64 n0, uint K, bool useTargetHaps) const {

    uint64 Nhaps = 2*((Nref==0 || useTargetHaps) ? N : Nref);
    if (K > Nhaps) K = Nhaps;
    vector <uint> bestHaps; bestHaps.reserve(K);
    if (K == Nhaps) {
      for (uint64 nHap = 0; nHap < Nhaps; nHap++)
	if (nHap/2 != n0)
	  bestHaps.push_back(nHap);
    }
    else {
      vector < pair <uint, uint> > hapErrInds(Nhaps);
      vector <uint64_masks> genoBitsT(Mseg64);
      for (uint64 m64 = 0; m64 < Mseg64; m64++) genoBitsT[m64] = genoBits[m64*N + n0];
      for (uint64 nHap = 0; nHap < Nhaps; nHap++) {
	uint numErrs = (nHap/2 == n0 ? 1000000 : 0);
	for (uint64 m64 = 0; m64 < Mseg64; m64++) {
	  uint64 is1 = haploBitsT[nHap*Mseg64 + m64];
	  uint64 wrongBits = (genoBitsT[m64].is0 & is1) | (genoBitsT[m64].is2 & ~is1);
	  numErrs += popcount64_012(wrongBits);
	}
	hapErrInds[nHap].first = numErrs;
	hapErrInds[nHap].second = nHap;
      }
      std::sort(hapErrInds.begin(), hapErrInds.end());
      for (uint k = 0; k < K; k++)
	if (hapErrInds[k].second/2 != n0)
	  bestHaps.push_back(hapErrInds[k].second);
    }
    return bestHaps;
  }

  float Eagle::runPBWT(uint64 n0, uint64 nF1, uint64 nF2, int Kpbwt, bool runReverse,
		       bool useTargetHaps) {

    const int SPEED_FACTOR = 1; const float lnPerr = logf(powf(10.0f, logPerr));

    bool print = (int) nF1 != -1;
    
    /***** PROCESS TARGET GENOTYPES *****/
    
    // create vectors of genos, genoBits, and hets
    vector <uchar> genos64j(Mseg64*64);
    vector <uint64> hets64j;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      genos64j[m64j] = getGeno0123(m64j, n0);
      if (maskSnps64j[m64j] && genos64j[m64j]==1)
	hets64j.push_back(m64j);
    }
    vector <uint64_masks> tgtGenoBits(Mseg64);
    for (uint64 m64 = 0; m64 < Mseg64; m64++)
      tgtGenoBits[m64] = genoBits[m64*N + n0];
    
    vector <uint64> pbwtBitsFine(Mseg64);
    float conf = 0;

    // check for 0 or 1 het (no need to phase)
    if (hets64j.size() <= 1) {
      cerr << "WARNING: Sample " << n0-Nref << " has a het count of " << hets64j.size() << endl;
    }
    else {

    // find split sites (for PBWT HapTree starts): hets and occasionally inter-het sites
    vector < pair <int, int> > tCallLocs; vector <int> tHomLocs;
    vector <uint64> splits64j;
    splits64j.push_back(hets64j[0]);
    const double cMmaxSplit = 0.5;
    for (uint64 h = 1; h < hets64j.size(); h++) {
      int lastCallLoc = splits64j.size(); // old het ind + 1

      for (uint64 m64j = hets64j[h-1]+1; m64j <= hets64j[h]; m64j++)
	if (maskSnps64j[m64j] && genos64j[m64j] <= 2)
	  if (m64j == hets64j[h] || cMs64j[m64j] > cMs64j[splits64j.back()] + cMmaxSplit) {
	    splits64j.push_back(m64j);
	    if (m64j < hets64j[h])
	      tHomLocs.push_back(splits64j.size()); // hom ind + 1
	  }
      int nextCallLoc = splits64j.size(); // new het ind + 1
      tCallLocs.push_back(make_pair(lastCallLoc, nextCallLoc));
    }
    if (print)
      cout << "num hets: " << hets64j.size() << " num splits: " << splits64j.size() << endl;
    
    // create vector of genos at split sites (padded on left and right to match hapBitsT)
    vector <uchar> splitGenos;
    splitGenos.push_back(0); // pad on left with 0 (to match hapBitsT)
    for (uint64 s = 0; s < splits64j.size(); s++)
      splitGenos.push_back(genos64j[splits64j[s]]);
    splitGenos.push_back(0); // pad on right with 0 (to match hapBitsT)

    // compute recombination probabilities
    vector <float> cMcoords(splits64j.size()+2);
    for (uint64 s = 0; s <= splits64j.size(); s++) {
      uint64 splitStart = (s == 0 ? 0 : splits64j[s-1]);
      uint64 splitStop = (s == splits64j.size() ? Mseg64*64 : splits64j[s]);
      cMcoords[s] = cMs64j[splitStart]; cMcoords[s+1] = cMs64j[splitStop];
      int homs = 0;
      for (uint64 m64j = splitStart+1; m64j < splitStop; m64j++)
	if (genos64j[m64j] == 0 || genos64j[m64j] == 2)
	  homs++;
    }

    Timer timer;
    
    /**** BUILD PBWT DATA STRUCTURE *****/

    // select top ref haps
    if (print) cout << "selecting " << Kpbwt << " ref haps...    " << std::flush;
    vector <uint> bestHaps = findMinErrDipHap(n0, Kpbwt, useTargetHaps);
    if (print) cout << " done " << timer.update_time() << endl;

    // create HapBitsT encoding of ref hets and hom errs
    if (print) cout << "making HapBitsT...             " << std::flush;
    HapBitsT hapBitsT(haploBitsT, Mseg64, splits64j, splitGenos, tgtGenoBits, bestHaps);
    if (print) cout << " done " << timer.update_time() << endl;

    // create HapHedge PBWT data structure
    if (print) cout << "making HapHedge...             " << std::flush;
    HapHedgeErr *hapHedgePtr = new HapHedgeErr(hapBitsT);
    if (print) cout << " done " << timer.update_time() << endl;


    /***** RUN COARSE (UNCONSTRAINED) DIPTREE SEARCH *****/

    // initialize DipTree object
    if (print) cout << "making DipTree (unconstr)...   " << std::flush;
    vector <char> constraints(splitGenos.size(), NO_CONSTRAINT);
    const int histLengthFast = 30, pbwtBeamWidthFast = 30/SPEED_FACTOR;
    DipTree dipTreeFast(*hapHedgePtr, splitGenos, &constraints[0], cMcoords, histLengthFast,
			pbwtBeamWidthFast, lnPerr, 0);
    if (print) cout << " done " << timer.update_time() << endl;

    // explore search space; make phase calls
    if (print) cout << "making phase calls (uncon)...  " << std::flush;
    const int callLengthFast = 10;
    const float minFix = 0.5f, maxFix = 0.9f, fixThresh = 0.01f;
    vector <ProbInd> probInds;
    vector <uint64> pbwtBitsFast(Mseg64);
    uint64 lastBit = 0;
    for (uint64 i = 0; i < tCallLocs.size(); i++) {
      float probAA = dipTreeFast.callProbAA(tCallLocs[i].first, tCallLocs[i].second,
					    callLengthFast);
      ProbInd probInd(probAA, tCallLocs[i].first, tCallLocs[i].second);
      if (probAA < 0.5f)
	lastBit = !lastBit;
      uint64 m64j = hets64j[i+1];
      pbwtBitsFast[m64j/64] |= lastBit<<(m64j&63);
      if (i > 0) { // try calling rel phase vs. 2 hets back (in case prev het is err)
	float probAA2 = dipTreeFast.callProbAA(tCallLocs[i-1].first, tCallLocs[i].second,
					       callLengthFast);
	ProbInd probInd2(probAA2, tCallLocs[i-1].first, tCallLocs[i].second);
	if (probInd2 < probInd)
	  probInd = probInd2;
      }
      probInds.push_back(probInd);
    }
    int T = splitGenos.size(); // splits64j.size()+2
    vector <char> revConstraints(T, NO_CONSTRAINT);
    // set relative phase constraints for most confident hets
    std::sort(probInds.begin(), probInds.end());
    //float fracFixed = 0;
    for (int f = 0; f < minFix*probInds.size() ||
	   (f < maxFix*probInds.size()
	    && (probInds[f].prob < fixThresh || probInds[f].prob > 1-fixThresh)) ||
	   (f < (int) probInds.size() && (probInds[f].prob==0 || probInds[f].prob==1)); f++) {
      //fracFixed = (f+1.0f) / probInds.size();
      constraints[probInds[f].ind2] = // het constraint: fix relative phase of ind2 wrt ind1
	revConstraints[T-1-probInds[f].ind1] =
	((probInds[f].ind2 - probInds[f].ind1)<<1) | (probInds[f].prob < 0.5f);
      if (constraints[probInds[f].ind1] == NO_CONSTRAINT)
	constraints[probInds[f].ind1] = OPP_CONSTRAINT; // set to -2 = "start of het block"
      revConstraints[T-1-probInds[f].ind2] = OPP_CONSTRAINT;
    }
    // make hom dosage calls (and divide by 2 to sort properly: uncertainty = dist from 0 or 1)
    probInds.clear();
    for (uint64 i = 0; i < tHomLocs.size(); i++)
      probInds.push_back(ProbInd(dipTreeFast.callDosage(tHomLocs[i], callLengthFast) / 2,
				 tHomLocs[i]));
    std::sort(probInds.begin(), probInds.end());
    for (int f = 0; f < minFix*probInds.size() ||
	   (f < maxFix*probInds.size()
	    && (probInds[f].prob < fixThresh || probInds[f].prob > 1-fixThresh)) ||
	   (f < (int) probInds.size() && (probInds[f].prob==0 || probInds[f].prob==1)); f++)
      constraints[probInds[f].ind1] = revConstraints[T-1-probInds[f].ind1] =
	(probInds[f].prob >= 0.5f); // hom constraint: no err allowed
    if (print) cout << " done " << timer.update_time() << endl;
    //cout << "frac fixed: " << fracFixed << endl;


    /***** RUN FINE (CONSTRAINED) DIPTREE SEARCH *****/

    // initialize DipTree object
    if (print) cout << "making DipTree (constrained)..." << std::flush;
    const int histLengthFine = 100, pbwtBeamWidthFine = 50/SPEED_FACTOR;
    DipTree dipTreeFine(*hapHedgePtr, splitGenos, &constraints[0], cMcoords, histLengthFine,
			pbwtBeamWidthFine, lnPerr, 0);
    if (print) cout << " done " << timer.update_time() << endl;

    // explore search space; make phase calls
    if (print) cout << "making phase calls (constr)... " << std::flush;
    const int callLengthFine = 20;
    vector <float> probAAsCur;
    lastBit = 0;
    for (uint64 i = 0; i < tCallLocs.size(); i++) {
      float probAA = dipTreeFine.callProbAA(tCallLocs[i].first, tCallLocs[i].second,
					    callLengthFine);
      probAAsCur.push_back(probAA);
      conf += std::max(probAA, 1-probAA);
      if (probAA < 0.5f)
	lastBit = !lastBit;
      uint64 m64j = hets64j[i+1];
      pbwtBitsFine[m64j/64] |= lastBit<<(m64j&63);
    }
    conf /= tCallLocs.size();
    if (print) cout << " done " << timer.update_time() << endl;

    delete hapHedgePtr;


    if (runReverse) {
      // create HapBitsT encoding of ref hets and hom errs
      if (print) cout << "making revHapBitsT...          " << std::flush;
      HapBitsT revHapBitsT(hapBitsT, -2);
      if (print) cout << " done " << timer.update_time() << endl;

      // create HapHedge PBWT data structure
      if (print) cout << "making revHapHedge...          " << std::flush;
      HapHedgeErr revHapHedge(revHapBitsT);
      if (print) cout << " done " << timer.update_time() << endl;

      vector <uchar> revSplitGenos(splitGenos);
      std::reverse(revSplitGenos.begin(), revSplitGenos.end());
      vector <float> revcMcoords(T);
      for (int t = 0; t < T; t++) revcMcoords[t] = cMcoords[T-1] - cMcoords[T-1-t];

      // initialize DipTree object
      if (print) cout << "making revDipTree (constr)...  " << std::flush;
      DipTree revDipTreeFine(revHapHedge, revSplitGenos, &revConstraints[0], revcMcoords,
			     histLengthFine, pbwtBeamWidthFine, lnPerr, 0);
      if (print) cout << " done " << timer.update_time() << endl;

      // explore search space; make phase calls
      if (print) cout << "making rev phase calls (con)..." << std::flush;
      vector <uint64> revPbwtBitsFine(Mseg64);
      lastBit = 0;
      for (int i = tCallLocs.size()-1; i >= 0; i--) {
	float probAA = revDipTreeFine.callProbAA(T-1-tCallLocs[i].second, T-1-tCallLocs[i].first,
						 callLengthFine);
	if (probAA + probAAsCur[i] < 1)
	  lastBit = !lastBit;
	uint64 m64j = hets64j[i];
	revPbwtBitsFine[m64j/64] |= lastBit<<(m64j&63);
      }
      if (print) cout << " done " << timer.update_time() << endl;
      pbwtBitsFine = revPbwtBitsFine;
    }

    // trio analysis output
    if (print) {
      /*
      for (uint64 m64 = 0; m64 < Mseg64; m64++) {
	checkHaploBits(n0, nF1, nF2, pbwtBitsFast[m64], m64, 25);
	checkHaploBits(n0, nF1, nF2, pbwtBitsFine[m64], m64, 25);
	cout << endl;
      }
      */
      vector <bool> phaseVec;
      for (uint64 m64 = 0; m64 < Mseg64; m64++) {
	vector <bool> phaseSeg = checkHaploBits(n0, nF1, nF2, pbwtBitsFast[m64], m64, -1);
	phaseVec.insert(phaseVec.end(), phaseSeg.begin(), phaseSeg.end());
      }
      printf("FAST# major SE: %2d   # tot SE: %2d / %d\n", countMajorSE(phaseVec),
	     countSE(phaseVec), (int) phaseVec.size()-1);

      phaseVec.clear();
      for (uint64 m64 = 0; m64 < Mseg64; m64++) {
	vector <bool> phaseSeg = checkHaploBits(n0, nF1, nF2, pbwtBitsFine[m64], m64, -1);
	phaseVec.insert(phaseVec.end(), phaseSeg.begin(), phaseSeg.end());
      }
      printf("FINE# major SE: %2d   # tot SE: %2d / %d\n", countMajorSE(phaseVec),
	     countSE(phaseVec), (int) phaseVec.size()-1);


      // check accuracy of fast calls
      vector <int> trioRelPhaseVec = trioRelPhase(n0, nF1, nF2);
      const int NUM_CALL_LENGTHS = 2;
      int callLengths[NUM_CALL_LENGTHS] = {10, 20/*, 50, 100*/};
      for (int l = 0; l < NUM_CALL_LENGTHS; l++) {
	cout << "callLength = " << callLengths[l] << endl;
	vector <float> probAAs;
	for (uint64 i = 0; i < tCallLocs.size(); i++)
	  probAAs.push_back(dipTreeFast.callProbAA(tCallLocs[i].first, tCallLocs[i].second,
						   callLengths[l]));
	vector < vector <float> > probCorFastSlow;
	for (uint i = 0; i < probAAs.size(); i++)
	  if (trioRelPhaseVec[i] >= 0) {
	    vector <float> tmp(4);
	    tmp[0] = std::min(probAAs[i], 1-probAAs[i]);
	    tmp[1] = (probAAs[i] > 0.5f) == !trioRelPhaseVec[i];
	    tmp[2] = (probAAsCur[i] > 0.5f) == !trioRelPhaseVec[i];
	    probCorFastSlow.push_back(tmp);
	  }
	std::sort(probCorFastSlow.begin(), probCorFastSlow.end());
	const int NUM_PCTS = 7;
	int pcts[NUM_PCTS] = {50, 80, 90, 95, 98, 99, 100};
	for (int p = 0; p < NUM_PCTS; p++) {
	  int pct = pcts[p];
	  uint iCut = probCorFastSlow.size()*pct/100;
	  int numErrs = 0, numErrsCur = 0;
	  for (uint i = 0; i < iCut; i++) {
	    if (probCorFastSlow[i][1] == 0)
	      numErrs++;
	    if (probCorFastSlow[i][2] == 0)
	      numErrsCur++;
	  }
	  printf("  len=%d,cut=%d%%: p opp = %f, %d errs, %d cur / %d calls\n",
		 callLengths[l], pct, probCorFastSlow[iCut-1][0], numErrs, numErrsCur, iCut);
	}
      }
    }

    }

    // write phase calls
    uint64 nTargetHap = 2*(n0-Nref);
    uint64 nTargetOpp = 2*(n0-Nref) + 1;
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      tmpHaploBitsT[nTargetHap*Mseg64 + m64] = pbwtBitsFine[m64];
      tmpHaploBitsT[nTargetOpp*Mseg64 + m64] = ~pbwtBitsFine[m64];
      for (uint64 j = 0; j < 64ULL; j++) {
	uint64 m64j = m64*64+j;
	if (maskSnps64j[m64j]) {
	  if (genos64j[m64j] == 0) {
	    tmpHaploBitsT[nTargetHap*Mseg64 + m64] &= ~(1ULL<<j);
	    tmpHaploBitsT[nTargetOpp*Mseg64 + m64] &= ~(1ULL<<j);
	  }
	  else if (genos64j[m64j] == 2) {
	    tmpHaploBitsT[nTargetHap*Mseg64 + m64] |= 1ULL<<j;
	    tmpHaploBitsT[nTargetOpp*Mseg64 + m64] |= 1ULL<<j;
	  }
	}
      }
    }
    return conf;
  }

}
