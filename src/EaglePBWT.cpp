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

  vector <bool> Eagle::computeRefIsMono(const vector <uint> &bestHaps) const {
    vector <bool> refIsMono(Mseg64*64, true);
    vector <uint64> anyIs0(Mseg64), anyIs1(Mseg64);
    for (uint i = 0; i < bestHaps.size(); i++) {
      uint64 nHap = bestHaps[i];
      for (uint64 m64 = 0; m64 < Mseg64; m64++) {
	uint64 is1 = haploBitsT[nHap*Mseg64 + m64];
	anyIs0[m64] |= ~is1;
	anyIs1[m64] |= is1;
      }
    }
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      if (maskSnps64j[m64j]) {
	uint64 m64 = m64j/64, j = m64j&63;
	refIsMono[m64j] = !((anyIs0[m64]>>j)&1) || !((anyIs1[m64]>>j)&1);
      }
    }
    return refIsMono;
  }

  float Eagle::runPBWT(uint64 n0, uint64 nF1, uint64 nF2, int Kpbwt, double cMexpect,
		       double histFactor, bool runReverse, bool useTargetHaps, bool impMissing,
		       bool isChrX) {
    vector < pair <int, int> > noConPS;
    return runPBWT(n0, nF1, nF2, Kpbwt, cMexpect, histFactor, runReverse, useTargetHaps,
		   impMissing, 0, noConPS, isChrX);
  }

  float Eagle::runPBWT(uint64 n0, uint64 nF1, uint64 nF2, int Kpbwt, double cMexpect,
		       double histFactor, bool runReverse, bool useTargetHaps, bool impMissing,
		       int usePS, const vector < pair <int, int> > &conPS, bool isChrX) {
    Timer timer;
    
    vector <uint> m64jInds(Mseg64*64+1);

    const int SPEED_FACTOR = 1; const float lnPerr = logf(powf(10.0f, logPerr));
    const int CALL_LENGTH_FACTOR = 1;

    bool print = (int) nF1 != -1;
    

    /***** SELECT BEST REFERENCE HAPLOTYPES *****/

    if (print) cout << "selecting " << Kpbwt << " ref haps...    " << std::flush;
    vector <uint> bestHaps = findMinErrDipHap(n0, Kpbwt, useTargetHaps);
    // find sites at which only one allele is represented in bestHaps => can't be used as split
    vector <bool> refIsMono = computeRefIsMono(bestHaps); // size = Mseg64*64
    if (print) cout << " done " << timer.update_time() << endl;


    /***** PROCESS TARGET GENOTYPES *****/
    
    // create vectors of genos, genoBits, and hets
    vector <uchar> genos64j(Mseg64*64);
    vector <uint64> hets64j, refMonoHets64j;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      genos64j[m64j] = getGeno0123(m64j, n0);
      if (maskSnps64j[m64j] && genos64j[m64j]==1) {
	if (!refIsMono[m64j])
	  hets64j.push_back(m64j);
	else
	  refMonoHets64j.push_back(m64j);
      }
    }
    vector <uint64_masks> tgtGenoBits(Mseg64);
    for (uint64 m64 = 0; m64 < Mseg64; m64++)
      tgtGenoBits[m64] = genoBits[m64*N + n0];
    
    vector <uint64> pbwtBitsFine(Mseg64);
    float conf = 0;

    // find split sites (for PBWT HapTree starts): hets and occasionally inter-het sites
    vector < pair <int, int> > tCallLocs; vector <int> tHomLocs;
    vector <uint64> splits64j;
    const double cMmaxSplit = 0.5;
    if (!hets64j.empty()) {
      splits64j.push_back(hets64j[0]);
      for (uint64 h = 1; h < hets64j.size(); h++) {
	int lastCallLoc = splits64j.size(); // old het ind + 1: tree indices are split indices + 1

	for (uint64 m64j = hets64j[h-1]+1; m64j <= hets64j[h]; m64j++)
	  if (maskSnps64j[m64j] && !refIsMono[m64j] && genos64j[m64j] <= 2)
	    if (m64j == hets64j[h] || cMs64j[m64j] > cMs64j[splits64j.back()] + cMmaxSplit) {
	      splits64j.push_back(m64j);
	      if (m64j < hets64j[h])
		tHomLocs.push_back(splits64j.size()); // hom ind + 1
	    }
	int nextCallLoc = splits64j.size(); // new het ind + 1
	tCallLocs.push_back(make_pair(lastCallLoc, nextCallLoc));
      }
    }
    else { // all hom or missing or mono in ref; put in splits as necessary
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
	if (maskSnps64j[m64j] && !refIsMono[m64j] && genos64j[m64j] <= 2)
	  if (splits64j.empty() || cMs64j[m64j] > cMs64j[splits64j.back()] + cMmaxSplit) {
	    splits64j.push_back(m64j);
	    tHomLocs.push_back(splits64j.size()); // hom ind + 1
	  }
    }
    if (print) {
      cout << "num hets (poly in best haps): " << hets64j.size();
      cout << " num hets (mono in best haps): " << refMonoHets64j.size();
      cout << " num splits: " << splits64j.size() << endl;
    }
    
    // allocate storage for reference haplotype samples (to use for impMissing and singleton hets)
    const int samples = 10;
    vector < vector <HapPair> > refSamples[2];
    refSamples[0].resize(splits64j.size()+1);
    refSamples[1].resize(splits64j.size()+1);
    // store which inter-split chunks contain at least one missing site
    vector <bool> tHasMissing(splits64j.size()+1); // "missing" = missing or singleton het
    for (int t = 0; t <= (int) splits64j.size(); t++) {
      uint64 m64jPrev = t==0 ? -1ULL : splits64j[t-1];
      uint64 m64jNext = t==(int) splits64j.size() ? Mseg64*64 : splits64j[t];
      for (uint64 m64j = m64jPrev+1; m64j < m64jNext; m64j++)
	if (maskSnps64j[m64j] && ((genos64j[m64j] == 3 && impMissing) || // missing
				  (genos64j[m64j] == 1 && refIsMono[m64j]))) // singleton het
	  tHasMissing[t] = true;
    }

    // create vector of genos at split sites (padded on left and right to match hapBitsT)
    vector <uchar> splitGenos;
    splitGenos.push_back(0); // pad on left with 0 (to match hapBitsT)
    for (uint64 s = 0; s < splits64j.size(); s++)
      splitGenos.push_back(genos64j[splits64j[s]]);
    splitGenos.push_back(0); // pad on right with 0 (to match hapBitsT)

    // check for 0 or 1 het (warn)
    if (!isChrX && hets64j.size() <= 1) {
      cerr << "WARNING: Sample " << n0-Nref+1 << " (1-indexed) has a het count of "
	   << hets64j.size() << endl;
    }
    
    // compute recombination probabilities
    vector <double> cMcoords(splits64j.size()+2);
    for (uint64 s = 0; s <= splits64j.size(); s++) {
      uint64 splitStart = (s == 0 ? 0 : splits64j[s-1]);
      uint64 splitStop = (s == splits64j.size() ? Mseg64*64 : splits64j[s]);
      cMcoords[s] = cMs64j[splitStart]; cMcoords[s+1] = cMs64j[splitStop];
      int homs = 0;
      for (uint64 m64j = splitStart+1; m64j < splitStop; m64j++)
	if (genos64j[m64j] == 0 || genos64j[m64j] == 2)
	  homs++;
    }

    /**** BUILD PBWT DATA STRUCTURE *****/

    // create HapBitsT encoding of ref hets and hom errs
    if (print) cout << "making HapBitsT...             " << std::flush;
    HapBitsT hapBitsT(haploBitsT, Mseg64, splits64j, splitGenos, tgtGenoBits, bestHaps);
    if (print) cout << " done " << timer.update_time() << endl;

    // create HapHedge PBWT data structure
    if (print) cout << "making HapHedge...             " << std::flush;
    HapHedgeErr *hapHedgePtr = new HapHedgeErr(hapBitsT);
    if (print) cout << " done " << timer.update_time() << endl;

    //hapHedgePtr->printTree(0);

    /***** RUN COARSE (UNCONSTRAINED) DIPTREE SEARCH *****/

    // initialize DipTree object
    if (print) cout << "making DipTree (unconstr)...   " << std::flush;
    vector <char> constraints(splitGenos.size(), NO_CONSTRAINT);
    vector <int> splitInds(Mseg64*64+1); // index map for FORMAT:PS constraints
    if (usePS) {
      // populate splitInds: 1-based indices t=1..T-2 in splits64j[t-1] of 1-based SNPs m+1
      for (uint64 m64j = 0, m = 0, t = 0; m64j < Mseg64*64; m64j++)
	if (maskSnps64j[m64j]) {
	  m++;
	  m64jInds[m] = m64j;
	  if (t < splits64j.size() && splits64j[t]==m64j) {
	    t++;
	    splitInds[m] = t;
	  }
	  else
	    splitInds[m] = 0;
	}

      // set constraints for fast search
      for (uint c = 0; c < conPS.size(); c++)
	if (splitInds[conPS[c].first] && splitInds[abs(conPS[c].second)])
	  constraints[splitInds[conPS[c].first]] =
	    ((splitInds[conPS[c].first]-splitInds[abs(conPS[c].second)])<<1)|(conPS[c].second<0);
    }

    const int histLengthFast = 30*histFactor, pbwtBeamWidthFast = 30/SPEED_FACTOR;
    DipTree dipTreeFast(*hapHedgePtr, splitGenos, &constraints[0], cMcoords, cMexpect,
			histLengthFast, pbwtBeamWidthFast, lnPerr, 0);
    if (print) cout << " done " << timer.update_time() << endl;

    // explore search space; make phase calls
    if (print) cout << "making phase calls (uncon)...  " << std::flush;
    const int callLengthFast = 10*CALL_LENGTH_FACTOR;
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
    int T = splitGenos.size(); // splits64j.size()+2; tree indices are split indices + 1
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

    // set constraints for fine search
    if (usePS == 2) {
      for (uint c = 0; c < conPS.size(); c++)
	if (splitInds[conPS[c].first] && splitInds[abs(conPS[c].second)])
	  constraints[splitInds[conPS[c].first]] =
	    revConstraints[T-1-splitInds[abs(conPS[c].second)]] =
	    ((splitInds[conPS[c].first]-splitInds[abs(conPS[c].second)])<<1)|(conPS[c].second<0);
    }

    /***** RUN FINE (CONSTRAINED) DIPTREE SEARCH *****/

    // initialize DipTree object
    if (print) cout << "making DipTree (constrained)..." << std::flush;
    const int histLengthFine = 100*histFactor, pbwtBeamWidthFine = 50/SPEED_FACTOR;
    DipTree dipTreeFine(*hapHedgePtr, splitGenos, &constraints[0], cMcoords, cMexpect,
			histLengthFine, pbwtBeamWidthFine, lnPerr, 0);
    if (print) cout << " done " << timer.update_time() << endl;

    // explore search space; make phase calls
    if (print) cout << "making phase calls (constr)... " << std::flush;
    const int callLengthFine = 20*CALL_LENGTH_FACTOR;
    const int callLengthSample = 20;

    // sample refs (BEFORE callProbAA: sampleRefs needs recent history that gets overwritten)
    for (int t = 0; t < T-1; t++)
      if (tHasMissing[t])
	refSamples[0][t] = dipTreeFine.sampleRefs(t, callLengthSample, samples, bestHaps, true);

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
      vector <double> revcMcoords(T);
      for (int t = 0; t < T; t++) revcMcoords[t] = cMcoords[T-1] - cMcoords[T-1-t];

      // initialize DipTree object
      if (print) cout << "making revDipTree (constr)...  " << std::flush;
      DipTree revDipTreeFine(revHapHedge, revSplitGenos, &revConstraints[0], revcMcoords, cMexpect,
			     histLengthFine, pbwtBeamWidthFine, lnPerr, 0);
      if (print) cout << " done " << timer.update_time() << endl;

      // explore search space; make phase calls
      if (print) cout << "making rev phase calls (con)..." << std::flush;

      // sample refs (BEFORE callProbAA: sampleRefs needs recent history that gets overwritten)
      for (int t = T-2; t >= 0; t--)
	if (tHasMissing[t])
	  refSamples[1][t] =
	    revDipTreeFine.sampleRefs(T-2-t, callLengthSample, samples, bestHaps, false);

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

    // write phase calls
    uint64 nTargetHap = 2*(n0-Nref);
    uint64 nTargetOpp = 2*(n0-Nref) + 1;
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      tmpHaploBitsT[nTargetHap*Mseg64 + m64] = pbwtBitsFine[m64];
      tmpHaploBitsT[nTargetOpp*Mseg64 + m64] = ~pbwtBitsFine[m64];
    }
    
    // set hom bits
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
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

    // impute missing genos and phase hets monomorphic in bestHaps
    for (int t = 0; t <= (int) splits64j.size(); t++) {
      if (!tHasMissing[t]) continue;

      if (!runReverse) {
	// no reverse samples available; just use fwd samples and don't bother with ends
	for (int s = 0; s < samples; s++)
	  for (int h = 0; h < 2; h++)
	    refSamples[0][t][s].haps[h].isEnd = false;
	refSamples[1][t] = refSamples[0][t];
      }
	
      // identify missing sites
      const uint64 m64jPrev = t==0 ? -1ULL : splits64j[t-1];
      const uint64 m64jNext = t==(int) splits64j.size() ? Mseg64*64 : splits64j[t];
      vector <uint64> miss64j;
      for (uint64 m64j = m64jPrev+1; m64j < m64jNext; m64j++)
	if (maskSnps64j[m64j] && genos64j[m64j] == 3 && impMissing) // missing
	  miss64j.push_back(m64j);

      // orient each hap pair wrt called phase
      if (!hets64j.empty())
	for (int fb = 0; fb < 2; fb++) {
	  for (int s = 0; s < samples; s++) {
	    double cMmid;
	    if (t==0) cMmid = cMs64j[splits64j[t]];
	    else if (t==(int) splits64j.size()) cMmid = cMs64j[splits64j[t-1]];
	    else cMmid = (cMs64j[splits64j[t-1]] + cMs64j[splits64j[t]]) / 2;

	    vector < pair <double, int> > cMdistSigns;
	    for (int i = 15; i >= 0; i--)
	      if ((((refSamples[fb][t][s].haps[0].tMaskRev ^
		     refSamples[fb][t][s].haps[1].tMaskRev)>>i)&1)
		  && splitGenos[t-i]==1) {
		uint64 m64j = splits64j[t-i-1];
		cMdistSigns.push_back(make_pair(fabs(cMs64j[m64j]-cMmid), ((refSamples[fb][t][s].haps[0].tMaskRev>>i)&1) == ((pbwtBitsFine[m64j/64]>>(m64j&63))&1)));
	      }
	    for (int i = 0; i < 16; i++)
	      if ((((refSamples[fb][t][s].haps[0].tMaskFwd ^
		     refSamples[fb][t][s].haps[1].tMaskFwd)>>i)&1)
		  && splitGenos[t+i+1]==1) {
		uint64 m64j = splits64j[t+i];
		cMdistSigns.push_back(make_pair(fabs(cMs64j[m64j]-cMmid), ((refSamples[fb][t][s].haps[0].tMaskFwd>>i)&1) == ((pbwtBitsFine[m64j/64]>>(m64j&63))&1)));
	      }
	    if (!cMdistSigns.empty()) {
	      sort(cMdistSigns.begin(), cMdistSigns.end());
	      if (!cMdistSigns[0].second)
		std::swap(refSamples[fb][t][s].haps[0], refSamples[fb][t][s].haps[1]);
	    }
	  }
	}

      // for each haplotype in turn, call missing sites (and save mean hap length)
      double hMeanLens[2];
      for (int h = 0; h < 2; h++) {
	vector <int> endRefs[2];
	vector <int> nonEndRefs[2], nonEndLens[2];
	// split sampled haplotypes into buckets: those that end in (t,t+1) and those that don't
	for (int fb = 0; fb < 2; fb++)
	  for (int s = 0; s < samples; s++) {
	    if (refSamples[fb][t][s].haps[h].isEnd)
	      endRefs[fb].push_back(refSamples[fb][t][s].haps[h].refSeq);
	    else {
	      nonEndRefs[fb].push_back(refSamples[fb][t][s].haps[h].refSeq);
	      nonEndLens[fb].push_back(refSamples[fb][t][s].haps[h].tLength);
	    }
	  }

	// initialize allele dosages
	int nMiss = miss64j.size();
	struct DosePair { double d[2]; };
	vector <DosePair> alleleDoses(nMiss);
	for (int m = 0; m < nMiss; m++) alleleDoses[m].d[0] = alleleDoses[m].d[1] = 0;

	// process non-ends: call using longer of fwd, rev ref samples
	double meanLens[2] = {0, 0};
	for (int fb = 0; fb < 2; fb++)
	  if (!nonEndLens[fb].empty())
	    meanLens[fb] = std::accumulate(nonEndLens[fb].begin(), nonEndLens[fb].end(), 0) /
	      (double) nonEndLens[fb].size();
	int fbLong = meanLens[0] > meanLens[1] ? 0 : 1;
	for (uint k = 0; k < nonEndRefs[fbLong].size(); k++) {
	  int refSeq = nonEndRefs[fbLong][k];
	  for (int m = 0; m < nMiss; m++)
	    alleleDoses[m].d[(haploBitsT[refSeq*Mseg64+miss64j[m]/64]>>(miss64j[m]&63))&1] += 1;
	}
	  
	// compute mean lengths including ends=0 for phasing singletons later
	for (int fb = 0; fb < 2; fb++)
	  meanLens[fb] = std::accumulate(nonEndLens[fb].begin(), nonEndLens[fb].end(), 0) /
	    (double) samples;
	hMeanLens[h] = std::max(meanLens[0], meanLens[1]); // set to max of fwd, rev

	// process ends: find most likely recombination points at which fwd and rev haps meet
	if (!endRefs[0].empty() && !endRefs[1].empty()) {
	  for (uint k = 0; k < endRefs[0].size() && k < endRefs[1].size(); k++) {
	    int refSeqFwd = endRefs[0][k], refSeqRev = endRefs[1][k];
	    int errFwd = 0, errRev = 0;
	    for (uint64 m64j = m64jPrev+1; m64j < m64jNext; m64j++) {
	      if ((genos64j[m64j] == 0 || genos64j[m64j] == 2) &&
		  (((haploBitsT[refSeqRev*Mseg64+m64j/64]>>(m64j&63))&1) != genos64j[m64j]/2))
		errRev++;
	    }
	      
	    // find recombination points that minimize errors (usually 0 errors)
	    int minErr = 1<<30;
	    vector <double> cMdiffs; vector <uint64> revStarts;

	    for (uint64 m64j = m64jPrev; m64j == m64jPrev || m64j < m64jNext; m64j++) {
	      if (m64j != m64jPrev) { // update err counts
		if ((genos64j[m64j] == 0 || genos64j[m64j] == 2) &&
		    (((haploBitsT[refSeqFwd*Mseg64+m64j/64]>>(m64j&63))&1) != genos64j[m64j]/2))
		  errFwd++;
		if ((genos64j[m64j] == 0 || genos64j[m64j] == 2) &&
		    (((haploBitsT[refSeqRev*Mseg64+m64j/64]>>(m64j&63))&1) != genos64j[m64j]/2))
		  errRev--;
	      }
		
	      // rev starts at m64j+1
	      double cMdiff = cMs64j[m64j+1] - (m64j==-1ULL ? 0 : cMs64j[m64j]) + 1e-9;

	      if (errFwd+errRev < minErr) {
		minErr = errFwd+errRev;
		cMdiffs.clear();
		revStarts.clear();
	      }
	      if (errFwd+errRev == minErr) {
		cMdiffs.push_back(cMdiff);
		revStarts.push_back(m64j+1);
	      }
	    }

	    // augment dosages proportionally to cMdiffs (btwn consecutive SNPs) at recomb points
	    double cMtot = std::accumulate(cMdiffs.begin(), cMdiffs.end(), 0.0);
	    for (int m = 0; m < nMiss; m++) {
	      double cMcum = 0;
	      for (uint x = 0; x < revStarts.size(); x++) {
		if (miss64j[m] < revStarts[x])
		  cMcum += cMdiffs[x];
		else
		  break;
	      }
	      alleleDoses[m].d[(haploBitsT[refSeqFwd*Mseg64+miss64j[m]/64]>>(miss64j[m]&63))&1] +=
		cMcum / cMtot;
	      alleleDoses[m].d[(haploBitsT[refSeqRev*Mseg64+miss64j[m]/64]>>(miss64j[m]&63))&1] +=
		1 - cMcum / cMtot;
	    }
	  }
	}

	// make final calls
	for (int m = 0; m < nMiss; m++) {
	  uint64 m64 = miss64j[m]/64, j = miss64j[m]&63;
	  if (alleleDoses[m].d[0] >= alleleDoses[m].d[1])
	    tmpHaploBitsT[(nTargetHap+h)*Mseg64 + m64] &= ~(1ULL<<j);
	  else
	    tmpHaploBitsT[(nTargetHap+h)*Mseg64 + m64] |= 1ULL<<j;
	}
      }

      // call phase at "singleton" hets monomorphic among bestHaps
      for (uint64 m64j = m64jPrev+1; m64j < m64jNext; m64j++)
	if (genos64j[m64j] == 1 && refIsMono[m64j]) { // "singleton" het
	  uint64 m64 = m64j/64, j = m64j&63;
	  uint64 commonBit = haploBitsT[bestHaps[0]*Mseg64 + m64] & (1ULL<<j);
	  uint64 rareBit = commonBit ^ (1ULL<<j);

	  int hShorter = hMeanLens[0] < hMeanLens[1] ? 0 : 1; // put rare allele on shorter hap

	  for (int h = 0; h < 2; h++)
	    tmpHaploBitsT[(nTargetHap+h)*Mseg64 + m64] &= ~(1ULL<<j); // clear bit	    
	  tmpHaploBitsT[(nTargetHap+hShorter)*Mseg64 + m64] |= rareBit;
	  tmpHaploBitsT[(nTargetHap+!hShorter)*Mseg64 + m64] |= commonBit;
	  /*
	  cout << "common bit: " << commonBit << endl;
	  cout << "rare bit: " << rareBit << endl;
	  cout << "hMeanLens[0]: " << hMeanLens[0] << endl;
	  cout << "hMeanLens[1]: " << hMeanLens[1] << endl;
	  cout << "hShorter: " << hShorter << endl;
	  */
	}
    }
    /*
    for (int t = 0; t < 500; t += 100) {
      cout << "==== t: " << t << " ====" << endl;
      for (int fb = 0; fb < 2; fb++) {
	cout << "--- fb: " << fb << " ---" << endl;
	for (int s = 0; s < 3; s++) {
	  cout << ".. sample: " << s << " ..     " << refSamples[fb][t][s].haps[0].tLength << "," << refSamples[fb][t][s].haps[1].tLength << endl;
	  for (int h = 0; h < 2; h++) {
	    for (int i = 15; i >= 0; i--) {
	      if ((((refSamples[fb][t][s].haps[h].tMaskRev ^ refSamples[fb][t][s].haps[1-h].tMaskRev)>>i)&1) && splitGenos[t-i]==1) {
		uint64 m64j = splits64j[t-i-1];
		cout << (((refSamples[fb][t][s].haps[h].tMaskRev>>i)&1) == ((pbwtBitsFine[m64j/64]>>(m64j&63))&1) ? "+" : "-");
	      }
	      else
		cout << ".";
	    }
	    cout << "|";
	    for (int i = 0; i < 16; i++) {
	      if ((((refSamples[fb][t][s].haps[h].tMaskFwd ^ refSamples[fb][t][s].haps[1-h].tMaskFwd)>>i)&1) && splitGenos[t+i+1]==1) {
		uint64 m64j = splits64j[t+i];
		cout << (((refSamples[fb][t][s].haps[h].tMaskFwd>>i)&1) == ((pbwtBitsFine[m64j/64]>>(m64j&63))&1) ? "+" : "-");

	      }
	      else
		cout << ".";
	    }
	    cout << "   ";
	  }
	  cout << endl;
	}
      }
    }
    */

    if (print && usePS) {
      int correct = 0;
      for (uint c = 0; c < conPS.size(); c++) {
	int m1 = conPS[c].first, m2 = abs(conPS[c].second), isOpp = conPS[c].second<0;
	uint m64j1 = m64jInds[m1], m64j2 = m64jInds[m2];
	assert(genos64j[m64j1]==1 && genos64j[m64j2]==1);
	correct += (uint) isOpp == (((pbwtBitsFine[m64j1/64]>>(m64j1&63))&1) ^
				    ((pbwtBitsFine[m64j2/64]>>(m64j2&63))&1));
      }
      cout << "constraints respected: " << correct << " / " << conPS.size() << endl;
    }

    return conf;
  }

}
