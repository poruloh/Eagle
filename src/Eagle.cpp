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
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <utility>
#include <numeric>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "omp.h"

#include <htslib/thread_pool.h>
#include <htslib/vcf.h>

#include "Types.hpp"
#include "FileUtils.hpp"
#include "MemoryUtils.hpp"
#include "NumericUtils.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"
#include "HapHedge.hpp"
#include "Version.hpp"
#include "Eagle.hpp"

//#define DETAILS

namespace EAGLE {

  using std::vector;
  using std::string;
  using std::pair;
  using std::make_pair;
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::max;
  using std::min;

  const double MEMO_UNSET = -1000;
  const char noTrioInfo = '`';
  const string trio1 = "\033[1;36mo\033[0m";
  const string trio2 = "\033[1;31mx\033[0m";
  const char IBDx2char = '_';
  const char ROHchar = '=';
  const char wrongChar = '@';
  const char conflictChar = '#';

  inline uint popcount64(uint64 i) {
    i = i - ((i >> 1) & 0x5555555555555555);
    i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
    i = (i + (i >> 4)) & 0xF0F0F0F0F0F0F0F;
    return (i * 0x101010101010101) >> 56;
  }

  void Eagle::init() {
    totTicks = 0; extTicks = 0; diphapTicks = 0; lshTicks = 0; lshCheckTicks = 0;
    dpTicks = 0; dpStaticTicks = 0; dpSwitchTicks = 0; dpUpdateTicks = 0; dpSortTicks = 0;
    dpUpdateCalls = 0; blipFixTicks = 0; blipPopTicks = 0; blipVoteTicks = 0; blipLshTicks = 0;

    maskSnps64j = ALIGNED_MALLOC_UCHARS(Mseg64*64);
    cMs64j = ALIGNED_MALLOC_DOUBLES(Mseg64*64+1);
    double cMlast = 0;
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      for (uint64 j = 0; j < seg64cMvecs[m64].size(); j++) {
	maskSnps64j[m64*64+j] = 1;
	cMs64j[m64*64+j] = cMlast = seg64cMvecs[m64][j];
      }
      for (uint64 j = seg64cMvecs[m64].size(); j < 64; j++) {
	maskSnps64j[m64*64+j] = 0;
	cMs64j[m64*64+j] = cMlast;
      }
    }
    cMs64j[Mseg64*64] = cMlast;

    haploBits = ALIGNED_MALLOC_UINT64S(Mseg64*2*N);
    haploBitsT = ALIGNED_MALLOC_UINT64S(2*N*Mseg64);
    segConfs = ALIGNED_MALLOC_UCHARS(2*N*Mseg64);

    maskIndivs = vector <uchar> (N, 1);

    for (uint wrongBitsA = 0; wrongBitsA < (1U<<switchScoreLutBits); wrongBitsA++)
      for (uint wrongBitsB = 0; wrongBitsB < (1U<<switchScoreLutBits); wrongBitsB++)
	for (uint hetBits = 0; hetBits < (1U<<switchScoreLutBits); hetBits++) {
	  uint wrongHomBitsA = wrongBitsA & ~hetBits;
	  uint wrongHetBitsA = wrongBitsA & hetBits;
	  uint wrongHomBitsB = wrongBitsB & ~hetBits;
	  uint wrongHetBitsB = wrongBitsB & hetBits;
	  uint lutInd = (wrongBitsA<<(switchScoreLutBits+switchScoreLutBits))
	    | (wrongBitsB<<(switchScoreLutBits))
	    | hetBits;
	  char &minDiff = switchScoreLut[lutInd][0]; minDiff = 0;
	  char &cumDiff = switchScoreLut[lutInd][1]; cumDiff = 0;
	  for (uint k = 0; k < switchScoreLutBits; k++) {
	    cumDiff += ((wrongHomBitsA>>k)&1)*homErrCost + ((wrongHetBitsA>>k)&1)*hetErrCost
	      - ((wrongHomBitsB>>k)&1)*homErrCost - ((wrongHetBitsB>>k)&1)*hetErrCost;
	    if (cumDiff < minDiff) minDiff = cumDiff;
	  }
	}
  }

  Eagle::Eagle(uint64 _N, uint64 _Mseg64, const uint64_masks *_genoBits,
	       vector < vector <double> > _seg64cMvecs, const AlleleFreqs *_seg64logPs,
	       vector <double> _invLD64j, const vector <IndivInfoX> &_indivs,
	       const vector <SnpInfoX> &_snps, const string &maskFile,
	       const vector <bool> &_isFlipped64j, double _pErr, int runStep2) :
    N(_N), Nref(0), Mseg64(_Mseg64), genoBits(_genoBits), seg64cMvecs(_seg64cMvecs),
    seg64logPs(_seg64logPs), invLD64j(_invLD64j), indivs(_indivs), snps(_snps),
    isFlipped64j(_isFlipped64j), logPerr(log10(_pErr)) {

    init();

    if (runStep2) {
      phaseConfs = ALIGNED_MALLOC_UCHARS(2*N*Mseg64*64);
      phaseConfs2 = ALIGNED_MALLOC_UCHARS(2*N*Mseg64*64);
      tmpHaploBitsT = NULL;
    }
    else {
      phaseConfs = phaseConfs2 = NULL;
      tmpHaploBitsT = ALIGNED_MALLOC_UINT64S(2*N*Mseg64);
      memset(tmpHaploBitsT, 0, 2*N*Mseg64*sizeof(tmpHaploBitsT[0]));
    }

    if (!maskFile.empty()) {
      int masked = 0;
      vector < pair <string, string> > maskFidIids = FileUtils::readFidIids(maskFile);
      std::set < pair < string, string> > maskSet(maskFidIids.begin(), maskFidIids.end());
      for (uint64 n = 0; n < N; n++)
	if (maskSet.count(make_pair(indivs[n].famID, indivs[n].indivID))) {
	  maskIndivs[n] = 0;
	  masked++;
	}
      cout << "Number of indivs masked: " << masked << endl;
    }
  }

  // constructor for ref-mode
  Eagle::Eagle(uint64 _Nref, uint64 _Ntarget, uint64 _Mseg64, const uint64_masks *_genoBits,
	       vector < vector <double> > _seg64cMvecs, double _pErr) :
    N(_Nref+_Ntarget), Nref(_Nref), Mseg64(_Mseg64), genoBits(_genoBits),
    seg64cMvecs(_seg64cMvecs), logPerr(log10(_pErr)) {

    init();
    isFlipped64j = vector <bool> (Mseg64*64); // no flipping in ref mode

    phaseConfs = phaseConfs2 = NULL;
    tmpHaploBitsT = ALIGNED_MALLOC_UINT64S(2*(N-Nref)*Mseg64);

    memset(segConfs, 0, 2*N*Mseg64*sizeof(segConfs[0]));
    for (uint64 nRef = 0; nRef < Nref; nRef++)
      for (uint64 m64 = 0; m64 < Mseg64; m64++)	{ // copy ref haploBits stored in genoBits
	haploBits[m64*2*N + 2*nRef] = genoBits[m64*N + nRef].is0;
	haploBits[m64*2*N + 2*nRef+1] = genoBits[m64*N + nRef].is2;
	for (uint64 nHap = 2*nRef; nHap <= 2*nRef+1; nHap++)
	  haploBitsT[nHap*Mseg64 + m64] = haploBits[m64*2*N + nHap];
      }
  }

  void Eagle::reallocLRPtoPBWT(void) { // non-ref mode transition: LRP iters 1-2 -> PBWT iters 3+
    assert(phaseConfs != NULL);
    ALIGNED_FREE(phaseConfs2); phaseConfs2 = NULL;
    ALIGNED_FREE(phaseConfs); phaseConfs = NULL;

    assert(tmpHaploBitsT == NULL);
    tmpHaploBitsT = ALIGNED_MALLOC_UINT64S(2*N*Mseg64);
  }

  Eagle::~Eagle() {
    ALIGNED_FREE(segConfs);
    ALIGNED_FREE(haploBitsT);
    ALIGNED_FREE(haploBits);
    if (phaseConfs != NULL) {
      ALIGNED_FREE(phaseConfs2);
      ALIGNED_FREE(phaseConfs);
    }
    if (tmpHaploBitsT != NULL) {
      ALIGNED_FREE(tmpHaploBitsT); // allocated only in ref-mode
    }
    ALIGNED_FREE(cMs64j);
    ALIGNED_FREE(maskSnps64j);
  }

  inline uint getNonMissingGeno(const uint64_masks &bits, uint64 j) {
    if (bits.is0 & (1ULL<<j)) return 0;
    if (bits.is2 & (1ULL<<j)) return 2;
    return 1; // assumed to be non-missing
  }

  inline uint bgetGeno0123(const uint64_masks &bits, uint64 j) {
    if (bits.is0 & (1ULL<<j)) return 0;
    if (bits.is2 & (1ULL<<j)) return 2;
    if (bits.is9 & (1ULL<<j)) return 3;
    return 1;
  }

  uint Eagle::getGeno0123(uint64 m64j, uint64 n) const {
    return bgetGeno0123(genoBits[m64j/64*N + n], m64j&63);
  }

  void Eagle::retractMatch(uint n0, Match &match, double memoLogBF[][4]) const {
    for (int dir = 0; dir < 2; dir++) {
      double cumLogBF = 0;
      while (cumLogBF < log10(4)) {
	uint m64j;
	if (dir == 0) m64j = match.m64jStart++;
	else m64j = match.m64jEnd--;
	cumLogBF += memoLogBF[m64j][getGeno0123(m64j, match.n)];
      }
    }
    match.m64jStart--;
    match.m64jEnd++;
  }

  Match Eagle::computeDuoLogBF(double memoLogBF[][4], double workLogBF[], uint64 n0, uint64 n1, uint64 m64cur) const {
    //double snpsChecked = 1;
    Match match(n1, m64cur*64, m64cur*64, 0);
    workLogBF[m64cur*64] = 0;
    for (int dir = 0; dir < 2; dir++) {
      uint64 inc, m64start, m64end, jStart, jEnd;
      if (dir == 0) {
	inc = 1; m64start = m64cur; m64end = Mseg64; jStart = 0; jEnd = 64;
      }
      else {
	inc = -1ULL; m64start = m64cur-1; m64end = -1ULL; jStart = 63; jEnd = -1ULL;
      }
      double maxLogBF = 0, curLogBF = 0;
      for (uint64 m64 = m64start; m64 != m64end; m64 += inc) {
	const uint64_masks &bits0 = genoBits[m64*N + n0], &bits1 = genoBits[m64*N + n1];
	uint64 wrongBits = (bits0.is0 & bits1.is2) | (bits0.is2 & bits1.is0);
	if (wrongBits & (wrongBits-1)) // 2+ wrong => fail
	  break;
	uint64 missMask = bits0.is9 | bits1.is9;
	for (uint64 j = jStart; j != jEnd; j += inc) {
	  uint geno1 = bgetGeno0123(bits1, j);
	  double logBFj = 0;
	  if (!(missMask & (1ULL<<j))) {
	    double &memoLogBFj = memoLogBF[m64*64+j][geno1];
	    if (memoLogBFj == MEMO_UNSET) {
	      uint geno0 = getNonMissingGeno(bits0, j);
	      double logP_geno1_null = seg64logPs[m64*64+j].cond[geno1][3];
	      double logP_geno1_duo = seg64logPs[m64*64+j].cond[geno1][geno0];
	      memoLogBFj = min(max((logP_geno1_duo - logP_geno1_null) * invLD64j[m64*64+j],
				   logPerr), -logPerr);
	    }
	    logBFj = memoLogBFj;
	  }
	  workLogBF[m64*64+j] = logBFj;
	  curLogBF += logBFj;
	  //if (logBFj != 0) snpsChecked += invLD64j[m64*64+j];
	  if (curLogBF > maxLogBF) {
	    maxLogBF = curLogBF;
	    if (inc == 1) match.m64jEnd = m64*64+j;
	    else match.m64jStart = m64*64+j;
	  }
	}
      }
      match.logBF += maxLogBF;
    }
    double minLogBF = 0, curLogBF = 0;
    for (uint64 m64j = match.m64jStart; m64j <= match.m64jEnd; m64j++) {
      curLogBF += workLogBF[m64j];
      if (curLogBF < minLogBF) {
	minLogBF = curLogBF;
	match.m64jStart = m64j+1;
	while (!maskSnps64j[match.m64jStart]) match.m64jStart++;
      }
    }
    match.logBF -= minLogBF;
    //match.logBF -= log10(snpsChecked);
    match.cMlenInit = cMs64j[match.m64jEnd] - cMs64j[match.m64jStart];
    return match;
  }

  void Eagle::trim(Match &match, const Match &ref, uint64 n0, int orientation, uint64 trimStart,
		   int inc, double workLogBF[]) const {

    uint64 n1 = match.n, n2 = ref.n;

    // find IBDx2; store IBDx2 status in workLogBF (0 or 1) to compute probabilities accordingly
    double IBDx2logBF = 0; uint64 IBDx2start = trimStart;
    for (uint64 m64j = trimStart; m64j+1!=match.m64jStart && m64j!=match.m64jEnd+1; m64j += inc)
      workLogBF[m64j] = 0; // initialize workLogBF to 0 (not IBDx2) in *MATCH* (n1)
    // check for IBDx2 in *REF* (n2)
    uint64 m64jLast; // last SNP to check: go beyond end of match to ensure detection of ref IBDx2
    if (inc == 1)
      m64jLast = min(match.m64jEnd + 50ULL, Mseg64*64-1);
    else
      m64jLast = max((int) match.m64jStart - 50, 0);
    for (uint64 m64j = trimStart; m64j!=m64jLast+inc; m64j += inc) {
      uint g0 = getGeno0123(m64j, n0);
      uint g2 = getGeno0123(m64j, n2);
      bool mismatch = false;
      if (g0 != 3 && g2 != 3) {
	if (g0 == g2)
	  IBDx2logBF += seg64logPs[m64j].cond[g2][g0] * invLD64j[m64j];
	else
	  mismatch = true;
      }
      if (mismatch || m64j==m64jLast) { // end of IBDx2 segment
	if (IBDx2logBF < -1) { // 10:1 IBDx2
#ifdef VERBOSE
	  printf("IBDx2 detected in %d: %.1f-%.1f (%d SNPs)\n", (int) n2, cMs64j[IBDx2start], cMs64j[m64j], (int) (m64j-IBDx2start));
#endif
	  for (uint64 m64j2 = IBDx2start; m64j2 != m64j; m64j2 += inc)
	    workLogBF[m64j2] = 1;
	}
	IBDx2logBF = 0; // reset
	IBDx2start = m64j+inc;
      }
    }

    double maxLogBF = 0, curLogBF = 0; uint64 m64jBest = trimStart;
    for (uint64 m64j = trimStart; m64j+1!=match.m64jStart && m64j!=match.m64jEnd+1; m64j += inc) {
      double logBF = 0;
      uint g0 = getGeno0123(m64j, n0);
      uint g1 = getGeno0123(m64j, n1);
      uint g2 = getGeno0123(m64j, n2);
      if (g0 != 3 && g1 != 3) {
	uint g0eff = g0;
	if (g0 == 1 && ref.m64jStart <= m64j && m64j <= ref.m64jEnd && g2 != 3) {
	  // n0 het and n2 (ref) match info available
	  if (g2 != 1) { // n2 hom
	    if (orientation == 1)
	      g0eff = g2; // treat g0 as n2 hom
	    else
	      g0eff = 2-g2; // treat g0 as opp n2 hom
	  }
	  else if (workLogBF[m64j] == 0) { // n0 and n2 both hets and not IBDx2
	    if (orientation == 1)
	      g0eff = 4; // same orientation as het-het => p(hap=1) = 1-p
	    else
	      g0eff = 5; // opp orientation to het-het => p(hap=1) = p
	  }
	}
	double logP_geno1_null = seg64logPs[m64j].cond[g1][3];
	double logP_geno1_duo = seg64logPs[m64j].cond[g1][g0eff];
	logBF = min(max((logP_geno1_duo - logP_geno1_null) * invLD64j[m64j],
			logPerr), -logPerr);
      }
      curLogBF += logBF;
      if (curLogBF > maxLogBF) {
	maxLogBF = curLogBF;
	m64jBest = m64j;
      }
      workLogBF[m64j] = curLogBF;
    }

    uint64 m64jTrim = m64jBest;
    if (inc == 1) match.m64jEnd = m64jTrim;
    else match.m64jStart = m64jTrim;

    // conservative trimming: backtrack to 10x higher prob
    uint m64jTrimCons = trimStart;
    for (uint64 m64j = trimStart; m64j != m64jBest; m64j += inc)
      if (workLogBF[m64j] < maxLogBF - log10(10))
	m64jTrimCons = m64j;
    if (inc == 1) match.m64jEndCons = std::min(match.m64jEndCons, m64jTrimCons);
    else match.m64jStartCons = std::max(match.m64jStartCons, m64jTrimCons);
  }

  vector <int> searchSigns(const vector <Match> &matches, const vector < vector <uint> > &sameEdges, const vector < vector <uint> > &oppEdges, const vector <bool> &kept) {
    // process left to right so that when sign choice is arbitrary, adjacent matches have same sign
    uint V = matches.size();
    vector <int> signs(V);
    vector < pair <uint, uint> > order(V);
    for (uint v = 0; v < V; v++)
      order[v] = make_pair(matches[v].m64jStart, v);
    sort(order.begin(), order.end());
    uint lastEnd = 0; int lastSign = 1; // sign of farthest-right match seen so far
    std::queue <uint> q;
    for (uint i = 0; i < V; i++) {
      uint v = order[i].second;
      if (!kept[v]) continue; // not used
      if (signs[v]) continue; // already visited
      signs[v] = lastSign;
      q.push(v);
      while (!q.empty()) {
	uint u = q.front(); q.pop();
	for (uint i = 0; i < sameEdges[u].size(); i++) {
	  uint w = sameEdges[u][i];
	  if (!kept[w]) continue;
	  if (signs[w]) {
	    if (signs[w] != signs[u])
	      return vector <int> ();
	  }
	  else {
	    signs[w] = signs[u];
	    q.push(w);
	    if (matches[w].m64jEnd > lastEnd) {
	      lastEnd = matches[w].m64jEnd;
	      lastSign = signs[w];
	    }
	  }
	}
	for (uint i = 0; i < oppEdges[u].size(); i++) {
	  uint w = oppEdges[u][i];
	  if (!kept[w]) continue;
	  if (signs[w]) {
	    if (signs[w] == signs[u])
	      return vector <int> ();
	  }
	  else {
	    signs[w] = -signs[u];
	    q.push(w);
	    if (matches[w].m64jEnd > lastEnd) {
	      lastEnd = matches[w].m64jEnd;
	      lastSign = signs[w];
	    }
	  }
	}
      }
    }
    return signs;
  }

  void updateVote(int &votesCur, int votesThresh, int vote) {
    if (abs(votesCur) >= votesThresh) return;
    votesCur += vote;
  }

  void Eagle::computePhaseConfs(uint64 n0, const vector <Match> &matches,
				const vector <int> &signs, bool cons) {

    vector < vector <int> > votes(2, vector <int> (Mseg64*64));
    vector <int> votesThresh(Mseg64*64);
    const int votesMax = 1000000;

    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      if (maskSnps64j[m64j]) {
        votesThresh[m64j] = // 2 / log10((1-p)/p) = number of votes needed to get >100:1 odds
	  (int) (2 / fabs(seg64logPs[m64j].cond[0][4] - seg64logPs[m64j].cond[0][5])) + 1;
	if (!(votesThresh[m64j] < votesMax)) votesThresh[m64j] = votesMax;
	uint g0 = getGeno0123(m64j, n0);
	if (g0 == 0) { votes[0][m64j] = votes[1][m64j] = -votesMax; }
	else if (g0 == 2) { votes[0][m64j] = votes[1][m64j] = votesMax; }
	else if (g0 == 3) { votes[0][m64j] = votes[1][m64j] = -1; } // missing: default P(1) = p
      }
    }

    vector <uchar> isIBDx2(Mseg64*64);

    for (uint i = 0; i < matches.size(); i++) {
      if (!signs[i]) continue;
      uint64 start, end;
      if (cons) {
	start = std::max(matches[i].m64jStartCons, matches[i].m64jStart);
	end = std::min(matches[i].m64jEndCons, matches[i].m64jEnd);
      }
      else {
	start = matches[i].m64jStart;
	end = matches[i].m64jEnd;
      }
#ifdef VERBOSE
      printf("match %d (%.1f-%.1f)\n", (int) i, cMs64j[start], cMs64j[end]);
#endif

      // find IBDx2 regions
      vector < pair <uint64, uint64> > IBDx2regions;
      uint64 m64jFirst = max((int) start - 50, 0); // go beyond ends to detect overhanging IBDx2
      uint64 m64jLast = min(end + 50ULL, Mseg64*64-1);
      double IBDx2logBF = 0; uint64 IBDx2start = m64jFirst;
      for (uint64 m64j = m64jFirst; m64j <= m64jLast; m64j++) {
	uint g0 = getGeno0123(m64j, n0);
	uint g1 = getGeno0123(m64j, matches[i].n);
	bool mismatch = false;
	if (g0 != 3 && g1 != 3) {
	  if (g0 == g1)
	    IBDx2logBF += seg64logPs[m64j].cond[g1][g0] * invLD64j[m64j];
	  else
	    mismatch = true;
	}
	if (mismatch || m64j==m64jLast) { // end of IBDx2 segment
	  if (IBDx2logBF < -1) { // 10:1 IBDx2
	    IBDx2regions.push_back(make_pair(IBDx2start, m64j));
#ifdef VERBOSE
	    printf("IBDx2 detected in %d: %.1f-%.1f (%d SNPs)\n", (int) matches[i].n, cMs64j[IBDx2start], cMs64j[m64j], (int) (m64j-IBDx2start));
#endif
	  }
	  IBDx2logBF = 0; // reset
	  IBDx2start = m64j+1;
	}
      }

      for (uint r = 0; r < IBDx2regions.size(); r++) // set IBDx2 region flags
	memset(&isIBDx2[IBDx2regions[r].first], 1, IBDx2regions[r].second-IBDx2regions[r].first);

      // accumulate votes
      for (uint64 m64j = start; m64j <= end; m64j++) {
	if (!maskSnps64j[m64j]) continue;
	uint g0 = getGeno0123(m64j, n0);
	uint g1 = getGeno0123(m64j, matches[i].n);

	int vote = 0;
	if (g1 == 0 || g1 == 2) // n1 hom: IBDx2 status irrelevant; phase determined (votesMax)
	  vote = (g1-1)*votesMax*2; // super strong vote (overrides any previous small votes)
	else if (g1 == 1 && !isIBDx2[m64j]) // n1 het and not IBDx2; weak phase info
	  vote = 1;

	if (vote) {
	  int q = (signs[i] == 1);
	  updateVote(votes[q][m64j], votesThresh[m64j], vote);
	  if (g0 == 1) // n0 het: pass info to opp chromosome
	    updateVote(votes[!q][m64j], votesThresh[m64j], -vote);
	}
      }

      for (uint r = 0; r < IBDx2regions.size(); r++) // unset IBDx2 region flags
	memset(&isIBDx2[IBDx2regions[r].first], 0, IBDx2regions[r].second-IBDx2regions[r].first);
    }

    // fast rng: last 16 bits of Marsaglia's MWC
    uint w = 521288629;
    if (phaseConfs != NULL) { // need to make hard calls here
      for (uint i = 0; i < (n0 & 0xff); i++)
	w=18000*(w&65535)+(w>>16);
    }

    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      if (maskSnps64j[m64j]) {
	for (uint64 q = 0; q <= 1; q++) {
	  double phaseConf;
	  if (votes[q][m64j] >= votesMax)
	    phaseConf = 1;
	  else if (votes[q][m64j] <= -votesMax)
	    phaseConf = 0;
	  else {
	    double OR = pow(10.0, fabs(seg64logPs[m64j].cond[0][4] - seg64logPs[m64j].cond[0][5])
			    * votes[q][m64j]);
	    phaseConf = OR / (1 + OR);
	  }
	  if (phaseConfs != NULL)
	    phaseConfs[(2*n0+q)*Mseg64*64 + m64j] = (uchar) (phaseConf * 255);
	  else {
	    uchar uPhaseConf = (uchar) (phaseConf * 255);
	    if (uPhaseConf == (uchar) 255 || ((w=18000*(w&65535)+(w>>16))&255) < uPhaseConf)
	      tmpHaploBitsT[(2*n0+q)*Mseg64 + (m64j/64)] |= 1ULL<<(m64j&63);
	  }
	}
      }
      else {
	if (phaseConfs != NULL)
	  phaseConfs[2*n0*Mseg64*64 + m64j] = phaseConfs[(2*n0+1)*Mseg64*64 + m64j] = 0;
      }
    }
  }

  vector <int> Eagle::trioRelPhase(uint64 n0, uint64 nF1, uint64 nF2) const {

    bool isParent = false;
    if (((int) nF1) < 0) { // nF1 is the child; n0 is a parent
      isParent = true;
      nF1 = -nF1;
    }

    vector <int> trioPhaseVec;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      int trioPhase = 0;
      if (!isParent) {
	if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
	if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      }
      else { // n0 is a parent; nF1 is the child
	int g0 = getGeno0123(m64j, nF1); // child
	int g2 = getGeno0123(m64j, nF2); // other parent
	if (g0+g2 != 2) { // not Mendel error or triple het
	  if (g0 == 0) trioPhase = -1;
	  if (g0 == 2) trioPhase = 1;
	  if (g0 == 1) { // child is a het
	    if (g2 == 0) trioPhase = 1;
	    if (g2 == 2) trioPhase = -1;
	  }
	}
      }
      trioPhaseVec.push_back(trioPhase); // 0 => unknown; +/-1 => pat/mat
      if (trioPhase == 0) continue;
    }
    vector <int> trioRelPhaseVec(trioPhaseVec.size()-1);
    for (uint i = 1; i < trioPhaseVec.size(); i++)
      trioRelPhaseVec[i-1] = (trioPhaseVec[i-1]==0 || trioPhaseVec[i]==0) ? -1 :
	(trioPhaseVec[i-1]==trioPhaseVec[i] ? 0 : 1); // -1 => unknown; 0 => same; 1 => opp
    return trioRelPhaseVec;
  }

  void Eagle::checkPhase(uint64 n0, uint64 nF1, uint64 nF2, double thresh) const {
    cout << "checking at thresh=" << thresh << ": ";

    double lastPhased = cMs64j[0];
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      if ((bits0.is0|bits0.is2)&(1ULL<<j)) continue; // hom
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      double phaseConf = phaseConfs[2*n0*Mseg64*64 + m64j] / 255.0;
      if (std::min(phaseConf, 1-phaseConf) <= thresh) {
	double cM = cMs64j[m64j];
	for (int tick = (int) (10*lastPhased) + 1; tick < 10*cM; tick++) {
	  if (tick % 10 == 0) cout << StringUtils::itos(tick/10); //(char) ('0' + (tick/10)%10);
	  else cout << '-';
	}
	if ((phaseConf < 0.5) == (trioPhase == 1))
	  cout << trio1;
	else
	  cout << trio2;
	lastPhased = cM;
      }
      else
	cout << '?';
    }
    cout << endl;
  }

  vector <bool> Eagle::checkPhaseConfsPhase(uint64 n0, uint64 nF1, uint64 nF2) const {
    vector <bool> ret;
    int lastPhased64j = -1;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // hom
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      int hapBit = (int) phaseConfs[2*n0*Mseg64*64 + m64j] >= 128;
      bool phase = hapBit == (trioPhase == 1);
      if (lastPhased64j != -1 && ret.back() != phase)
	printf(" %.2f", (cMs64j[lastPhased64j] + cMs64j[m64j]) / 2);
      lastPhased64j = m64j;
      ret.push_back(phase);
    }
    cout << endl;
    return ret;
  }

  void Eagle::checkHapPhase(uint64 n0, uint64 nF1, uint64 nF2, const uint64 curHaploBitsT[],
			    uint64 m64, uint64 side, vector < vector <int> > votes) const {
    if ((int) nF1 == -1) return;
    for (uint64 m64j = (m64-side)*64; m64j < (m64+side+1)*64; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      //if (((haploBits[m64cur*2*N + n1hap]>>j)&1) == (trioPhase == 1))
      //if (((haploBitsT[n1hap*Mseg64 + m64cur]>>j)&1) == (trioPhase == 1))
      if (((curHaploBitsT[m64cur]>>j)&1) == (trioPhase == 1))
	cout << trio1;
      else
	cout << trio2;
      if (!votes.empty())
	cout << "[" << votes[j][(curHaploBitsT[m64cur]>>j)&1] << "|" << votes[j][!((curHaploBitsT[m64cur]>>j)&1)] << ";" << votes[j][((curHaploBitsT[m64cur]>>j)&1)+2] << "|" << votes[j][!((curHaploBitsT[m64cur]>>j)&1)+2] << "]";
    }
    cout << endl;
  }

  vector <bool> Eagle::checkHapPhase1(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap,
				      uint64 m64start, uint64 m64end, int sign) const {
    vector <bool> ret;
    if ((int) nF1 == -1) return ret;
    cout << "n1hap = " << n1hap << "; m64 = [" << m64start << "," << m64end << "): ";
    for (uint64 m64j = m64start*64; m64j < m64end*64; m64j++) {
      if (m64j != m64start*64 && (m64j&63)==0) cout << "|";
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      if (((haploBits[m64cur*2*N + n1hap]>>j)&1) == (trioPhase == sign)) {
	cout << trio1;
	ret.push_back(0);
      }
      else {
	cout << trio2;
	ret.push_back(1);
      }
    }
    cout << endl;
    return ret;
  }

  vector <bool> Eagle::checkHapPhase1j(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap,
				       uint64 m64jStart, uint64 m64jEnd, int sign) const {
    vector <bool> ret;
    if ((int) nF1 == -1) return ret;
    //cout << "n1hap = " << n1hap << "; m64 = [" << m64start << "," << m64end << "): ";
    for (uint64 m64j = m64jStart; m64j < m64jEnd; m64j++) {
      if (m64j != m64jStart && (m64j&63)==0) cout << m64j/64;//"|";
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      if (((haploBits[m64cur*2*N + n1hap]>>j)&1) == (trioPhase == sign)) {
	cout << trio1;
	ret.push_back(0);
      }
      else {
	cout << trio2;
	ret.push_back(1);
      }
    }
    //cout << endl;
    return ret;
  }

  vector <bool> Eagle::checkHapPhase1jCall(uint64 n0, uint64 nF1, uint64 nF2, uint64 callBitsT[],
					   uint64 m64jStart, uint64 m64jEnd, bool print, int sign) const {
    vector <bool> ret;
    if ((int) nF1 == -1) return ret;
    for (uint64 m64j = m64jStart; m64j < m64jEnd; m64j++) {
      if (m64j != m64jStart && (m64j&63)==0)
	if (print) cout << m64j/64;//"|";
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      if (((callBitsT[m64cur]>>j)&1) == (trioPhase == sign)) {
	if (print) cout << trio1;
	ret.push_back(0);
      }
      else {
	if (print) cout << trio2;
	ret.push_back(1);
      }
    }
    if (print) cout << endl;
    return ret;
  }

  int Eagle::checkHapPhase2(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap,
				      uint64 n2hapA, uint64 n2hapB, uint64 m64, int sign) const {
    vector <bool> ret;
    if ((int) nF1 == -1) return 0;/*ret*/;

    uint64 n1is1 = haploBitsT[n1hap*Mseg64 + m64];
    uint64 n2is1A = haploBitsT[n2hapA*Mseg64 + m64];
    uint64 n2is1B = haploBitsT[n2hapB*Mseg64 + m64];
    const uint64_masks &bits0 = genoBits[m64*N + n0];
    uint64 wrongHomBitsA = (bits0.is0 & (n1is1 | n2is1A)) | (bits0.is2 & ~(n1is1 & n2is1A));
    uint64 wrongHetBitsA = (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1A));
    uint64 wrongHomBitsB = (bits0.is0 & (n1is1 | n2is1B)) | (bits0.is2 & ~(n1is1 & n2is1B));
    uint64 wrongHetBitsB = (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1B));
    uint score = popcount64(wrongHomBitsB)*homErrCost
      + popcount64(wrongHetBitsB)*hetErrCost;
    uint minScore = score; uint64 kSwitch = 0;
    for (uint64 k = 0; k < 64; k++) {
      score += ((wrongHomBitsA>>k)&1)*homErrCost + ((wrongHetBitsA>>k)&1)*hetErrCost
	- ((wrongHomBitsB>>k)&1)*homErrCost - ((wrongHetBitsB>>k)&1)*hetErrCost;
      if (score < minScore) {
	minScore = score;
	kSwitch = k+1;
      }
    }

    cout << "m64 = " << m64 << ": (" << n1hap << "," << n2hapA;
    if (n2hapA != n2hapB) cout << "-" << n2hapB;
    cout << ") score = " << minScore << " ";
    cout << " conf = " << (int) segConfs[n1hap*Mseg64+m64] << "," << (int) segConfs[n2hapA*Mseg64+m64];
    if (n2hapA != n2hapB) cout << "-" << (int) segConfs[n2hapB*Mseg64+m64];
    cout << " ";

    for (uint64 j = 0; j < 64; j++) {
      uint64 m64j = m64*64+j;
      if (!maskSnps64j[m64j]) continue;
      const uint64_masks &bits0 = genoBits[m64*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      bool phase;
      if (((haploBits[m64*2*N + n1hap]>>j)&1) == (trioPhase == sign))
	phase = 0;
      else
	phase = 1;
      bool hetErr =
	((haploBits[m64*2*N + n1hap]>>j)&1) ==
	((haploBits[m64*2*N + (j<kSwitch?n2hapA:n2hapB)]>>j)&1);
      uchar conf1 = 0, conf2 = 0;
      if (hetErr) {
	conf1 = phaseConfs[n1hap*Mseg64*64 + m64j];
	conf2 = phaseConfs[(j<kSwitch?n2hapA:n2hapB)*Mseg64*64 + m64j];
	if (min((int) conf2, 255-conf2) < min((int) conf1, 255-conf1))
	  phase = !phase;
      }
      cout << (phase==0?trio1:trio2);
      if (hetErr) cout << "?" << "[" << (int) conf1 << "|" << (int) conf2 << "]";
      ret.push_back(phase);
    }
    //cout << endl;
    return minScore/*ret*/;
  }

  vector <bool> Eagle::checkHaploBits(uint64 n0, uint64 nF1, uint64 nF2, uint64 hapBits,
				      uint64 m64, int pad) const {
    vector <bool> ret;
    if ((int) nF1 == -1) return ret;

    bool isParent = false;
    if (((int) nF1) < 0) { // nF1 is the child; n0 is a parent
      isParent = true;
      nF1 = -nF1;
    }

    int printed = 0;
    for (uint64 j = 0; j < 64; j++) {
      uint64 m64j = m64*64+j;
      if (!maskSnps64j[m64j]) continue;
      const uint64_masks &bits0 = genoBits[m64*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64*N + nF2];
      int trioPhase = 0;
      if (!isParent) {
	if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
	if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      }
      else { // n0 is a parent; nF1 is the child
	int g0 = getGeno0123(m64j, nF1); // child
	int g2 = getGeno0123(m64j, nF2); // other parent
	if (g0+g2 != 2) { // not Mendel error or triple het
	  if (g0 == 0) trioPhase = -1;
	  if (g0 == 2) trioPhase = 1;
	  if (g0 == 1) { // child is a het
	    if (g2 == 0) trioPhase = 1;
	    if (g2 == 2) trioPhase = -1;
	  }
	}
      }
      if (!trioPhase) continue;
      if (pad >= 0) {
	if (((hapBits>>j)&1) == (trioPhase == 1))
	  cout << trio1;
	else
	  cout << trio2;
	printed++;
      }
      ret.push_back(((hapBits>>j)&1) == (trioPhase == 1));
    }
    while (printed < pad) { cout << " "; printed++; }
    return ret;
  }

  pair <uint64, uint64> Eagle::phaseSegHMM(uint64 n0, uint64 n1hap, uint64 n2hapA, uint64 n2hapB,
					   uint64 m64, uint64 &hetErrMask) const {
    uint64 n1is1 = haploBitsT[n1hap*Mseg64 + m64];
    uint64 n2is1A = haploBitsT[n2hapA*Mseg64 + m64];
    uint64 n2is1B = haploBitsT[n2hapB*Mseg64 + m64];
    const uint64_masks &bits0 = genoBits[m64*N + n0];
    uint64 wrongHomBitsA = (bits0.is0 & (n1is1 | n2is1A)) | (bits0.is2 & ~(n1is1 & n2is1A));
    uint64 wrongHetBitsA = (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1A));
    uint64 wrongHomBitsB = (bits0.is0 & (n1is1 | n2is1B)) | (bits0.is2 & ~(n1is1 & n2is1B));
    uint64 wrongHetBitsB = (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1B));
    uint score = popcount64(wrongHomBitsB)*homErrCost
      + popcount64(wrongHetBitsB)*hetErrCost;
    uint minScore = score;
    double cMdiffMax = m64==0 ? 0.0 : cMs64j[m64*64] - cMs64j[m64*64-1];
    uint64 kSwitch = 0;
    uint64 kSeg = seg64cMvecs[m64].size();
    for (uint64 k = 0; k < kSeg; k++) {
      score += ((wrongHomBitsA>>k)&1)*homErrCost + ((wrongHetBitsA>>k)&1)*hetErrCost
	- ((wrongHomBitsB>>k)&1)*homErrCost - ((wrongHetBitsB>>k)&1)*hetErrCost;
      double cMdiff = (k+1==kSeg ? cMs64j[(m64+1)*64] : cMs64j[m64*64+k+1]) - cMs64j[m64*64+k];
      if (score < minScore || (score == minScore && cMdiff > cMdiffMax)) {
	minScore = score;
	cMdiffMax = cMdiff;
	kSwitch = k+1;
      }
    }

    uint64 phaseBits1 = 0, phaseBits2 = 0; hetErrMask = 0;
    for (uint64 j = 0; j < 64; j++) {
      uint64 m64j = m64*64+j;
      if (!maskSnps64j[m64j]) continue;
      const uint64_masks &bits0 = genoBits[m64*N + n0];
      if (bits0.is0&(1ULL<<j)) // dip = 0: hap1 = hap2 = 0
	;
      else if (bits0.is2&(1ULL<<j)) { // dip = 2: hap1 = hap2 = 1
	phaseBits1 |= 1ULL<<j;
	phaseBits2 |= 1ULL<<j;
      }
      else {
	uint64 phase1 = (haploBits[m64*2*N + n1hap]>>j)&1;
        uint64 phase2 = ((haploBits[m64*2*N + (j<kSwitch?n2hapA:n2hapB)]>>j)&1);
	if (bits0.is9&(1ULL<<j)) { // missing
	  phaseBits1 |= phase1<<j;
	  phaseBits2 |= phase2<<j;
	}
	else { // het
	  bool phase = phase1;
	  bool hetErr = phase1 == phase2;
	  if (hetErr) {
	    hetErrMask |= 1ULL<<j;
	    if (Nref == 0) {
	      uchar conf1 = phaseConfs[n1hap*Mseg64*64 + m64j];
	      uchar conf2 = phaseConfs[(j<kSwitch?n2hapA:n2hapB)*Mseg64*64 + m64j];
	      if (min((int) conf2, 255-conf2) < min((int) conf1, 255-conf1))
		phase = !phase;
	    }
	  }
	  phaseBits1 |= ((uint64) phase)<<j;
	  phaseBits2 |= ((uint64) !phase)<<j;
	}
      }
    }
    return make_pair(phaseBits1, phaseBits2);
  }

  vector <bool> Eagle::checkSegPhase(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap, uint64 n2hap,
				     int sign, uint64 m64) const {
    vector <bool> ret;
    for (uint64 j = 0; j < 64; j++) {
      uint64 m64j = m64*64+j;
      if (!maskSnps64j[m64j]) continue;
      const uint64_masks &bits0 = genoBits[m64*N + n0];
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
      const uint64_masks &bitsF1 = genoBits[m64*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64*N + nF2];
      int trioPhase = 0;
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase++;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase--;
      if (!trioPhase) continue;
      /*
      int hapBit1 = (haploBits[m64*2*N + n1hap]>>j)&1;
      if (sign == -1) hapBit1 = 1-hapBit1;
      int hapBit2 = (haploBits[m64*2*N + n2hap]>>j)&1;
      if (sign == 1) hapBit2 = 1-hapBit2;
      uchar conf1 = phaseConfs[n1hap*Mseg64*64 + m64j];
      uchar conf2 = phaseConfs[n2hap*Mseg64*64 + m64j];
      //if (hapBit1 != hapBit2) cout << "?[" << (uint) conf1 << "|" << (uint) conf2 << "]";

      int hapBitFinal = hapBit1;
      if (!(conf1 == 0 || conf1 == 255) && (conf2 == 0 || conf2 == 255)) hapBitFinal = hapBit2;
      */
      int hapBitFinal = (int) phaseConfs2[2*n0*Mseg64*64 + m64j] >= 128;
      if (hapBitFinal == (trioPhase == 1)) {
	cout << trio1;
	ret.push_back(0);
      }
      else {
	cout << trio2;
	ret.push_back(1);
      }

    }
    cout << endl;
    return ret;
  }

  void Eagle::computeSegPhaseConfs(uint64 n0, uint64 n1hap, uint64 n2hap, int sign, uint64 m64,
				   int err) {
    const double maxOffsetFrac = 0.25, offsetFrac = ((n0&63)-31.5)/31.5 * maxOffsetFrac;
    uint64 m64jStart, m64jEnd;
    if (offsetFrac >= 0) {
      m64jStart = m64==0 ? 0 : m64*64 + (uint) (offsetFrac * seg64cMvecs[m64].size());
      m64jEnd = (m64+1)*64 + (m64+1==Mseg64 ? 0 : (uint) (offsetFrac * seg64cMvecs[m64+1].size()));
    }
    else {
      m64jStart = m64==0 ? 0 : (m64-1)*64 + (uint) ((1+offsetFrac) * seg64cMvecs[m64-1].size());
      m64jEnd = m64*64 + (m64+1==Mseg64 ? 64 : (uint) ((1+offsetFrac) * seg64cMvecs[m64].size()));
    }
    /*
    cout << "m64 = " << m64 << ": " << m64jStart/64 << "." << (m64jStart&63) << " - "
	 << m64jEnd/64 << "." << (m64jEnd&63) << endl;
    */
    int cropErr = max(err-1, 0);
    uchar hetConfs[2]; hetConfs[0] = (uchar) cropErr; hetConfs[1] = (uchar) (255-cropErr);
    for (uint64 m64j = m64jStart; m64j < m64jEnd; m64j++) {
      if (maskSnps64j[m64j]) {
	uint64 m64cur = m64j/64; uint64 j = m64j&63;
	uint g0 = getGeno0123(m64j, n0);
	if (g0 == 0)
	  phaseConfs2[2*n0*Mseg64*64 + m64j] = phaseConfs2[(2*n0+1)*Mseg64*64 + m64j] = 0;
	else if (g0 == 2)
	  phaseConfs2[2*n0*Mseg64*64 + m64j] = phaseConfs2[(2*n0+1)*Mseg64*64 + m64j] = 255;
	else {
	  int hapBit1 = (haploBits[m64cur*2*N + n1hap]>>j)&1;
	  int hapBit2 = (haploBits[m64cur*2*N + n2hap]>>j)&1;
	  if (g0 == 1) { // het (not missing)
	    uchar conf1 = phaseConfs[n1hap*Mseg64*64 + m64j];
	    uchar conf2 = phaseConfs[n2hap*Mseg64*64 + m64j];
	    if (!(conf1 == 0 || conf1 == 255) && (conf2 == 0 || conf2 == 255))
	      hapBit1 = 1-hapBit2; // only use n2hap if n1hap conf <100% and n2hap conf = 100%
	    else
	      hapBit2 = 1-hapBit1; // default: go with n1hap
	  }
	  if (sign == -1) std::swap(hapBit1, hapBit2);
	  phaseConfs2[2*n0*Mseg64*64 + m64j] = hetConfs[hapBit1];
	  phaseConfs2[(2*n0+1)*Mseg64*64 + m64j] = hetConfs[hapBit2];
	}
      }
      else
	phaseConfs2[2*n0*Mseg64*64 + m64j] = phaseConfs2[(2*n0+1)*Mseg64*64 + m64j] = 0;
    }
  }

  string Eagle::computePhaseString(uint64 n0, uint64 nF1, uint64 nF2,
				   const vector <Match> &matches, const vector <int> &signs,
				   uint64 start, double cMend, bool cons)
    const {

    if ((int) nF1 == -1 || (int) nF2 == -1) return "";
    vector <uint64> starts(matches.size()), ends(matches.size());
    for (uint i = 0; i < matches.size(); i++) {
      if (cons) {
	starts[i] = std::max(matches[i].m64jStartCons, matches[i].m64jStart);
	ends[i] = std::min(matches[i].m64jEndCons, matches[i].m64jEnd);
      }
      else {
	starts[i] = matches[i].m64jStart;
	ends[i] = matches[i].m64jEnd;
      }
    }
    string phase;
    double lastPhased = cMs64j[start], lastHet = lastPhased;
    int hetCount = 0, snpCount = 0;
    for (uint64 m64j = start; m64j < Mseg64*64 && cMs64j[m64j] < cMend; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];

      int trioPhase = 1; // default
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase = 1;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase = -1;
      int votes1 = 0, votes2 = 0;
      snpCount++;
      if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) {
	bool wrong = false;
	for (uint i = 0; i < matches.size(); i++)
	  if (signs[i] && starts[i] <= m64j && m64j <= ends[i]) {
	    const uint64_masks &bits1 = genoBits[m64cur*N + matches[i].n];
	    if (((bits0.is0&bits1.is2)|(bits0.is2&bits1.is0))&(1ULL<<j)) wrong = true;
	  }
	if (wrong)
	  phase += wrongChar;
	continue; // discard non-hets
      }
      for (uint i = 0; i < matches.size(); i++)
	if (signs[i] && starts[i] <= m64j && m64j <= ends[i]) {
	  const uint64_masks &bits1 = genoBits[m64cur*N + matches[i].n];
	  int phase1 = 0;
	  if (bits1.is0&(1ULL<<j)) phase1 = 1 * signs[i] * trioPhase;
	  if (bits1.is2&(1ULL<<j)) phase1 = -1 * signs[i] * trioPhase;
	  if (phase1 == 1) votes1++;
	  else if (phase1 == -1) votes2++;
	}
      trioPhase = 0; // change default
      if ((bitsF1.is0|bitsF2.is2)&(1ULL<<j)) trioPhase = 1;
      if ((bitsF1.is2|bitsF2.is0)&(1ULL<<j)) trioPhase = -1;

      if (votes1+votes2) {
	double cM = cMs64j[m64j];
	for (int tick = (int) (10*lastPhased) + 1; tick < 10*cM; tick++) {
	  if (tick % 10 == 0) phase += StringUtils::itos(tick/10); //(char) ('0' + (tick/10)%10);
	  else if (lastHet <= lastPhased) phase += ROHchar;
	  else phase += IBDx2char;
	}
	if (cM-lastPhased > 0.5) {
	  char buf[20]; sprintf(buf, "[%.1fcM:%d/%d]", cM-lastPhased, hetCount, snpCount);
	  phase += string(buf);
	}
	if (votes1&votes2)
	  phase += conflictChar;
	else if (trioPhase) {
	  if (votes1)
	    phase += trio1;
	  else
	    phase += trio2;
	}
	else
	  phase += noTrioInfo;

	hetCount = 0; snpCount = 0;
	lastPhased = cM;
      }
      else
	phase += '?';

      lastHet = cMs64j[m64j];
      hetCount++;
    }

    if (cMend >= cMs64j[Mseg64*64-1]) {
      double cM = cMs64j[Mseg64*64-1];
      for (int tick = (int) (10*lastPhased) + 1; tick < 10*cM; tick++) {
	if (tick % 10 == 0) phase += (char) ('0' + (tick/10)%10);
	else if (lastHet <= lastPhased) phase += ROHchar;
	else phase += IBDx2char;
      }
      if (cM-lastPhased > 0.5) {
	char buf[20]; sprintf(buf, "[%.1fcM:%d/%d]", cM-lastPhased, hetCount, snpCount);
	phase += string(buf);
      }
    }
    return phase;
  }

  void Eagle::printMatch(uint64 n0, uint64 nF1, uint64 nF2, const Match &duoMatch,
			 double memoLogBF[][4]) const {
    if ((int) nF1 == -1 || (int) nF2 == -1) return;
    uint64 n1 = duoMatch.n;
    uint64 x = duoMatch.m64jStart, y = duoMatch.m64jEnd;
    double logBF = duoMatch.logBF;
    int same = 0, opp = 0; string phase;
    double lastPhased = cMs64j[x], maxPhasedGap = 0, lastHet = cMs64j[x], maxHetGap = 0;
    int hetCount = 0, snpCount = 0, numErr = 0;
    for (uint m64j = duoMatch.m64jStart; m64j <= duoMatch.m64jEnd; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      const uint64_masks &bits1 = genoBits[m64cur*N + n1];
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      uint64 wrongBits = (bits0.is0 & bits1.is2) | (bits0.is2 & bits1.is0);
      if (wrongBits&(1ULL<<j)) {
	phase += wrongChar;
	numErr++;
      }

      uint64 trioPhased = ~(bits0.is0|bits0.is2|bits0.is9) &
	((bitsF1.is0^bitsF2.is0) | (bitsF1.is2^bitsF2.is2));
      uint64 phased1 = ~(bits0.is0|bits0.is2|bits0.is9) & (bits1.is0|bits1.is2);

      if (phased1&(1ULL<<j)) {
	double cM = cMs64j[m64j];
	for (int tick = (int) (10*lastPhased) + 1; tick < 10*cM; tick++) {
	  if (tick % 10 == 0) phase += (char) ('0' + (tick/10)%10);
	  else if (lastHet <= lastPhased) phase += ROHchar;
	  else phase += IBDx2char;
	}
	if (cM-lastPhased > 0.5) {
	  char buf[20]; sprintf(buf, "[%.1fcM:%d/%d]", cM-lastPhased, hetCount, snpCount);
	  phase += string(buf);
	}
	if (trioPhased&(1ULL<<j)) {
	  if ((bits1.is0&(1ULL<<j)) == ((bitsF1.is0|bitsF2.is2)&(1ULL<<j))) {
	    same++;
	    phase += trio1;
	  }
	  else {
	    opp++;
	    phase += trio2;
	  }
	}
	else
	  phase += noTrioInfo;

	hetCount = 0; snpCount = 0;
	if (cM - lastPhased > maxPhasedGap)
	  maxPhasedGap = cM - lastPhased;
	lastPhased = cM;
      }
      else
	snpCount++;

      if ((~(bits0.is0|bits0.is2|bits0.is9))&(1ULL<<j)) {
	double cM = cMs64j[m64j];
	if (cM - lastHet > maxHetGap)
	  maxHetGap = cM - lastHet;
	lastHet = cM;
	hetCount++;
      }
    }
    double cM = cMs64j[y];
    if (cM - lastPhased > maxPhasedGap)
      maxPhasedGap = cM - lastPhased;
    lastPhased = cM;
    if (cM - lastHet > maxHetGap)
      maxHetGap = cM - lastHet;
    lastHet = cM;

    const uint WINDOW = 20; double minLogBFwindow = 0, minLoc = 0;
    for (uint64 wStart = x; wStart+WINDOW <= y; wStart++) {
      if (!maskSnps64j[wStart]) continue;
      double logBFwindow = 0;
      for (uint64 m64j = wStart; m64j < wStart+WINDOW; m64j++) {
	const uint64_masks &bits1 = genoBits[m64j/64*N + n1];
	uint geno1 = bgetGeno0123(bits1, m64j&63);
	if (memoLogBF[m64j][geno1] != MEMO_UNSET) logBFwindow += memoLogBF[m64j][geno1];
      }
      if (logBFwindow < minLogBFwindow) {
	minLogBFwindow = logBFwindow;
	minLoc = cMs64j[wStart];
      }
    }

    printf("n0=%-5d n1=%-5d (%s) BF= %.1f cM= %.1f (%.1f-%.1f):", (int) n0, (int) n1,
	   Nref==0 ? indivs[n1].indivID.c_str() : "", logBF, cMs64j[y]-cMs64j[x], cMs64j[x],
	   cMs64j[y]);
    cout << " (" << numErr << " errs) ";
    cout << same << "|" << opp;
    cout << " " << phase;
    cout << " max gap: " << maxPhasedGap;
    cout << " max ROH: " << maxHetGap;
    printf(" min window logBF: %.1f minLoc: %.1fcM", minLogBFwindow, minLoc);
    cout << endl;
  }

  void Eagle::checkTrioErrorRate(uint64 n0, uint64 nF1, uint64 nF2) const {

    if ((int) nF1 == -1) return;
    int numOppHom1 = 0, numOppHom2 = 0, numWrongHet = 0, snpCount = 0;
    for (uint m64j = 0; m64j < Mseg64*64; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      snpCount++;
      uint64 m64cur = m64j/64; uint64 j = m64j&63;
      const uint64_masks &bits0 = genoBits[m64cur*N + n0];
      const uint64_masks &bitsF1 = genoBits[m64cur*N + nF1];
      const uint64_masks &bitsF2 = genoBits[m64cur*N + nF2];
      uint64 oppHom1 = (bits0.is0 & bitsF1.is2) | (bits0.is2 & bitsF1.is0);
      uint64 oppHom2 = (bits0.is0 & bitsF2.is2) | (bits0.is2 & bitsF2.is0);
      uint64 wrongHet = ~(bits0.is0|bits0.is2|bits0.is9) &
	((bitsF1.is0&bitsF2.is0) | (bitsF1.is2&bitsF2.is2));
      if (oppHom1&(1ULL<<j)) numOppHom1++;
      if (oppHom2&(1ULL<<j)) numOppHom2++;
      if (wrongHet&(1ULL<<j)) numWrongHet++;
    }
    cout << "oppHom1: " << numOppHom1 << " oppHom2: " << numOppHom2 << " wrongHet: " << numWrongHet
	 << " / " << snpCount << endl;
  }

  void Eagle::findLongHalfIBD(uint64 n0, vector <uint> topInds[2], vector <uint> topIndsLens[2],
			      uint K) const {
    /*
    for (uint e = 0; e < 2; e++) {
      topInds[e] = vector <uint> (Mseg64 * K);
      topIndsLens[e] = vector <uint> (Mseg64);
    }
    uint *runStarts[2][2]; // [even/odd][max err]; lengths are N
    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	runStarts[p][e] = ALIGNED_MALLOC_UINTS(N);
    uint *runStartFreqs[2][2]; // [even/odd][max err]; lengths are Mseg64+1 (b/c m64+1 can go over)
    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	runStartFreqs[p][e] = ALIGNED_MALLOC_UINTS(Mseg64+1);

    // initialize "prev" p=1
    for (uint e = 0; e < 2; e++) {
      const uint p = 1;
      memset(runStarts[p][e], 0, N*sizeof(runStarts[p][e][0])); // all runs start at 0
      runStartFreqs[p][e][0] = N;
      memset(&runStartFreqs[p][e][1], 0, Mseg64*sizeof(runStarts[p][e][1]));
    }

    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      // initialize cur [p=(m64&1)] arrays
      const uint p = m64&1;
      for (uint e = 0; e < 2; e++) {
	memcpy(runStarts[p][e], runStarts[!p][e], N*sizeof(runStarts[p][e][0]));
	memset(runStartFreqs[p][e], 0, (Mseg64+1)*sizeof(runStartFreqs[p][e][0]));
      }

      uint64_masks bits0 = genoBits[m64*N + n0];
      for (uint64 n1 = 0; n1 < N; n1++) {
	const uint64_masks &bits1 = genoBits[m64*N + n1];
	uint64 wrongBits = (bits0.is0 & bits1.is2) | (bits0.is2 & bits1.is0);
	// update cur runStarts based on bit matches
	if (wrongBits) {
	  if (wrongBits & (wrongBits-1)) // 2+ wrong => fail; move both starts forward
	    runStarts[p][0][n1] = runStarts[p][1][n1] = m64+1;
	  else { // 1 wrong => move 1-err start to 0-err start; fail 0-err
	    runStarts[p][1][n1] = runStarts[p][0][n1];
	    runStarts[p][0][n1] = m64+1;
	  }
	}
	// compute cur runStartFreqs based on (updated) cur runStarts
	for (uint e = 0; e < 2; e++)
	  runStartFreqs[p][e][runStarts[p][e][n1]]++;
      }

      for (uint e = 0; e < 2; e++)
	for (uint s64 = 0; s64 <= m64; s64++)
	  if ((runStartFreqs[p][e][s64] >= K && m64+1 == Mseg64) || // reached end
	      (runStartFreqs[p][e][s64] < K && runStartFreqs[!p][e][s64] >= K)) {
	    // all n1 s.t. runStarts[p][e][n1] == s64 are in top K starting from s64
	    // some n1 s.t. runStarts[!p][e][n1] == s64 are in top K starting from s64
	    uint allPos = 0, somePos = runStartFreqs[p][e][s64];
	    uint *topIndsChunk = &topInds[e][s64*K];
	    for (uint64 n1 = 0; n1 < N; n1++) {
	      if (allPos < K && runStarts[p][e][n1] == s64)
		topIndsChunk[allPos++] = n1;
	      else if (somePos < K && runStarts[!p][e][n1] == s64)
		topIndsChunk[somePos++] = n1;
	    }
	    topIndsLens[e][s64] = min(somePos, K);
	  }
    }

    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	ALIGNED_FREE(runStartFreqs[p][e]);
    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	ALIGNED_FREE(runStarts[p][e]);
    */
    for (uint e = 0; e < 2; e++) {
      topInds[e] = vector <uint> (Mseg64 * K);
      topIndsLens[e] = vector <uint> (Mseg64);
    }
    uint *runStarts[2]; // [even/odd]; lengths are N
    for (uint p = 0; p < 2; p++)
      runStarts[p] = ALIGNED_MALLOC_UINTS(N);
    uint *runStartFreqs[2]; // [even/odd]; lengths are Mseg64+1 (b/c m64+1 can go over)
    for (uint p = 0; p < 2; p++)
      runStartFreqs[p] = ALIGNED_MALLOC_UINTS(Mseg64+1);

    // initialize "prev" p=1
    {
      const uint p = 1;
      memset(runStarts[p], 0, N*sizeof(runStarts[p][0])); // all runs start at 0
      runStartFreqs[p][0] = N;
      memset(&runStartFreqs[p][1], 0, Mseg64*sizeof(runStarts[p][1]));
    }

    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      // initialize cur [p=(m64&1)] arrays
      const uint p = m64&1;
      {
	memcpy(runStarts[p], runStarts[!p], N*sizeof(runStarts[p][0]));
	memset(runStartFreqs[p], 0, (Mseg64+1)*sizeof(runStartFreqs[p][0]));
      }

      uint64_masks bits0 = genoBits[m64*N + n0];
      for (uint64 n1 = 0; n1 < N; n1++) {
	const uint64_masks &bits1 = genoBits[m64*N + n1];
	uint64 wrongBits = (bits0.is0 & bits1.is2) | (bits0.is2 & bits1.is0);
	// update cur runStarts based on bit matches
	if (wrongBits) {
	  runStarts[p][n1] = m64+1;
	}
	// compute cur runStartFreqs based on (updated) cur runStarts
	runStartFreqs[p][runStarts[p][n1]]++;
      }

	for (uint s64 = 0; s64 <= m64; s64++)
	  if ((runStartFreqs[p][s64] >= K && m64+1 == Mseg64) || // reached end
	      (runStartFreqs[p][s64] < K && runStartFreqs[!p][s64] >= K)) {
	    // all n1 s.t. runStarts[p][n1] == s64 are in top K starting from s64
	    // some n1 s.t. runStarts[!p][n1] == s64 are in top K starting from s64
	    uint allPos = 0, somePos = runStartFreqs[p][s64];
	    uint *topIndsChunk = &topInds[1][s64*K];
	    for (uint64 n1 = 0; n1 < N; n1++) {
	      if (allPos < K && runStarts[p][n1] == s64)
		topIndsChunk[allPos++] = n1;
	      else if (somePos < K && runStarts[!p][n1] == s64)
		topIndsChunk[somePos++] = n1;
	    }
	    topIndsLens[1][s64] = min(somePos, K);
	  }
    }

    for (uint p = 0; p < 2; p++)
      ALIGNED_FREE(runStartFreqs[p]);
    for (uint p = 0; p < 2; p++)
      ALIGNED_FREE(runStarts[p]);
  }

  void Eagle::findLongDipHap(uint64 n0, vector <uint> topInds[2], vector <uint> topIndsLens[2],
			     uint K, uint errStart=1) const {

    uint64 Nhaps = 2*(Nref==0 ? N : Nref);
    for (uint e = 0; e < 2; e++) {
      topInds[e] = vector <uint> (Mseg64 * K);
      topIndsLens[e] = vector <uint> (Mseg64);
    }
    uint *runStarts[2][2]; // [even/odd][max err]; lengths are Nhaps
    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	runStarts[p][e] = ALIGNED_MALLOC_UINTS(Nhaps);
    uint *runStartFreqs[2][2]; // [even/odd][max err]; lengths are Mseg64+1 (b/c m64+1 can go over)
    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	runStartFreqs[p][e] = ALIGNED_MALLOC_UINTS(Mseg64+1);

    // initialize "prev" p=1
    for (uint e = 0; e < 2; e++) {
      const uint p = 1;
      memset(runStarts[p][e], 0, Nhaps*sizeof(runStarts[p][e][0])); // all runs start at 0
      runStartFreqs[p][e][0] = Nhaps;
      memset(&runStartFreqs[p][e][1], 0, Mseg64*sizeof(runStarts[p][e][1]));
    }

    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      // initialize cur [p=(m64&1)] arrays
      const uint p = m64&1;
      for (uint e = 0; e < 2; e++) {
	memcpy(runStarts[p][e], runStarts[!p][e], Nhaps*sizeof(runStarts[p][e][0]));
	memset(runStartFreqs[p][e], 0, (Mseg64+1)*sizeof(runStartFreqs[p][e][0]));
      }

      uint64_masks bits0 = genoBits[m64*N + n0];
      //uint *runStarts_p_0 = runStarts[p][0], *runStarts_p_1 = runStarts[p][1];
      //uint *runStartFreqs_p_0 = runStartFreqs[p][0], *runStartFreqs_p_1 = runStartFreqs[p][1];
      for (uint64 n1 = 0; n1 < Nhaps; n1++) {
	uint64 is1 = haploBits[m64*2*N + n1];
	uint64 wrongBits = (bits0.is0 & is1) | (bits0.is2 & ~is1);

	// update cur runStarts based on bit matches
	if (wrongBits) {
	  if (wrongBits & (wrongBits-1)) // 2+ wrong => fail; move both starts forward
	    runStarts[p][0][n1] = runStarts[p][1][n1] = m64+1;
	  else { // 1 wrong => move 1-err start to 0-err start; fail 0-err
	    runStarts[p][1][n1] = runStarts[p][0][n1];
	    runStarts[p][0][n1] = m64+1;
	  }
	}
	// compute cur runStartFreqs based on (updated) cur runStarts
	for (uint e = errStart; e < 2; e++)
	  runStartFreqs[p][e][runStarts[p][e][n1]]++;
	/*
	if (wrongBits) {
	  if (wrongBits & (wrongBits-1)) // 2+ wrong => fail; move both starts forward
	    runStarts_p_0[n1] = runStarts_p_1[n1] = m64+1;
	  else { // 1 wrong => move 1-err start to 0-err start; fail 0-err
	    runStarts_p_1[n1] = runStarts_p_0[n1];
	    runStarts_p_0[n1] = m64+1;
	  }
	}
	runStartFreqs_p_0[runStarts_p_0[n1]]++;
	runStartFreqs_p_1[runStarts_p_1[n1]]++;
	*/
      }

      for (uint e = errStart; e < 2; e++)
	for (uint s64 = 0; s64 <= m64; s64++)
	  if ((runStartFreqs[p][e][s64] >= K && m64+1 == Mseg64) || // reached end
	      (runStartFreqs[p][e][s64] < K && runStartFreqs[!p][e][s64] >= K)) {
	    // all n1 s.t. runStarts[p][e][n1] == s64 are in top K starting from s64
	    // some n1 s.t. runStarts[!p][e][n1] == s64 are in top K starting from s64
	    uint allPos = 0, somePos = runStartFreqs[p][e][s64];
	    uint *topIndsChunk = &topInds[e][s64*K];
	    for (uint64 n1 = 0; n1 < Nhaps; n1++) {
	      if (allPos < K && runStarts[p][e][n1] == s64)
		topIndsChunk[allPos++] = n1;
	      else if (somePos < K && runStarts[!p][e][n1] == s64)
		topIndsChunk[somePos++] = n1;
	    }
	    topIndsLens[e][s64] = min(somePos, K);
	  }
    }

    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	ALIGNED_FREE(runStartFreqs[p][e]);
    for (uint p = 0; p < 2; p++)
      for (uint e = 0; e < 2; e++)
	ALIGNED_FREE(runStarts[p][e]);
  }

  void Eagle::randomlyPhaseTmpHaploBitsT(uint64 n0) {
    // fast rng: last 16 bits of Marsaglia's MWC
    uint w = 521288629;
    for (uint i = 0; i < (n0 & 0xff); i++)
      w=18000*(w&65535)+(w>>16);

    uint64 n1 = (n0+1) % N;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
      if (maskSnps64j[m64j]) {
	uint g = getGeno0123(m64j, n0);
	if (g == 3) g = getGeno0123(m64j, n1); // if missing, try filling with geno of next sample
	if (g == 3) g = 1; // if still missing, set to het

	if (g == 0) // nothing to do; tmpHaploBitsT is already cleared in init()
	  ;
	else if (g == 2) { // set bit in both parental haplotypes
	  for (uint64 q = 0; q <= 1ULL; q++)
	    tmpHaploBitsT[(2*n0+q)*Mseg64 + (m64j/64)] |= 1ULL<<(m64j&63);
	}
	else { // set bit in random parental haplotype
	  uint64 q = (w=18000*(w&65535)+(w>>16))&1;
	  tmpHaploBitsT[(2*n0+q)*Mseg64 + (m64j/64)] |= 1ULL<<(m64j&63);
	}
      }
  }

  pair <double, vector <double> > Eagle::findLongDipMatches(uint64 n0, uint64 nF1, uint64 nF2) {

    if (!maskIndivs[n0]) return make_pair(0.0, vector <double> ());

    vector <uint> topInds[2]; // [max err]; lengths are Mseg64 * K
    vector <uint> topIndsLens[2]; // [max err]; lengths are Mseg64

    const uint K = 10;
    Timer timer;
    findLongHalfIBD(n0, topInds, topIndsLens, K);
    double halfIBDtime = timer.update_time();

    const double duoMatchThresh = log10(N*10), longMatchMin = 4.0;


    /***** COMPUTE BAYES FACTORS FOR LONG HALF-IBD REGIONS TO SELECT LONG MATCHES *****/

    vector <Match> cumTopMatches;
    vector <Match> longMatches;
    double *workLogBF = ALIGNED_MALLOC_DOUBLES(Mseg64*64);
    double (*memoLogBF)[4] = (double (*)[4]) ALIGNED_MALLOC(Mseg64*64*4*sizeof(double));
    // initialize memo lookup table of logBFs (n1 given n0)
    for (uint m64j = 0; m64j < Mseg64*64; m64j++)
      for (uint g = 0; g < 4; g++)
	memoLogBF[m64j][g] = MEMO_UNSET;

    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      std::set <uint> n1s; n1s.insert(n0);
      // delete previous matches that have ended: sort by greaterEnd and erase range at end
      sort(cumTopMatches.begin(), cumTopMatches.end(), Match::greaterEnd);
      for (uint k = 0; k < cumTopMatches.size(); k++) {
	if (m64 >= cumTopMatches[k].m64jEnd/64) {
	  cumTopMatches.resize(k);
	  break;
	}
	else
	  n1s.insert(cumTopMatches[k].n);
      }

      // incorporate new matches
      for (uint e = 0; e < 2; e++) {
	for (uint k = 0; k < topIndsLens[e][m64]; k++) {
	  uint n1 = topInds[e][m64*K + k];
	  if (maskIndivs[n1] && !n1s.count(n1)) {
	    n1s.insert(n1);
	    Match duoMatch = computeDuoLogBF(memoLogBF, workLogBF, n0, n1, m64);
	    assert(maskSnps64j[duoMatch.m64jStart]);
	    assert(maskSnps64j[duoMatch.m64jEnd]);

	    if (duoMatch.logBF > duoMatchThresh && duoMatch.cMlenInit >= longMatchMin) {
	      longMatches.push_back(duoMatch);
#ifdef VERBOSE
	      printMatch(n0, nF1, nF2, duoMatch, memoLogBF);
#endif
	      //retractMatch(n0, duoMatch, memoLogBF);
	      cumTopMatches.push_back(duoMatch);
	    }
	  }
	}
      }
    }
#ifdef VERBOSE
    cout << "num longMatches: " << longMatches.size() << endl;
#endif

    /***** TRIM MATCHES UNTIL CONSISTENT; DETERMINE SAME/OPP ORIENTATION OF PAIRS *****/

    const double longMatchMinTrim = 3.0, minSameOppDiff = 0.5, minMatchLenDiff = 0.5;

    sort(longMatches.begin(), longMatches.end(), Match::greaterLen);
    vector < vector <uint> > sameEdges(longMatches.size()), oppEdges(longMatches.size());
    for (uint i = 0; i < longMatches.size(); i++)
      for (uint j = 0; j < i; j++) { // trim against previous (longer) matches
	double iLen = cMs64j[longMatches[i].m64jEnd] - cMs64j[longMatches[i].m64jStart];
	double jLen = cMs64j[longMatches[j].m64jEnd] - cMs64j[longMatches[j].m64jStart];
	// check if too short (after previous trimming)
	if (iLen < longMatchMinTrim || jLen < longMatchMinTrim) continue;

	Match &d1 = longMatches[iLen<jLen?i:j]; // shorter
	Match &d2 = longMatches[iLen<jLen?j:i]; // longer

	if (d1.n == nF1 || d1.n == nF2 || d2.n == nF1 || d2.n == nF2) continue;
	if (d1.m64jEnd < d2.m64jStart || d2.m64jEnd < d1.m64jStart)
	  continue; // no overlap
#ifdef VERBOSE
	printf("  cM= %.1f (%.1f-%.1f), cM= %.1f (%.1f-%.1f):",
	       cMs64j[d1.m64jEnd]-cMs64j[d1.m64jStart], cMs64j[d1.m64jStart], cMs64j[d1.m64jEnd],
	       cMs64j[d2.m64jEnd]-cMs64j[d2.m64jStart], cMs64j[d2.m64jStart], cMs64j[d2.m64jEnd]);
#endif
	uint64 maxStart = max(d1.m64jStart, d2.m64jStart);
	uint64 minEnd = min(d1.m64jEnd, d2.m64jEnd);
	uint64 lastSame = maxStart, lastOpp = lastSame;
	double longestSame = 0, longestOpp = 0;
	uint64 sameStart = 0, sameEnd = 0, oppStart = 0, oppEnd = 0;
	for (uint64 m64j = maxStart; m64j <= minEnd; m64j++) {
	  bool same = false, opp = false;
	  if (m64j == minEnd) {
	    same = true; opp = true;
	  }
	  else if (getGeno0123(m64j, n0) == 1) {
	    uint g1 = getGeno0123(m64j, d1.n);
	    uint g2 = getGeno0123(m64j, d2.n);
	    if ((g1==0||g1==2) && (g2==0||g2==2)) {
	      same = g1==g2;
	      opp = g1!=g2;
	    }
	  }
	  if (same) {
	    double oppLen = cMs64j[m64j] - cMs64j[lastSame];
	    if (longestOpp < oppLen) {
	      longestOpp = oppLen;
	      oppStart = lastSame; oppEnd = m64j;
	    }
	    lastSame = m64j;
	  }
	  if (opp) {
	    double sameLen = cMs64j[m64j] - cMs64j[lastOpp];
	    if (longestSame < sameLen) {
	      longestSame = sameLen;
	      sameStart = lastOpp; sameEnd = m64j;
	    }
	    lastOpp = m64j;
	  }
	}
	double fullOverlap = cMs64j[minEnd] - cMs64j[maxStart];
#ifdef VERBOSE
	printf(" %.1f|%.1f (%.1f-%.1f)|(%.1f-%.1f)\n", longestSame, longestOpp,
	       cMs64j[sameStart], cMs64j[sameEnd], cMs64j[oppStart], cMs64j[oppEnd]);
#endif

	uint64 trimStart = 0; int orientation = 0;
	if (longestSame > longestOpp + minSameOppDiff
	    || (longestSame == fullOverlap && longestSame > longestOpp)) {
	  sameEdges[i].push_back(j);
	  sameEdges[j].push_back(i);
	  orientation = 1;
	  trimStart = (sameStart+sameEnd)/2;
	}
	else if (longestOpp > longestSame + minSameOppDiff
		 || (longestOpp == fullOverlap && longestOpp > longestSame)) {
	  oppEdges[i].push_back(j);
	  oppEdges[j].push_back(i);
	  orientation = -1;
	  trimStart = (oppStart+oppEnd)/2;
	}

	if (orientation != 0) { // clearly same or opp
	  // left end
	  if (cMs64j[d1.m64jStart] < cMs64j[d2.m64jStart] - minMatchLenDiff) // d1 extends further
	    trim(d2, d1, n0, orientation, trimStart, -1, workLogBF); // => trim d2
	  else
	    trim(d1, d2, n0, orientation, trimStart, -1, workLogBF);
	  // right end
	  if (cMs64j[d1.m64jEnd] > cMs64j[d2.m64jEnd] + minMatchLenDiff) // d1 extends further
	    trim(d2, d1, n0, orientation, trimStart, 1, workLogBF); // => trim d2
	  else
	    trim(d1, d2, n0, orientation, trimStart, 1, workLogBF);
#ifdef VERBOSE
	  printf("  --> %.1f (%.1f-%.1f)  --> %.1f (%.1f-%.1f)\n",
		 cMs64j[d1.m64jEnd]-cMs64j[d1.m64jStart], cMs64j[d1.m64jStart], cMs64j[d1.m64jEnd],
		 cMs64j[d2.m64jEnd]-cMs64j[d2.m64jStart], cMs64j[d2.m64jStart], cMs64j[d2.m64jEnd]);
#endif
	}
	else { // can't determine same vs. opp => chop shorter so that matches no longer overlap
	  if (d1.m64jStart < d2.m64jStart) d1.m64jEnd = d2.m64jStart-1;
	  else d1.m64jStart = min(d2.m64jEnd+1ULL, Mseg64*64-1);
#ifdef VERBOSE
	  printf("  trimmed first to cM= %.1f (%.1f-%.1f)\n",
		 cMs64j[d1.m64jEnd]-cMs64j[d1.m64jStart], cMs64j[d1.m64jStart], cMs64j[d1.m64jEnd]);
#endif
	}
      }


    /***** REDUCE TO SET OF TRIMMED MATCHES WITH CONSISTENT SIGNS *****/

    vector <bool> kept(longMatches.size()), checked(longMatches.size());
    for (uint t = 0; t < longMatches.size(); t++) {
      // find longest remaining trimmed match
      uint i = 0; double longest = 0;
      for (uint iTest = 0; iTest < longMatches.size(); iTest++)
	if (!checked[iTest]) {
	  const Match &d1 = longMatches[iTest];
	  if (cMs64j[d1.m64jEnd]-cMs64j[d1.m64jStart] > longest) {
	    longest = cMs64j[d1.m64jEnd]-cMs64j[d1.m64jStart];
	    i = iTest;
	  }
	}
      checked[i] = true;
      if (longest < longMatchMinTrim) break;
      if (longMatches[i].n == nF1 || longMatches[i].n == nF2) continue;

      kept[i] = true;
#ifdef VERBOSE
      const Match &d1 = longMatches[i];
      printf("cM= %.1f (%.1f-%.1f)\n",
	     cMs64j[d1.m64jEnd]-cMs64j[d1.m64jStart], cMs64j[d1.m64jStart], cMs64j[d1.m64jEnd]);
#endif
      vector <int> signs = searchSigns(longMatches, sameEdges, oppEdges, kept);
      if (signs.empty()) { // inconsistent signs
	kept[i] = false;
#ifdef VERBOSE
	cout << "  WARNING: sign inconsistency: eliminating" << endl;
#endif
      }
    }

    // compute final signs
    vector <int> signs = searchSigns(longMatches, sameEdges, oppEdges, kept);

#ifdef VERBOSE
    for (uint i = 0; i < longMatches.size(); i++)
      if (kept[i]) {
	const Match &d1 = longMatches[i];
	printf("cM= %.1f (%.1f-%.1f): ",
	vector <Match> curMatch(1, d1); vector <int> curSign(1, signs[i]);
	       cMs64j[d1.m64jEnd]-cMs64j[d1.m64jStart], cMs64j[d1.m64jStart], cMs64j[d1.m64jEnd]);
	cout << computePhaseString(n0, nF1, nF2, curMatch, curSign, d1.m64jStart,
				   cMs64j[d1.m64jEnd]+1e-9, false) << endl;
	cout << computePhaseString(n0, nF1, nF2, curMatch, curSign, d1.m64jStart,
				   cMs64j[d1.m64jEnd]+1e-9, true) << endl;
      }
    cout << endl << endl
	 << "phase: " << computePhaseString(n0, nF1, nF2, longMatches, signs, 0, 1e100, false) << endl;
    cout << endl << endl
	 << "phase: " << computePhaseString(n0, nF1, nF2, longMatches, signs, 0, 1e100, true) << endl;
#endif
    computePhaseConfs(n0, longMatches, signs, true);
    if ((int) nF1 != -1 && (int) nF2 != -1) {
      checkPhase(n0, nF1, nF2, 0.1);
      checkPhase(n0, nF1, nF2, 0.5);
    }

    ALIGNED_FREE(workLogBF);
    ALIGNED_FREE(memoLogBF);

    // record longest match length per seg64 for output
    vector <double> cMs(Mseg64);
    bool cons = true;
    for (uint i = 0; i < longMatches.size(); i++) {
      if (!signs[i]) continue;
      uint64 start, end;
      if (cons) {
	start = std::max(longMatches[i].m64jStartCons, longMatches[i].m64jStart);
	end = std::min(longMatches[i].m64jEndCons, longMatches[i].m64jEnd);
      }
      else {
	start = longMatches[i].m64jStart;
	end = longMatches[i].m64jEnd;
      }
      for (uint64 m64 = (start+63)/64; m64 < end/64; m64++)
	cMs[m64] = std::max(cMs[m64], cMs64j[end]-cMs64j[start]);
    }

    return make_pair(halfIBDtime, cMs);
  }

  int Eagle::numDipHapWrongBits(uint64 m64, uint64 n0, uint64 n1hap) const {
    uint64 is1 = haploBits[m64*2*N + n1hap];
    const uint64_masks &bits0 = genoBits[m64*N + n0];
    uint64 wrongBits = (bits0.is0 & is1) | (bits0.is2 & ~is1);
    return popcount64(wrongBits);
  }

  int Eagle::firstDipHapGoodBit(uint64 m64, uint64 n0, uint64 n1hap) const {
    uint64 is1 = haploBits[m64*2*N + n1hap];
    const uint64_masks &bits0 = genoBits[m64*N + n0];
    uint64 wrongBits = (bits0.is0 & is1) | (bits0.is2 & ~is1);
    return wrongBits ? 64 - __builtin_clzll(wrongBits) : 0; // MSB
  }

  int Eagle::firstDipHapWrongBit(uint64 m64, uint64 n0, uint64 n1hap) const {
    if (m64 >= Mseg64) return 0;
    uint64 is1 = haploBits[m64*2*N + n1hap];
    const uint64_masks &bits0 = genoBits[m64*N + n0];
    uint64 wrongBits = (bits0.is0 & is1) | (bits0.is2 & ~is1);
    return popcount64((wrongBits & ~(wrongBits-1))-1); // LSB
  }

    struct DipHapSeg {
      uint n, start, end;
      DipHapSeg(uint _n, uint _start, uint _end) : n(_n), start(_start), end(_end) {}
      bool operator < (const DipHapSeg &seg2) const {
	return end-start > seg2.end-seg2.start
	  || (end-start == seg2.end-seg2.start && n < seg2.n);
      }
    };
    struct DipHapSegFarther {
      bool operator() (const DipHapSeg &seg1, const DipHapSeg &seg2) const {
	return seg1.end > seg2.end || (seg1.end == seg2.end && seg1.n < seg2.n);
      }
    };

  int Eagle::countSE(const vector <bool> &phaseVec) {
    int ans = 0;
    for (uint h = 1; h < phaseVec.size(); h++)
      ans += (phaseVec[h] != phaseVec[h-1]);
    return ans;
  }

  int Eagle::countMajorSE(const vector <bool> &phaseVec) {
    vector <bool> phaseVec7;
    for (uint h7 = 0; h7+7 < phaseVec.size(); h7 += 7) {
      int votes = 0;
      for (uint h = h7; h < h7+7; h++)
	votes += phaseVec[h];
      phaseVec7.push_back(votes >= 4);
    }
    return countSE(phaseVec7);
  }

  double Eagle::findLongHapMatches(uint64 n0, uint64 nF1, uint64 nF2, int iter) {

    if (!maskIndivs[n0]) return 0;

    if (Mseg64 < 3U) {
      cerr << "Too few SNP segments for analysis" << endl;
      exit(1);
    }

    uint seed = n0; // for rand_r()

    vector <uint> topInds[2]; // [max err]; lengths are Mseg64 * K
    vector <uint> topIndsLens[2]; // [max err]; lengths are Mseg64

    const uint K = 20;
    Timer timer;
    findLongDipHap(n0, topInds, topIndsLens, K);
    double hapTime = timer.update_time();
    /*
    // VISUALIZE RESULTS
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      cout << "m64 = " << m64 << endl;
      int m64xMin = std::max(0, (int) m64-5);
      int m64xMax = std::min(m64+11, Mseg64);
      for (int m64x = m64xMin; m64x < m64xMax; m64x++) {
	int numSnps = 0;
	for (int j = 0; j < 64; j++)
	  numSnps += maskSnps64j[m64x*64+j];
	printf("%2d ", numSnps);
      }
      cout << endl;
      for (uint e = 0; e < 2; e++) {
	cout << "e = " << e << endl;
	for (uint k = 0; k < topIndsLens[e][m64]; k++) {
	  uint64 n1hap = topInds[e][m64*K + k];
	  for (int m64x = m64xMin; m64x < m64xMax; m64x++) {
	    printf("%2d ", numDipHapWrongBits(m64x, n0, n1hap));
	  }
	  cout << "   (" << n1hap << ")" << endl;
	}
      }
    }
    */
    const int maxWrongBits = 3;
    std::set <DipHapSeg> longDipHapSegs[Mseg64];
    const uint maxSegs = (iter == 2 ? 10 : 20); // TODO: increase accuracy with 20?
    vector <uint> lastEnd64(2*N);
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      //if ((int) nF1 != -1) cout << "m64 = " << m64 << endl;
      for (uint e = 0; e < 2; e++) {
	//if ((int) nF1 != -1) cout << "e = " << e << endl;
	for (uint k = 0; k < topIndsLens[e][m64]; k++) {
	  uint64 n1hap = topInds[e][m64*K + k];
	  if (!maskIndivs[n1hap/2]) continue;
	  if (n1hap/2 == n0) continue;
	  if (lastEnd64[n1hap] > m64) continue;

	  // find start
	  uint segStart = m64;
	  while ((int) segStart >= 0 && numDipHapWrongBits(segStart, n0, n1hap) <= maxWrongBits)
	    segStart--;
	  segStart++;
	  // find end
	  uint segEnd = m64;
	  while (segEnd < Mseg64 && numDipHapWrongBits(segEnd, n0, n1hap) <= maxWrongBits)
	    segEnd++;

	  //checkHapPhase1(n0, nF1, nF2, n1hap, segStart, segEnd);

	  lastEnd64[n1hap] = segEnd;
	  for (uint64 m64x = segStart; m64x < segEnd; m64x++) {
	    longDipHapSegs[m64x].insert(DipHapSeg(n1hap, segStart, segEnd));
	    if (longDipHapSegs[m64x].size() > maxSegs)
	      longDipHapSegs[m64x].erase(--longDipHapSegs[m64x].end());
	  }
	}
      }
    }
    /*
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      cout << "m64 = " << m64 << endl;
      for (std::set <DipHapSeg>::iterator it = longDipHapSegs[m64].begin();
	   it != longDipHapSegs[m64].end(); it++) {
	printf("%2d: %2d-%2d   (%d)\n", it->end-it->start, it->start, it->end, it->n);
      }
    }
    */
    uint64 curHaploBitsT[Mseg64];
    const uint64 side = 1;
    //const uint numHapHaps = 5;
    vector <uint> errBests(Mseg64);
    vector < pair <uint, uint > > n12hapBests(Mseg64);
    for (uint64 m64 = 0+side; m64+side < Mseg64; m64++) {
      /*
      if ((int) nF1 != -1) {
	cout << "m64 = " << m64 << endl;
	for (uint i = 0; i < numHapHaps; i++) printf("   ");
      }
      checkHapPhase(n0, nF1, nF2, haploBitsT + 2*n0*Mseg64, m64, side);
      */
      vector < pair <uint, uint> > bestHitPairs; uint minWrongBits = 99;
      for (std::set <DipHapSeg>::iterator it = longDipHapSegs[m64].begin();
	   it != longDipHapSegs[m64].end(); it++) { // for each dip-hap long match
	uint64 n1hap = it->n;
	// ---------- SLOW ----------
	/*
	std::set < pair <uint, uint> > wrongBitsHaps; uint worstInSet = 1<<30;
	for (uint64 n2hap = 0; n2hap < 2*N; n2hap++) {
	  if (!maskIndivs[n2hap/2]) continue;
	  if (n2hap/2 == n0) continue;
	  if (n2hap/2 == n1hap/2) continue; // opp haps of same indiv
	  uint numWrongBits = 0;
	  for (uint64 m64x = m64-side; m64x <= m64+side; m64x++) {
	    uint64 n1is1 = haploBits[m64x*2*N + n1hap];
	    uint64 n2is1 = haploBits[m64x*2*N + n2hap];
	    const uint64_masks &bits0 = genoBits[m64x*N + n0];
	    uint64 wrongBits = (bits0.is0 & (n1is1 | n2is1)) | (bits0.is2 & ~(n1is1 & n2is1))
	      | (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1));
	    numWrongBits += popcount64(wrongBits);
	  }
	  if (numWrongBits < worstInSet) {
	    wrongBitsHaps.insert(make_pair(numWrongBits, n2hap));
	    if (wrongBitsHaps.size() > numHapHaps) {
	      wrongBitsHaps.erase(--wrongBitsHaps.end());
	      worstInSet = (--wrongBitsHaps.end())->first;
	    }
	  }
	}
	for (std::set < pair <uint, uint> >::iterator itH = wrongBitsHaps.begin();
	     itH != wrongBitsHaps.end(); itH++)
	  printf("%2d ", itH->first);
	  //printf("%2d   (%d)\n", itH->first, itH->second);
	*/
	// ---------- FAST ----------
	if ((int) nF1 != -1)
	  printf("| ");
	for (uint64 m64x = m64-side; m64x <= m64+side; m64x++) {
	  curHaploBitsT[m64x] = 0;
	  for (uint64 j = 0; j < 64; j++) {
	    uint64 m64j = m64x*64+j, bit = 0;
	    if (maskSnps64j[m64j]) {
	      uint g0 = getGeno0123(m64j, n0);
	      if (g0 == 0 || g0 == 3) bit = 0;
	      else if (g0 == 2) bit = 1;
	      else bit = 1-((haploBitsT[n1hap*Mseg64 + m64x]>>j)&1);
	    }
	    curHaploBitsT[m64x] |= bit<<j;
	  }
	}
	std::ostringstream oss;
	for (uint h = 0; h < hashLookups[m64].size(); h++) { // for each hashing
	  uint numHits;
	  const uint *lenHapInds =
	    hashLookups[m64][h].query(computeHash(curHaploBitsT, hashBits[m64][h]));
	  if (lenHapInds == NULL)
	    numHits = 0;
	  else
	    numHits = lenHapInds[0];
	  //uint best = 99;
	  for (uint k = 1; k <= numHits; k++) {
	    uint64 n2hap = lenHapInds[k];
	    if (!maskIndivs[n2hap/2]) continue;
	    if (n2hap/2 == n0) continue;
	    if (n2hap/2 == n1hap/2) continue; // opp haps of same indiv
	    uint numWrongBits = 0;
	    for (uint64 m64x = m64-side; m64x <= m64+side; m64x++) {
	      uint64 n1is1 = haploBitsT[n1hap*Mseg64 + m64x];
	      uint64 n2is1 = haploBitsT[n2hap*Mseg64 + m64x];
	      const uint64_masks &bits0 = genoBits[m64x*N + n0];
	      uint64 wrongBits = (bits0.is0 & (n1is1 | n2is1)) | (bits0.is2 & ~(n1is1 & n2is1))
		| (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1));
	      uint popcnt = popcount64(wrongBits);
	      numWrongBits += popcnt;
	    }
	    if (numWrongBits < minWrongBits) {
	      minWrongBits = numWrongBits;
	      bestHitPairs.clear();
	      bestHitPairs.push_back(make_pair((uint) n1hap, (uint) n2hap));
	    }
	    else if (numWrongBits == minWrongBits)
	      bestHitPairs.push_back(make_pair((uint) n1hap, (uint) n2hap));
	    /*
	    if (numWrongBits < best)
	      best = numWrongBits;
	    */
	  }
	  /*
	  if ((int) nF1 != -1) {
	    printf("%2d ", best);
	    char buf[10]; sprintf(buf, "%2d ", numHits);
	    oss << string(buf);
	  }
	  */
	}
	if (/*(int) nF1 != -1*/false) {
	  cout << ": " << oss.str();
	  checkHapPhase(n0, nF1, nF2, haploBitsT + /*wrongBitsHaps.begin()->second*/n1hap*Mseg64, m64, side); // n1hap is better!
	}
      }
      if (bestHitPairs.empty()) { // keep current phasing
	errBests[m64] = 99;
	n12hapBests[m64] = make_pair((uint) n0*2, (uint) n0*2+1);
      }
      else {
	errBests[m64] = minWrongBits;
	n12hapBests[m64] = bestHitPairs[rand_r(&seed) % bestHitPairs.size()];
      }
    }

    if ((int) nF1 != -1)
      cout << endl << "2nd-iter phase:" << endl;
    int sign = 1; vector <int> signs(Mseg64);
    uint64 prevInd = side;
    vector <bool> phaseVec;
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      uint64 curInd = min(max(m64, side), Mseg64-1-side); /*int minWrongBits = 99;
      for (int diff = - (int) side; diff <= (int) side; diff++) {
	uint64 m64d = m64+diff;
	if (m64d >= Mseg64) continue;
	uint64 n1is1 = haploBitsT[n12hapBests[m64d].first*Mseg64 + m64];
	uint64 n2is1 = haploBitsT[n12hapBests[m64d].second*Mseg64 + m64];
	const uint64_masks &bits0 = genoBits[m64*N + n0];
	uint64 wrongBits = (bits0.is0 & (n1is1 | n2is1)) | (bits0.is2 & ~(n1is1 & n2is1))
	  | (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1));
	int numWrongBits = popcount64(wrongBits);
	if (numWrongBits < minWrongBits) {
	  minWrongBits = numWrongBits;
	  curInd = m64d;
	}
      }
							  */ // best of 3 neighbors doesn't help?
      vector < pair <int, int> > offsetMults;
      if (m64 > side && m64 < Mseg64-side) { // update sign
	uint64 n1prev = n12hapBests[prevInd].first, n2prev = n12hapBests[prevInd].second;
	uint64 n1cur = n12hapBests[curInd].first, n2cur = n12hapBests[curInd].second;
	for (uint64 m64j = (m64-side)*64; m64j < (m64+side)*64; m64j++) {
	  uint h1prev = (haploBitsT[n1prev*Mseg64 + (m64j/64)]>>(m64j&63))&1;
	  uint h2prev = (haploBitsT[n2prev*Mseg64 + (m64j/64)]>>(m64j&63))&1;
	  uint h1cur = (haploBitsT[n1cur*Mseg64 + (m64j/64)]>>(m64j&63))&1;
	  uint h2cur = (haploBitsT[n2cur*Mseg64 + (m64j/64)]>>(m64j&63))&1;
	  if (h1prev + h2prev == 1 && h1cur + h2cur == 1) {
	    int offset = abs((int) (m64j - m64*64));
	    offsetMults.push_back(make_pair(offset, h1prev == h1cur ? 1 : -1));
	  }
	}
      }
      if (!offsetMults.empty()) {
	sort(offsetMults.begin(), offsetMults.end());
	int totVotes = min((int) offsetMults.size(), 5), sameVotes = 0;
	for (int k = 0; k < totVotes; k++)
	  if (offsetMults[k].second == 1)
	    sameVotes++;
	if (sameVotes < (totVotes+1)/2) sign *= -1; // swap sign
      }
      signs[m64] = sign;
      uint64 n1hap = n12hapBests[curInd].first, n2hap = n12hapBests[curInd].second;
      computeSegPhaseConfs(n0, n1hap, n2hap, sign, m64, errBests[curInd]);
      prevInd = curInd;
    }

    if ((int) nF1 != -1) {
      for (uint64 m64 = 0; m64 < Mseg64; m64++) {
	uint64 curInd = min(max(m64, side), Mseg64-1-side);
	char usedInd = ' ';
	if (curInd > m64) usedInd = 'v';
	else if (curInd < m64) usedInd = '^';
	uint64 n1hap = n12hapBests[curInd].first, n2hap = n12hapBests[curInd].second;
	int sign = signs[m64];
	printf("m64 = %2d   err: %2d %c (%6d,%6d)   ",
	       (int) m64, (int) errBests[curInd], usedInd, (int) (sign==1 ? n1hap : n2hap),
	       (int) (sign==1 ? n2hap : n1hap));
	//cout << closestOffset << endl;
	vector <bool> phaseSeg = checkSegPhase(n0, nF1, nF2, n1hap, n2hap, sign, m64);
	phaseVec.insert(phaseVec.end(), phaseSeg.begin(), phaseSeg.end());
      }

      printf("# major SE: %2d   # tot SE: %2d / %d\n", countMajorSE(phaseVec), countSE(phaseVec),
	     (int) phaseVec.size()-1);
      fflush(stdout);
    }
    /*
    for (uint64 n1hap = 2*n0; n1hap <= 2*n0+1; n1hap++) {
      cout << "n0 = " << n0 << "; n1hap = " << n1hap << ": ";
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
	cout << (int) phaseConfs2[n1hap*Mseg64*64 + m64j] << " ";
      cout << endl;
    }
    */
    return hapTime;
  }

  inline uint64 pairToULL(pair <uint, uint> p) { return ((uint64) p.first<<32ULL)|p.second; }
  inline pair <uint, uint> ullToPair(uint64 ull) {
    return make_pair((uint) (ull>>32ULL), (uint) ull);
  }

  bool Eagle::updateHelper(std::unordered_map <uint64, DPState> &dpTab, uint &dpBestScore,
			   pair <uint, uint> cur, pair <uint, uint> next, uint score) const {
#ifdef RDTSC_TIMING
    uint64 tscStart = Timer::rdtsc();
#endif
    /* ---- SLOW ----
    if (dpTab.find(pairToULL(next)) == dpTab.end() || dpTab[pairToULL(next)].score > score)
      dpTab[pairToULL(next)] = DPState(score, cur);
    */
    if (score > dpBestScore + 2*switchCost) return false;
    DPState &nextState = dpTab[pairToULL(next)];
    if (nextState.score == 0 || nextState.score > score) {
      nextState.score = score;
      nextState.from = cur;
    }
    if (score < dpBestScore) dpBestScore = score;
#ifdef RDTSC_TIMING
    dpUpdateTicks += Timer::rdtsc() - tscStart;
    dpUpdateCalls++;
#endif
    return true;
  }

  uint Eagle::computeStaticScore(uint n0, uint n1hap, uint n2hap, uint64 m64) const {
#ifdef RDTSC_TIMING
    uint64 tscStart = Timer::rdtsc();
#endif
    uint64 n1is1 = haploBitsT[n1hap*Mseg64 + m64];
    uint64 n2is1 = haploBitsT[n2hap*Mseg64 + m64];
    const uint64_masks &bits0 = genoBits[m64*N + n0];
    uint64 wrongHomBits = (bits0.is0 & (n1is1 | n2is1)) | (bits0.is2 & ~(n1is1 & n2is1));
    uint64 wrongHetBits = (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1));
    uint score = popcount64(wrongHomBits)*homErrCost
      + popcount64(wrongHetBits)*hetErrCost
      + segConfs[n1hap*Mseg64+m64] / 5 + segConfs[n2hap*Mseg64+m64] / 5;
#ifdef RDTSC_TIMING
    dpStaticTicks += Timer::rdtsc() - tscStart;
#endif
    return score;
  }

  uint Eagle::computeSwitchScore(uint n0, uint n1hap, uint n2hapA, uint n2hapB, uint64 m64) const {
#ifdef RDTSC_TIMING
    uint64 tscStart = Timer::rdtsc();
#endif
    uint64 n1is1 = haploBitsT[n1hap*Mseg64 + m64];
    uint64 n2is1A = haploBitsT[n2hapA*Mseg64 + m64];
    uint64 n2is1B = haploBitsT[n2hapB*Mseg64 + m64];
    const uint64_masks &bits0 = genoBits[m64*N + n0];

    uint64 wrongHomBitsA = (bits0.is0 & (n1is1 | n2is1A)) | (bits0.is2 & ~(n1is1 & n2is1A));
    uint64 wrongHetBitsA = (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1A));
    uint64 wrongHomBitsB = (bits0.is0 & (n1is1 | n2is1B)) | (bits0.is2 & ~(n1is1 & n2is1B));
    uint64 wrongHetBitsB = (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1B));
    /* ---- SLOW ----
    uint score = popcount64(wrongHomBitsB)*homErrCost
      + popcount64(wrongHetBitsB)*hetErrCost;
    uint minScore = score;
    for (uint64 k = 0; k < 64; k++) {
      score += ((wrongHomBitsA>>k)&1)*homErrCost + ((wrongHetBitsA>>k)&1)*hetErrCost
	- ((wrongHomBitsB>>k)&1)*homErrCost - ((wrongHetBitsB>>k)&1)*hetErrCost;
      if (score < minScore) minScore = score;
    }
    assert(score == popcount64(wrongHomBitsA)*homErrCost
	   + popcount64(wrongHetBitsA)*hetErrCost);
    */
    uint64 wrongBitsA = wrongHomBitsA | wrongHetBitsA;
    uint64 wrongBitsB = wrongHomBitsB | wrongHetBitsB;
    uint64 hetBits = ~(bits0.is0|bits0.is2|bits0.is9);
    uint mask = (1U<<switchScoreLutBits)-1;
    int curScore = popcount64(wrongHomBitsB)*homErrCost
      + popcount64(wrongHetBitsB)*hetErrCost;
    int bestScore = curScore;
    for (uint64 b = 0; b < 64; b += switchScoreLutBits) {
      uint lutInd = ((((uint) (wrongBitsA>>b))&mask)<<(switchScoreLutBits+switchScoreLutBits))
	| ((((uint) (wrongBitsB>>b))&mask)<<switchScoreLutBits)
	| ((((uint) (hetBits>>b))&mask));
      bestScore = min(bestScore, curScore + switchScoreLut[lutInd][0]);
      curScore += switchScoreLut[lutInd][1];
    }
    curScore = bestScore + segConfs[n1hap*Mseg64+m64] / 5 + segConfs[n2hapA*Mseg64+m64] / 5
      + segConfs[n2hapB*Mseg64+m64] / 5;

    // additional penalty for errors in (n1hap,n2hapB) at end of m64
    curScore += ((wrongBitsB & 0xf000000000000000) != 0) + ((wrongBitsB & 0xff00000000000000) != 0)
      + ((wrongBitsB & 0xffff000000000000) != 0);

    if (m64+1 < Mseg64) { // additional penalty for errors in (n1hap,n2hapB) at start of m64+1
      n1is1 = haploBitsT[n1hap*Mseg64 + m64+1];
      n2is1B = haploBitsT[n2hapB*Mseg64 + m64+1];
      const uint64_masks &bits0next = genoBits[(m64+1)*N + n0];
      wrongBitsB = (bits0next.is0 & (n1is1 | n2is1B)) | (bits0next.is2 & ~(n1is1 & n2is1B))
	| (~(bits0next.is0|bits0next.is2|bits0next.is9) & ~(n1is1 ^ n2is1B));
      curScore += ((wrongBitsB & 0xf) != 0) + ((wrongBitsB & 0xff) != 0)
	+ ((wrongBitsB & 0xffff) != 0);
    }

#ifdef RDTSC_TIMING
    dpSwitchTicks += Timer::rdtsc() - tscStart;
#endif
    return curScore;
  }

  void Eagle::updateTable(std::unordered_map <uint64, DPState> dpTable[], uint dpBestScores[],
			  uint64 m64, uint64 dist, uint n0, uint n1hapA, uint n2hapA, uint n1hapB,
			  uint n2hapB, uint score) const {
    if (n1hapB/2 == n2hapB/2) return; // disallow copying both haps from an indiv
    if ((n1hapA == n1hapB && n2hapA == n2hapB) || (n1hapA != n1hapB && n2hapA != n2hapB))
      score += computeStaticScore(n0, n1hapB, n2hapB, m64);
    else {
      if (n1hapA == n1hapB)
	score += computeSwitchScore(n0, n1hapA, n2hapA, n2hapB, m64);
      else /* (n2hapA == n2hapB) */
	score += computeSwitchScore(n0, n2hapA, n1hapA, n1hapB, m64);
    }
    if (!updateHelper(dpTable[m64], dpBestScores[m64], make_pair(n1hapA, n2hapA),
		      make_pair(n1hapB, n2hapB), score))
      return;
    for (uint64 m64x = m64+1; m64x < m64+dist && m64x < Mseg64; m64x++) {
      score += computeStaticScore(n0, n1hapB, n2hapB, m64x);
      if (!updateHelper(dpTable[m64x], dpBestScores[m64x], make_pair(n1hapB, n2hapB),
			make_pair(n1hapB, n2hapB), score))
	return;
    }
  }

  void updateErrHits(vector <uint> &hitVec, uint64 &bestErrLoc, uint64 errLoc, uint n2hap) {
    if (errLoc > bestErrLoc) {
      bestErrLoc = errLoc;
      hitVec.resize(1, n2hap);
    }
    else if (errLoc == bestErrLoc)
      hitVec.push_back(n2hap);
  }

  void Eagle::safeInsert(std::set <uint> &refHapSet, uint n1hap, uint n0) const {
    if (!maskIndivs[n1hap/2]) return;
    if (n1hap/2 == n0) return;
    refHapSet.insert(n1hap);
  }

  vector < pair <uint64, uint64> > Eagle::findGoodSegs(uint64 n0, uint64 nF1, uint64 nF2, uint64 n1hap) const {
    vector < pair <uint64, uint64> > goodSegs;
    uint64 firstGood = 0;
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      uint64_masks bits0 = genoBits[m64*N + n0];
      uint64_masks bitsF1 = genoBits[m64*N + nF1];
      uint64_masks bitsF2 = genoBits[m64*N + nF2];
      uint64 is0 = bits0.is0 | (~(bits0.is0|bits0.is2|bits0.is9) &
				((bitsF1.is0 & ~bitsF2.is0) | (~bitsF1.is2 & bitsF2.is2)));
      uint64 is2 = bits0.is2 | (~(bits0.is0|bits0.is2|bits0.is9) &
				((bitsF1.is2 & ~bitsF2.is2) | (~bitsF1.is0 & bitsF2.is0)));
      uint64 is1 = haploBits[m64*2*N + n1hap];
      uint64 wrongBits = (is0 & is1) | (is2 & ~is1);
      if (wrongBits) {
	uint64 firstWrong = m64*64 + popcount64((wrongBits & ~(wrongBits-1))-1);
	if (firstWrong - firstGood >= 64)
	  goodSegs.push_back(make_pair(firstGood, firstWrong));
	firstGood = m64*64 + (wrongBits ? 64 - __builtin_clzll(wrongBits) : 0);
      }
    }
    goodSegs.push_back(make_pair(firstGood, Mseg64*64));
    return goodSegs;
  }

  void Eagle::updateFarHaps(vector < pair <uint, uint> > &farHaps, uint n1hap, uint m64jStart, uint m64jEnd) const {
    const double cMminLen = 1.0;
    if (cMs64j[m64jEnd] - cMs64j[m64jStart] < cMminLen) return;
    //cout << m64jStart/64 << "." << (m64jStart&63) << "-" << m64jEnd/64 << "." << (m64jEnd&63) << " ";
    if (m64jEnd > farHaps[m64jStart].first) {
      farHaps[m64jStart].first = m64jEnd;
      farHaps[m64jStart].second = n1hap;
    }
  }

  double Eagle::runHMM(uint64 n0, uint64 nF1, uint64 nF2, int iter, uint beamWidth,
		       uint maxHapStates) {

    if (!maskIndivs[n0]) return 0;

    if (Mseg64 < 3U) {
      cerr << "Too few SNP segments for analysis" << endl;
      exit(1);
    }
    if (Mseg64 > 50000U) {
      cerr << "Too many SNP segments for analysis: " << Mseg64 << " (max = 50000)" << endl;
      exit(1);
    }

    uint seed = n0; // for rand_r()

    /***** FIND LONGEST DIP-HAP MATCHES *****/

#ifdef RDTSC_TIMING
    uint64 tscStart = Timer::rdtsc();
#endif
    vector <uint> topInds[2]; // [max err]; lengths are Mseg64 * K
    vector <uint> topIndsLens[2]; // [max err]; lengths are Mseg64

    const uint K = 100;
    Timer timer;
    findLongDipHap(n0, topInds, topIndsLens, K, iter >= 4 ? 0 : 1);
    double hapTime = timer.update_time();
#ifdef RDTSC_TIMING
    diphapTicks += Timer::rdtsc() - tscStart;
#endif
    /*
    // VISUALIZE RESULTS
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      cout << "m64 = " << m64 << endl;
      int m64xMin = std::max(0, (int) m64-5);
      int m64xMax = std::min(m64+11, Mseg64);
      for (int m64x = m64xMin; m64x < m64xMax; m64x++) {
	int numSnps = 0;
	for (int j = 0; j < 64; j++)
	  numSnps += maskSnps64j[m64x*64+j];
	printf("%2d ", numSnps);
      }
      cout << endl;
      for (uint e = 0; e < 2; e++) {
	cout << "e = " << e << endl;
	for (uint k = 0; k < topIndsLens[e][m64]; k++) {
	  uint64 n1hap = topInds[e][m64*K + k];
	  for (int m64x = m64xMin; m64x < m64xMax; m64x++) {
	    printf("%2d ", numDipHapWrongBits(m64x, n0, n1hap));
	  }
	  cout << "   (" << n1hap << ")" << endl;
	}
      }
    }
    */

    /***** EXTEND DIP-HAP MATCHES TO OBTAIN LONGEST MATCHES COVERING EACH SEG *****/

#ifdef RDTSC_TIMING
    uint64 tscExtStart = Timer::rdtsc();
#endif
    const int maxWrongBits = 1;
    std::set <DipHapSeg> longDipHapSegsForward[Mseg64];
    const uint maxSegs = std::max(50U, maxHapStates/4);
    vector <uint> lastEnd64(2*N);
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
#ifdef DETAILS
      if ((int) nF1 != -1) cout << "m64 = " << m64 << endl;
#endif
      for (uint e = 0; e < 2; e++) {
#ifdef DETAILS
	if ((int) nF1 != -1) cout << "e = " << e << endl;
#endif
	for (uint k = 0; k < topIndsLens[e][m64]; k++) {
	  uint64 n1hap = topInds[e][m64*K + k];
	  if (!maskIndivs[n1hap/2]) continue;
	  if (n1hap/2 == n0) continue;
	  if (lastEnd64[n1hap] > m64) continue;

	  // find start
	  uint segStart = m64;
	  while ((int) segStart >= 0 && numDipHapWrongBits(segStart, n0, n1hap) <= maxWrongBits)
	    segStart--;
	  if ((int) segStart < 0) segStart = 0; //segStart++; NEW: start one chunk before!
	  uint segStart64j = segStart*64 + firstDipHapGoodBit(segStart, n0, n1hap);

	  // find end
	  uint segEnd = m64+1;
	  while (segEnd < Mseg64 && numDipHapWrongBits(segEnd, n0, n1hap) <= maxWrongBits)
	    segEnd++;
	  uint segEnd64j = segEnd*64 + firstDipHapWrongBit(segEnd, n0, n1hap);

#ifdef DETAILS
	  checkHapPhase1(n0, nF1, nF2, n1hap, segStart, segEnd);
#endif
	  lastEnd64[n1hap] = segEnd;
	  for (uint64 m64x = segStart; m64x+1 /* extend forward */ < segEnd; m64x++) {
	    longDipHapSegsForward[m64x].insert(DipHapSeg(n1hap, segStart64j, segEnd64j));
	    if (longDipHapSegsForward[m64x].size() > maxSegs)
	      longDipHapSegsForward[m64x].erase(--longDipHapSegsForward[m64x].end());
	  }
	}
      }
    }
#ifdef RDTSC_TIMING
    extTicks += Timer::rdtsc() - tscExtStart;
#endif

    /***** COMPILE SET OF REFERENCE HAPLOTYPES FOR EACH TRANSITIONS AT EACH SEG *****/

    std::set <uint> refHapSets[Mseg64];
    vector <uint> refHapVecs[Mseg64];
    vector < pair <uint, uint> > refHapOppPairs[Mseg64];
    const uint numTopShort = maxHapStates/4, numTopLong = maxHapStates/4,
      minRefHaps = maxHapStates/4;
    const uint64 side = 1;
    uint64 curHaploBitsT[Mseg64];
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      uint cumTop = 0;
      // dip-hap hits starting at m64 and m64+1
      for (uint64 m64x = m64; m64x < m64+2 && m64x < Mseg64; m64x++)
	for (uint e = 1; e != -1U; e--) {
	  cumTop += numTopShort / 4;
	  for (uint k = 0; k < topIndsLens[e][m64] && refHapSets[m64].size() < cumTop; k++) {
	    uint n1hap = topInds[e][m64*K + k];
	    if (segConfs[n1hap*Mseg64+m64] < 5 &&
		(m64+1==Mseg64 || segConfs[n1hap*Mseg64+m64+1] < 5))
	      safeInsert(refHapSets[m64], n1hap, n0);
	  }
	}
      // long dip hap covering m64, m64+1
      cumTop += numTopLong;
      for (std::set <DipHapSeg>::iterator it = longDipHapSegsForward[m64].begin();
	   it != longDipHapSegsForward[m64].end() && refHapSets[m64].size() < cumTop; it++) {
	uint n1hap = it->n;
	if (segConfs[n1hap*Mseg64+m64] < 5 &&
	    (m64+1==Mseg64 || segConfs[n1hap*Mseg64+m64+1] < 5))
	  safeInsert(refHapSets[m64], n1hap, n0);
      }
#ifdef DETAILS
      cout << "m64 = " << m64 << ": " << refHapSets[m64].size() << endl;
#endif
      /***** AUGMENT REFERENCE HAPLOTYPE SET WITH COMPLEMENTS *****/
      vector <uint> refHapVec(refHapSets[m64].begin(), refHapSets[m64].end());
      if (m64+1 >= side && m64+1+side < Mseg64) { // look up m64+1 in LSH
	for (uint k1 = 0; k1 < refHapVec.size(); k1++) {
	  uint64 n1hap = refHapVec[k1];
	  // require no err on right part of m64
	  uint64 m64errMaskMax = 1ULL<<((uint64) rand_r(&seed)&63);

	  for (uint64 m64x = m64+1-side; m64x <= m64+1+side; m64x++) {
	    curHaploBitsT[m64x] = 0;
	    for (uint64 j = 0; j < 64; j++) {
	      uint64 m64j = m64x*64+j, bit = 0;
	      if (maskSnps64j[m64j]) {
		uint g0 = getGeno0123(m64j, n0);
		if (g0 == 0 || g0 == 3) bit = 0;
		else if (g0 == 2) bit = 1;
		else bit = 1-((haploBitsT[n1hap*Mseg64 + m64x]>>j)&1);
	      }
	      curHaploBitsT[m64x] |= bit<<j;
	    }
	  }
	  vector <uint> hits[3]; // n2haps with 1st err at [1]: m64+1; [2]: m64+2 (allow 1@m64+1)
	  vector <uint64> bestErrLocs(3); // farthest error locations seen
	  for (uint h = 0; h < hashLookups[m64+1].size(); h++) { // for each hashing
	    uint numHits;
#ifdef RDTSC_TIMING
	    uint64 tscLshStart = Timer::rdtsc();
#endif
	    const uint *lenHapInds =
	      hashLookups[m64+1][h].query(computeHash(curHaploBitsT, hashBits[m64+1][h]));
#ifdef RDTSC_TIMING
	    lshTicks += Timer::rdtsc() - tscLshStart;
	    uint64 tscLshCheckStart = Timer::rdtsc();
#endif
	    if (lenHapInds == NULL)
	      numHits = 0;
	    else
	      numHits = lenHapInds[0];
	    for (uint k = 1; k <= numHits; k++) {
	      uint64 n2hap = lenHapInds[k];
	      if (!maskIndivs[n2hap/2]) continue;
	      if (n2hap/2 == n0) continue;
	      if (n2hap/2 == n1hap/2) continue; // opp haps of same indiv
	      int errFail = 0, err1 = 0; // err1: 1 err in m64+1
	      uint64 err1Loc = 0, err2Loc = 0;
	      for (uint64 m64x = m64; m64x <= m64+2; m64x++) {
		uint64 n1is1 = haploBitsT[n1hap*Mseg64 + m64x];
		uint64 n2is1 = haploBitsT[n2hap*Mseg64 + m64x];
		const uint64_masks &bits0 = genoBits[m64x*N + n0];
		uint64 wrongBits = (bits0.is0 & (n1is1 | n2is1)) | (bits0.is2 & ~(n1is1 & n2is1))
		  | (~(bits0.is0|bits0.is2|bits0.is9) & ~(n1is1 ^ n2is1));
		if (m64x == m64) {
		  if (wrongBits >= m64errMaskMax) { // not perfect in right part of m64
		    errFail = 2; // fail
		    break;
		  }
		}
		else if (m64x == m64+1) {
		  if (wrongBits & (wrongBits-1)) { // 2+ err
		    errFail = 1; // fail
		    err1Loc = (wrongBits & ~(wrongBits-1))-1; // LSB-1: big is good
		    break;
		  }
		  else if (wrongBits != 0)
		    err1 = 1;
		}
		else {
		  err2Loc = (wrongBits & ~(wrongBits-1))-1; // LSB-1: big is good
		}
	      }
	      if (!errFail) { // did not fail
		if (err2Loc != 0) err2Loc -= err1; // slightly worse to have 1 error in m64+1
		updateErrHits(hits[2], bestErrLocs[2], err2Loc, n2hap);
	      }
	      else if (errFail == 1)
		updateErrHits(hits[1], bestErrLocs[1], err1Loc, n2hap);
	    }
#ifdef RDTSC_TIMING
	    lshCheckTicks += Timer::rdtsc() - tscLshCheckStart;
#endif
	    if (bestErrLocs[2] == -1ULL) break; // early exit if perfect match found
	  }

	  //checkHapPhase1(n0, nF1, nF2, refHapVec[k1], m64, min(m64+3, Mseg64));
	  uint64 n2hap = -1ULL;
	  for (int xLoc = 2; xLoc >= 0; xLoc--)
	    if (!hits[xLoc].empty()) {
	      n2hap = hits[xLoc][rand_r(&seed) % hits[xLoc].size()];
	      //cout << "2." << popcount64(bestErrLocs[2]) << ": ";
	      break;
	    }
	  if (n2hap != -1ULL) {
#ifdef DETAILS
	    checkHapPhase1(n0, nF1, nF2, n2hap, m64, min(m64+3, Mseg64));
#endif
	    safeInsert(refHapSets[m64], n2hap, n0);
	    refHapOppPairs[m64].push_back(make_pair((uint) n1hap, (uint) n2hap));
	  }
	}
      }

      // make sure at least some ref haps are chosen: relax earlier conf filter
      for (std::set <DipHapSeg>::iterator it = longDipHapSegsForward[m64].begin();
	   it != longDipHapSegsForward[m64].end() && refHapSets[m64].size() < minRefHaps/2; it++)
	safeInsert(refHapSets[m64], it->n, n0);
      for (uint64 m64x = m64; m64x < m64+2 && m64x < Mseg64; m64x++)
	for (uint e = 0; e < 2; e++)
	  for (uint k = 0; k < topIndsLens[e][m64] && refHapSets[m64].size() < minRefHaps; k++)
	    safeInsert(refHapSets[m64], topInds[e][m64*K + k], n0);

      refHapVecs[m64] = vector <uint> (refHapSets[m64].begin(), refHapSets[m64].end());
      //cout << "m64 = " << m64 << ": " << refHapSets[m64].size() << endl;
    }

#ifdef DETAILS
    // TEMPORARY: CHECKING TRIO HAP - REF HAP MATCHES
    std::set <DipHapSeg> longTrioHapHapSegsForward[2][Mseg64];
    uint longestTrioHaps[2][Mseg64]; memset(&longestTrioHaps[0][0], 0, 2*Mseg64*sizeof(longestTrioHaps[0][0]));
    for (int p = 0; p < 2; p++) {
      cout << "Checking trio hap - ref hap matches" << endl;

      for (uint64 n1hap = 0; n1hap < 2*N; n1hap++) {
	if (n1hap % 10000 == 0) cout << "at " << n1hap << endl;
	if (!maskIndivs[n1hap/2]) continue;
	if (n1hap/2 == n0) continue;
	vector < pair <uint64, uint64> > goodSegs
	  = findGoodSegs(n0, p==0?nF1:nF2, p==0?nF2:nF1, n1hap);
	for (uint s = 0; s < goodSegs.size(); s++)
	  for (uint64 m64x = (goodSegs[s].first + 32) / 64; m64x < (goodSegs[s].second + 32) / 64;
	       m64x++) {
	    if (goodSegs[s].second - goodSegs[s].first > longestTrioHaps[p][m64x]) {
	      longTrioHapHapSegsForward[p][m64x].insert(DipHapSeg(n1hap, goodSegs[s].first,
								  goodSegs[s].second));
	      if (longTrioHapHapSegsForward[p][m64x].size() > 10) {
		longTrioHapHapSegsForward[p][m64x].erase(--longTrioHapHapSegsForward[p][m64x].end());
		longestTrioHaps[p][m64x] = ((--longTrioHapHapSegsForward[p][m64x].end())->end)
		  - ((--longTrioHapHapSegsForward[p][m64x].end())->start);
	      }
	    }
	  }
      }
    }
#endif
    /***** RUN HMM (BEAM SEARCH) *****/

    std::unordered_map <uint64, DPState> dpTable[Mseg64];
    uint dpBestScores[Mseg64]; for (uint64 m64 = 0; m64 < Mseg64; m64++) dpBestScores[m64] = 1<<30;
    pair <uint, uint> finalHapPairs[Mseg64];
    vector <int> bestPathScores(Mseg64);

#ifdef RDTSC_TIMING
    uint64 tscDpStart = Timer::rdtsc();
#endif
    for (uint64 m64 = 0; m64 <= Mseg64; m64++) {
      // find best from prev
#ifdef RDTSC_TIMING
      uint64 tscDpSortStart = Timer::rdtsc();
#endif
      vector <DPState> curStates;
      if (m64 == 0) curStates.push_back(DPState(0, make_pair(-1U, -1U)));
      if (m64 > 0) {
	uint maxDiffScore = 2*switchCost;
	uint bestScore = dpBestScores[m64-1];

	// prune best states to beamWidth (via bucket sort)
	vector < vector <uint64> > curStateBuckets(maxDiffScore+1);
	for (std::unordered_map <uint64, DPState>::iterator it = dpTable[m64-1].begin();
	     it != dpTable[m64-1].end(); it++)
	  if (it->second.score <= bestScore + maxDiffScore)
	    curStateBuckets[it->second.score - bestScore].push_back(it->first);
	for (uint dScore = 0; dScore <= maxDiffScore && curStates.size() < beamWidth; dScore++) {
	  //if (curStates.size() + curStateBuckets[dScore].size() > beamWidth)
	  //sort(curStateBuckets[dScore].begin(), curStateBuckets[dScore].end());
	  for (uint k = 0; k < curStateBuckets[dScore].size() && curStates.size() < beamWidth; k++)
	    curStates.push_back(DPState(bestScore+dScore, ullToPair(curStateBuckets[dScore][k])));
	}
	/* ---- SLOW ----
	for (std::unordered_map <uint64, DPState>::iterator it = dpTable[m64-1].begin();
	     it != dpTable[m64-1].end(); it++)
	  if (it->second.score <= bestScore + maxDiffScore)
	    curStates.push_back(DPState(it->second.score, ullToPair(it->first))); // .from = prev
        sort(curStates.begin(), curStates.end());
	*/
      }
#ifdef RDTSC_TIMING
      dpSortTicks += Timer::rdtsc() - tscDpSortStart;
#endif
      if (m64 < Mseg64) {
	const vector <uint> &refHapVec = refHapVecs[m64]; uint Kref = refHapVec.size();

	// iterate through best in beam
	if (m64 > 0) {
	  for (uint s = 0; s < curStates.size() && s < beamWidth; s++) {
	    uint n1hap = curStates[s].from.first, n2hap = curStates[s].from.second;
	    uint score = curStates[s].score;
	    // continue
	    updateTable(dpTable, dpBestScores, m64, 1, n0, n1hap, n2hap, n1hap, n2hap, score);
	  }
	  for (uint s = 0; s < curStates.size() && s < beamWidth; s++) {
	    uint n1hap = curStates[s].from.first, n2hap = curStates[s].from.second;
	    uint score = curStates[s].score;
	    // switch n1hap or n2hap
	    for (uint k = 0; k < Kref; k++) {
	      updateTable(dpTable, dpBestScores, m64, 3, n0, n1hap, n2hap, refHapVec[k], n2hap,
			  score+switchCost);
	      updateTable(dpTable, dpBestScores, m64, 3, n0, n1hap, n2hap, n1hap, refHapVec[k],
			  score+switchCost);
	    }
	  }
	}

	// clean start (from curStates[0])
	if (m64 == 0) {
	  for (uint k1 = 0; k1 < Kref; k1++)
	    for (uint k2 = k1+1; k2 < Kref; k2++)
	      updateTable(dpTable, dpBestScores, m64, 3, n0, curStates[0].from.first,
			  curStates[0].from.second, refHapVec[k1], refHapVec[k2],
			  curStates[0].score+2*switchCost);
	}
	else { // only try clean start from opp pairs
	  for (uint k12 = 0; k12 < refHapOppPairs[m64].size(); k12++)
	    updateTable(dpTable, dpBestScores, m64, 3, n0, curStates[0].from.first,
			curStates[0].from.second, refHapOppPairs[m64][k12].first,
			refHapOppPairs[m64][k12].second, curStates[0].score+3*switchCost);
	}
      }
      else { // finished; backtrack through DP table
	finalHapPairs[Mseg64-1] = curStates[0].from;
	bestPathScores[Mseg64-1] = curStates[0].score;
	for (int m64x = Mseg64-2; m64x >= 0; m64x--) {
	  finalHapPairs[m64x] = dpTable[m64x+1][pairToULL(finalHapPairs[m64x+1])].from;
	  //cout << "m64 = " << m64x << ": n1hap = " << finalHapPairs[m64x].first << ", n2hap = " << finalHapPairs[m64x].second << "   score = " << dpTable[m64x+1][pairToULL(finalHapPairs[m64x+1])].score << endl;
	  bestPathScores[m64x] = dpTable[m64x+1][pairToULL(finalHapPairs[m64x+1])].score;
	}
      }
    }
#ifdef RDTSC_TIMING
    dpTicks += Timer::rdtsc() - tscDpStart;
#endif

    /***** FIND TRANSITIONS WITHIN 64-SNP CHUNKS *****/

    uint64 hmmHaploBitsT[2][Mseg64], fixHaploBitsT[2][Mseg64], hetErrMasks[Mseg64],
      uncertainMasks[Mseg64];
    vector <bool> phaseVec;
    vector <int> n1haps, n2haps, n3haps, signs;
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      uint64 n1hap, n2hap, n3hap; int sign;
      if (m64 == 0) {
	n1hap = finalHapPairs[m64].first;
	n2hap = n3hap = finalHapPairs[m64].second;
	sign = 1;
      }
      else if (finalHapPairs[m64].first == finalHapPairs[m64-1].first) {
	n1hap = finalHapPairs[m64].first;
	n2hap = finalHapPairs[m64-1].second;
	n3hap = finalHapPairs[m64].second;
	sign = 1;
      }
      else if (finalHapPairs[m64].second == finalHapPairs[m64-1].second) {
	n1hap = finalHapPairs[m64].second;
	n2hap = finalHapPairs[m64-1].first;
	n3hap = finalHapPairs[m64].first;
	sign = -1;
      }
      else { // restart
	n1hap = finalHapPairs[m64].first;
	n2hap = n3hap = finalHapPairs[m64].second;
	sign = 1;
      }
      n1haps.push_back(n1hap); n2haps.push_back(n2hap); n3haps.push_back(n3hap); signs.push_back(sign);
      pair <uint64, uint64> phaseBits
	= phaseSegHMM(n0, n1hap, n2hap, n3hap, m64, hetErrMasks[m64]);
      hmmHaploBitsT[0][m64] = phaseBits.first;
      hmmHaploBitsT[1][m64] = phaseBits.second;
      if (sign == -1) std::swap(hmmHaploBitsT[0][m64], hmmHaploBitsT[1][m64]);
      //checkHapPhase(n0, nF1, nF2, hmmHaploBitsT[0], m64, 0);

      //vector <bool> phaseSeg = checkHapPhase2(n0, nF1, nF2, n1hap, n2hap, n3hap, m64, sign);
      //cout << endl;
      //phaseVec.insert(phaseVec.end(), phaseSeg.begin(), phaseSeg.end());
    }

    /***** DETECT AND FIX BLIPS BASED ON HAPLOTYPE FREQUENCIES *****/

    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
#ifdef RDTSC_TIMING
      uint64 tscBlipFixStart = Timer::rdtsc();
#endif
      vector < vector <int> > votes(64, vector <int> (4));
      uint64 m64mid = min(max(m64, side), Mseg64-1-side);

      std::unordered_set <uint> seen;
      //cout << "m64 = " << m64 << endl;
      for (int opp = 0; opp < 2; opp++)
	for (uint h = 0; h < min((uint) hashLookups[m64mid].size(), 10U); h++) { // for each hash
	  uint numHits;
#ifdef RDTSC_TIMING
	  uint64 tscBlipLshStart = Timer::rdtsc();
#endif
	  const uint *lenHapInds =
	    hashLookups[m64mid][h].query(computeHash(hmmHaploBitsT[opp], hashBits[m64mid][h]));
#ifdef RDTSC_TIMING
	  blipLshTicks += Timer::rdtsc() - tscBlipLshStart;
#endif
	  if (lenHapInds == NULL)
	    numHits = 0;
	  else
	    numHits = lenHapInds[0];
	  for (uint k = 1; k <= min(numHits, 25U); k++) {
	    uint64 nXhap = lenHapInds[k];
	    if (!maskIndivs[nXhap/2]) continue;
	    if (nXhap/2 == n0) continue;
	    if (seen.count(nXhap)) continue;
	    seen.insert(nXhap);

#ifdef RDTSC_TIMING
	    uint64 tscBlipPopStart = Timer::rdtsc();
#endif
	    int errs = 0;
	    for (uint64 m64x = m64mid-side; m64x <= m64mid+side; m64x++) {
	      const uint64_masks &bits0 = genoBits[m64x*N + n0];
	      errs += popcount64((hmmHaploBitsT[opp][m64x] ^ haploBitsT[nXhap*Mseg64+m64x])
				 & ~bits0.is9);
	    }
#ifdef RDTSC_TIMING
	    blipPopTicks += Timer::rdtsc() - tscBlipPopStart;
	    uint64 tscBlipVoteStart = Timer::rdtsc();
#endif
	    if (errs <= 2) {
	      //checkHapPhase(n0, nF1, nF2, haploBitsT + nXhap*Mseg64, m64mid, side);
	      for (uint64 j = 0; j < 64; j++)
		votes[j][(((haploBitsT[nXhap*Mseg64+m64]>>j)&1)^opp) + 2*opp]++;
	    }
#ifdef RDTSC_TIMING
	    blipVoteTicks += Timer::rdtsc() - tscBlipVoteStart;
#endif
	  }
	}

      fixHaploBitsT[0][m64] = hmmHaploBitsT[0][m64]; fixHaploBitsT[1][m64] = hmmHaploBitsT[1][m64];
      const uint64_masks &bits0 = genoBits[m64*N + n0];
      for (uint64 j = 0; j < 64; j++) {
	if ((bits0.is0|bits0.is2|bits0.is9)&(1ULL<<j)) continue; // not het
	double relVoteDiff = (hetErrMasks[m64]&(1ULL<<j)) ? 2 : 10; // weak if uncertain phase call
	double eps = 0.5; // pseudocount
	double ratioOR = (votes[j][0]+eps)*(votes[j][2]+eps)/(votes[j][1]+eps)/(votes[j][3]+eps);
	if (ratioOR > relVoteDiff) {
	  fixHaploBitsT[0][m64] &= ~(1ULL<<j);
	  fixHaploBitsT[1][m64] |= 1ULL<<j;
	}
	if (ratioOR < 1.0/relVoteDiff) {
	  fixHaploBitsT[0][m64] |= 1ULL<<j;
	  fixHaploBitsT[1][m64] &= ~(1ULL<<j);
	}
      }
#ifdef RDTSC_TIMING
      blipFixTicks += Timer::rdtsc() - tscBlipFixStart;
#endif
      /*
      if ((int) nF1 != -1) { // output blip fix + vote info
	if (fixHaploBitsT[0][m64] != hmmHaploBitsT[0][m64])
	  cout << "*** ";
	cout << "m64 = " << m64 << ": ";
	checkHaploBits(n0, nF1, nF2, hmmHaploBitsT[0][m64], m64);
	if (fixHaploBitsT[0][m64] != hmmHaploBitsT[0][m64]) {
	  cout << " -> ";
	  checkHapPhase(n0, nF1, nF2, fixHaploBitsT[0], m64, 0);
	}
	cout << " ";
	checkHapPhase(n0, nF1, nF2, hmmHaploBitsT[0], m64, 0, votes);
      }
      */
    }
    // set phaseConfs2: 0|255 unless hetErr or flip or miss; 1|254 at those sites
    uchar *phaseConfsHap0, *phaseConfsHap1;
    if (Nref) { // allocate temp arrays for phaseConfs2 haplotype confidences
      phaseConfsHap0 = ALIGNED_MALLOC_UCHARS(Mseg64*64);
      phaseConfsHap1 = ALIGNED_MALLOC_UCHARS(Mseg64*64);
    }
    else {
      phaseConfsHap0 = phaseConfs2 + 2*n0*Mseg64*64;
      phaseConfsHap1 = phaseConfs2 + (2*n0+1)*Mseg64*64;
    }
    uchar hapConfs[2][2];
    hapConfs[0][0] = 0; hapConfs[0][1] = 255; // [0][*]: not uncertain
    hapConfs[1][0] = 1; hapConfs[1][1] = 254; // [1][*]: uncertain
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      uncertainMasks[m64] = hetErrMasks[m64] | (fixHaploBitsT[0][m64] ^ hmmHaploBitsT[0][m64])
	| genoBits[m64*N + n0].is9;
      for (uint64 j = 0; j < 64; j++) {
	uint64 m64j = m64*64 + j;
	bool uncertain = (uncertainMasks[m64]>>j)&1;
	phaseConfsHap0[m64j] = hapConfs[uncertain][(fixHaploBitsT[0][m64]>>j)&1];
	phaseConfsHap1[m64j] = hapConfs[uncertain][(fixHaploBitsT[1][m64]>>j)&1];
      }
    }

#ifdef CHECK_TRUE_DIP_HAP
    /***** CHECK TRUE BEST DIP-HAP: ARE THEY RESPECTED BY FINAL PHASE? *****/
    cout << "CHECK TRUE BEST DIP-HAP: ARE THEY RESPECTED BY FINAL PHASE?" << endl;
    std::set <DipHapSeg> segsOutput;
    for (uint64 m64 = 0; m64 < Mseg64; m64++)
      for (int p = 0; p < 2; p++) {
	//cout << "m64 = " << m64 << " true longest cover, parent " << p+1 << endl;
	for (std::set <DipHapSeg>::iterator it = longTrioHapHapSegsForward[p][m64].begin();
	     it != longTrioHapHapSegsForward[p][m64].end(); it++) {
	  if (!segsOutput.count(*it)) {
          segsOutput.insert(*it);
	  /*
	  cout << "n1hap = " << it->n << "; m64 = [" << it->start/64 << "." << (it->start&63) << "," << it->end/64 << "." << (it->end&63) << "): ";
	  cout << ((std::find(refHapVecs[m64].begin(), refHapVecs[m64].end(), it->n)
		    != refHapVecs[m64].end()) ? "YES" : "\033[1;33mNO\033[0m") << endl;//" ";
	  */
	  //checkHapPhase1j(n0, nF1, nF2, it->n, it->start, it->end); cout << endl;
	  vector <bool> ret = checkHapPhase1jCall(n0, nF1, nF2, fixHaploBitsT[p], it->start, it->end, false);
	  if (find(ret.begin(), ret.end(), 0) != ret.end() &&
	      find(ret.begin(), ret.end(), 1) != ret.end())
	    checkHapPhase1jCall(n0, nF1, nF2, fixHaploBitsT[p], it->start, it->end, true);
	  }
	  break;
	}
      }
#endif

    /***** POST-PROCESS: FIND LONG DIP-HAP MATCHES ALMOST CONSISTENT WITH PHASING *****/
    vector < pair <uint, uint> > farHaps(Mseg64*64);
    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      //cout << "m64 = " << m64 << ":" << endl;
      for (uint e = 0; e < 2; e++) {
	//cout << "e = " << e << ":" << endl;
	for (uint k = 0; k < topIndsLens[e][m64] && k < 10U /*TODO*/; k++) {
	  uint64 n1hap = topInds[e][m64*K + k];
	  // TODO: save this info from earlier
	  if (!maskIndivs[n1hap/2]) continue;
	  if (n1hap/2 == n0) continue;
	  // find start
	  uint segStart = m64;
	  while ((int) segStart >= 0 && numDipHapWrongBits(segStart, n0, n1hap) <= maxWrongBits)
	    segStart--;
	  if ((int) segStart < 0) segStart = 0;
	  uint segStart64j = segStart*64 + firstDipHapGoodBit(segStart, n0, n1hap);
	  // find end
	  uint segEnd = m64+1;
	  while (segEnd < Mseg64 && numDipHapWrongBits(segEnd, n0, n1hap) <= maxWrongBits)
	    segEnd++;
	  uint segEnd64j = segEnd*64 + firstDipHapWrongBit(segEnd, n0, n1hap);

	  vector <uint> err64j(1, segStart64j);
	  bool foundHet = false, prevPhase = false;
	  for (uint64 m64j = segStart64j; m64j < segEnd64j; m64j++)
	    if (maskSnps64j[m64j]) {
	      uint64 m64 = m64j/64, j = m64j&63;
	      uint g0 = getGeno0123(m64j, n0); // TODO: speed up
	      bool relPhase = ((haploBitsT[n1hap*Mseg64+m64] ^ fixHaploBitsT[0][m64])>>j)&1;
	      if (g0 == 1) {
		if (!foundHet)
		  foundHet = true;
		else if (relPhase != prevPhase)
		  err64j.push_back(m64j);
		prevPhase = relPhase;
	      }
	      else if (g0 == 0 || g0 == 2) {
		if (relPhase) { // dip-hap err
		  err64j.push_back(m64j);
		  err64j.push_back(m64j);
		}
	      }
	  }
	  err64j.push_back(segEnd64j);
	  /*
          // output post-process debug info
          cout << "true: "; checkHapPhase1j(n0, nF1, nF2, n1hap, segStart64j, segEnd64j); cout << endl;
	  cout << "call: "; checkHapPhase1jCall(n0, nF1, nF2, fixHaploBitsT[0], segStart64j, segEnd64j, true);
	  cout << "err vs. call:";
	  for (uint i = 0; i < err64j.size(); i++) {
	    cout << " " << err64j[i]/64 << "." << (err64j[i]&63);
	    printf("(%.2f)", cMs64j[err64j[i]]);
	  }
	  cout << endl;
	  */

	  const double cMconsecMin = 1.5, cMendMin = 1.5;
	  uint iStart = 0;
	  for (uint i = 1; i < err64j.size(); i++) {
	    double cMseg = cMs64j[err64j[i]] - cMs64j[err64j[i-1]];
	    if (i == iStart+1) { // piece ending at err64j[i] is first in new chunk
	      if (cMseg < cMendMin) // no good; can't start yet
		iStart = i;
	    }
	    else {
	      double cMprev = cMs64j[err64j[i-1]] - cMs64j[err64j[i-2]];
	      if (cMprev < cMconsecMin && cMseg < cMconsecMin) { // consec short => split
		// deal with chunk that just ended at either err64j[i-2] or err64j[i-1]
		if (cMprev < cMendMin) // last is too short
		  updateFarHaps(farHaps, n1hap, err64j[iStart], err64j[i-2]);
		else
		  updateFarHaps(farHaps, n1hap, err64j[iStart], err64j[i-1]);
		// deal with beginning of next chunk
		if (cMseg < cMendMin)
		  iStart = i;
		else
		  iStart = i-1;
	      }
	    }
	    if (i+1 == err64j.size()) {
	      if (cMseg < cMendMin)
		updateFarHaps(farHaps, n1hap, err64j[iStart], err64j[i-1]);
	      else
		updateFarHaps(farHaps, n1hap, err64j[iStart], err64j[i]);
	    }
	  }
	  //cout << endl;
	}
      }
    }
    uint farEnd = 0;
    vector <DipHapSeg> hapSegs;
    for (uint64 segStart64j = 0; segStart64j < Mseg64*64; segStart64j++)
      if (farHaps[segStart64j].first > farEnd) {
	farEnd = farHaps[segStart64j].first;
	hapSegs.push_back(DipHapSeg(farHaps[segStart64j].second, segStart64j, farEnd));
      }
    std::set <uint> postSwitches;
    for (uint i = 0; i < hapSegs.size(); i++) {
      uint64 n1hap = hapSegs[i].n, segStart64j = hapSegs[i].start, segEnd64j = hapSegs[i].end;
      uint64 useStart64j = (i==0 || hapSegs[i-1].end < segStart64j) ?
	segStart64j : (hapSegs[i-1].end + segStart64j)/2;
      uint64 useEnd64j = (i+1==hapSegs.size() || segEnd64j < hapSegs[i+1].start) ?
	segEnd64j : (segEnd64j + hapSegs[i+1].start)/2;
      /*
      // output post-process debug info
      double cMstart = cMs64j[segStart64j], cMend = cMs64j[segEnd64j];
      if ((int) nF1 != -1) {
	updateFarHaps(farHaps, n1hap, segStart64j, segEnd64j); // just to print
	cout << endl;
	cout << "true: "; checkHapPhase1j(n0, nF1, nF2, n1hap, segStart64j, segEnd64j); cout << endl;
	cout << "call: "; checkHapPhase1jCall(n0, nF1, nF2, fixHaploBitsT[0], segStart64j, segEnd64j, true);
	printf("%.2f-%.2f cM; ", cMstart, cMend);
	cout << "use " << useStart64j/64 << "." << (useStart64j&63) << "-"
	     << useEnd64j/64 << "." << (useEnd64j&63) << endl;
      }
      */

      vector <int> numSinceSwitches; vector <uint64> switchLocs;
      bool foundHet = false, prevPhase = false; int numSinceSwitch = 0;
      for (uint64 m64j = segStart64j; m64j < segEnd64j; m64j++)
	if (maskSnps64j[m64j]) {
	  uint64 m64 = m64j/64, j = m64j&63;
	  uint g0 = getGeno0123(m64j, n0); // TODO: speed up
	  if (g0 == 1) {
	    numSinceSwitch++;
	    bool relPhase = ((haploBitsT[n1hap*Mseg64+m64] ^ fixHaploBitsT[0][m64])>>j)&1;
	    if (!foundHet)
	      foundHet = true;
	    else if (relPhase != prevPhase) {
	      numSinceSwitches.push_back(numSinceSwitch);
	      switchLocs.push_back(m64j);
	      numSinceSwitch = 0;
	    }
	    prevPhase = relPhase;
	  }
	}
      numSinceSwitches.push_back(numSinceSwitch);
      for (uint s = 0; s < switchLocs.size(); s++)
	if (numSinceSwitches[s] > 2 && numSinceSwitches[s+1] > 2 &&
	    useStart64j <= switchLocs[s] && switchLocs[s] < useEnd64j) {
	  uint64 m64j = switchLocs[s];
	  /* output post-process debug info
	  if ((int) nF1 != -1) {
	    cout << "SUGGEST SWITCH: " << m64j/64 << "." << (m64j&63);
	    printf("(%.2f; split %.2f=%.2f+%.2f)", cMs64j[m64j], cMend-cMstart, cMs64j[m64j]-cMstart, cMend-cMs64j[m64j]);
	  }
	  */
	  postSwitches.insert(m64j);
	}
    }

    // apply switches
    int sign = 1;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      uint64 m64 = m64j/64, j = m64j&63;
      if (postSwitches.count(m64j))
	sign = -sign;
      if (sign == -1) {
	uint64 tmp = (fixHaploBitsT[0][m64] ^ fixHaploBitsT[1][m64]) & (1ULL<<j);
	fixHaploBitsT[0][m64] ^= tmp; fixHaploBitsT[1][m64] ^= tmp;
	tmp = (hmmHaploBitsT[0][m64] ^ hmmHaploBitsT[1][m64]) & (1ULL<<j);
	hmmHaploBitsT[0][m64] ^= tmp; hmmHaploBitsT[1][m64] ^= tmp;
	std::swap(phaseConfsHap0[m64j], phaseConfsHap1[m64j]);
      }
    }

    if (Nref) { // store phased haploBits for target sample
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
	if (!maskSnps64j[m64j]) continue;
	uint64 hapBits[2];
	hapBits[0] = (int) phaseConfsHap0[m64j] >= 128;
	hapBits[1] = (int) phaseConfsHap1[m64j] >= 128;
	for (uint64 h = 0; h <= 1ULL; h++) {
	  uint64 nTargetHap = 2*(n0-Nref) + h;
	  tmpHaploBitsT[nTargetHap*Mseg64 + (m64j/64)] |= hapBits[h]<<(m64j&63);
	}
      }
      ALIGNED_FREE(phaseConfsHap1);
      ALIGNED_FREE(phaseConfsHap0);
    }

    if ((int) nF1 != -1) { // final output
      cout << "3rd-iter phase:" << endl;
      for (uint64 m64 = 0; m64 < Mseg64; m64++) {
	/*
	checkHaploBits(n0, nF1, nF2, hmmHaploBitsT[0][m64], m64);
	if (fixHaploBitsT[0][m64] != hmmHaploBitsT[0][m64]) {
	  cout << " -> ";
	  checkHaploBits(n0, nF1, nF2, fixHaploBitsT[0][m64], m64);
	}
	*/
	vector <bool> phaseSeg = checkHaploBits(n0, nF1, nF2, fixHaploBitsT[0][m64], m64, -1);
	//cout << endl;
	phaseVec.insert(phaseVec.end(), phaseSeg.begin(), phaseSeg.end());
      }
      printf("# major SE: %2d   # tot SE: %2d / %d\n", countMajorSE(phaseVec), countSE(phaseVec),
	     (int) phaseVec.size()-1);
      fflush(stdout);
    }

#ifdef RDTSC_TIMING
    totTicks += Timer::rdtsc() - tscStart;
#endif
    return hapTime;
  }

  void Eagle::writePhaseConfs(const string &tmpPhaseFile) const {
    FILE *fout = fopen(tmpPhaseFile.c_str(), "wb");
    fwrite(phaseConfs, 1, 2*N*Mseg64*64, fout);
    fclose(fout);
  }

  void Eagle::readPhaseConfs(const string &tmpPhaseFile) {
    FILE *fout = fopen(tmpPhaseFile.c_str(), "rb");
    fread(phaseConfs, 1, 2*N*Mseg64*64, fout);
    fclose(fout);
  }

  void Eagle::cpPhaseConfs(uint64 n0start, uint64 n0end) {
    memcpy(phaseConfs + 2*n0start*Mseg64*64, phaseConfs2 + 2*n0start*Mseg64*64,
	   2*(n0end-n0start)*Mseg64*64);
  }

  void Eagle::cpTmpHaploBitsT(uint64 n0start, uint64 n0end) {
    memcpy(haploBitsT + 2*n0start*Mseg64, tmpHaploBitsT + 2*n0start*Mseg64,
	   2*(n0end-n0start)*Mseg64 * sizeof(haploBitsT[0]));
    for (uint64 nHap = 2*n0start; nHap < 2*n0end; nHap++)
      for (uint64 m64 = 0; m64 < Mseg64; m64++)
	haploBits[m64*2*N + nHap] = haploBitsT[nHap*Mseg64 + m64];
  }

  void Eagle::outputSE(const vector <uint> &children, const vector <uint> &nF1s,
		       const vector <uint> &nF2s, int step) const {

    if (children.empty()) return;
    vector <int> majorSEs, totSEs; vector <double> majorSErates, totSErates;
    for (uint att = 0; att < children.size(); att++) {
      printf("Switch error locations (step %d):", step);
      vector <bool> phaseVec = checkPhaseConfsPhase(children[att], nF1s[att], nF2s[att]);
      majorSEs.push_back(countMajorSE(phaseVec));
      totSEs.push_back(countSE(phaseVec));
      majorSErates.push_back(majorSEs.back() * 100.0 / (phaseVec.size()-1));
      totSErates.push_back(totSEs.back() * 100.0 / (phaseVec.size()-1));
      printf("# major SE: %2d   # tot SE: %2d / %d   (step %d)\n", majorSEs.back(),
	     totSEs.back(), (int) phaseVec.size()-1, step);
    }
    sort(majorSEs.begin(), majorSEs.end());
    sort(totSEs.begin(), totSEs.end());
    int useTrios = 70;
    printf("%d-trio avg # major SE: %.2f   avg # tot SE: %.2f   (step %d)\n", useTrios,
	   std::accumulate(majorSEs.begin(), majorSEs.begin()+useTrios, 0) / (double) useTrios,
	   std::accumulate(totSEs.begin(), totSEs.begin()+useTrios, 0) / (double) useTrios, step);
    sort(majorSErates.begin(), majorSErates.end());
    sort(totSErates.begin(), totSErates.end());
    printf("median major SE rate: %.2f%%  median tot SE rate: %.2f%%   (step %d)\n",
	   majorSErates[(nF1s.size()-1)/2], totSErates[(nF1s.size()-1)/2], step);
    fflush(stdout);
  }

  void Eagle::writeHapsGzSample(const string &prefix) const {
    FileUtils::AutoGzOfstream hapsGzOut; hapsGzOut.openOrExit(prefix + ".haps.gz");
    uint64 m = 0; // index in snps vector
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
      if (!maskSnps64j[m64j]) continue;
      hapsGzOut << snps[m].chrom << " " << snps[m].ID << " " << snps[m].physpos
		<< " " << snps[m].allele1 << " " << snps[m].allele2;
      for (uint64 n0 = 0; n0 < N; n0++) {
	int hapBit1, hapBit2;
	if (phaseConfs != NULL) {
	  hapBit1 = (int) phaseConfs[2*n0*Mseg64*64 + m64j] < 128;
	  hapBit2 = (int) phaseConfs[(2*n0+1)*Mseg64*64 + m64j] < 128;
	}
	else {
	  hapBit1 = !((haploBits[(m64j/64)*2*N + 2*n0]>>(m64j&63))&1);
	  hapBit2 = !((haploBits[(m64j/64)*2*N + 2*n0+1]>>(m64j&63))&1);
	}
	if (isFlipped64j[m64j]) {
	  hapBit1 = !hapBit1;
	  hapBit2 = !hapBit2;
	}
	hapsGzOut << " " << hapBit1 << " " << hapBit2;
      }
      hapsGzOut << endl;
      m++;
    }
    hapsGzOut.close();

    FileUtils::AutoGzOfstream sampleOut; sampleOut.openOrExit(prefix + ".sample");
    sampleOut << std::setprecision(3);
    sampleOut << "ID_1 ID_2 missing" << endl;
    sampleOut << "0 0 0" << endl;
    for (uint64 n0 = 0; n0 < N; n0++) {
      /*
      int ctrMiss = 0, ctrTot = 0;
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
	if (maskSnps64j[m64j]) {
	  ctrTot++;
	  ctrMiss += getGeno0123(m64j, n0)==3;
	}
      */
      double missing = 0;//ctrMiss / (double) ctrTot;
      sampleOut << indivs[n0].famID << " " << indivs[n0].indivID << " " << missing << endl;
    }
    sampleOut.close();
  }

  void bcf_hdr_append_eagle_version(bcf_hdr_t *hdr, int argc, char **argv)
  {
    kstring_t str = {0,0,0};
    const char cmd[] = "eagle";
    ksprintf(&str,"##%sVersion=%s+htslib-%s\n", cmd, EAGLE_VERSION, hts_version());
    bcf_hdr_append(hdr,str.s);

    str.l = 0;
    ksprintf(&str,"##%sCommand=%s", cmd, "eagle");
    int i;
    for (i=1; i<argc; i++)
      {
        if ( strchr(argv[i],' ') )
	  ksprintf(&str, " '%s'", argv[i]);
        else
	  ksprintf(&str, " %s", argv[i]);
      }
    kputc('\n', &str);
    bcf_hdr_append(hdr,str.s);
    free(str.s);

    bcf_hdr_sync(hdr);
  }

  void Eagle::writeVcf(const string &tmpFile, const vector <bool> &isTmpPhased,
		       const string &outFile, int chromX, double bpStart, double bpEnd,
		       const string &writeMode, bool noImpMissing, bool keepMissingPloidyX,
		       int argc, char **argv) const {

    htsFile *htsTmp = hts_open(tmpFile.c_str(), "r");
    if (htsTmp == NULL) {
      cerr << "ERROR: Could not open temporary file " << tmpFile << endl;
      exit(1);
    }
    htsFile *out = hts_open(outFile.c_str(), writeMode.c_str());
    if (out == NULL) {
      cerr << "ERROR: Could not write to file " << outFile << endl;
      exit(1);
    }
    // Get a threadpool for writing equal to
    // the maximum number of threads that have been configured.
    htsThreadPool p = {hts_tpool_init(omp_get_max_threads()),
                       0};
    hts_set_thread_pool(out, &p);
    hts_set_thread_pool(htsTmp, &p);

    bcf_hdr_t *hdr = bcf_hdr_read(htsTmp);
    bcf_hdr_append_eagle_version(hdr, argc, argv);
    bcf_hdr_write(out, hdr);

    bcf1_t *rec = bcf_init1();
    int mtgt_gt = 0, *tgt_gt = NULL;

    uint64 m64j = 0; // SNP index; update to correspond to current record

    vector <int> mostRecentPloidy(N-Nref, 2);

    int tmpLineNum = -1; // index in tmp file (corresponding to isTmpPhased)
    while (bcf_read(htsTmp, hdr, rec) >= 0) {
      tmpLineNum++;
      if (!isTmpPhased[tmpLineNum]) { // site was not phased; remove phase information and output
	int bp = rec->pos+1;
	if (bpStart <= bp && bp <= bpEnd) { // check if within output region
	  int ntgt_gt = bcf_get_genotypes(hdr, rec, &tgt_gt, &mtgt_gt);
	  for (int k = 0; k < ntgt_gt; k++)
	    if (tgt_gt[k] != bcf_int32_vector_end && !bcf_gt_is_missing(tgt_gt[k])) {
	      int idx = bcf_gt_allele(tgt_gt[k]); // allele index
	      tgt_gt[k] = bcf_gt_unphased(idx); // convert allele index to bcf value (unphased)
	    }
	  bcf_update_genotypes(hdr, rec, tgt_gt, ntgt_gt);
	  bcf_write(out, hdr, rec);
	}
	continue;
      }

      int chrom = StringUtils::bcfNameToChrom(bcf_hdr_id2name(hdr, rec->rid), 1, chromX);

      int ntgt_gt = bcf_get_genotypes(hdr, rec, &tgt_gt, &mtgt_gt);

      int bp = rec->pos+1;
      if (bpStart <= bp && bp <= bpEnd) { // check if within output region
	if (bcf_hdr_nsamples(hdr) == ntgt_gt) { // site is all-haploid
	  bcf_write(out, hdr, rec);
	}
	else {
	  for (int i = 0; i < (int) (N-Nref); i++) {
	    int ploidy = 2;
	    int *ptr = tgt_gt + i*ploidy;

	    // --keepMissingPloidyX: assume missing genotypes in target VCF have correct ploidy
	    if (chrom == chromX && keepMissingPloidyX && bcf_gt_is_missing(ptr[0]))
	      mostRecentPloidy[i] = (ptr[1] == bcf_int32_vector_end ? 1 : 2);

	    if (chrom != chromX || (bcf_gt_is_missing(ptr[0]) && mostRecentPloidy[i] == 2)
		|| ptr[1] != bcf_int32_vector_end) { // diploid... be careful about missing '.'
	      mostRecentPloidy[i] = 2;
	      bool missing = false;
	      int minIdx = 1000, maxIdx = 0;
	      for (int j = 0; j < ploidy; j++) {
		if ( bcf_gt_is_missing(ptr[j]) ) { // missing allele
		  missing = true;
		}
		else {
		  int idx = bcf_gt_allele(ptr[j]); // allele index
		  minIdx = std::min(minIdx, idx);
		  maxIdx = std::max(maxIdx, idx);
		}
	      }

	      if (!missing && minIdx == maxIdx) { // hom => same allele
		ptr[0] = ptr[1] = bcf_gt_phased(minIdx);
	      }
	      else if (!missing && minIdx > 0) { // ALT1/ALT2 het => don't phase
		ptr[0] = ptr[1] = bcf_gt_missing;
	      }
	      else { // REF/ALT* het => phase as called by Eagle
		if (missing && noImpMissing) { // don't call alleles
		  ptr[0] = ptr[1] = bcf_gt_missing;
		}
		else {
		  for (int j = 0; j < ploidy; j++) {
		    uint64 nTargetHap = 2*i + j;
		    int altIdx = missing ? 1 : maxIdx;
		    int hapBit = (tmpHaploBitsT[nTargetHap*Mseg64+(m64j/64)]>>(m64j&63))&1;
		    if (isFlipped64j[m64j]) hapBit = !hapBit;
		    int idx = hapBit ? altIdx : 0;
		    ptr[j] = bcf_gt_phased(idx); // convert allele index to bcf value (phased)
		  }
		}
	      }
	    }
	    else { // haploid
	      mostRecentPloidy[i] = 1;
	      if ( bcf_gt_is_missing(ptr[0]) && !noImpMissing ) { // missing allele
		int j = 0;
		uint64 nTargetHap = 2*i + j;
		int altIdx = 1;
		int hapBit = (tmpHaploBitsT[nTargetHap*Mseg64+(m64j/64)]>>(m64j&63))&1;
		if (isFlipped64j[m64j]) hapBit = !hapBit;
		int idx = hapBit ? altIdx : 0;
		ptr[j] = bcf_gt_phased(idx); // convert allele index to bcf value (phased)
	      }
	    }
	  }

	  bcf_update_genotypes(hdr, rec, tgt_gt, ntgt_gt);

	  bcf_write(out, hdr, rec);
	}
      }

      m64j++;
      if (m64j != Mseg64*64 && (m64j&63) == seg64cMvecs[m64j/64].size())
	m64j = (m64j + 64ULL) & ~63ULL; // move to next segment
    }

    assert(m64j == Mseg64*64);

    free(tgt_gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(out);
    hts_close(htsTmp);
    hts_tpool_destroy(p.pool);
    remove(tmpFile.c_str());
  }

  // write phased output in non-ref mode
  // differences from the above (ref-mode) are as follows:
  // - checks chrom
  // - does not increment m64j outside the bpStart-bpEnd region
  // - does not delete tmpFile (now vcfFile with original input)
  void Eagle::writeVcfNonRef(const string &vcfFile, const string &outFile, int inputChrom,
			     int chromX, double bpStart, double bpEnd, const string &writeMode,
			     bool noImpMissing, bool keepMissingPloidyX, int argc, char **argv)
    const {

    htsFile *htsIn = hts_open(vcfFile.c_str(), "r");
    htsFile *out = hts_open(outFile.c_str(), writeMode.c_str());
    // Get a threadpool for writing equal to
    // the maximum number of threads that have been configured.
    htsThreadPool p = {hts_tpool_init(omp_get_max_threads()),
                       0};
    hts_set_thread_pool(out, &p);
    hts_set_thread_pool(htsIn, &p);

    bcf_hdr_t *hdr = bcf_hdr_read(htsIn);
    bcf_hdr_append_eagle_version(hdr, argc, argv);
    bcf_hdr_write(out, hdr);

    bcf1_t *rec = bcf_init1();
    int mtgt_gt = 0, *tgt_gt = NULL;

    uint64 m64j = 0; // SNP index; update to correspond to current record

    vector <int> mostRecentPloidy(N-Nref, 2);

    while (bcf_read(htsIn, hdr, rec) >= 0) {
      // check CHROM
      int chrom = StringUtils::bcfNameToChrom(bcf_hdr_id2name(hdr, rec->rid), 1, chromX);
      if (inputChrom != 0) {
	if (chrom != inputChrom) {
	  continue;
	}
      }

      int bp = rec->pos+1;
      if (bpStart <= bp && bp <= bpEnd) { // check if within output region
	int ntgt_gt = bcf_get_genotypes(hdr, rec, &tgt_gt, &mtgt_gt);

	for (int i = 0; i < (int) (N-Nref); i++) {
	  int ploidy = 2;
	  int *ptr = tgt_gt + i*ploidy;

	  // --keepMissingPloidyX: assume missing genotypes in target VCF have correct ploidy
	  if (chrom == chromX && keepMissingPloidyX && bcf_gt_is_missing(ptr[0]))
	    mostRecentPloidy[i] = (ptr[1] == bcf_int32_vector_end ? 1 : 2);

	  if (chrom != chromX || (bcf_gt_is_missing(ptr[0]) && mostRecentPloidy[i] == 2)
	      || ptr[1] != bcf_int32_vector_end) { // diploid... be careful about missing '.'
	    mostRecentPloidy[i] = 2;
	    bool missing = false;
	    int minIdx = 1000, maxIdx = 0; // (shouldn't matter; SNPs should be biallelic)
	    for (int j = 0; j < ploidy; j++) {
	      if ( bcf_gt_is_missing(ptr[j]) ) { // missing allele
		missing = true;
	      }
	      else {
		int idx = bcf_gt_allele(ptr[j]); // allele index
		minIdx = std::min(minIdx, idx);
		maxIdx = std::max(maxIdx, idx);
	      }
	    }

	    if (!missing && minIdx == maxIdx) { // hom => same allele
	      ptr[0] = ptr[1] = bcf_gt_phased(minIdx);
	    }
	    else if (!missing && minIdx > 0) { // ALT1/ALT2 het => don't phase (shouldn't happen)
	      ptr[0] = ptr[1] = bcf_gt_missing;
	    }
	    else if (!missing || !noImpMissing) { // miss||REF/ALT* het => phase as called by Eagle
	      for (int j = 0; j < ploidy; j++) {
		uint64 nTargetHap = 2*i + j;
		int altIdx = missing ? 1 : maxIdx;
		int hapBit = (tmpHaploBitsT[nTargetHap*Mseg64+(m64j/64)]>>(m64j&63))&1;
		if (isFlipped64j[m64j]) hapBit = !hapBit;
		int idx = hapBit ? altIdx : 0;
		ptr[j] = bcf_gt_phased(idx); // convert allele index to bcf value (phased)
	      }
	    }
	  }
	  else { // haploid
	    mostRecentPloidy[i] = 1;
	    if ( bcf_gt_is_missing(ptr[0]) && !noImpMissing ) { // missing allele
	      int j = 0;
	      uint64 nTargetHap = 2*i + j;
	      int altIdx = 1;
	      int hapBit = (tmpHaploBitsT[nTargetHap*Mseg64+(m64j/64)]>>(m64j&63))&1;
	      if (isFlipped64j[m64j]) hapBit = !hapBit;
	      int idx = hapBit ? altIdx : 0;
	      ptr[j] = bcf_gt_phased(idx); // convert allele index to bcf value (phased)
	    }
	  }
	}

	bcf_update_genotypes(hdr, rec, tgt_gt, ntgt_gt);

	bcf_write(out, hdr, rec);

	m64j++;
	if ((m64j&63) == seg64cMvecs[m64j/64].size())
	  m64j = (m64j + 64ULL) & ~63ULL; // move to next segment
      }
    }

    assert(m64j == Mseg64*64);

    free(tgt_gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(out);
    hts_close(htsIn);
    hts_tpool_destroy(p.pool);
  }

  void Eagle::makeHardCalls(uint64 n0start, uint64 n0end, uint seed) {
    // fast rng: last 16 bits of Marsaglia's MWC
    uint w = 521288629;
    for (uint i = 0; i < seed % 12345; i++)
      w=18000*(w&65535)+(w>>16);
    //memset(haploBits, 0, 2*N*Mseg64*sizeof(haploBits[0]));
    memset(segConfs + 2*n0start*Mseg64, 0, 2*(n0end-n0start)*Mseg64*sizeof(segConfs[0]));
    for (uint64 nHap = 2*n0start; nHap < 2*n0end; nHap++)
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++) {
	if ((m64j&63)==0)
	  haploBits[(m64j/64)*2*N + nHap] = 0;
	uchar phaseConf = phaseConfs[nHap*Mseg64*64 + m64j];
	segConfs[nHap*Mseg64+m64j/64] = max(segConfs[nHap*Mseg64+m64j/64],
					    min(phaseConf, (uchar) (255-phaseConf)));
	if (phaseConf == (uchar) 255 || ((w=18000*(w&65535)+(w>>16))&255) < phaseConf)
	  haploBits[(m64j/64)*2*N + nHap] |= 1ULL<<(m64j&63);
	if ((m64j&63)==63)
	  haploBitsT[nHap*Mseg64 + (m64j/64)] = haploBits[(m64j/64)*2*N + nHap];
      }
  }

  uint Eagle::computeHash(const uint64 *curHaploBitsT, const uint64 *curHashBits, uint B) const {
    uint hash = 0;
    for (uint b = 0; b < B; b++)
      hash |= ((curHaploBitsT[curHashBits[b]>>6]>>(curHashBits[b]&63))&1)<<b;
    return hash;
  }

  uint Eagle::computeHash(const uint64 *curHaploBitsT, const vector <uint64> &curHashBits) const {
    return computeHash(curHaploBitsT, &curHashBits[0], curHashBits.size());
  }

  double Eagle::computeLogHetP(uint64 m64j) const {
    assert(Nref!=0); // only call this function when in ref-mode
    int sumHaps = 0;
    for (uint64 nHap = 0; nHap < 2*Nref; nHap++)
      sumHaps += (haploBits[(m64j/64)*2*N + nHap]>>(m64j&63))&1;
    double p = sumHaps / (2.0 * Nref);
    p = std::min(p, 1-p);
    return log10(p);
  }

  void Eagle::initRefIter(int refIter) {
    uint64 Ntarget = N - Nref;
    if (refIter > 1) { // copy tmpHaploBitsT from previous iter -> haploBits, haploBitsT
      memcpy(haploBitsT + 2*Nref*Mseg64, tmpHaploBitsT, 2*Ntarget*Mseg64*sizeof(tmpHaploBitsT[0]));
      for (uint64 nHap = 2*Nref; nHap < 2*N; nHap++) // copy transpose
	for (uint64 m64 = 0; m64 < Mseg64; m64++)
	  haploBits[m64*2*N + nHap] = haploBitsT[nHap*Mseg64 + m64];
    }
    // clear tmpHaploBitsT (temp storage of phased target haplotypes)
    memset(tmpHaploBitsT, 0, 2*Ntarget*Mseg64*sizeof(tmpHaploBitsT[0]));
  }

  // input arg iter = non-ref mode iter (ref mode iter 1 = non-ref mode iter 3)
  void Eagle::buildHashTables(int iter, int batch, int seed) {

    std::srand(1000000*seed + 1000*iter + batch); // seed random_shuffle

    const uint maxValuesPerKey = 99;
    const uint baseLSH = 10, bonusLSH = iter > 2 ? 4 : 0;
    const uint numLSH = baseLSH + bonusLSH;
    const uint maxBits = 32;
    const double minLogHetP = log10(0.02);
    hashLookups = vector < vector <StaticMultimap> > (Mseg64, vector <StaticMultimap> (numLSH));
    hashBits = vector < vector < vector <uint64> > > (Mseg64, vector < vector <uint64> > (numLSH));

    const double reduction = (iter == 2 ? 0 : 0.05);

    const uint64 side = 1;
    for (uint64 m64 = 0+side; m64+side < Mseg64; m64++) {
      for (uint h = 0; h < numLSH; h++) {
	vector <uint64> m64js;

	if (h < baseLSH) { // standard hash regions: 3x down to 2x m64
	  for (uint64 m64j = (uint64) ((m64-side+reduction*h)*64);
	       m64j < (uint64) ((m64+side+1-reduction*h)*64); m64j++)
	    if (maskSnps64j[m64j] &&
		(Nref==0 ? seg64logPs[m64j].cond[1][3] : computeLogHetP(m64j)) > minLogHetP)
	      m64js.push_back(m64j);
	}
	else { // small hash regions
	  int offStart = 0, offEnd = 0;
	  switch (h-baseLSH) {
	  case 0: offStart = -32; offEnd = 32; break;
	  case 1: offStart = 0; offEnd = 64; break;
	  case 2: offStart = -32; offEnd = 0; break;
	  case 3: offStart = 0; offEnd = 32; break;
	  }
	  for (uint64 m64j = (uint64) (m64*64 + offStart); m64j < (uint64) (m64*64 + offEnd);
	       m64j++)
	    if (maskSnps64j[m64j])
	      m64js.push_back(m64j);
	}

	if (m64js.empty())
	  for (uint64 m64j = (m64-side)*64; m64j < (m64+side+1)*64; m64j++)
	    if (maskSnps64j[m64j])
	      m64js.push_back(m64j);

	uint bitsInHash = (h < baseLSH ? maxBits-h : 24);

	// randomly select SNPs m64j to use in hash
	uint m64jInd = m64js.size();
	for (uint b = 0; b < bitsInHash; b++) {
	  // choose next SNP (in random order); if at end, reshuffle
	  if (m64jInd == m64js.size()) {
	    std::random_shuffle(m64js.begin(), m64js.end());
	    m64jInd = 0;
	  }
	  hashBits[m64][h].push_back(m64js[m64jInd++]);
	}
      }
    }

    uint64 nRefHaps = 2*((Nref!=0 && iter==3) ? Nref : N); // ref-mode iter 1 -> iter 3
    vector < vector <uint> > keyVecs(omp_get_max_threads(), vector <uint> (nRefHaps));
#pragma omp parallel for
    for (uint64 m64 = 0+side; m64 < Mseg64-side; m64++) {
      cout << "." << std::flush;
      for (uint h = 0; h < numLSH; h++) {
	// compute hashes
	vector <uint> &keyVec = keyVecs[omp_get_thread_num()];
	for (uint64 nHap = 0; nHap < nRefHaps; nHap++) // in ref-mode, only use ref
	  keyVec[nHap] = maskIndivs[nHap/2] ?
	    computeHash(haploBitsT + nHap*Mseg64, &hashBits[m64][h][0], hashBits[m64][h].size())
	    : -1U;
	hashLookups[m64][h].init(keyVec, maxValuesPerKey);
      }
    }
  }

  const uint64 *Eagle::getHaploBitsT(void) const { return haploBitsT; }
  uint64 Eagle::getNlib(int iter) const { return ((Nref!=0 && iter==3) ? Nref : N); }
  uint64 Eagle::getMseg64(void) const { return Mseg64; }
  const uchar *Eagle::getMaskSnps64j(void) const { return maskSnps64j; }

  double Eagle::computeHetRate(void) const {
    uint64 homCtr = 0, totCtr = 0;
    for (uint64 m64 = 0; m64 < Mseg64; m64++)
      for (uint64 n = Nref; n < N; n++) { // ref genoBits aren't initialized!
	const uint64_masks &bits = genoBits[m64*N + n];
	homCtr += popcount64(bits.is0 | bits.is2);
	totCtr += popcount64(~bits.is9);
      }
    return 1 - homCtr / (double) totCtr;
  }

}
