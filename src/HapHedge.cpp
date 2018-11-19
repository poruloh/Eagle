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
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "MemoryUtils.hpp"
#include "Types.hpp"
#include "HapHedge.hpp"

namespace EAGLE {

  using std::vector;
  using std::string;
  using std::cout;
  using std::endl;

  void HapBitsT::setBit(uint64 n, uint64 m) { haploBitsT[n*M64 + (m>>6)] |= 1ULL<<(m&63); }

#ifdef TEST_VS
  HapBitsT::HapBitsT(const vector <string> &_hapBitsVS)
    : hapBitsVS(_hapBitsVS), Nhaps(hapBitsVS.size()), M(hapBitsVS[0].length()) {}
  int HapBitsT::getBit(uint64 n, uint64 m) const { return hapBitsVS[n][m]-'0'; }
  HapBitsT::~HapBitsT(void) { }
#else
  HapBitsT::HapBitsT(const uint64 *_haploBitsT, uint64 _Nhaps, uint64 Mseg64,
		     const uchar maskSnps64j[]) {
    Nhaps = _Nhaps;
    M = 0;
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
      if (maskSnps64j[m64j])
	M++;
    M64 = (M+63)/64;
    haploBitsT = ALIGNED_MALLOC_UINT64S(Nhaps * M64);
    memset(haploBitsT, 0, Nhaps * M64 * sizeof(haploBitsT[0]));
    for (uint64 n = 0; n < Nhaps; n++) {
      uint64 mCur = 0;
      for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
	if (maskSnps64j[m64j]) {
	  if ((_haploBitsT[n*Mseg64 + (m64j>>6)]>>(m64j&63))&1)
	    setBit(n, mCur);
	  mCur++;
	}
    }
  }
  inline int popcount64_01(uint64 i) {
    return i!=0;
  }
  inline int popcount64_012(uint64 i) {
    if (i == 0) return 0;
    else if ((i & (i-1ULL)) == 0) return 1;
    else return 2;
  }
  HapBitsT::HapBitsT(const uint64 *inBitsT, uint64 Mseg64, const vector <uint64> &splits64j,
		     const vector <uchar> &splitGenos, const vector <uint64_masks> &tgtGenoBits,
		     const vector <uint> &bestHaps) {
    Nhaps = bestHaps.size();
    M = 2*(splits64j.size()+2);
    M64 = (M+63)/64;
    haploBitsT = ALIGNED_MALLOC_UINT64S(Nhaps * M64);
    memset(haploBitsT, 0, Nhaps * M64 * sizeof(haploBitsT[0]));

    for (uint k = 0; k < bestHaps.size(); k++) {
      uint64 n = bestHaps[k];
      if (k+1<bestHaps.size() && bestHaps[k+1]==n+1) {
	// process both haplotypes from this individual
	uint64 hLastDiff = 0; // for double-IBD
	uint64 hLastErrA = 0, hLastErrB = 0;
	for (uint64 h = 0; h <= splits64j.size(); h++) {
	  // set hom before splits64j[h]
	  int numErrsA = 0, numErrsB = 0;
	  uint64 homStart = (h == 0 ? 0 : splits64j[h-1]+1);
	  uint64 homStop = (h == splits64j.size() ? Mseg64*64 : std::max(splits64j[h], 1ULL)) - 1;
	  for (uint64 m64 = (homStart>>6); m64 <= (homStop>>6); m64++) {
	    uint64 mask = -1ULL;
	    if (m64 == (homStart>>6))
	      mask &= (-1ULL>>(homStart&63))<<(homStart&63);
	    if (m64 == (homStop>>6))
	      mask &= (-1ULL>>(63-(homStop&63)));

	    numErrsA += popcount64_01(mask&((tgtGenoBits[m64].is0 & inBitsT[n*Mseg64 + m64]) |
					    (tgtGenoBits[m64].is2 & ~inBitsT[n*Mseg64 + m64])));
	    numErrsB += popcount64_01(mask&((tgtGenoBits[m64].is0 & inBitsT[(n+1)*Mseg64 + m64]) |
					    (tgtGenoBits[m64].is2 & ~inBitsT[(n+1)*Mseg64 + m64])));
	    if (numErrsA && numErrsB)
	      break;
	  }

	  if (numErrsA)
	    setBit(k, 2*h+1);
	  if (numErrsB)
	    setBit(k+1, 2*h+1);

	  // set bits at split
	  int splitA = 0, splitB = 0, splitGeno = 0;
	  if (h < splits64j.size()) {
	    splitA = (inBitsT[n*Mseg64 + (splits64j[h]>>6)]>>(splits64j[h]&63))&1;
	    splitB = (inBitsT[(n+1)*Mseg64 + (splits64j[h]>>6)]>>(splits64j[h]&63))&1;
	    splitGeno = splitGenos[h];
	    if (splitA) setBit(k, 2*(h+1));
	    if (splitB) setBit(k+1, 2*(h+1));
	  }

	  // check for double-IBD
	  if (numErrsA || numErrsB || splitA+splitB != splitGeno || h==splits64j.size()) {
	    if (h > hLastDiff+20) {
	      for (uint h2 = hLastDiff+1; h2 < h; h2++) {
		setBit(k, 2*h2+1);
		setBit(k+1, 2*h2+1);
	      }
	    }
	    else if (h > hLastDiff+10) {
	      uint64 kMask = hLastErrA >= hLastErrB ? k : k+1;
	      for (uint h2 = hLastDiff+1; h2 < h; h2++)
		setBit(kMask, 2*h2+1);
	    }
	    hLastDiff = h;
	  }
	  if (numErrsA) hLastErrA = h;
	  if (numErrsB) hLastErrB = h;
	}
	k++;
      }
      else { // process lone reference haplotype
	for (uint64 h = 0; h <= splits64j.size(); h++) {
	  // set hom
	  uint64 homStart = (h == 0 ? 0 : splits64j[h-1]+1);
	  uint64 homStop = (h == splits64j.size() ? Mseg64*64 : std::max(splits64j[h], 1ULL)) - 1;
	  int numErrs = 0;
	  for (uint64 m64 = (homStart>>6); m64 <= (homStop>>6); m64++) {
	    uint64 mask = -1ULL;
	    if (m64 == (homStart>>6))
	      mask &= (-1ULL>>(homStart&63))<<(homStart&63);
	    if (m64 == (homStop>>6))
	      mask &= (-1ULL>>(63-(homStop&63)));
	    numErrs += popcount64_01(mask&((tgtGenoBits[m64].is0 & inBitsT[n*Mseg64 + m64]) |
					   (tgtGenoBits[m64].is2 & ~inBitsT[n*Mseg64 + m64])));
	    if (numErrs)
	      break;
	  }
	  if (numErrs)
	    setBit(k, 2*h+1);

	  // set bit at split
	  if (h < splits64j.size() && ((inBitsT[n*Mseg64+(splits64j[h]>>6)]>>(splits64j[h]&63))&1))
	    setBit(k, 2*(h+1));
	}
      }
    }
  }
  const unsigned char flipByte[] = {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
    0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
    0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
    0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
    0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
    0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
    0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
    0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
    0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
    0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
    0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
    0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
    0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
    0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
  };
  uint64 reverse64(uint64 v) {
    return
      ((uint64) flipByte[v & 0xff] << 56ULL) |
      ((uint64) flipByte[(v >> 8) & 0xff] << 48ULL) |
      ((uint64) flipByte[(v >> 16) & 0xff] << 40ULL) |
      ((uint64) flipByte[(v >> 24) & 0xff] << 32ULL) |
      ((uint64) flipByte[(v >> 32) & 0xff] << 24ULL) |
      ((uint64) flipByte[(v >> 40) & 0xff] << 16ULL) |
      ((uint64) flipByte[(v >> 48) & 0xff] << 8ULL) |
      flipByte[v >> 56];
  }
  HapBitsT::HapBitsT(const HapBitsT &hapBitsFwdT, int dir)
    : Nhaps(hapBitsFwdT.getNhaps()), M(hapBitsFwdT.getM()), M64((M+63)/64) {
    
    assert(dir<0); // reverse existing HapBitsT

    haploBitsT = ALIGNED_MALLOC_UINT64S(Nhaps * M64);
    memset(haploBitsT, 0, Nhaps * M64 * sizeof(haploBitsT[0]));

    if (dir == -1) { // simple reversal
      for (uint64 n = 0; n < Nhaps; n++)
	for (uint64 m = 0; m < M; m++)
	  if (hapBitsFwdT.getBit(n, m))
	    setBit(n, M-1-m);
    }
    else { // repad 0 xyz 0 0 -> 0 zyx 0 0
      uint64 Mtrim = M-1, Mtrim64 = (Mtrim+63)/64, offset = Mtrim&63;
      for (uint64 n = 0; n < Nhaps; n++) {
	/*
	for (uint64 m = 1; m < M-2; m++)
	  if (hapBitsFwdT.getBit(n, m))
	    setBit(n, M-2-m);
	*/
	for (uint64 m64 = 0; m64 < Mtrim64; m64++) {
	  uint64 fwdBits = offset==0 ? hapBitsFwdT.getBits64(n, Mtrim64-1-m64) :
	    ((m64<Mtrim64-1 ? (hapBitsFwdT.getBits64(n, Mtrim64-2-m64)>>offset) : 0) |
	     (hapBitsFwdT.getBits64(n, Mtrim64-1-m64)<<(64-offset)));
	  //assert(haploBitsT[n * M64 + m64] == reverse64(fwdBits));
	  haploBitsT[n * M64 + m64] = reverse64(fwdBits);
	}
      }
    }
  }
  int HapBitsT::getBit(uint64 n, uint64 m) const { return (haploBitsT[n*M64 + (m>>6)]>>(m&63))&1; }
  uint64 HapBitsT::getBits64(uint64 n, uint64 m64) const { return haploBitsT[n*M64 + m64]; }
  HapBitsT::~HapBitsT(void) { ALIGNED_FREE(haploBitsT); }
#endif
  int HapBitsT::getNhaps(void) const { return Nhaps; }
  int HapBitsT::getM(void) const { return M; }


  int d, a, count, up, left, right;
  WorkTreeNode::WorkTreeNode(int _d, int _a, int _count, int _up, int _left, int _right) :
    d(_d), a(_a), count(_count), up(_up), left(_left), right(_right) {}


  void dfsPreOrder(WorkTreeNode workNodes[], HapTreeNode *nodes, int &pos, int cur) {
    if (cur == -1) return;
    nodes[pos].mSplit = workNodes[cur].d;
    int left = workNodes[cur].left;
    nodes[pos].count0 = 1 + (left==-1 ? 0 : workNodes[left].count);
    nodes[pos].seq1 = workNodes[cur].a;
    pos++;
    dfsPreOrder(workNodes, nodes, pos, workNodes[cur].left);
    dfsPreOrder(workNodes, nodes, pos, workNodes[cur].right);
  }

  void dfsPreOrderMulti(WorkTreeNode workNodes[], HapTreeMultiNode *nodes, int &pos, int cur) {
    if (cur == -1) return;
    int curPos = pos;
    nodes[curPos].mSplit = workNodes[cur].d;
    nodes[curPos].count0 = workNodes[cur].count;
    nodes[curPos].node0 = curPos+1;
    nodes[curPos].seq1 = workNodes[cur].a;
    pos++;
    dfsPreOrderMulti(workNodes, nodes, pos, workNodes[cur].left);
    if (pos == curPos+1) { // no next node 0
      nodes[curPos].node0 = -1;
      /* MEMORY-SAVING ALTERNATIVE: use 1 bit (sign of state.count) to encode next node
      nodes[curPos].count0 = -nodes[curPos].count0;
      */
    }
    nodes[curPos].node1 = pos;
    dfsPreOrderMulti(workNodes, nodes, pos, workNodes[cur].right);
    if (pos == nodes[curPos].node1) // no next node 1
      nodes[curPos].node1 = -1;
  }

  HapTree::HapTree(const HapBitsT &_hapBitsT, int a[], int d[])
    : hapBitsT(_hapBitsT), Nhaps(hapBitsT.getNhaps()), invNhaps(1.0f/Nhaps) {

    seq0 = a[0];
    nodes = (HapTreeNode *) ALIGNED_MALLOC((Nhaps-1) * sizeof(nodes[0]));

    WorkTreeNode *workNodes = (WorkTreeNode *) ALIGNED_MALLOC(Nhaps * sizeof(workNodes[0]));

    // perform in-order traversal
    workNodes[0] = WorkTreeNode(-1, a[0], 0, -1, -1, -1);
    for (int n = 1; n < Nhaps; n++) {
      if (d[n] >= d[n-1]) {
	workNodes[n] = WorkTreeNode(d[n], a[n], 1, n-1, -1, -1);
	workNodes[n-1].right = n;
      }
      else {
	int cur = n-1, up = workNodes[cur].up;
	while (d[n] < workNodes[up].d) {
	  workNodes[up].count += workNodes[cur].count;
	  cur = up; up = workNodes[cur].up;
	}
	workNodes[n] = WorkTreeNode(d[n], a[n], workNodes[cur].count+1, up, cur, -1);
	workNodes[cur].up = n;
	workNodes[up].right = n;
      }
    }
    int cur = Nhaps-1, up = workNodes[cur].up;
    while (up != -1) {
      workNodes[up].count += workNodes[cur].count;
      cur = up; up = workNodes[cur].up;
    }

    // perform pre-order traversal
    int pos = 0;
    dfsPreOrder(workNodes, nodes, pos, workNodes[0].right);

    ALIGNED_FREE(workNodes);
  }
  HapTree::~HapTree(void) {
    ALIGNED_FREE(nodes);
  }
  float HapTree::getInvNhaps(void) const {
    return invNhaps;
  }
  HapTreeState HapTree::getRootState(void) const {
    HapTreeState h;
    h.seq = seq0; h.node = 0; h.count = Nhaps;
    return h;
  }
  bool HapTree::next(int m, HapTreeState &state, int nextBit) const {
    if (state.count == 1 || nodes[state.node].mSplit > m) {
      return hapBitsT.getBit(state.seq, m) == nextBit;
    }
    else {
      if (nextBit == 0) {
	state.count = nodes[state.node++].count0;
      }
      else {
	state.seq = nodes[state.node].seq1;
	int c0 = nodes[state.node].count0;
	state.node += c0;
	state.count -= c0;
      }
      return true;
    }
  }

  HapTreeMulti::HapTreeMulti(const HapBitsT &_hapBitsT, SortDiv ad[], int M,
			     WorkTreeNode workNodes[])
    : hapBitsT(_hapBitsT), Nhaps(hapBitsT.getNhaps()), invNhaps(1.0f/Nhaps) {

    // perform in-order traversal
    workNodes[0] = WorkTreeNode(-1, ad[0].a, 0, -1, -1, -1);
    int n = 0, nextMult = 1, uniqHaps = Nhaps;
    for (int nHap = 1; nHap < Nhaps; nHap++) {
      if (ad[nHap].d == M) { // sequence is identical to previous; merge
	nextMult++;
	uniqHaps--;
      }
      else {
	n++;
	if (ad[nHap].d >= workNodes[n-1].d) {
	  workNodes[n] = WorkTreeNode(ad[nHap].d, ad[nHap].a, nextMult, n-1, -1, -1);
	  workNodes[n-1].right = n;
	}
	else {
	  int cur = n-1, up = workNodes[cur].up;
	  int leftCount = nextMult;
	  while (ad[nHap].d < workNodes[up].d) {
	    leftCount += workNodes[cur].count;
	    cur = up; up = workNodes[cur].up;
	  }
	  leftCount += workNodes[cur].count;
	  workNodes[n] = WorkTreeNode(ad[nHap].d, ad[nHap].a, leftCount, up, cur, -1);
	  workNodes[cur].up = n;
	  workNodes[up].right = n;
	}
	nextMult = 1;
      }
    }

    nodes = (HapTreeMultiNode *) ALIGNED_MALLOC((uniqHaps-1) * sizeof(nodes[0]));

    rootState.seq = ad[0].a; rootState.node = workNodes[0].right==-1?-1:0; rootState.count = Nhaps;
    int pos = 0;
    dfsPreOrderMulti(workNodes, nodes, pos, workNodes[0].right);
  }
  HapTreeMulti::~HapTreeMulti(void) {
    ALIGNED_FREE(nodes);
  }
  float HapTreeMulti::getInvNhaps(void) const {
    return invNhaps;
  }
  HapTreeState HapTreeMulti::getRootState(void) const {
    return rootState;
  }
  bool HapTreeMulti::next(int m, HapTreeState &state, int nextBit) const {
    if (state.node == -1 || nodes[state.node].mSplit > m) {
      return hapBitsT.getBit(state.seq, m) == nextBit;
    }
    else {
      if (nextBit == 0) {
	state.count = nodes[state.node].count0;
	state.node = nodes[state.node].node0;
	/* MEMORY-SAVING ALTERNATIVE: use 1 bit (sign of state.count) to encode next node
	if (state.count > 0)
	  state.node++;
	else {
	  state.node = -1;
	  state.count = -state.count;
	}
	*/
      }
      else {
	state.count -= nodes[state.node].count0;
	/* MEMORY-SAVING ALTERNATIVE: use 1 bit (sign of state.count) to encode next node
	state.count -= abs(nodes[state.node].count0);
	*/
	state.seq = nodes[state.node].seq1;
	state.node = nodes[state.node].node1; // update last: overwrites state.node!
      }
      return true;
    }
  }
  void HapTreeMulti::nextAtFrac(int m, HapTreeState &state, double nextFrac) const {
    // see above for MEMORY-SAVING ALTERNATIVES
    if (state.node == -1 || nodes[state.node].mSplit > m) {
      return;
    }
    else {
      if (nodes[state.node].count0 >= nextFrac * state.count) { // nextBit = 0
	state.count = nodes[state.node].count0;
	state.node = nodes[state.node].node0;
      }
      else {
	state.count -= nodes[state.node].count0;
	state.seq = nodes[state.node].seq1;
	state.node = nodes[state.node].node1; // update last: overwrites state.node!
      }
    }
  }
  // for debugging
  void HapTreeMulti::dfsPrint(string curPrefix, int m, int M, const HapTreeState &state) const {
    if (m > M) return;
    cout << "m = " << m << ", prefix = " << curPrefix << ": count = " << state.count << endl;
    for (int b = 0; b < 2; b++) {
      HapTreeState nextState = state;
      if (next(m, nextState, b))
	dfsPrint(curPrefix + (char) ('0'+b), m+1, M, nextState);
    }
  }

  HapHedge::HapHedge(const HapBitsT &_hapBitsT, int _skip/*, const vector <int> &treeStarts*/) :
    hapBitsT(_hapBitsT), skip(_skip), T((hapBitsT.getM()+skip-1) / skip) {
      
    treePtrs = new HapTree *[T];

    int N = hapBitsT.getNhaps(), M = hapBitsT.getM();

    // initialize work arrays
    int *a1 = new int[N], *d1 = new int[N], *a = new int[N], *b = new int[N], *d = new int[N],
      *e = new int[N];
    for (int n = 0; n < N; n++) {
      a1[n] = n;
      d1[n] = M;
    }

    for (int m = M-1; m >= 0; m--) {
      // compute sort order and divergence array
      int u = 0, v = 0, p = m, q = m;
      for (int n = 0; n < N; n++) {
	if (d1[n] < p) p = d1[n];
	if (d1[n] < q) q = d1[n];
	if (hapBitsT.getBit(a1[n], m) == 0) {
	  a[u] = a1[n]; d[u] = p; u++; p = M;
	}
	else {
	  b[v] = a1[n]; e[v] = q; v++; q = M;
	}
      }
      memcpy(a1, a, u * sizeof(a1[0])); memcpy(a1+u, b, v * sizeof(a1[0]));
      memcpy(d1, d, u * sizeof(d1[0])); memcpy(d1+u, e, v * sizeof(d1[0]));
      
      // perform pre-order traversal and store to HapTree
      if (m % skip == 0) {
	treePtrs[m/skip] = new HapTree(hapBitsT, a1, d1);
      }
    }
    delete[] a1;
    delete[] d1;
    delete[] a;
    delete[] b;
    delete[] d;
    delete[] e;
  }
  HapHedge::~HapHedge(void) {
    for (int t = 0; t < T; t++)
      delete treePtrs[t];
    delete[] treePtrs;
  }
  const HapTree &HapHedge::getHapTree(int t) const {
    return *treePtrs[t];
  };
  int HapHedge::getM(void) const {
    return hapBitsT.getM();
  }
  int HapHedge::getSkip(void) const {
    return skip;
  }
  int HapHedge::getNumTrees(void) const {
    return T;
  }
  const HapBitsT &HapHedge::getHapBitsT(void) const {
    return hapBitsT;
  }

  HapHedgeErr::HapHedgeErr(const HapBitsT &_hapBitsT) :
    hapBitsT(_hapBitsT), T((hapBitsT.getM()+1) / 2) {
    
    treePtrs = new HapTreeMulti *[T];

    int N = hapBitsT.getNhaps(), M = hapBitsT.getM();

    // initialize work arrays
    SortDiv *ad = new SortDiv[N+1], *ad1 = new SortDiv[N+1]; // N+1 for convenience below
    WorkTreeNode *workNodes = (WorkTreeNode *) ALIGNED_MALLOC(N * sizeof(workNodes[0]));
    for (int n = 0; n <= N; n++) { // N+1 for convenience below
      ad[n].a = n;
      ad[n].d = M;
    }

    uint64 *curBits64 = new uint64[N];
    uchar *curBits8 = new uchar[N]; // N-byte buffer hopefully fits in L1 cache for random access

    for (int m = M-1; m >= 0; m--) {
      if (m == M-1 || (m&63) == 63) { // move current 64-bit m-block for each sample to curBits64
	for (int n = 0; n < N; n++)
	  curBits64[n] = hapBitsT.getBits64(n, m>>6);
      }
      //uint64 curBit64 = 1ULL<<(m&63);
      if (m == M-1 || (m&7) == 7) { // move current 8-bit m-block for each sample to curBits8
	uint64 shift = (m&63)&~7;
	for (int n = 0; n < N; n++)
	  curBits8[n] = (uchar) (curBits64[n]>>shift);
      }
      uchar curBit8 = 1<<(m&7);

      // compute sort order and divergence array
      int u = 0, v = 0, p = m, q = m;
      if (m % 2 != 1) {
	for (int n = 0; n < N; n++) {
	  int a_n = ad[n].a;
	  if (/*curBits64[a_n] & curBit64*/curBits8[a_n] & curBit8) {
	    ad1[v].a = a_n; ad1[v].d = q; v++; q = ad[n+1].d; if (q < p) p = q;
	  }
	  else {
	    ad[u].a = a_n; ad[u].d = p; u++; p = ad[n+1].d; if (p < q) q = p;
	  }
	}
      }
      else {
	for (int n = 0; n < N; n++) {
	  int a_n = ad[n].a;
	  if (/*curBits64[a_n] & curBit64*/curBits8[a_n] & curBit8) {
	    ad1[v].a = a_n; ad1[v].d = M; v++; q = ad[n+1].d; if (q < p) p = q;
	  }
	  else {
	    ad[u].a = a_n; ad[u].d = p; u++; p = ad[n+1].d; if (p < q) q = p;
	  }
	}
	ad1[0].d = m;
      }
      memcpy(ad+u, ad1, v * sizeof(ad[0]));
      
      // perform pre-order traversal and store to HapTree
      if (m % 2 == 0)
	treePtrs[m/2] = new HapTreeMulti(hapBitsT, ad, M, workNodes);
    }

    delete[] curBits8;
    delete[] curBits64;
    ALIGNED_FREE(workNodes);
    delete[] ad1;
    delete[] ad;
  }
  HapHedgeErr::~HapHedgeErr(void) {
    for (int t = 0; t < T; t++)
      delete treePtrs[t];
    delete[] treePtrs;
  }
  const HapTreeMulti &HapHedgeErr::getHapTreeMulti(int t) const {
    return *treePtrs[t];
  };
  int HapHedgeErr::getNumTrees(void) const {
    return T;
  }
  const HapBitsT &HapHedgeErr::getHapBitsT() const {
    return hapBitsT;
  }
  // for debugging
  void HapHedgeErr::printTree(int t) const {
    treePtrs[t]->dfsPrint("", 2*t, 2*T, treePtrs[t]->getRootState());
  }

}
