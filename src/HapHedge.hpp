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

#ifndef HAPHEDGE_HPP
#define HAPHEDGE_HPP

#include <vector>

#include "Types.hpp"

using namespace std;

namespace EAGLE {

  struct HapTreeState {
    int seq, node, count;
  };

  struct SortDiv {
    int a, d;
  };

  class HapBitsT {
#ifdef TEST_VS
    const std::vector <std::string> &hapBitsVS;
#endif
    uint64 *haploBitsT;
    uint64 Nhaps, M, M64;
    void setBit(uint64 n, uint64 m);
  public:
#ifdef TEST_VS
    HapBitsT(const std::vector <std::string> &_hapBitsVS);
#else
    HapBitsT(const uint64 *_haploBitsT, uint64 _Nhaps, uint64 Mseg64, const uchar maskSnps64j[]);
    HapBitsT(const uint64 *inBitsT, uint64 Mseg64, const std::vector <uint64> &splits64j,
	     const std::vector <uchar> &splitGenos, const std::vector <uint64_masks> &tgtGenoBits,
	     const std::vector <uint> &bestHaps);
    HapBitsT(const HapBitsT &hapBitsFwdT, int dir); // reverse constructor (dir must equal -1)
    uint64 getBits64(uint64 n, uint64 m64) const;
#endif
    ~HapBitsT(void);
    int getBit(uint64 n, uint64 m) const;
    int getNhaps(void) const;
    int getM(void) const;
  };

  struct HapTreeNode {
    int mSplit, count0, seq1;
  };

  struct HapTreeMultiNode {
    int mSplit, count0, node0, seq1, node1;
    // MEMORY-SAVING ALTERNATIVE: use 1 bit (sign of state.count) vs. node0 to encode next node
  };

  struct WorkTreeNode {
    int d, a, count, up, left, right;
    WorkTreeNode(int _d, int _a, int _count, int _up, int _left, int _right);
  };

  class HapTree {
    const HapBitsT &hapBitsT;
    const int Nhaps; float invNhaps;
    int seq0; // lexicographically first ref seq
    HapTreeNode *nodes; // [Nhaps-1]
  public:
    HapTree(const HapBitsT &_hapBitsT, int a[], int d[]);
    ~HapTree(void);
    float getInvNhaps(void) const;
    HapTreeState getRootState(void) const;
    bool next(int m, HapTreeState &state, int nextBit) const;
  };

  class HapTreeMulti {
    const HapBitsT &hapBitsT;
    const int Nhaps; float invNhaps;
    HapTreeState rootState;
    HapTreeMultiNode *nodes; // [Nhaps-1]
  public:
    HapTreeMulti(const HapBitsT &_hapBitsT, SortDiv ad[], int M, WorkTreeNode workNodes[]);
    ~HapTreeMulti(void);
    float getInvNhaps(void) const;
    HapTreeState getRootState(void) const;
    bool next(int m, HapTreeState &state, int nextBit) const;
    void nextAtFrac(int m, HapTreeState &state, double nextFrac) const;
    void dfsPrint(std::string curPrefix, int m, int M, const HapTreeState &state) const;
  };

  class HapHedge {
    const HapBitsT &hapBitsT;
    const int skip, T;
    HapTree **treePtrs;
  public:
    HapHedge(const HapBitsT &_hapBitsT, int skip/*, const std::vector <int> &treeStarts*/);
    ~HapHedge(void);
    const HapTree &getHapTree(int t) const;
    int getM(void) const;
    int getSkip(void) const;
    int getNumTrees(void) const;
    const HapBitsT &getHapBitsT(void) const;
  };

  class HapHedgeErr {
    const HapBitsT &hapBitsT;
    const int T;
    HapTreeMulti **treePtrs;
  public:
    HapHedgeErr(const HapBitsT &_hapBitsT);
    ~HapHedgeErr(void);
    const HapTreeMulti &getHapTreeMulti(int t) const;
    int getNumTrees(void) const;
    const HapBitsT &getHapBitsT() const;
    void printTree(int t) const;
  };

}

#endif
