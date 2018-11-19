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
#include <map>
#include <algorithm>

#include "MemoryUtils.hpp"
#include "Types.hpp"
#include "StaticMultimap.hpp"

namespace EAGLE {

  using std::vector;
  using std::map;
  using std::min;
  using std::swap;

  uint randMWC(uint &z, uint &w, uint mod) {
    z=36969*(z&65535)+(z>>16);
    w=18000*(w&65535)+(w>>16);
    return ((z<<16)+w) % mod;
  }

  StaticMultimap::StaticMultimap() : initialized(false) {}

  // for value=[0..len(keyVec)), keyVec[value] = key
  // ignore if key == -1
  // only store up to maxValuesPerKey (if more values, choose randomly)
  StaticMultimap::StaticMultimap(const vector <uint> &keyVec, uint maxValuesPerKey) {
    init(keyVec, maxValuesPerKey);
  }

  void StaticMultimap::init(const vector <uint> &keyVec, uint maxValuesPerKey) {
    // populate tmp map: non-ignored key -> indices i with keyVec[i]=key
    map < uint, vector <uint> > tmpMap;
    for (uint i = 0; i < keyVec.size(); i++)
      if (keyVec[i] != -1U)
	tmpMap[keyVec[i]].push_back(i);

    // iterate through tmp map to determine total number of values to keep (chop at max)
    nKeys = tmpMap.size();
    nValues = 0;
    for (map < uint, vector <uint> >::iterator it = tmpMap.begin(); it != tmpMap.end(); it++)
      nValues += min((uint) it->second.size(), maxValuesPerKey);
    
    // allocate arrays
    keys = ALIGNED_MALLOC_UINTS(nKeys);
    startInds = ALIGNED_MALLOC_UINTS(nKeys);
    lenValues = ALIGNED_MALLOC_UINTS(nKeys + nValues);

    // initialize Marsaglia's MWC
    uint z = 362436069, w = 521288629;    

    // iterate through tmp map to randomly select and store up to maxValuesPerKey
    uint keysPos = 0, lenValuesPos = 0;
    for (map < uint, vector <uint> >::iterator it = tmpMap.begin(); it != tmpMap.end(); it++) {
      keys[keysPos] = it->first; // store key
      startInds[keysPos] = lenValuesPos; // store address of record in storage array
      keysPos++;
      vector <uint> &values = it->second;
      if (values.size() > maxValuesPerKey) // randomly move maxValuesPerKey elements to front
	for (uint j = 0; j < maxValuesPerKey; j++)
	  swap(values[j], values[j + randMWC(z, w, values.size() - j)]);
      uint nValuesCurKey = min((uint) values.size(), maxValuesPerKey);
      lenValues[lenValuesPos++] = nValuesCurKey; // store length as first entry of record
      for (uint j = 0; j < nValuesCurKey; j++) // store list of values
	lenValues[lenValuesPos++] = values[j];
    }    

    initialized = true;
  }

  StaticMultimap::~StaticMultimap() {
    if (initialized) {
      ALIGNED_FREE(lenValues);
      ALIGNED_FREE(startInds);
      ALIGNED_FREE(keys);
    }
  }

  // returns pointer to record: list size followed by list
  // returns NULL if key not found
  const uint *StaticMultimap::query(uint key) const {
    uint lo = 0, hi = nKeys; // condition: keys[lo] <= key < keys[hi] (hi = nKeys ok)
    while (lo+1<hi) {
      uint mid = (lo+hi)/2;
      if (keys[mid] <= key)
	lo = mid;
      else
	hi = mid;
    }
    if (keys[lo] == key)
      return lenValues + startInds[lo];
    else
      return NULL;
  }

}
