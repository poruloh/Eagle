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

#ifndef STATICMULTIMAP_HPP
#define STATICMULTIMAP_HPP

#include <vector>

#include "Types.hpp"

namespace EAGLE {

  class StaticMultimap {

  private:

    uint nKeys, nValues;
    uint *keys; // [VECTOR]: nKeys
    uint *startInds; // [VECTOR]: nKeys
    uint *lenValues; // [VECTOR]: nKeys + nValues (each record: list size followed by list)
    bool initialized;

  public:
  
    StaticMultimap();
    // for value=[0..len(keyVec)), keyVec[value] = key
    // ignore if key == -1
    // only store up to maxValuesPerKey (if more values, choose randomly)
    StaticMultimap(const std::vector <uint> &keyVec, uint maxValuesPerKey);
    void init(const std::vector <uint> &keyVec, uint maxValuesPerKey);
    ~StaticMultimap();

    // returns pointer to record: list size followed by list
    // returns NULL if key not found
    const uint *query(uint key) const;

  };

}

#endif
