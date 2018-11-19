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

#include <cstdlib>
#include <iostream>

#include "MemoryUtils.hpp"
#include "Types.hpp"

void *ALIGNED_MALLOC(uint64 size) {
#ifdef USE_MKL_MALLOC
  void *p = mkl_malloc(size, MEM_ALIGNMENT);
#else
  void *p = _mm_malloc(size, MEM_ALIGNMENT);
#endif
  // TODO: change to assert() or dispense with altogether and change ALIGNED_MALLOC to macro?
  if (p == NULL) {
    std::cerr << "ERROR: Failed to allocate " << size << " bytes" << std::endl;
    exit(1);
  } else if ((uint64) p & 0xf) {
    std::cerr << "ERROR: Memory alignment of " << size << " bytes failed" << std::endl;
    exit(1);
  }
  return p;
}
