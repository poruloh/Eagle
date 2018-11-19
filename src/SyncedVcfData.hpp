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

#ifndef SYNCEDVCFDATA_HPP
#define SYNCEDVCFDATA_HPP

#include <vector>
#include <string>
#include <utility>
#include <boost/utility.hpp>

#include "Types.hpp"

namespace EAGLE {
  
  class SyncedVcfData : boost::noncopyable {

  private:
    uint64 Nref, Ntarget, M; // N = Nref + Ntarget
    uint64 Mseg64; // number of <=64-SNP chunks
    uint64_masks *genoBits; // [[MATRIX]]: M64 x N (is0, is2, is9 64-bit masks; 3 bits/base)
                            //             n = [0..Nref) contain ref haploBits in is0, is2
                            //             n = [Nref..Nref+Ntarget) contain target genoBits
    std::vector <std::vector <double> > seg64cMvecs;
    std::vector <std::string> targetIDs;

    std::vector < std::pair <int, int> > processVcfs
    (const std::string &vcfRef, const std::string &vcfTarget, const std::string &vcfExclude,
     bool allowRefAltSwap, int &chrom, int chromX, double bpStart, double bpEnd,
     std::vector <bool> &hapsRef, std::vector <uchar> &genosTarget, const std::string &tmpFile,
     const std::string &writeMode, int usePS,
     std::vector < std::vector < std::pair <int, int> > > &conPSall, bool outputUnphased,
     std::vector <bool> &isTmpPhased);
    std::vector <double> processMap(std::vector < std::pair <int, int> > &chrBps,
				    const std::string &geneticMapFile);
    void buildGenoBits(const std::vector <bool> &hapsRef, const std::vector <uchar> &genosTarget,
		       const std::vector <double> &cMs, double cMmax);

  public:
    /**    
     * reads ref+target vcf data
     * writes target[isec] to tmpFile
     * fills in cM coordinates and seg64cMvecs, genoBits
     */
    SyncedVcfData(const std::string &vcfRef, const std::string &vcfTarget,
		  const std::string &vcfExclude, bool allowRefAltSwap, int &chrom, int chromX,
		  double bpStart, double bpEnd, const std::string &geneticMapFile, double cMmax,
		  const std::string &tmpFile, const std::string &writeMode, int usePS,
		  std::vector < std::vector < std::pair <int, int> > > &conPSall, double &snpRate,
		  bool outputUnphased, std::vector <bool> &isTmpPhased);

    ~SyncedVcfData();

    uint64 getNref(void) const;
    uint64 getNtarget(void) const;
    uint64 getMseg64(void) const;
    const uint64_masks *getGenoBits(void) const;
    std::vector <std::vector <double> > getSeg64cMvecs(void) const;
    const std::string &getTargetID(int n) const;

  };
}

#endif
