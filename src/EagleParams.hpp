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

#ifndef EAGLEPARAMS_HPP
#define EAGLEPARAMS_HPP

#include <vector>
#include <string>

#include "Types.hpp"

namespace EAGLE {

  class EagleParams {
  public:

    // main input files
    std::string famFile, bimFile, bedFile, vcfFile, vcfRef, vcfTarget, vcfExclude;
    int chrom, chromX;

    // optional reference map file for filling in genpos
    std::string geneticMapFile;

    std::vector <std::string> removeFiles; // list(s) of indivs to remove
    std::vector <std::string> excludeFiles; // list(s) of SNPs to exclude
    double bpStart, bpEnd, bpFlanking;

    std::string outPrefix; // .haps.gz .sample
    std::string vcfOutFormat, vcfOutSuffix, vcfWriteMode;
    // outFormat b|u|z|v -> outSuffix bcf|bcf|vcf.gz|vcf, writeMode wb|wbu|wz|w
    bool noImpMissing;
    bool outputUnphased; // output unphased sites (target-only, multi-allelic, etc.)
    bool keepMissingPloidyX; // assume missing genotypes in VCF have correct ploidy (.=haploid, ./.=diploid)
    int usePS; // use FORMAT:PS phase constraints: 1=soft, 2=harder

    bool allowRefAltSwap; // in reference-based phasing mode
    bool usePBWT;
    bool pbwtOnly; // in non-ref mode, don't run Steps 1 or 2
    int runStep2; // in non-ref mode, do/don't run Step 2
    int pbwtIters;
    double expectIBDcM; // expected length of an IBD segment (for transition probabilities)
    double histFactor; // history length multiplier

    // QC params
    double maxMissingPerSnp, maxMissingPerIndiv;
    
    int numThreads;

    double cMmax;
    int beamWidth3, beamWidth4;
    int maxHapStates;
    double fracStep4;
    double pErr;
    uint seed;
    bool noMapCheck;

    int Kpbwt;

    // testing
    int iter;
    std::string tmpPhaseConfsPrefix;
    std::string maskFile; // list of indivs to mask (e.g., relatives)
    bool trioCheck;

    // populates members; error-checks
    bool processCommandLineArgs(int argc, char *argv[]);
  };
}

#endif
