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

#ifndef GENODATA_HPP
#define GENODATA_HPP

#include <vector>
#include <string>
#include <map>
#include <boost/utility.hpp>

#include "Types.hpp"
#include "FileUtils.hpp"

namespace EAGLE {
  
  struct SnpInfoX {
    int chrom;
    std::string ID;
    double genpos; // Morgans
    int physpos;
    std::string allele1, allele2;
    double MAF; // note: MAFs are computed on preQC indivs
    double miss;
    bool passQC;
  };

  struct IndivInfoX {
    std::string famID;
    std::string indivID;
    std::string paternalID;
    std::string maternalID;
    int sex; // (1=male; 2=female; other=unknown)
    double pheno;
    double miss;
    bool passQC;
  };

  struct AlleleFreqs {
    double cond[3][6];
    /*
    double dip[3];
    double hap[2];
    */
  };

  class GenoData : boost::noncopyable {

  private:
    uint64 NpreQC, MpreQC; // # of indivs/snps in fam/bim file not in --remove/exclude file(s)
    uint64 N, M; // post-QC

    std::vector <IndivInfoX> indivsPreQC, indivs; // [VECTOR]: NpreQC, N
    std::vector <SnpInfoX> snpsPreQC, snps; // [VECTOR]: MpreQC, M

    uint64 Mseg64; // number of <=64-SNP chunks
    uint64_masks *genoBits; // [[MATRIX]]: M64 x N (is0, is2, is9 64-bit masks; 3 bits/base)
    std::vector <std::vector <double> > seg64cMvecs;
    std::vector <std::vector <uint64> > seg64preQCsnpInds;
    AlleleFreqs *seg64logPs;
    std::vector <bool> isFlipped64j;

    std::vector <bool> processIndivs(const std::string &famFile,
				     const std::vector <std::string> &removeFiles);
    std::vector <bool> processSnps(const std::string &bimFile, int chrom, double bpStart,
				   double bpEnd, const std::vector <std::string> &excludeFiles);
    void processMap(std::vector <SnpInfoX> &snpsVec, const std::string &geneticMapFile,
		    bool noMapCheck);
    void buildGenoBits(uchar *genosPreQC, const std::vector <bool> &genos2bit, double cMmax);
    bool fillSnpSubRowNorm1(float x[], uint64 m64j, const std::vector <int> &indivInds) const;

  public:
    /**    
     * reads indiv info from fam file, snp info from bim file
     * allocates memory, reads genotypes, and does QC
     * assumes numbers of bim and bed files match
     */
    void initBed(const std::string &famFile, const std::string &bimFile,
		 const std::string &bedFile, int chrom, double bpStart, double bpEnd,
		 const std::string &geneticMapFile, const std::vector <std::string> &excludeFiles,
		 const std::vector <std::string> &removeFiles, double maxMissingPerSnp,
		 double maxMissingPerIndiv, bool noMapCheck, double cMmax);
    /**    
     * reads genotypes from VCF/BCF file
     * does not save indiv info (will be reread from VCF during output)
     * only saves chrom, physpos, genpos in snp info (rest will be reread from VCF during output)
     * allocates memory, reads genotypes, and restricts to region if specified; does not do QC
     */
    void initVcf(const std::string &vcfFile, const int inputChrom, const int chromX,
		 double bpStart, double bpEnd, const std::string &geneticMapFile, bool noMapCheck,
		 double cMmax);

    void printRange(void) const;

    ~GenoData();

    static int plinkChromCode(const std::string &chrom);
    static std::vector <SnpInfoX> readBimFile(const std::string &bimFile);
    /**
     * assumes Nbed = bedIndivRemoved.size()
     * reads (Nbed+3)>>2 bytes into bedLineIn
     * stores sum(!bedIndivRemoved) bytes into genoLine if storeGenoLine == true
     */
    static void readBedLine(FileUtils::AutoGzIfstream &fin, uchar bedLineIn[], uchar genoLine[], 
			    std::vector <bool> &bedIndivRemoved, bool storeGenoLine);
    static double computeAlleleFreq(const uchar genoLine[], uint64 genoN);
    static double computeMAF(const uchar genoLine[], uint64 genoN);
    static double computeSnpMissing(const uchar genoLine[], uint64 genoN);

    const std::vector <SnpInfoX> &getSnps(void) const;
    uint64 getN(void) const;
    uint64 getMseg64(void) const;
    const uint64_masks *getGenoBits(void) const;
    std::vector <std::vector <double> > getSeg64cMvecs(void) const;
    const AlleleFreqs *getSeg64logPs(void) const;
    IndivInfoX getIndiv(uint64 n) const;
    std::vector <double> computeInvLD64j(uint64 NsubMax) const;
    const std::vector <IndivInfoX> &getIndivs(void) const;
    const std::vector <bool> &getIsFlipped64j(void) const;
    double computeSnpRate(void) const;

  };
}

#endif
