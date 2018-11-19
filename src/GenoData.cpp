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

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <cstring>
#include <cmath>

#include <htslib/thread_pool.h>
#include <htslib/vcf.h>

#include "omp.h"

#include "Types.hpp"
#include "FileUtils.hpp"
#include "MemoryUtils.hpp"
#include "MapInterpolater.hpp"
#include "LapackConst.hpp"
#include "GenoData.hpp"

namespace EAGLE {

  using std::vector;
  using std::string;
  using std::cout;
  using std::cerr;
  using std::endl;
  using FileUtils::getline;

  int GenoData::plinkChromCode(const string &chrom) {
    if (isdigit(chrom[0])) return atoi(chrom.c_str());
    if (chrom == "X") return 23;
    if (chrom == "Y") return 24;
    if (chrom == "XY") return 25;
    if (chrom == "MT") return 26;
    return -1;
  }

  // set indivsPreQC, NpreQC
  // return bedIndivRemoved
  vector <bool> GenoData::processIndivs(const string &famFile,
					const vector <string> &removeFiles) {
    std::map <string, uint64> FID_IID_to_ind;
    string line;

    vector <IndivInfoX> bedIndivs;
    cout << "Reading fam file: " << famFile << endl;
    FileUtils::AutoGzIfstream fin; fin.openOrExit(famFile);
    while (getline(fin, line)) {
      std::istringstream iss(line);
      IndivInfoX indiv;
      if (!(iss >> indiv.famID >> indiv.indivID >> indiv.paternalID >> indiv.maternalID
	    >> indiv.sex >> indiv.pheno)) {
	cerr << "ERROR: Incorrectly formatted fam file: " << famFile << endl;
	cerr << "Line " << bedIndivs.size()+1 << ":" << endl;
	cerr << line << endl;
	cerr << "Unable to input 6 values (4 string, 1 int, 1 double)" << endl;
	exit(1);
      }
      string combined_ID = indiv.famID + " " + indiv.indivID;
      if (FID_IID_to_ind.find(combined_ID) != FID_IID_to_ind.end()) {
	cerr << "ERROR: Duplicate individual in fam file at line " << bedIndivs.size()+1 << endl;
	exit(1);
      }
      FID_IID_to_ind[combined_ID] = bedIndivs.size();
      bedIndivs.push_back(indiv);
    }
    fin.close();
    uint64 Nbed = bedIndivs.size();

    cout << "Total indivs in PLINK data: Nbed = " << Nbed << endl;

    // process individuals to remove
    vector <bool> bedIndivRemoved(Nbed);
    for (uint f = 0; f < removeFiles.size(); f++) {
      const string &removeFile = removeFiles[f];
      cout << "Reading remove file (indivs to remove): " << removeFile << endl;
      fin.openOrExit(removeFile);
      int lineCtr = 0;
      int numRemoved = 0;
      int numAbsent = 0;
      while (getline(fin, line)) {
	lineCtr++;
	std::istringstream iss(line);
	string FID, IID;
	if (!(iss >> FID >> IID)) {
	  cerr << "ERROR: Incorrectly formatted remove file: " << removeFile << endl;
	  cerr << "Line " << lineCtr << ":" << endl;
	  cerr << line << endl;
	  cerr << "Unable to input FID and IID" << endl;
	  exit(1);
	}
	string combined_ID = FID + " " + IID;
	if (FID_IID_to_ind.find(combined_ID) == FID_IID_to_ind.end()) {
	  if (numAbsent < 5)
	    cerr << "WARNING: Unable to find individual to remove: " << combined_ID << endl;
	  numAbsent++;
	}
	else if (!bedIndivRemoved[FID_IID_to_ind[combined_ID]]) {
	  bedIndivRemoved[FID_IID_to_ind[combined_ID]] = true;
	  numRemoved++;
	}
      }
      fin.close();
      cout << "Removed " << numRemoved << " individual(s)" << endl;
      if (numAbsent)
	cerr << "WARNING: " << numAbsent << " individual(s) not found in data set" << endl;
    }

    for (uint64 nbed = 0; nbed < Nbed; nbed++)
      if (!bedIndivRemoved[nbed])
	indivsPreQC.push_back(bedIndivs[nbed]);
    NpreQC = indivsPreQC.size();
    cout << "Total indivs stored in memory: NpreQC = " << NpreQC << endl;

    return bedIndivRemoved;
  }

  vector <SnpInfoX> GenoData::readBimFile(const string &bimFile) {
    vector <SnpInfoX> ret;
    string line;
    FileUtils::AutoGzIfstream fin; fin.openOrExit(bimFile);
    int numOutOfOrder = 0;
    while (getline(fin, line)) {
      std::istringstream iss(line);
      SnpInfoX snp; string chrom_str;
      if (!(iss >> chrom_str >> snp.ID >> snp.genpos >> snp.physpos >> snp.allele1 >> snp.allele2))
	{
	cerr << "ERROR: Incorrectly formatted bim file: " << bimFile << endl;
	cerr << "Line " << ret.size()+1 << ":" << endl;
	cerr << line << endl;
	cerr << "Unable to input 6 values (2 string, 1 double, 1 int, 2 string)" << endl;
	exit(1);
      }
      snp.chrom = plinkChromCode(chrom_str);
      if (snp.chrom == -1) {
	cerr << "ERROR: Unknown chromosome code in bim file: " << bimFile << endl;
	cerr << "Line " << ret.size()+1 << ":" << endl;
	cerr << line << endl;
	exit(1);
      }
      if (!ret.empty() &&
	  (snp.chrom < ret.back().chrom ||
	   (snp.chrom == ret.back().chrom && (snp.physpos < ret.back().physpos ||
					       snp.genpos < ret.back().genpos)))) {
	if (numOutOfOrder < 5) {
	  cerr << "WARNING: Out-of-order snp in bim file: " << bimFile << endl;
	  cerr << "Line " << ret.size()+1 << ":" << endl;
	  cerr << line << endl;
	}
	numOutOfOrder++;
	//exit(1);
      }
      ret.push_back(snp);
    }
    if (numOutOfOrder)
      cerr << "WARNING: Total number of out-of-order snps in bim file: " << numOutOfOrder << endl;
    fin.close();
    // TODO: exit with error or sort SNPs?
    return ret;
  }

  // set snpsPreQC, MpreQC
  // return bedSnpExcluded
  vector <bool> GenoData::processSnps(const string &bimFile, int chrom, double bpStart,
				      double bpEnd, const vector <string> &excludeFiles) {
    FileUtils::AutoGzIfstream fin;
    string line;

    // read bim file
    cout << "Reading bim file: " << bimFile << endl;
    vector <SnpInfoX> bedSnps = readBimFile(bimFile);

    uint64 Mbed = bedSnps.size();
    cout << "Total snps in PLINK data: Mbed = " << Mbed << endl;

    vector <bool> bedSnpExcluded(Mbed);

    if (chrom == 0) {
      if (bedSnps[0].chrom != bedSnps.back().chrom) {
	cerr << "ERROR: Only one chromosome may be analyzed at a time; use --chrom" << endl;
	exit(1);
      }
      else
	chrom = bedSnps[0].chrom;
    }

    uint64 MbedOnChrom = 0;
    for (uint64 mbed = 0; mbed < Mbed; mbed++) {
      if (bedSnps[mbed].chrom != chrom || bedSnps[mbed].physpos < bpStart ||
	  bedSnps[mbed].physpos > bpEnd)
	bedSnpExcluded[mbed] = true;
      else
	MbedOnChrom++;
    }
    if (MbedOnChrom < Mbed) {
      cout << "Restricting to " << MbedOnChrom << " SNPs on chrom " << chrom
	   << " in region [bpStart,bpEnd] = [" << bpStart << "," << bpEnd << "]" << endl;
    }

    // create dictionary rsID -> index in full bed snp list
    std::map <string, uint64> rsID_to_ind;
    for (uint64 mbed = 0; mbed < Mbed; mbed++) {
      if (rsID_to_ind.find(bedSnps[mbed].ID) != rsID_to_ind.end()) {
	cerr << "WARNING: Duplicate snp ID " << bedSnps[mbed].ID
	     << " -- masking duplicate" << endl;
	bedSnpExcluded[mbed] = true;
      }
      else
	rsID_to_ind[bedSnps[mbed].ID] = mbed;
    }

    // TODO: also limit SNP density (e.g., Omni 2.5M) and/or MAF?
    // process snps to exclude
    for (uint f = 0; f < excludeFiles.size(); f++) {
      const string &excludeFile = excludeFiles[f];
      cout << "Reading exclude file (SNPs to exclude): " << excludeFile << endl;
      fin.openOrExit(excludeFile);
      int numExcluded = 0;
      int numAbsent = 0;
      while (getline(fin, line)) {
	std::istringstream iss(line);
	string rsID; iss >> rsID;
	if (rsID_to_ind.find(rsID) == rsID_to_ind.end()) {
	  if (numAbsent < 5)
	    cerr << "WARNING: Unable to find SNP to exclude: " << rsID << endl;
	  numAbsent++;
	}
	else if (!bedSnpExcluded[rsID_to_ind[rsID]]) {
	  bedSnpExcluded[rsID_to_ind[rsID]] = true;
	  numExcluded++;
	}
      }
      fin.close();
      cout << "Excluded " << numExcluded << " SNP(s)" << endl;
      if (numAbsent)
	cerr << "WARNING: " << numAbsent << " SNP(s) not found in data set" << endl;
    }

    for (uint64 mbed = 0; mbed < Mbed; mbed++)
      if (!bedSnpExcluded[mbed])
	snpsPreQC.push_back(bedSnps[mbed]);
    MpreQC = snpsPreQC.size();
    cout << "Total SNPs stored in memory: MpreQC = " << MpreQC << endl;

    return bedSnpExcluded;
  }

  void GenoData::processMap(vector <SnpInfoX> &snpsVec, const string &geneticMapFile,
			    bool noMapCheck) {
    // fill in map if external file provided
    if (geneticMapFile != "USE_BIM") {
      cout << "Filling in genetic map coordinates using reference file:" << endl;
      cout << "  " << geneticMapFile << endl;
      Genetics::MapInterpolater mapInterpolater(geneticMapFile);
      for (uint64 m = 0; m < snpsVec.size(); m++)
	snpsVec[m].genpos = mapInterpolater.interp(snpsVec[m].chrom, snpsVec[m].physpos);
    }
    else {
      // check map and rescale if in cM units: calculate d(genpos)/d(physpos)
      double scale = (snpsVec.back().genpos - snpsVec[0].genpos)
	/ (snpsVec.back().physpos - snpsVec[0].physpos);
      if (0.5e-6 < scale && scale < 3e-6) {
	cerr << "WARNING: Genetic map appears to be in cM units; rescaling by 0.01" << endl;
	for (uint64 m = 0; m < snpsVec.size(); m++)
	  snpsVec[m].genpos *= 0.01;
      }
      else if (!(0.5e-8 < scale && scale < 3e-8)) {
	if (noMapCheck) {
	  cerr << "WARNING: Genetic map appears wrong based on overall cM/Mb" << endl;
	  cerr << "         Proceeding anyway because --noMapCheck is set" << endl;
	}
	else {
	  cerr << "ERROR: Genetic map appears wrong based on overall cM/Mb" << endl;
	  cerr << "       To proceed anyway, set --noMapCheck" << endl;
	  exit(1);
	}
      }
    }
  }

  inline double log10safe(double x) { return x > 0 ? log10(x) : -1000; }

  void GenoData::buildGenoBits(uchar *genosPreQC, const vector <bool> &genos2bit, double cMmax) {
    const uint segMin = 16;
    vector <uint64> preQCsnpInds; vector <double> cMvec;
    for (uint64 m = 0; m < MpreQC; m++)
      if (snpsPreQC[m].passQC) {
	if (preQCsnpInds.size() == 64 ||
	    (preQCsnpInds.size() >= segMin &&
	     snpsPreQC[m].genpos > snpsPreQC[preQCsnpInds[0]].genpos + cMmax/100)) {
	  seg64preQCsnpInds.push_back(preQCsnpInds); seg64cMvecs.push_back(cMvec);
	  preQCsnpInds.clear(); cMvec.clear();
	}
	preQCsnpInds.push_back(m); cMvec.push_back(100 * snpsPreQC[m].genpos);
      }
    seg64preQCsnpInds.push_back(preQCsnpInds); seg64cMvecs.push_back(cMvec);

    Mseg64 = seg64preQCsnpInds.size();
    cout << "Number of <=(64-SNP, " << cMmax << "cM) segments: " << Mseg64 << endl;
    cout << "Average # SNPs per segment: " << M / Mseg64 << endl;

    isFlipped64j = vector <bool> (Mseg64*64);
    genoBits = ALIGNED_MALLOC_UINT64_MASKS(Mseg64 * N);
    memset(genoBits, 0, Mseg64 * N * sizeof(genoBits[0]));
    seg64logPs = (AlleleFreqs *) ALIGNED_MALLOC(Mseg64*64 * sizeof(AlleleFreqs));
    for (uint64 m64j = 0; m64j < Mseg64*64; m64j++)
      for (uint g1 = 0; g1 <= 2; g1++)
	for (uint g0 = 0; g0 <= 3; g0++)
	  seg64logPs[m64j].cond[g1][g0] = NAN;

    for (uint64 m64 = 0; m64 < Mseg64; m64++) {
      for (uint64 j = 0; j < seg64preQCsnpInds[m64].size(); j++) {
        uint64 m = seg64preQCsnpInds[m64][j]; // m, n indices are preQC
	int genoCounts[3]; genoCounts[0] = genoCounts[1] = genoCounts[2] = 0;
	for (uint64 n = 0; n < NpreQC; n++)
	  if (indivsPreQC.empty() || indivsPreQC[n].passQC) {
	    uchar geno = genosPreQC != NULL ? genosPreQC[m * NpreQC + n] :
	      (genos2bit[2*(m * NpreQC + n)] + 2*genos2bit[2*(m * NpreQC + n)+1]);
	    if (geno <= 2) genoCounts[geno]++;
	  }
	uchar is0geno = 0, is2geno = 2;
	if (genoCounts[2] > genoCounts[0]) {
	  isFlipped64j[m64*64+j] = true;
	  is0geno = 2; is2geno = 0;
	  std::swap(genoCounts[0], genoCounts[2]);
	}
	uint64 nPostQC = 0;
	for (uint64 n = 0; n < NpreQC; n++)
	  if (indivsPreQC.empty() || indivsPreQC[n].passQC) {
	    uchar geno = genosPreQC != NULL ? genosPreQC[m * NpreQC + n] :
	      (genos2bit[2*(m * NpreQC + n)] + 2*genos2bit[2*(m * NpreQC + n)+1]);
	    genoBits[m64 * N + nPostQC].is0 |= ((uint64) (geno == is0geno))<<j;
	    genoBits[m64 * N + nPostQC].is2 |= ((uint64) (geno == is2geno))<<j;
	    genoBits[m64 * N + nPostQC].is9 |= ((uint64) (geno > 2))<<j;
	    nPostQC++;
	  }
	double tot = genoCounts[0] + genoCounts[1] + genoCounts[2];
	double p0 = genoCounts[0]/tot, p1half = 0.5*genoCounts[1]/tot, p2 = genoCounts[2]/tot;
	if (p1half == 0) p1half = 1e-9; // avoid division by 0
	double p = (genoCounts[1] + 2*genoCounts[2]) / (2*tot);
	AlleleFreqs &af = seg64logPs[m64*64+j];

	for (uint g1 = 0; g1 <= 2; g1++)
	  af.cond[g1][3] = log10safe(genoCounts[g1] / tot); // [3]: unconditioned

	af.cond[2][2] = log10safe(p2 / (p1half + p2));
	af.cond[1][2] = log10safe(p1half / (p1half + p2));
	af.cond[0][2] = log10safe(0);

	af.cond[2][1] = log10safe(0.5 * p2 / (p1half + p2));
	af.cond[1][1] = log10safe(0.5 * p1half / (p1half + p2) + 0.5 * p1half / (p0 + p1half));
	af.cond[0][1] = log10safe(0.5 * p0 / (p0 + p1half));

	af.cond[2][0] = log10safe(0);
	af.cond[1][0] = log10safe(p1half / (p0 + p1half));
	af.cond[0][0] = log10safe(p0 / (p0 + p1half));

	// same orientation as het-het => p(hap=1) = 1-p
	af.cond[2][4] = log10safe((1-p) * p2 / (p1half + p2));
	af.cond[1][4] = log10safe((1-p) * p1half / (p1half + p2) + p * p1half / (p0 + p1half));
	af.cond[0][4] = log10safe(p * p0 / (p0 + p1half));

	// opp orientation to het-het => p(hap=1) = p
	af.cond[2][5] = log10safe(p * p2 / (p1half + p2));
	af.cond[1][5] = log10safe(p * p1half / (p1half + p2) + (1-p) * p1half / (p0 + p1half));
	af.cond[0][5] = log10safe((1-p) * p0 / (p0 + p1half));

	if (p > 0.55) {
	  cerr << "INTERNAL ERROR: Minor/major allele coding bug" << endl;
	  exit(1);
	}
      }
      for (uint64 n = 0; n < N; n++)
	for (uint64 j = seg64preQCsnpInds[m64].size(); j < 64; j++)
	  genoBits[m64*N+n].is9 |= 1ULL<<j;
    }
  }

  /**
   * fills x[] with indivInds.size() elements corresponding to chosen subset of indivInds
   * replaces missing values with average; mean-centers and normalizes vector length to 1
   * if monomorphic among non-missing, fills with all-0s
   *
   * return: true if snp is polymorphic in indivInds; false if not
   */
  bool GenoData::fillSnpSubRowNorm1(float x[], uint64 m64j, const vector <int> &indivInds) const {
    uint64 m64 = m64j/64, j = m64j&63, jBit = 1ULL<<j;
    float sumPresent = 0; int numPresent = 0;
    for (uint64 i = 0; i < indivInds.size(); i++) {
      uint64 n = indivInds[i];
      if (genoBits[m64 * N + n].is0 & jBit) x[i] = 0.0f;
      else if (genoBits[m64 * N + n].is2 & jBit) x[i] = 2.0f;
      else if (genoBits[m64 * N + n].is9 & jBit) x[i] = 9.0f;
      else x[i] = 1.0f;
      if (x[i] != 9.0f) {
	sumPresent += x[i];
	numPresent++;
      }
    }
    float avg = sumPresent / numPresent;
    float sum2 = 0;
    for (uint64 i = 0; i < indivInds.size(); i++) {
      if (x[i] != 9.0f) { // non-missing; mean-center
	x[i] -= avg;
	sum2 += x[i]*x[i];
      }
      else // missing; replace with mean (centered to 0)
	x[i] = 0;
    }
    if (sum2 < 0.001) { // monomorphic among non-missing
      for (uint64 i = 0; i < indivInds.size(); i++) x[i] = 0; // set to 0
      return false;
    }
    else { // polymorphic
      float invNorm = 1.0f / sqrtf(sum2);
      for (uint64 i = 0; i < indivInds.size(); i++) x[i] *= invNorm; // normalize to vector len 1
      return true;
    }
  }

  float dotProdToAdjR2(float dotProd, int n) {
    float r2 = dotProd*dotProd;
    return r2 - (1-r2)/(n-2);
  }

  vector <double> GenoData::computeInvLD64j(uint64 NsubMax) const {
    uint64 Mchr = Mseg64*64;
    vector <double> chipLDscores(Mchr, 1.0);

    uint64 step = std::max(N / NsubMax, 1ULL);
    vector <int> indivInds;
    for (uint64 n = 0; n < N && indivInds.size() < NsubMax; n += step)
      indivInds.push_back(n);
    uint64 Nsub = indivInds.size();
    cout << "Estimating LD scores using " << Nsub << " indivs" << endl;

    // allocate memory
    uchar *chrMaskSnps = ALIGNED_MALLOC_UCHARS(Mchr);
    memset(chrMaskSnps, 0, Mchr * sizeof(chrMaskSnps[0]));
    float *chrNormalizedGenos = ALIGNED_MALLOC_FLOATS(Mchr * Nsub);
    memset(chrNormalizedGenos, 0, Mchr * Nsub * sizeof(chrNormalizedGenos[0]));
    const int mBlock = 64;
    float *dotProds = ALIGNED_MALLOC_FLOATS(Mchr * mBlock);

    // fill and normalize genotypes
    for (uint64 mchr = 0; mchr < Mchr; mchr++) {
      if ((mchr&63) < seg64cMvecs[mchr/64].size())
	chrMaskSnps[mchr] = fillSnpSubRowNorm1(chrNormalizedGenos + mchr*Nsub, mchr, indivInds);
      else
	chipLDscores[mchr] = 0;
    }

    uint64 mchrWindowStart = 0;
    for (uint64 mchr0 = 0; mchr0 < Mchr; mchr0 += mBlock) { // sgemm to compute r2s
      uint64 mBlockCrop = std::min(Mchr, mchr0+mBlock) - mchr0;
      while ((mchrWindowStart&63) >= seg64cMvecs[mchrWindowStart/64].size() ||
	     seg64cMvecs[mchrWindowStart/64][mchrWindowStart&63] + 1 <
	     seg64cMvecs[mchr0/64][mchr0&63])
	mchrWindowStart++;
      uint64 prevWindowSize = mchr0+mBlockCrop-1 - mchrWindowStart;

      // [mchrWindowStart..mchr0+mBlockCrop-1) x [mchr0..mchr0+mBlockCrop)
      {
	char TRANSA_ = 'T';
	char TRANSB_ = 'N';
	int M_ = prevWindowSize;
	int N_ = mBlockCrop;
	int K_ = Nsub;
	float ALPHA_ = 1;
	float *A_ = chrNormalizedGenos + mchrWindowStart*Nsub;
	int LDA_ = Nsub;
	float *B_ = chrNormalizedGenos + mchr0*Nsub;
	int LDB_ = Nsub;
	float BETA_ = 0;
	float *C_ = dotProds;
	int LDC_ = prevWindowSize;
	SGEMM_MACRO(&TRANSA_, &TRANSB_, &M_, &N_, &K_, &ALPHA_, A_, &LDA_, B_, &LDB_,
		    &BETA_, C_, &LDC_);
      }

      for (uint64 mPlus = 0; mPlus < mBlockCrop; mPlus++) {
	uint64 m = mchr0 + mPlus;
	if (!chrMaskSnps[m]) continue;
	for (uint64 mPlus2 = 0; mchrWindowStart+mPlus2 < mchr0+mPlus; mPlus2++) {
	  uint64 m2 = mchrWindowStart + mPlus2;
	  if (!chrMaskSnps[m2]) continue;
	  float adjR2 = dotProdToAdjR2(dotProds[mPlus2 + mPlus*prevWindowSize], Nsub);
	  chipLDscores[m] += adjR2;
	  chipLDscores[m2] += adjR2;
	}
      }
    }

    ALIGNED_FREE(dotProds);
    ALIGNED_FREE(chrNormalizedGenos);
    ALIGNED_FREE(chrMaskSnps);

    for (uint mchr = 0; mchr < Mchr; mchr++) chipLDscores[mchr] = 1/chipLDscores[mchr];
    return chipLDscores; // reciprocals taken above
  }

  void GenoData::printRange(void) const {

    int physRange = snps.back().physpos - snps[0].physpos;
    double cMrange = 100*(snps.back().genpos - snps[0].genpos);

    cout << "Physical distance range: " << physRange << " base pairs" << endl;
    cout << "Genetic distance range:  " << cMrange << " cM" << endl;
    cout << "Average # SNPs per cM:   " << (int) (M/cMrange+0.5) << endl;

    if (physRange == 0 || cMrange == 0) {
      cerr << "ERROR: Physical and genetic distance ranges must be positive" << endl;
      cerr << "       First SNP: chr=" << snps[0].chrom << " pos=" << snps[0].physpos
	   << " cM=" << 100*snps[0].genpos << endl;
      cerr << "       Last SNP:  chr=" << snps.back().chrom << " pos=" << snps.back().physpos
	   << " cM=" << 100*snps.back().genpos << endl;
      exit(1);
    }
  }

  double GenoData::computeSnpRate(void) const {
    double cMrange = 100*(snps.back().genpos - snps[0].genpos);
    return M/cMrange;
  }

  /**
   * reads indiv info from fam file, snp info from bim file
   * allocates memory, reads genotypes, and does QC
   */
  void GenoData::initBed(const string &famFile, const string &bimFile, const string &bedFile,
			 int chrom, double bpStart, double bpEnd, const string &geneticMapFile,
			 const vector <string> &excludeFiles, const vector <string> &removeFiles,
			 double maxMissingPerSnp, double maxMissingPerIndiv, bool noMapCheck,
			 double cMmax) {

    // indivsPreQC (without --remove indivs)
    vector <bool> bedIndivRemoved = processIndivs(famFile, removeFiles);
    // snpsPreQC (restricted to chrom:bpStart-bpEnd and without --exclude snps)
    vector <bool> bedSnpExcluded = processSnps(bimFile, chrom, bpStart, bpEnd, excludeFiles);
    processMap(snpsPreQC, geneticMapFile, noMapCheck); // modify snpsPreQC

    // allocate genotypes
    cout << "Allocating " << MpreQC << " x " << NpreQC << " bytes to temporarily store genotypes"
	 << endl;
    uchar *genosPreQC = ALIGNED_MALLOC_UCHARS(MpreQC * NpreQC); // temporary

    cout << "Reading genotypes and performing QC filtering on snps and indivs..." << endl;

    // open bed file
    const uint64 Nbed = bedIndivRemoved.size(), Mbed = bedSnpExcluded.size();
    cout << "Reading bed file: " << bedFile << endl;
    cout << "    Expecting " << Mbed * ((Nbed+3)>>2) << " (+3) bytes for "
	 << Nbed << " indivs, " << Mbed << " snps" << endl;
    FileUtils::AutoGzIfstream fin;
    fin.openOrExit(bedFile, std::ios::in | std::ios::binary);
    uchar header[3];
    fin.read((char *) header, 3);
    if (!fin || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
      cerr << "ERROR: Incorrect first three bytes of bed file: " << bedFile << endl;
      exit(1);
    }

    // read genos + QC snps (and record indiv miss rates)
    vector <int> numMissingPerIndiv(NpreQC);
    uchar *bedLineIn = ALIGNED_MALLOC_UCHARS((Nbed+3)>>2);
    int numSnpsFailedQC = 0;
    uint64 m = 0;
    for (uint64 mbed = 0; mbed < Mbed; mbed++) {
      uchar *genoLine = genosPreQC + m*NpreQC;
      readBedLine(fin, bedLineIn, genoLine, bedIndivRemoved, !bedSnpExcluded[mbed]);
      if (!bedSnpExcluded[mbed]) {
	snpsPreQC[m].MAF = computeMAF(genoLine, NpreQC);
	snpsPreQC[m].miss = computeSnpMissing(genoLine, NpreQC);
	snpsPreQC[m].passQC = snpsPreQC[m].miss <= maxMissingPerSnp;
	if (snpsPreQC[m].passQC) {
	  for (uint64 n = 0; n < NpreQC; n++)
	    numMissingPerIndiv[n] += genoLine[n] == 9;
	}
	else {
	  if (numSnpsFailedQC < 5)
	    cout << "Filtering snp " << snpsPreQC[m].ID << ": "
		 << snpsPreQC[m].miss << " missing" << endl;
	  numSnpsFailedQC++;
	}
	m++;
      }
    }
    ALIGNED_FREE(bedLineIn);

    if (numSnpsFailedQC)
      cout << "Filtered " << numSnpsFailedQC << " SNPs with > " << maxMissingPerSnp << " missing"
	   << endl;

    if (!fin || fin.get() != EOF) {
      cerr << "ERROR: Wrong file size or reading error for bed file: "
	   << bedFile << endl;
      exit(1);
    }
    fin.close();

    // select subset of snps passing QC
    for (uint64 m = 0; m < MpreQC; m++)
      if (snpsPreQC[m].passQC)
	snps.push_back(snpsPreQC[m]);
    M = snps.size();

    // QC indivs for missingness
    int numIndivsFailedQC = 0;
    for (uint64 n = 0; n < NpreQC; n++) {
      indivsPreQC[n].miss = numMissingPerIndiv[n] / (double) M;
      indivsPreQC[n].passQC = indivsPreQC[n].miss <= maxMissingPerIndiv;
      if (!indivsPreQC[n].passQC) {
	if (numIndivsFailedQC < 5)
	  cout << "Filtering indiv " << indivsPreQC[n].famID << " " << indivsPreQC[n].indivID
	       << ": " << numMissingPerIndiv[n] << "/" << M << " missing" << endl;
	numIndivsFailedQC++;
      }
    }
    if (numIndivsFailedQC)
      cout << "Filtered " << numIndivsFailedQC << " indivs with > " << maxMissingPerIndiv
	   << " missing" << endl;

    // select subset of indivs passing QC
    for (uint64 n = 0; n < NpreQC; n++)
      if (indivsPreQC[n].passQC)
	indivs.push_back(indivsPreQC[n]);
    N = indivs.size();

    cout << endl;
    cout << "Total post-QC indivs: N = " << N << endl;
    cout << "Total post-QC SNPs: M = " << M << endl;
    if (N < 2) {
      cerr << "ERROR: N>=2 individuals are required" << endl;
      exit(1);
    }
    if (M == 0) {
      cerr << "ERROR: At least 1 SNP is required" << endl;
      exit(1);
    }

    cout << "MAF spectrum: " << endl;
    const double mafBounds6[7] = {0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.500001};
    vector <int> mafBinCounts(6);
    for (uint64 m = 0; m < M; m++)
      for (int b = 0; b < 6; b++)
	if (mafBounds6[b] <= snps[m].MAF && snps[m].MAF < mafBounds6[b+1])
	  mafBinCounts[b]++;
    for (int b = 0; b < 6; b++)
      printf("    %2.0f-%2.0f%%: %7d\n", 100*mafBounds6[b], 100*mafBounds6[b+1], mafBinCounts[b]);

    printRange();

    if (cMmax == 0) {
      cMmax = std::min(1.0, std::max(N / 1e5, 0.25));
      cout << "Auto-selecting --maxBlockLen: " << cMmax << " cM" << endl;
    }

    vector <bool> nullVec;
    buildGenoBits(genosPreQC, nullVec, cMmax);

    ALIGNED_FREE(genosPreQC);
  }

  /**
   * reads genotypes from VCF/BCF file
   * does not save indiv info (will be reread from VCF during output)
   * only saves chrom, physpos, genpos in snp info (rest will be reread from VCF during output)
   * allocates memory, reads genotypes, and restricts to region if specified; does not do QC
   */
  void GenoData::initVcf(const string &vcfFile, const int inputChrom, const int chromX,
			 double bpStart, double bpEnd, const string &geneticMapFile,
			 bool noMapCheck, double cMmax) {

    htsFile *fin = hts_open(vcfFile.c_str(), "r");

    if (fin == NULL) {
      cerr << "ERROR: Could not open " << vcfFile << " for reading" << endl;
      exit(1);
    }
    htsThreadPool p = {hts_tpool_init(omp_get_max_threads()),
                       0};
    hts_set_thread_pool(fin, &p);

    bcf_hdr_t *hdr = bcf_hdr_read(fin);
    bcf1_t *rec = bcf_init1();
    int mgt = 0, *gt = NULL;

    NpreQC = bcf_hdr_nsamples(hdr);
    cout << "Reading genotypes for N = " << NpreQC << " samples" << endl;
    vector <bool> genos2bit;

    int wantChrom = inputChrom; // might be 0; if so, update
    // read genos; save chrom and physpos for each SNP
    while (bcf_read(fin, hdr, rec) >= 0) {
      // check CHROM
      int chrom = StringUtils::bcfNameToChrom(bcf_hdr_id2name(hdr, rec->rid), 1, chromX);
      if (wantChrom == 0) wantChrom = chrom; // if --chrom was not specified, set to first
      if (chrom != wantChrom) { // only allow multi-chrom file if --chrom has been specified
	if (inputChrom == 0) {
	  cerr << "ERROR: File contains data for >1 chromosome; specify one with --chrom" << endl;
	  exit(1);
	}
	else
	  continue;
      }

      // check if POS is within selected region
      int bp = rec->pos+1;
      if (!(bpStart <= bp && bp <= bpEnd)) continue;

      // check for multi-allelics (TODO: ignore with warning and don't phase in output)
      if (rec->n_allele > 2) {
	cerr << "ERROR: Multi-allelic site found (i.e., ALT contains multiple alleles)" << endl;
	cerr << "       Either drop or split (bcftools norm -m) multi-allelic variants" << endl;
	exit(1);
      }

      // add chrom and bp to SNP list
      SnpInfoX snp; snp.chrom = chrom; snp.physpos = bp;
      snpsPreQC.push_back(snp);

      // read genotypes
      int ngt = bcf_get_genotypes(hdr, rec, &gt, &mgt);
      if (ngt != 2 * (int) NpreQC) {
	cerr << "ERROR: Samples are not diploid" << endl;
	exit(1);
      }
      for (int i = 0; i < (int) NpreQC; i++) {
	int ploidy = 2;
	int *ptr = gt + i*ploidy;

	uchar geno = 0;
	bool missing = false;
	for (int j = 0; j < ploidy; j++) {
	  if ( ptr[j]==bcf_int32_vector_end ) {
	    if (j == 0) {
	      cerr << "ERROR: ptr[0]==bcf_int32_vector_end... zero ploidy?" << endl;
	      exit(1);
	    }
	    else { // 2nd of ploidy==2 genotypes is set to bcf_int32_vector_end => haploid
	      if ( missing ) continue;  // missing diploid genotype can be written in VCF as "."
	      else if (wantChrom == chromX) // X chromosome => haploid ok
		geno *= 2; // encode as diploid homozygote
	      else {
		cerr << "ERROR: Haploid genotype found" << endl;
		exit(1);
	      }
	    }
	  }
	  else {
	    if ( bcf_gt_is_missing(ptr[j]) ) // missing allele
	      missing = true;
	    else
	      geno += bcf_gt_allele(ptr[j]); // 0=REF, 1=ALT (multi-allelics prohibited)
	  }
	}

	if (missing) geno = 3;
	genos2bit.push_back(geno&1);
	genos2bit.push_back(geno>>1);
      }
    }

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fin);
    hts_tpool_destroy(p.pool);

    cout << "Read M = " << snpsPreQC.size() << " variants" << endl;
    processMap(snpsPreQC, geneticMapFile, noMapCheck); // modify snpsPreQC

    // don't perform QC; use all SNPs
    MpreQC = snpsPreQC.size();
    for (uint64 m = 0; m < MpreQC; m++)
      snpsPreQC[m].passQC = true;
    snps = snpsPreQC;
    M = snps.size();

    // don't perform QC; use all samples
    N = NpreQC;

    printRange();

    if (cMmax == 0) {
      cMmax = std::min(1.0, std::max(N / 1e5, 0.25));
      cout << "Auto-selecting --maxBlockLen: " << cMmax << " cM" << endl;
    }

    buildGenoBits(NULL, genos2bit, cMmax);
  }

  GenoData::~GenoData() {
    ALIGNED_FREE(seg64logPs);
    ALIGNED_FREE(genoBits);
  }

  /**
   * assumes Nbed = bedIndivRemoved.size()
   * reads (Nbed+3)>>2 bytes into bedLineIn
   * stores sum(!bedIndivRemoved) bytes into genoLine if loadGenoLine == true
   */
  void GenoData::readBedLine(FileUtils::AutoGzIfstream &fin, uchar bedLineIn[], uchar genoLine[],
			     vector <bool> &bedIndivRemoved, bool storeGenoLine) {
    uint64 Nbed = bedIndivRemoved.size();
    fin.read((char *) bedLineIn, (Nbed+3)>>2);
    if (storeGenoLine) {
      const uchar bedToGeno[4] = {2, 9, 1, 0};
      uint64 n = 0;
      for (uint64 nbed = 0; nbed < Nbed; nbed++)
	if (!bedIndivRemoved[nbed])
	  genoLine[n++] = bedToGeno[(bedLineIn[nbed>>2]>>((nbed&3)<<1))&3];
    }
  }
  double GenoData::computeAlleleFreq(const uchar genoLine[], uint64 genoN) {
    double sum = 0; int num = 0;
    for (uint64 n = 0; n < genoN; n++)
      if (genoLine[n] != 9) {
	sum += genoLine[n];
	num++;
      }
    return 0.5 * sum / num;
  }
  double GenoData::computeMAF(const uchar genoLine[], uint64 genoN) {
    double alleleFreq = computeAlleleFreq(genoLine, genoN);
    return std::min(alleleFreq, 1.0-alleleFreq);
  }
  double GenoData::computeSnpMissing(const uchar genoLine[], uint64 genoN) {
    double sum = 0; int num = 0;
    for (uint64 n = 0; n < genoN; n++) {
      sum += (genoLine[n] == 9);
      num++;
    }
    return sum / num;
  }

  const vector <SnpInfoX> &GenoData::getSnps(void) const { return snps; }
  uint64 GenoData::getN(void) const { return N; }
  uint64 GenoData::getMseg64(void) const { return Mseg64; }
  const uint64_masks *GenoData::getGenoBits(void) const { return genoBits; }
  vector <vector <double> > GenoData::getSeg64cMvecs(void) const { return seg64cMvecs; }
  const AlleleFreqs *GenoData::getSeg64logPs(void) const { return seg64logPs; }
  IndivInfoX GenoData::getIndiv(uint64 n) const { return indivs[n]; }
  const vector <IndivInfoX> &GenoData::getIndivs(void) const { return indivs; }
  const vector <bool> &GenoData::getIsFlipped64j(void) const { return isFlipped64j; }

};
