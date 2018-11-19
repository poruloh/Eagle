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
#include <fstream>
#include <map>
#include <set>

#include "omp.h"

#include <boost/version.hpp>

#include "Eagle.hpp"
#include "EagleParams.hpp"
#include "GenoData.hpp"
#include "HapHedge.hpp"
#include "SyncedVcfData.hpp"
#include "Timer.hpp"
#include "StringUtils.hpp"
#include "Version.hpp"

using namespace EAGLE;
using namespace std;

//#define LOCAL_TEST
//#define OLD_IMP_MISSING

void adjustHistFactor(double &histFactor, double hetRate, double snpRate) {
  if (histFactor == 0) {
    if (snpRate * hetRate == 0)      
      histFactor = 1;
    else {
      // compute how far (in cM) default 100 hets will typically span
      double cMdefaultHist = 100 / snpRate / hetRate;
      printf("Typical span of default 100-het history length: %.2f cM\n", cMdefaultHist);
      const double cMminHist = 1.0;
      histFactor = max(1.0, min(20.0, cMminHist / cMdefaultHist));
      printf("Setting --histFactor=%.2f\n", histFactor);
      if (histFactor != 1)
	printf("Typical span of %d-het history length: %.2f cM\n", (int) (100*histFactor),
	       cMdefaultHist*histFactor);
      cout << flush;
    }
  }
}

void phaseWithRef(EagleParams &params, Timer &timer, double t0, int argc, char **argv) {

  string tmpFile = params.outPrefix + ".unphased." + params.vcfOutSuffix;
  string outFile = params.outPrefix + "." + params.vcfOutSuffix;
  vector < vector < pair <int, int> > > conPSall; double snpRate;
  vector <bool> isTmpPhased;
  SyncedVcfData vcfData(params.vcfRef, params.vcfTarget, params.vcfExclude, params.allowRefAltSwap,
			params.chrom, params.chromX, params.bpStart-params.bpFlanking,
			params.bpEnd+params.bpFlanking, params.geneticMapFile,
			params.cMmax==0 ? 1 : params.cMmax, tmpFile, params.vcfWriteMode,
			params.usePS, conPSall, snpRate, params.outputUnphased, isTmpPhased);
  cout << endl << "Time for reading input: " << timer.update_time() << " sec" << endl << endl;

  Eagle eagle(vcfData.getNref(), vcfData.getNtarget(), vcfData.getMseg64(),
	      vcfData.getGenoBits(), vcfData.getSeg64cMvecs(), params.pErr);

  double hetRate = eagle.computeHetRate();
  cout << "Fraction of heterozygous genotypes: " << hetRate << endl;
  if (params.usePBWT) adjustHistFactor(params.histFactor, hetRate, snpRate);
  
  uint64 Nref = vcfData.getNref(), Ntarget = vcfData.getNtarget();
  int iters = params.pbwtIters;
  if (iters == 0) {
    if (Ntarget < Nref/2)
      iters = 1;
    else if (Ntarget < 2*Nref)
      iters = 2;
    else
      iters = 3;
    cout << endl << "Auto-selecting number of phasing iterations: setting --pbwtIters to "
	 << iters << endl << endl;
  }
  bool internalImpMissing = !params.noImpMissing || iters > 1; // internally impute missing

  cout << endl << "BEGINNING PHASING" << endl;

  vector <float> confs(Ntarget);

  for (int iter = 1; iter <= iters; iter++) {
    double t23 = timer.get_time(); timer.update_time();
    double timeMN2 = 0, timeImpMissing = 0;
    cout << endl << "PHASING ITER " << iter << " OF " << iters << endl << endl;
    eagle.initRefIter(iter);

    if (params.usePBWT) { // run PBWT algorithm
      HapBitsT *hapBitsTptr = NULL;
      HapHedge *hapHedgePtr = NULL;
#ifdef OLD_IMP_MISSING
      if (internalImpMissing) {
	cout << "Making HapHedge" << endl;
	hapBitsTptr = new HapBitsT(eagle.getHaploBitsT(), 2*eagle.getNlib(2+iter),
				   eagle.getMseg64(), eagle.getMaskSnps64j());
	int skip = 25;
	hapHedgePtr = new HapHedge(*hapBitsTptr, skip);
	cout << "Built PBWT on " << hapBitsTptr->getNhaps() << " haplotypes" << endl;
	cout << "Time for HapHedge: " << timer.update_time() << endl;
      }
#endif
      cout << endl << "Phasing target samples" << endl;
      int numPhased = 0; const int dots = 80;
#pragma omp parallel for reduction(+:timeImpMissing) schedule(dynamic, 4)
      for (uint i = Nref; i < Nref+Ntarget; i++) {
	int nF1 = -1, nF2 = -1;
	if (params.trioCheck) {
	  if ((i-Nref)%3 == 0) { // child
	    nF1 = i+1; nF2 = i+2;
	  }
	  else if ((i-Nref)%3 == 1) { // parent 1
	    nF1 = -(i-1); nF2 = i+1;
	  }
	  else { // parent 2
	    nF1 = -(i-2); nF2 = i-1;
	  }
	}
	confs[i-Nref] = eagle.runPBWT
	  (i, nF1, nF2, params.Kpbwt/(iter<iters?2:1), params.expectIBDcM, params.histFactor,
	   iter==iters, iter>1, internalImpMissing, params.usePS, conPSall[i-Nref],
	   params.chrom==params.chromX);
#ifdef OLD_IMP_MISSING
	if (internalImpMissing) {
	  Timer tim;
	  eagle.imputeMissing(*hapHedgePtr, i);
	  timeImpMissing += tim.update_time();
	}
#endif
#pragma omp critical(NUM_PHASED)
	{
	  numPhased++;
	  int newDots = numPhased * dots / Ntarget - (numPhased-1) * dots / Ntarget;
	  if (newDots) cout << string(newDots, '.') << flush;
	}
      }
#ifdef OLD_IMP_MISSING
      if (internalImpMissing) {
	delete hapHedgePtr;
	delete hapBitsTptr;
      }
#endif
    }
    else { // run LRP algorithm
      cout << "Building hash tables" << endl;
      eagle.buildHashTables(2+iter, 0, params.seed); // in ref mode, first iter is 3
      cout << " (time: " << timer.update_time() << ")" << endl;

      cout << endl << "Phasing target samples" << endl;
#pragma omp parallel for reduction(+:timeMN2) schedule(dynamic, 4)
      for (uint i = Nref; i < Nref+Ntarget; i++)
	timeMN2 += eagle.runHMM(i, -1, -1, 3, params.beamWidth4, params.maxHapStates);
    }

    cout << endl << "Time for phasing iter " << iter << ": " << (timer.get_time()-t23) << endl;
    if (!params.usePBWT)
      cout << "Time for phasing iter " << iter << " MN^2: " << timeMN2 / params.numThreads
	   << endl;
#ifdef OLD_IMP_MISSING
    else if (internalImpMissing)
      cout << "Time for phasing iter " << iter << " impMissing: "
	   << timeImpMissing / params.numThreads << endl;
#endif
  }

  /***** FINAL OUTPUT *****/

  timer.update_time();
  cout << "Writing " << params.vcfOutSuffix << " output to " << outFile << endl;
  eagle.writeVcf(tmpFile, isTmpPhased, outFile, params.chromX, params.bpStart, params.bpEnd,
		 params.vcfWriteMode, params.noImpMissing, params.keepMissingPloidyX, argc, argv);
  cout << "Time for writing output: " << timer.update_time() << endl;

  cout << "Total elapsed time for analysis = " << (timer.get_time()-t0) << " sec" << endl;

  cout << endl;
  cout << "Mean phase confidence of each target individual:" << endl;
  cout << "ID" << "\t" << "PHASE_CONFIDENCE" << endl;
  for (uint i = Nref; i < Nref+Ntarget; i++) {
    cout << vcfData.getTargetID(i-Nref) << "\t" << confs[i-Nref] << endl;
  }
}

int main(int argc, char *argv[]) {

  Timer timer; double t0 = timer.get_time();

  cout << "                      +-----------------------------+" << endl;
  cout << "                      |                             |" << endl;
  cout << "                      |   Eagle v";
  printf("%-19s|\n", EAGLE_VERSION);
  cout << "                      |   ";
  printf("%-26s|\n", EAGLE_VERSION_DATE);
  cout << "                      |   Po-Ru Loh                 |" << endl;
  cout << "                      |                             |" << endl;
  cout << "                      +-----------------------------+" << endl;
  cout << endl;

  cout << "Copyright (C) 2015-2018 Harvard University." << endl;
  cout << "Distributed under the GNU GPLv3+ open source license." << endl << endl;

  //cout << "Boost version: " << BOOST_LIB_VERSION << endl;
  //cout << endl;

  printf("Command line options:\n\n");
  printf("%s ", argv[0]);
  for (int i = 1; i < argc; i++) {
    if (strlen(argv[i]) >= 2 && argv[i][0] == '-' && argv[i][1] == '-')
      printf("\\\n    ");
    bool hasSpace = false;
    for (uint j = 0; j < strlen(argv[i]); j++)
      if (isspace(argv[i][j]))
	hasSpace = true;
    if (hasSpace) {
      if (argv[i][0] == '-') {
	bool foundEquals = false;
	for (uint j = 0; j < strlen(argv[i]); j++) {
	  printf("%c", argv[i][j]);
	  if (argv[i][j] == '=' && !foundEquals) {
	    printf("\"");
	    foundEquals = true;
	  }
	}
	printf("\" ");
      }
      else
	printf("\"%s\" ", argv[i]);
    }
    else
      printf("%s ", argv[i]);
  }
  cout << endl << endl;

  EagleParams params;
  if (!params.processCommandLineArgs(argc, argv)) {
    cerr << "Aborting due to error processing command line arguments" << endl;
    exit(1);
  }

  cout << "Setting number of threads to " << params.numThreads << endl;
  omp_set_num_threads(params.numThreads);

  if (!params.vcfRef.empty()) { // use reference haplotypes
    phaseWithRef(params, timer, t0, argc, argv);
    return 0;
  }
  
  cout << endl << "=== Reading genotype data ===" << endl << endl;

  GenoData genoData;
  if (!params.vcfFile.empty())
    genoData.initVcf(params.vcfFile, params.chrom, params.chromX, params.bpStart, params.bpEnd,
		     params.geneticMapFile, params.noMapCheck, params.cMmax);
  else 
    genoData.initBed(params.famFile, params.bimFile, params.bedFile, params.chrom, params.bpStart,
		     params.bpEnd, params.geneticMapFile, params.excludeFiles, params.removeFiles,
		     params.maxMissingPerSnp, params.maxMissingPerIndiv, params.noMapCheck,
		     params.cMmax);

  if (!params.usePBWT) { // Eagle v1 algorithm
    params.pbwtOnly = false; // should already be false
    params.runStep2 = 1;
  }
  else { // PBWT algorithm
    // if SNP density is >1000 SNPs/Mb, don't run Steps 1+2 (even if --pbwtOnly is not set)
    int bpSpan = genoData.getSnps().back().physpos - genoData.getSnps()[0].physpos;
    double snpsPerMb = genoData.getSnps().size() / (bpSpan*1e-6);
    if (snpsPerMb > 1000) params.pbwtOnly = true;

    if (genoData.getN() > 200000U) params.pbwtOnly = true;

    if (params.pbwtOnly) params.runStep2 = 0;

    // if --runStep2 hasn't yet been set, SNP density must be low; run Step 2 unless too few chunks
    if (params.runStep2 != 0 && params.runStep2 != 1)
      params.runStep2 = (genoData.getMseg64() >= 3U); // can't run Step 2 with < 3 SNP segments
  }

  vector <double> invLD64j;
  if (!params.pbwtOnly) invLD64j = genoData.computeInvLD64j(1000);
  
  Eagle eagle(genoData.getN(), genoData.getMseg64(), genoData.getGenoBits(),
	      genoData.getSeg64cMvecs(), genoData.getSeg64logPs(), invLD64j, genoData.getIndivs(),
	      genoData.getSnps(), params.maskFile, genoData.getIsFlipped64j(), params.pErr,
	      params.runStep2);
  
  double hetRate = eagle.computeHetRate();
  cout << "Fraction of heterozygous genotypes: " << hetRate << endl;
  if (params.usePBWT) adjustHistFactor(params.histFactor, hetRate, genoData.computeSnpRate());
  
  map <string, pair <string, string> > trioIIDs;
  vector <uint> children, nF1s, nF2s;
  uint N = genoData.getN();
  double timeMN2 = 0;

#ifdef LOCAL_TEST
  {
    ifstream
      //finTrios("/groups/price/poru/HSPH_SVN/data/GERA/phasing/eur.CEU_gt_0.9.trios_indep.fam");
      finTrios("/groups/price/UKBiobank/download/ukb4777_trios.fam");
    if (finTrios) {
      string FID, IID, s1, s2, s3, s4;
      while (finTrios >> FID >> IID >> s1 >> s2 >> s3 >> s4)
	trioIIDs[IID] = make_pair(s1, s2);

      for (uint i = 0; i < N; i++) {
	int nF1 = -1, nF2 = -1;
	if (trioIIDs.find(genoData.getIndiv(i).indivID) == trioIIDs.end()) continue;
	pair <string, string> parents = trioIIDs[genoData.getIndiv(i).indivID];
	for (uint iF = 0; iF < N; iF++)
	  if (genoData.getIndiv(iF).indivID == parents.first ||
	      genoData.getIndiv(iF).indivID == parents.second) {
	    if (nF1 == -1) nF1 = iF;
	    else nF2 = iF;
	  }
	if (nF1 != -1 && nF2 != -1) {
	  children.push_back(i); nF1s.push_back(nF1); nF2s.push_back(nF2);
	}
      }
      cout << "Identified " << children.size() << " trio children" << endl;
    }
  }
#endif

  if (params.iter == 0) {
    if (params.outPrefix == "") {
      cerr << "ERROR: --outPrefix must be specified" << endl;
      exit(1);
    }

    /***** RUN STEP 1 *****/

    if (!params.pbwtOnly) {
      cout << endl << "BEGINNING STEP 1" << endl << endl;
      double t1 = timer.get_time(); timer.update_time(); timeMN2 = 0;

      for (uint att = 0; att < min(9U, (uint) children.size()); att++) // run on trio children
	eagle.findLongDipMatches(children[att], nF1s[att], nF2s[att]);
#pragma omp parallel for reduction(+:timeMN2) schedule(dynamic, 4)
      for (uint i = 0; i < N; i++) {	
	//cout << StringUtils::itos(i)+"\n" << flush;
	timeMN2 += eagle.findLongDipMatches(i, -1, -1).first;
	//cout << StringUtils::itos(-i)+"\n" << flush;
      }

      if (!params.tmpPhaseConfsPrefix.empty())
	eagle.writePhaseConfs(params.tmpPhaseConfsPrefix+".step1.bin");
      cout << "Time for step 1: " << (timer.get_time()-t1) << endl;
      cout << "Time for step 1 MN^2: " << timeMN2 / params.numThreads << endl;
      eagle.outputSE(children, nF1s, nF2s, 1);
    }
    else { // use PBWT only => phase randomly
      cout << endl << "SKIPPED STEP 1" << endl;
#pragma omp parallel for reduction(+:timeMN2) schedule(dynamic, 4)
      for (uint i = 0; i < N; i++)
	eagle.randomlyPhaseTmpHaploBitsT(i);
    }

    if (params.runStep2) { // running step 2 => Step 1 phase confs were saved to phaseConfs
      cout << endl << "Making hard calls" << flush; timer.update_time();
      eagle.makeHardCalls(0, N, params.seed);
      cout << " (time: " << timer.update_time() << ")" << endl << endl;
    }
    else { // not running Step 2 => Step 1 phase calls were saved to tmpHaploBitsT
      eagle.cpTmpHaploBitsT(0, N);
    }

    /***** RUN STEPS 2-4 (STEP 4 = STEP 3a in paper) *****/

    for (int step = 2; step <= (params.usePBWT ? (1+params.runStep2) : 4); step++) {
      cout << endl << "BEGINNING STEP " << step << endl << endl;
      double t23 = timer.get_time(); timer.update_time(); timeMN2 = 0;

      const uint64 numBatches = step == 2 ? 1 : 10;
      const uint64 runBatches = step <= 3 ? numBatches : (uint64) (numBatches * params.fracStep4);
      for (uint64 b = 1; b <= runBatches; b++) {
	cout << "BATCH " << b << " OF " << runBatches << endl;
	cout << "Building hash tables" << endl;
	eagle.buildHashTables(step, b, params.seed);
	cout << " (time: " << timer.update_time() << ")" << endl;

	if (b == 1)
	  for (uint att = 0; att < min(9U, (uint) children.size()); att++) // run on trio children
	    step == 2 ? eagle.findLongHapMatches(children[att], nF1s[att], nF2s[att], step)
	      : eagle.runHMM(children[att], nF1s[att], nF2s[att], step,
			     step==3 ? params.beamWidth3 : params.beamWidth4, params.maxHapStates);
	
	uint iStart = (b-1)*N/numBatches, iEnd = b*N/numBatches;
	cout << endl << "Phasing samples " << (iStart+1) << "-" << iEnd << endl;
#pragma omp parallel for reduction(+:timeMN2) schedule(dynamic, 4)
	for (uint i = iStart; i < iEnd; i++) {	
	  //if (step == 3) cout << StringUtils::itos(i)+"\n" << flush;
	  timeMN2 += step == 2 ? eagle.findLongHapMatches(i, -1, -1, step)
	    : eagle.runHMM(i, -1, -1, step, step==3 ? params.beamWidth3 : params.beamWidth4,
			   params.maxHapStates);
	  //if (step == 3) cout << StringUtils::itos(-i)+"\n" << flush;
	}

	eagle.cpPhaseConfs(iStart, iEnd);
	cout << "Time for phasing batch: " << timer.update_time() << endl;

	cout << endl << "Making hard calls" << flush;
	eagle.makeHardCalls(iStart, iEnd, params.seed + step);
	cout << " (time: " << timer.update_time() << ")" << endl << endl;
      }

      if (!params.tmpPhaseConfsPrefix.empty())
	eagle.writePhaseConfs(params.tmpPhaseConfsPrefix+".step"+StringUtils::itos(step)+".bin");
      cout << "Time for step " << step << ": " << (timer.get_time()-t23) << endl;
      cout << "Time for step " << step << " MN^2: " << timeMN2 / params.numThreads << endl;
      eagle.outputSE(children, nF1s, nF2s, step);
    }

    if (params.usePBWT) { // run PBWT iters
      if (!params.runStep2) // didn't run Step 2 => didn't allocate phaseConfs
	cout << endl << "SKIPPED STEP 2" << endl;
      else // ran Step 2 => allocated phaseConfs; didn't allocate tmpHaploBitsT
	eagle.reallocLRPtoPBWT();

      cout << endl << endl << "BEGINNING STEP 3 (PBWT ITERS)" << endl << endl;
      int iters = params.pbwtIters;
      if (iters == 0) {
	iters = 2 + params.pbwtOnly;
	cout << "Auto-selecting number of PBWT iterations: setting --pbwtIters to "
	     << iters << endl << endl;
      }
      for (int iter = 1; iter <= iters; iter++) {
	cout << endl << "BEGINNING PBWT ITER " << iter << endl << endl;
	double t23 = timer.get_time(); timer.update_time(); double timeImpMissing = 0;

	int skip = 16; int Kpbwt = params.Kpbwt; bool runReverse = true;
	if (iter < iters) { // run rougher computation
	  Kpbwt >>= (iters-iter);
	  skip *= 2;
	  runReverse = false;
	}

	const uint64 numBatches = 10;
	for (uint64 b = 1; b <= numBatches; b++) {
	  cout << "BATCH " << b << " OF " << numBatches << endl;
#ifdef OLD_IMP_MISSING
	  cout << endl << "Making HapHedge" << endl;
	  HapBitsT hapBitsT(eagle.getHaploBitsT(), 2*eagle.getNlib(2+iter),
			    eagle.getMseg64(), eagle.getMaskSnps64j());
	  HapHedge hapHedge(hapBitsT, skip);
	  cout << "Built PBWT on " << hapBitsT.getNhaps() << " haplotypes" << endl;
	  cout << "Time for HapHedge: " << timer.update_time() << endl;
#endif
	  if (b == 1)
	    for (uint att = 0; att < min(9U, (uint) children.size()); att++) // run on trios
	      eagle.runPBWT(children[att], nF1s[att], nF2s[att], Kpbwt, params.expectIBDcM,
			    params.histFactor, runReverse, true, false, params.chrom==params.chromX);

	  uint iStart = (b-1)*N/numBatches, iEnd = b*N/numBatches;
	  cout << endl << "Phasing samples " << (iStart+1) << "-" << iEnd << endl;
#pragma omp parallel for reduction(+:timeImpMissing) schedule(dynamic, 4)
	  for (uint i = iStart; i < iEnd; i++) {	
	    eagle.runPBWT(i, -1, -1, Kpbwt, params.expectIBDcM, params.histFactor, runReverse,
			  true, true, params.chrom==params.chromX);
#ifdef OLD_IMP_MISSING
	    Timer tim;
	    eagle.imputeMissing(hapHedge, i);
	    timeImpMissing += tim.update_time();
#endif
	  }

	  eagle.cpTmpHaploBitsT(iStart, iEnd);
	  cout << "Time for phasing batch: " << timer.update_time() << endl << endl;
	}

	cout << "Time for PBWT iter " << iter << ": " << (timer.get_time()-t23) << endl;
#ifdef OLD_IMP_MISSING
	cout << "Time for PBWT iter " << iter << " impMissing: "
	     << timeImpMissing / params.numThreads << endl;
#endif
	//eagle.outputSE(children, nF1s, nF2s, step); // currently requires phaseConfs
      }
    }

    /***** FINAL OUTPUT *****/

    if (!params.vcfFile.empty()) {
      string outFile = params.outPrefix + "." + params.vcfOutSuffix;
      cout << "Writing " << params.vcfOutSuffix << " output to " << outFile << endl;
      eagle.writeVcfNonRef(params.vcfFile, outFile, params.chrom, params.chromX, params.bpStart,
			   params.bpEnd, params.vcfWriteMode, params.noImpMissing,
			   params.keepMissingPloidyX, argc, argv);
    }
    else {
      cout << "Writing .haps.gz and .sample output" << endl; timer.update_time();
      eagle.writeHapsGzSample(params.outPrefix);
    }
    cout << "Time for writing output: " << timer.update_time() << endl;
  }
  else if (params.iter == 1) { // PERFORM 1ST-ITER PHASING FOR A SMALL SUBSET ONLY (FOR TESTING)
    double t1 = timer.get_time(); timer.update_time(); timeMN2 = 0;
    map <int, int> longestFreqs; int tot = 0;
    int att = 0, maxAtt = 10;
    for (uint i = 0; i < N; i++) {
      int nF1 = -1, nF2 = -1;

      if (trioIIDs.find(genoData.getIndiv(i).indivID) == trioIIDs.end()) continue;
      cout << "Testing n0 = " << i << ": " << genoData.getIndiv(i).famID << " "
	   << genoData.getIndiv(i).indivID << endl;
      pair <string, string> parents = trioIIDs[genoData.getIndiv(i).indivID];
      for (uint iF = 0; iF < N; iF++)
	if (genoData.getIndiv(iF).indivID == parents.first ||
	    genoData.getIndiv(iF).indivID == parents.second) {
	  cout << "Parent n1: " << iF << endl;
	  if (nF1 == -1) nF1 = iF;
	  else nF2 = iF;
	}
      eagle.checkTrioErrorRate(i, nF1, nF2);

      pair <double, vector <double> > ret = eagle.findLongDipMatches(i, nF1, nF2);
      timeMN2 += ret.first;
      for (uint j = 0; j < ret.second.size(); j++) {
	longestFreqs[(int) ret.second[j]]++;
	tot++;
      }
      att++;
      if (att == maxAtt) break;
    }
    int cum = 0;
    for (map <int, int>::iterator it = longestFreqs.begin(); it != longestFreqs.end(); it++) {
      cum += it->second;
      cout << it->first << " cM: " << it->second << "   cum: " << (double) cum/tot << endl;
    }
    cout << "Time for step 1: " << (timer.get_time()-t1) << endl;
    cout << "Time for step 1 MN^2: " << timeMN2 / params.numThreads << endl;
  }
  else { // iter > 1 (FOR TESTING)
    cout << "Reading phase confidences" << endl;
    eagle.readPhaseConfs(/*"test_UKBiobank/eagle_305_chr10_small_iter2_bsub.tmpPhaseConfs.bin"*/params.tmpPhaseConfsPrefix+".step"+StringUtils::itos(params.iter-1)+".bin");
    timer.update_time();
    cout << "Making hard calls" << endl;
    eagle.makeHardCalls(0, N, params.seed);
    cout << "Time for hard calls: " << timer.update_time() << endl;

    /* global PBWT
    cout << "Making forward HapHedge" << endl;
    HapBitsT hapBitsFwdT(eagle.getHaploBitsT(), 2*eagle.getNlib(params.iter), eagle.getMseg64(),
			 eagle.getMaskSnps64j());
    int skip = 1;
    HapHedge hapHedgeFwd(hapBitsFwdT, skip);
    cout << "Making backward HapHedge" << endl;
    HapBitsT hapBitsBwdT(hapBitsFwdT, -1);
    HapHedge hapHedgeBwd(hapBitsBwdT, skip);
    cout << "Time for HapHedge: " << timer.update_time() << endl;
    */
#define USE_PBWT
#ifndef USE_PBWT
    cout << "Building hash tables" << endl;
    eagle.buildHashTables(params.iter, 1, params.seed);
    cout << endl << "Time for hash tables: " << timer.update_time() << endl;
#else
    eagle.reallocLRPtoPBWT();
#endif
    int att = 0, maxAtt = N;
    for (uint i = 0; i < N; i++) {
      int nF1 = -1, nF2 = -1;

      if (trioIIDs.find(genoData.getIndiv(i).indivID) == trioIIDs.end()) continue;
      cout << "Testing n0 = " << i << ": " << genoData.getIndiv(i).famID << " "
	   << genoData.getIndiv(i).indivID << endl;
      pair <string, string> parents = trioIIDs[genoData.getIndiv(i).indivID];
      for (uint iF = 0; iF < N; iF++)
	if (genoData.getIndiv(iF).indivID == parents.first ||
	    genoData.getIndiv(iF).indivID == parents.second) {
	  cout << "Parent n1: " << iF << endl;
	  if (nF1 == -1) nF1 = iF;
	  else nF2 = iF;
	}
      eagle.checkTrioErrorRate(i, nF1, nF2);
#ifdef USE_PBWT
      eagle.runPBWT(i, nF1, nF2, params.Kpbwt, params.expectIBDcM, params.histFactor, true, false,
		    !params.noImpMissing, params.chrom==params.chromX);
#else
      timeMN2 += params.iter == 2 ? eagle.findLongHapMatches(i, nF1, nF2, params.iter)
	: eagle.runHMM(i, nF1, nF2, params.iter,
		       params.iter==3 ? params.beamWidth3 : params.beamWidth4,
		       params.maxHapStates/*, &hapHedgeFwd, &hapHedgeBwd*/);
#endif
      att++; if (att == maxAtt) break;
    }

    cout << "Time for step minus init: " << timer.update_time() << endl;
    cout << "Time for MN^2: " << timeMN2 / (params.iter == 0 ? params.numThreads : 1) << endl;
#ifdef RDTSC_TIMING
    printf("%.1f%% of time in dip-hap\n", 100*eagle.diphapTicks / (double) eagle.totTicks);
    printf("%.1f%% of time in ext\n", 100*eagle.extTicks / (double) eagle.totTicks);
    printf("%.1f%% of time in LSH\n", 100*eagle.lshTicks / (double) eagle.totTicks);
    printf("%.1f%% of time in LSH hit checks\n", 100*eagle.lshCheckTicks / (double) eagle.totTicks);
    printf("%.1f%% of time in DP\n", 100*eagle.dpTicks / (double) eagle.totTicks);
    printf("  rel %.1f%% in sort\n", 100*eagle.dpSortTicks / (double) eagle.dpTicks);
    printf("  rel %.1f%% in update\n", 100*eagle.dpUpdateTicks / (double) eagle.dpTicks);
    printf("  rel %.1f%% in computeStatic\n", 100*eagle.dpStaticTicks / (double) eagle.dpTicks);
    printf("  rel %.1f%% in computeSwitch\n", 100*eagle.dpSwitchTicks / (double) eagle.dpTicks);
    printf("%.1f%% of time in blip fix\n", 100*eagle.blipFixTicks / (double) eagle.totTicks);
    printf("  rel %.1f%% in LSH\n", 100*eagle.blipLshTicks / (double) eagle.blipFixTicks);
    printf("  rel %.1f%% in popcount\n", 100*eagle.blipPopTicks / (double) eagle.blipFixTicks);
    printf("  rel %.1f%% in vote update\n", 100*eagle.blipVoteTicks / (double) eagle.blipFixTicks);
    cout << "Number of update calls: " << eagle.dpUpdateCalls << endl;
#endif
  }

  cout << "Total elapsed time for analysis = " << (timer.get_time()-t0) << " sec" << endl;
}
