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
#include <cstdlib>

#include "StringUtils.hpp"
#include "FileUtils.hpp"
#include "EagleParams.hpp"

#include <boost/program_options.hpp>

namespace EAGLE {

  using std::vector;
  using std::string;
  using std::cout;
  using std::cerr;
  using std::endl;

  // populates members; error-checks
  bool EagleParams::processCommandLineArgs(int argc, char *argv[]) {

    vector <string> removeFileTemplates, excludeFileTemplates;
    string chromStr; // allow "X"

    namespace po = boost::program_options;

    po::options_description commonOptions;
    commonOptions.add_options()
      ("geneticMapFile", po::value<string>(&geneticMapFile)->required(),
       "HapMap genetic map provided with download: tables/genetic_map_hg##.txt.gz")
      // "chr pos rate(cM/Mb) map(cM)"
      ("outPrefix", po::value<string>(&outPrefix)->required(), "prefix for output files")
      ("numThreads", po::value<int>(&numThreads)->default_value(1),
       "number of computational threads")
      ;

    po::options_description nonRefMode
	("Input options for phasing without a reference");
    nonRefMode.add_options()
      // genotype data parameters
      ("bfile", po::value<string>(), "prefix of PLINK .fam, .bim, .bed files")
      ("bfilegz", po::value<string>(), "prefix of PLINK .fam.gz, .bim.gz, .bed.gz files")
      ("fam", po::value<string>(&famFile),
       "PLINK .fam file (note: file names ending in .gz are auto-decompressed)")
      ("bim", po::value<string>(&bimFile), "PLINK .bim file")
      ("bed", po::value<string>(&bedFile), "PLINK .bed file")
      ("vcf", po::value<string>(&vcfFile), "[compressed] VCF/BCF file containing input genotypes")
      ("remove", po::value< vector <string> >(&removeFileTemplates),
       "file(s) listing individuals to ignore (no header; FID IID must be first two columns)")
      ("exclude", po::value< vector <string> >(&excludeFileTemplates),
       "file(s) listing SNPs to ignore (no header; SNP ID must be first column)")
      ("maxMissingPerSnp", po::value<double>(&maxMissingPerSnp)->default_value(0.1, "0.1"),
       "QC filter: max missing rate per SNP")
      ("maxMissingPerIndiv", po::value<double>(&maxMissingPerIndiv)->default_value(0.1, "0.1"),
       "QC filter: max missing rate per person")      
      ;

    po::options_description refMode
      ("Input/output options for phasing using a reference panel");
    refMode.add_options()
      ("vcfRef", po::value<string>(&vcfRef),
       "tabix-indexed [compressed] VCF/BCF file for reference haplotypes")
      ("vcfTarget", po::value<string>(&vcfTarget),
       "tabix-indexed [compressed] VCF/BCF file for target genotypes")
      ("vcfOutFormat", po::value<string>(&vcfOutFormat)->default_value("."),
       "b|u|z|v: compressed BCF (b), uncomp BCF (u), compressed VCF (z), uncomp VCF (v)")
      ("noImpMissing", "disable imputation of missing target genotypes (. or ./.)")
      ("allowRefAltSwap", "allow swapping of REF/ALT in target vs. ref VCF")
      ("outputUnphased", "output unphased sites (target-only, multi-allelic, etc.)")
      ("keepMissingPloidyX",
       "assume missing genotypes have correct ploidy (.=haploid, ./.=diploid)")
      ("vcfExclude", po::value<string>(&vcfExclude),
       "tabix-indexed [compressed] VCF/BCF file containing variants to exclude from phasing")
      ;

    po::options_description bothModes("Region selection options");
    bothModes.add_options()
      ("chrom", po::value<string>(&chromStr)->default_value("0"),
       "chromosome to analyze (if input has many)")
      ("bpStart", po::value<double>(&bpStart)->default_value(0),
       "minimum base pair position to analyze")
      ("bpEnd", po::value<double>(&bpEnd)->default_value(1e9, "1e9"),
       "maximum base pair position to analyze")
      ("bpFlanking", po::value<double>(&bpFlanking)->default_value(0),
       "(ref-mode only) flanking region to use during phasing but discard in output")
      ;

    po::options_description algOptions("Algorithm options");
    algOptions.add_options()
      ("Kpbwt", po::value<int>(&Kpbwt)->default_value(10000),
       "number of conditioning haplotypes") // TODO: throw error if set in --v1 mode
      ("pbwtIters", po::value<int>(&pbwtIters)->default_value(0),
       "number of PBWT phasing iterations (0=auto)")
      ("expectIBDcM", po::value<double>(&expectIBDcM)->default_value(2, "2.0"),
       "expected length of haplotype copying (cM)")
      ("histFactor", po::value<double>(&histFactor)->default_value(0, "0"),
       "history length multiplier (0=auto)")
      ("genoErrProb", po::value<double>(&pErr)->default_value(0.003, "0.003"),
       "estimated genotype error probability")
      ("pbwtOnly", "in non-ref mode, use only PBWT iters (automatic for sequence data)")
      ("v1", "use Eagle1 phasing algorithm (instead of default Eagle2 algorithm)")
      ;

    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("help,h", "print help message with typical options")
      
      // experimental options
      ("usePS", po::value<int>(&usePS)->default_value(0),
       "use FORMAT:PS phase constraints in target VCF: 1=soft, 2=harder")
      ("runStep2", po::value<int>(&runStep2)->default_value(-1),
       "enable/disable Step 2 of non-ref algorithm (-1=auto)")
      ("chromX", po::value<int>(&chromX)->default_value(23), "maximum chromosome number (chrX)")

      // Eagle1 advanced options
      ("v1fast", "Eagle1 fast mode: --maxBlockLen=0.3, --maxStatePairsStep4=100, --fracStep4=0.5")
      ("maxBlockLen", po::value<double>(&cMmax)->default_value(0),
       "max length (in cM units) of a SNP block; increase to trade accuracy for speed (0=auto)")
      ("maxStatePairsStep3", po::value<int>(&beamWidth3)->default_value(100),
       "maximum state pairs per position in dynamic programming (HMM-like) search (step 3)")
      ("maxStatePairsStep4", po::value<int>(&beamWidth4)->default_value(200),
       "maximum state pairs per position in dynamic programming (HMM-like) search (step 4)")
      ("fracStep4", po::value<double>(&fracStep4)->default_value(1),
       "fraction of samples to re-phase in 4th step")
      ("seed", po::value<uint>(&seed)->default_value(0), "random seed (ignored in ref-mode)")

      // error-checking
      ("noMapCheck", "disable automatic check of genetic map scale")

      // testing options
      ("iter", po::value<int>(&iter)->default_value(0), "iter to run")
      ("maskFile", po::value<string>(&maskFile), "indivs to mask (e.g., relatives)")
      ("tmpPhaseConfsPrefix", po::value<string>(&tmpPhaseConfsPrefix),
       "prefix for tmp files of phase confidences")
      ("maxHapStates", po::value<int>(&maxHapStates)->default_value(80),
       "maximum copying haplotype states per position in dynamic programming search")
      ("trioCheck", "flag to output trio check; assumes target samples are in child,mat,pat order")
      ;

    po::options_description visible("Options");
    visible.add(commonOptions).add(nonRefMode).add(refMode).add(bothModes).add(algOptions);

    po::options_description all("All options");
    all.add(commonOptions).add(nonRefMode).add(refMode).add(bothModes).add(algOptions).add(hidden);
    all.add_options()
      ("bad-args", po::value< vector <string> >(), "bad args")
      ;
    po::positional_options_description positional_desc;
    positional_desc.add("bad-args", -1); // for error-checking command line
    
    po::variables_map vm;
    po::command_line_parser cmd_line(argc, argv);
    cmd_line.options(all);
    cmd_line.style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing);
    cmd_line.positional(positional_desc);
    try {
      po::store(cmd_line.run(), vm);

      if (vm.count("help")) {
	cout << endl;
	cout << visible << endl;
	exit(0);
      }
      
      po::notify(vm); // throws an error if there are any problems

      usePBWT = !vm.count("v1");
      pbwtOnly = vm.count("pbwtOnly");
      if (pbwtOnly && !usePBWT) {
	cerr << "ERROR: --pbwtOnly cannot be specified if using the --v1 algorithm" << endl;
	return false;
      }
      trioCheck = vm.count("trioCheck");

      if (vm.count("bfile") +
	  vm.count("bfilegz") +
	  (vm.count("fam") || vm.count("bim") || vm.count("bed")) +
	  vm.count("vcf") +
	  (vm.count("vcfRef") || vm.count("vcfTarget")) != 1) {
	cerr << "ERROR: Use exactly one of the --bfile, --bfilegz, --fam,bim,bed, --vcf, or"
	     << endl << "       --vcfRef,vcfTarget input formats" << endl;
	return false;
      }

      if (vm.count("vcf") + vm.count("vcfRef") + vm.count("vcfTarget") == 0 &&
	  vcfOutFormat != ".") {
	cerr << "ERROR: --vcfOutFormat can only be used with vcf input" << endl;
	return false;
      }

      if (vm.count("bfile")) {
	string bfile = vm["bfile"].as<string>();
	famFile = bfile + ".fam";
	bimFile = bfile + ".bim";
	bedFile = bfile + ".bed";
      }

      if (vm.count("bfilegz")) {
	string bfile = vm["bfilegz"].as<string>();
	famFile = bfile + ".fam.gz";
	bimFile = bfile + ".bim.gz";
	bedFile = bfile + ".bed.gz";
      }

      if (vm.count("bad-args")) {
	cerr << "ERROR: Unknown options:";
	vector <string> bad_args = vm["bad-args"].as< vector <string> >();
	for (uint i = 0; i < bad_args.size(); i++) cerr << " " << bad_args[i];
	cerr << endl;
	return false;
      }

      noMapCheck = vm.count("noMapCheck");
      noImpMissing = vm.count("noImpMissing");
      allowRefAltSwap = vm.count("allowRefAltSwap");
      outputUnphased = vm.count("outputUnphased");
      keepMissingPloidyX = vm.count("keepMissingPloidyX");

      if (vm.count("vcfRef") || vm.count("vcfTarget") || vm.count("vcf")) { // VCF mode
	if (vm.count("vcf")) { // non-ref mode
	  if (bpFlanking != 0) {
	    cerr << "ERROR: --bpFlanking is only supported in ref-mode" << endl;
	    return false;
	  }
	}
	else { // ref-mode
	  if (vcfRef.empty()) {
	    cerr << "ERROR: --vcfRef must be specified in reference-based phasing mode" << endl;
	    return false;
	  }
	  if (vcfTarget.empty()) {
	    cerr << "ERROR: --vcfTarget must be specified in reference-based phasing mode" << endl;
	    return false;
	  }
	  if (vcfRef.substr(vcfRef.length()-4) != string(".bcf")) {
	    cerr << "WARNING: --vcfRef does not end in '.bcf'; BCF input is fastest" << endl;
	  }
	  if (chromStr=="0" && (bpStart>0 || bpEnd<1e9)) {
	    cerr << "ERROR: --chrom is required when specifying --bpStart or --bpEnd in ref mode"
		 << endl;
	    return false;
	  }
	}

	// vcf input checks for both ref and non-ref mode
	if (geneticMapFile == "USE_BIM") {
	  cerr << "ERROR: --geneticMapFile must be specified when using VCF/BCF input"
	       << endl;
	  return false;
	}
	if (!removeFileTemplates.empty() || !excludeFileTemplates.empty() ||
	    maxMissingPerSnp != 0.1 || maxMissingPerIndiv != 0.1) {
	  cerr << "ERROR: --remove, --exclude, --maxMissingPerSnp, --maxMissingPerIndiv"
	       << "       are not supported for VCF/BCF input or in reference mode" << endl;
	  return false;
	}
	if (vcfOutFormat == "b") { vcfOutSuffix = "bcf"; vcfWriteMode = "wb"; }
	else if (vcfOutFormat == "u") { vcfOutSuffix = "bcf"; vcfWriteMode = "wbu"; }
	else if (vcfOutFormat == "z" || vcfOutFormat == ".") {
	  vcfOutSuffix = "vcf.gz"; vcfWriteMode = "wz";
	}
	else if (vcfOutFormat == "v") { vcfOutSuffix = "vcf"; vcfWriteMode = "w"; }
	else {
	  cerr << "ERROR: --vcfOutFormat must be one of {b,u,z,v}" << endl;
	  return false;
	}
	if (bpFlanking < 0) {
	  cerr << "ERROR: --bpFlanking cannot be negative" << endl;
	  return false;
	}
      }
      else { // non-ref mode; plink input
	if (famFile.empty()) {
	  cerr << "ERROR: fam file must be specified either using --fam or --bfile"
	       << endl;
	  return false;
	}
	if (bimFile.empty()) {
	  cerr << "ERROR: bim file must be specified either using --bim or --bfile"
	       << endl;
	  return false;
	}
	if (bedFile.empty()) {
	  cerr << "ERROR: bed file must be specified either using --bed or --bfile"
	       << endl;
	  return false;
	}
	if (noImpMissing) {
	  cerr << "ERROR: --noImpMissing is only supported with vcf input" << endl;
	  return false;
	}
	if (bpFlanking != 0) {
	  cerr << "ERROR: --bpFlanking is only supported in ref-mode" << endl;
	  return false;
	}
      }

      removeFiles = StringUtils::expandRangeTemplates(removeFileTemplates);
      excludeFiles = StringUtils::expandRangeTemplates(excludeFileTemplates);

      if (!(0 <= maxMissingPerSnp && maxMissingPerSnp <= 1)) {
	cerr << "ERROR: --maxMissingPerSnp must be between 0 and 1" << endl;
	return false;	
      }
      if (!(0 <= maxMissingPerIndiv && maxMissingPerIndiv <= 1)) {
	cerr << "ERROR: --maxMissingPerIndiv must be between 0 and 1" << endl;
	return false;
      }
      chrom = StringUtils::bcfNameToChrom(chromStr.c_str(), 0, chromX); // checks for range

      if (pbwtIters < 0 || pbwtIters > 10) {
	cerr << "ERROR: --pbwtIters must be either 0=auto or <=10" << endl;
	return false;
      }

      // check advanced options
      if (vm.count("v1fast")) {
	cMmax = 0.5;
	beamWidth3 = 100;
	beamWidth4 = 100;
	fracStep4 = 0.5;
      }
      if ((cMmax != 0 && cMmax < 0.1) || cMmax > 1.0) {
	cerr << "ERROR: --maxBlockLen must be 0=auto or between 0.1 and 1 cM" << endl;
	return false;
      }
      if (beamWidth3 < 10 || beamWidth3 > 1000) {
	cerr << "ERROR: --maxStatePairsStep3 must be between 10 and 1000" << endl;
	return false;
      }
      if (beamWidth4 < 10 || beamWidth4 > 1000) {
	cerr << "ERROR: --maxStatePairsStep4 must be between 10 and 1000" << endl;
	return false;
      }
      if (maxHapStates < 40 || maxHapStates > 1000) {
	cerr << "ERROR: --maxHapStates must be between 40 and 1000" << endl;
	return false;
      }
      if (fracStep4 < 0.0 || fracStep4 > 1.0) {
	cerr << "ERROR: --fracStep4 must be between 0.0 and 1.0" << endl;
	return false;
      }
      if (pErr < 1e-6 || pErr > 0.1) {
	cerr << "ERROR: --genoErrProb must be between 0.000001 and 0.1" << endl;
	return false;
      }

      // check that all files specified are readable/writeable
      FileUtils::requireEmptyOrReadable(famFile);
      FileUtils::requireEmptyOrReadable(bimFile);
      FileUtils::requireEmptyOrReadable(bedFile);
      FileUtils::requireEmptyOrReadable(vcfRef);
      FileUtils::requireEmptyOrReadable(vcfTarget);
      FileUtils::requireEmptyOrReadable(vcfExclude);
      if (geneticMapFile != "USE_BIM") {
	vector <string> reqHeader;
	reqHeader.push_back("chr"); reqHeader.push_back("position");
	reqHeader.push_back("COMBINED_rate(cM/Mb)"); reqHeader.push_back("Genetic_Map(cM)");
	if (FileUtils::parseHeader(geneticMapFile, " \t") != reqHeader) {
	  cerr << "ERROR: --geneticMapFile must have four columns with names:" << endl
	       << "       chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)" << endl;
	  return false;
	}
      }
      FileUtils::requireEachEmptyOrReadable(removeFiles);
      FileUtils::requireEachEmptyOrReadable(excludeFiles);
    }
    catch (po::error &e) {
      cerr << "ERROR: " << e.what() << endl << endl;
      cerr << visible << endl;
      return false;
    }
    return true;
  }
}

