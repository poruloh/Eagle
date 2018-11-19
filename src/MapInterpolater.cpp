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
#include <string>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>

#include "MapInterpolater.hpp"
#include "FileUtils.hpp"

namespace Genetics {

  using std::vector;
  using std::string;
  using std::pair;
  using std::make_pair;
  using std::map;
  using std::cout;
  using std::cerr;
  using std::endl;
  using FileUtils::getline;

  const string MapInterpolater::MAP_FILE_HEADER =
    "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)";

  MapInterpolater::MapInterpolater(const string &geneticMapFile) {
    if (geneticMapFile.empty()) return;

    chrBpToRateGen[make_pair(0, 0)] = make_pair(0.0, 0.0); // sentinel at beginning
    string line;
    FileUtils::AutoGzIfstream fin; fin.openOrExit(geneticMapFile);
    getline(fin, line);
    if (line != MAP_FILE_HEADER) {
      cerr << "ERROR: Wrong format of reference map " << geneticMapFile << endl;
      cerr << "       Expecting header: " << MAP_FILE_HEADER << endl;
      exit(1);
    }
    int chr0 = 0, bp0 = 0; double gen0 = 0;
    int chr, bp; double rate, gen;
    while (fin >> chr >> bp >> rate >> gen) {
      if (chr == chr0) {
	if (gen<gen0 || bp<bp0) {
	  cerr << "ERROR: Genetic map contains out-of-order row:" << endl;
	  cerr << "       " << chr << " " << bp << " " << rate << " " << gen << endl;
	  exit(1);
	}
	if (bp<=0) {
	  cerr << "ERROR: Genetic map positions must be positive:" << endl;
	  cerr << "       " << chr << " " << bp << " " << rate << " " << gen << endl;
	  exit(1);
	}
	if (bp==bp0) {
	  cerr << "ERROR: Genetic map contains duplicate position:" << endl;
	  cerr << "       " << chr << " " << bp << " " << rate << " " << gen << endl;
	  exit(1);
	}
	chrBpToRateGen[make_pair(chr, bp)] = make_pair((gen-gen0)/(1e-6*(bp-bp0)), gen);
      }
      chr0 = chr; bp0 = bp; gen0 = gen;
    }
    chrBpToRateGen[make_pair(chr, bp+1e9)] = make_pair(1.0, gen+1e3); // sentinel at end
  }
  
  // returns interpolated genetic position in Morgans
  double MapInterpolater::interp(int chr, int bp) const {
    if (chrBpToRateGen.empty()) return 0;
    map < pair <int, int>, pair <double, double> >::const_iterator ubIter =
      chrBpToRateGen.upper_bound(make_pair(chr, bp)); // map record > (chr, bp)
    if (ubIter == chrBpToRateGen.end()) {
      cerr << "ERROR: Chromosome " << chr << " is not in genetic map" << endl;
      exit(1);
    }
    int ubChr = ubIter->first.first;
    int ubBp = ubIter->first.second;
    double ubRate = ubIter->second.first;
    double ubGen = ubIter->second.second;

    if (chr == ubChr) return 0.01 * (ubGen + 1e-6 * ubRate * (bp-ubBp)); // interpolate interval
    else return 0.01 * (--ubIter)->second.second; // end of previous chromosome
  }
}
