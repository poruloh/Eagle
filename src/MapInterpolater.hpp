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

#ifndef MAPINTERPOLATER_HPP
#define MAPINTERPOLATER_HPP

#include <string>
#include <map>
#include <utility>

namespace Genetics {

  class MapInterpolater {
    std::map < std::pair <int, int>, std::pair <double, double> > chrBpToRateGen;
    static const std::string MAP_FILE_HEADER;
  public:
    // input file format: chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    // (Oxford map format preceded by chr column)
    MapInterpolater(const std::string &geneticMapFile);
    // returns interpolated genetic position in Morgans
    double interp(int chr, int bp) const;
  };

}
#endif
