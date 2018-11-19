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

#ifndef NUMERICUTILS_HPP
#define NUMERICUTILS_HPP

#include <cstdlib>
#include <vector>
#include <utility>

#include "Types.hpp"

namespace NumericUtils {

  inline double sq(double x) { return x*x; }
  double sum(const double x[], uint64 N);
  double mean(const std::vector <double> &x);

  // takes into account that some 0 values may indicate missing/ignored: divide out by Nused, not N
  double mean(const double x[], uint64 N, uint64 Nused);

  // regress y on x, assuming both have been 0-centered (so 0-filled missing values ok)
  double regCoeff(const double y[], const double x[], uint64 N);

  double dot(const double x[], const double y[], uint64 N);
  double norm2(const double x[], uint64 N);
  void normalize(double x[], uint64 N);

  void logSumExp(float &x, float y);

  std::pair <double, double> meanStdDev(const double x[], uint64 N);
  std::pair <double, double> meanStdErr(const double x[], uint64 N);
  std::pair <double, double> meanStdDev(const std::vector <double> &x);
  std::pair <double, double> meanStdErr(const std::vector <double> &x);
}

#endif
