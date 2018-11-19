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

#include <cmath>
#include <cstdlib>
#include <utility>

#include "NumericUtils.hpp"

namespace NumericUtils {
  double sum(const double x[], uint64 N) {
    double ans = 0;
    for (uint64 n = 0; n < N; n++)
      ans += x[n];
    return ans;
  }
  double mean(const std::vector <double> &x) {
    uint64 N = x.size(); return sum(&x[0], N) / N;
  }
  // takes into account that some 0 values may indicate missing/ignored: divide out by Nused, not N
  double mean(const double x[], uint64 N, uint64 Nused) {
    return sum(x, N) / Nused;
  }
  // regress y on x, assuming both have been 0-centered (so 0-filled missing values ok)
  double regCoeff(const double y[], const double x[], uint64 N) {
    /* WRONG! if not mean-centered already, need to mask missing indivs in loop
    double xbar = mean(x, N, Nused);
    double ybar = mean(y, N, Nused);
    cout << "xbar: " << xbar << " ybar: " << ybar << endl;
    double numer = 0, denom = 0;
    for (uint64 n = 0; n < N; n++) {
      numer += (x[n]-xbar) * (y[n]-ybar);
      denom += sq(x[n]-xbar);
    }
    */
    double numer = 0, denom = 0;
    for (uint64 n = 0; n < N; n++) {
      numer += x[n] * y[n];
      denom += sq(x[n]);
    }
    return numer / denom;
  }
  double dot(const double x[], const double y[], uint64 N) {
    double ans = 0;
    for (uint64 n = 0; n < N; n++)
      ans += x[n] * y[n];
    return ans;
  }
  double norm2(const double x[], uint64 N) {
    double ans = 0;
    for (uint64 n = 0; n < N; n++)
      ans += sq(x[n]);
    return ans;
  }
  void normalize(double x[], uint64 N) {
    double scale = 1.0 / sqrt(norm2(x, N));
    for (uint64 n = 0; n < N; n++)
      x[n] *= scale;
  }

  std::pair <double, double> meanStdDev(const double x[], uint64 N) {
    double mu = 0, s2 = 0;
    for (uint64 n = 0; n < N; n++) mu += x[n];
    mu /= N;
    for (uint64 n = 0; n < N; n++) s2 += sq(x[n]-mu);
    s2 /= (N-1);
    return std::make_pair(mu, sqrt(s2));
  }
  std::pair <double, double> meanStdErr(const double x[], uint64 N) {
    std::pair <double, double> ret = meanStdDev(x, N);
    ret.second /= sqrt((double) N);
    return ret;
  }
  std::pair <double, double> meanStdDev(const std::vector <double> &x) {
    return meanStdDev(&x[0], x.size());
  }
  std::pair <double, double> meanStdErr(const std::vector <double> &x) {
    return meanStdErr(&x[0], x.size());
  }
  void logSumExp(float &x, float y) {
    float big, diff;
    if (x > y) {
      big = x; diff = y-x;
    }
    else {
      big = y;
      diff = x-y;
    }
    if (diff < -10) x = big; // a < 1e-4 * b => ignore a
    else if (diff < -5) x = big + expf(diff); // a < 1e-2 * b => use 1st-order Taylor expansion
    else x = big + logf(1.0f + expf(diff));
  }
}
