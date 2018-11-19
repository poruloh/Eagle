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

#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include <vector>
#include <string>

namespace StringUtils {

  const std::string RANGE_DELIMS = "{:}"; // must have 3 chars

  int stoi(const std::string &s);
  double stod(const std::string &s);
  std::string itos(int i);
  std::string findDelimiters(const std::string &s, const std::string &c);

  // will not return blanks
  std::vector <std::string> tokenizeMultipleDelimiters(const std::string &s, const std::string &c);
  
  // basic range template: expand "{start:end}" to vector <string> with one entry per range element
  // if end==start-1, will return empty
  std::vector <std::string> expandRangeTemplate(const std::string &str);
  std::vector <std::string> expandRangeTemplates(const std::vector <std::string> &rangeTemplates);

  int bcfNameToChrom(const char *nameBuf, int chromMin, int chromX);
}

#endif
