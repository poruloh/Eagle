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
#include <cstdio>
#include <cstdlib>

#include "StringUtils.hpp"
#include "FileUtils.hpp"
#include "Types.hpp"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace FileUtils {

  using std::string;
  using std::vector;
  using std::cerr;
  using std::endl;

  void openOrExit(std::ifstream &stream, const string &file,
		  std::ios_base::openmode mode) {
    stream.open(file.c_str(), mode);
    if (!stream) {
      cerr << "ERROR: Unable to open file: " << file << endl;
      exit(1);
    }
  }
  void openWritingOrExit(std::ofstream &stream, const string &file,
			 std::ios_base::openmode mode) {
    stream.open(file.c_str(), mode);
    if (!stream) {
      cerr << "ERROR: Unable to open file for writing: " << file << endl;
      exit(1);
    }
  }
  void requireEmptyOrReadable(const std::string &file) {
    if (file.empty()) return;
    std::ifstream fin;
    fin.open(file.c_str());
    if (!fin) {
      cerr << "ERROR: Unable to open file: " << file << endl;
      exit(1);
    }
    fin.close();
  }  
  void requireEachEmptyOrReadable(const std::vector <std::string> &fileList) {
    for (uint i = 0; i < fileList.size(); i++)
      requireEmptyOrReadable(fileList[i]);
  }
  void requireEmptyOrWriteable(const std::string &file) {
    if (file.empty()) return;
    std::ofstream fout;
    fout.open(file.c_str(), std::ios::out|std::ios::app);
    if (!fout) {
      cerr << "ERROR: Output file is not writeable: " << file << endl;
      exit(1);
    }
    fout.close();
  }
  vector <string> parseHeader(const string &fileName, const string &delimiters) {
    AutoGzIfstream fin; fin.openOrExit(fileName);
    string header;
    getline(fin, header);
    vector <string> split = StringUtils::tokenizeMultipleDelimiters(header, delimiters);
    fin.close();
    return split;
  }
  int lookupColumnInd(const string &fileName, const string &delimiters, const string &columnName) {
    vector <string> headers = parseHeader(fileName, delimiters);
    int columnInd = -1;
    for (uint c = 0; c < headers.size(); c++)
      if (headers[c] == columnName)
	columnInd = c; // first column is snp ID, treated separately
    if (columnInd == -1) {
      cerr << "WARNING: Column " << columnName << " not found in headers of " << fileName << endl;
      //exit(1);
    }
    return columnInd;
  }
  double readDoubleNanInf(std::istream &stream) {
    string str;
    stream >> str;
    double x;
    sscanf(str.c_str(), "%lf", &x);
    return x;
  }

  vector < std::pair <string, string> > readFidIids(const string &file) {
    vector < std::pair <string, string> > ret;
    AutoGzIfstream fin;
    fin.openOrExit(file);
    string FID, IID, line;
    while (fin >> FID >> IID) {
      if (FID.empty() || IID.empty()) {
	cerr << "ERROR: In file " << file << endl;
	cerr << "       unable to read FID and IID; check format" << endl;
	exit(1);
      }
      ret.push_back(make_pair(FID, IID));
      getline(fin, line);
    }
    fin.close();
    return ret;
  }


  int AutoGzIfstream::lineCount(const std::string &file) {
    AutoGzIfstream fin; fin.openOrExit(file);
    int ctr = 0; string line;
    while (getline(fin, line))
      ctr++;
    return ctr;
  }

  /***** AutoGzIfstream class implementation *****/

  void AutoGzIfstream::openOrExit(const std::string &file, std::ios_base::openmode mode) {
    fin.open(file.c_str(), mode);
    if (!fin) {
      cerr << "ERROR: Unable to open file: " << file << endl;
      exit(1);
    }
    if ((int) file.length() > 3 && file.substr(file.length()-3) == ".gz")
      boost_in.push(boost::iostreams::gzip_decompressor());
    boost_in.push(fin);
  }

  void AutoGzIfstream::close() {
    fin.close();
    boost_in.reset();
  }

  AutoGzIfstream::operator bool() const {
    return !boost_in.fail();
  }

  AutoGzIfstream& AutoGzIfstream::read(char *s, std::streamsize n) {
    boost_in.read(s, n);
    return *this;
  }

  int AutoGzIfstream::get() {
    return boost_in.get();
  }

  double AutoGzIfstream::readDoubleNanInf() {
    return FileUtils::readDoubleNanInf(boost_in);
  }

  void AutoGzIfstream::clear() {
    boost_in.clear();
  }

  AutoGzIfstream& AutoGzIfstream::seekg(std::streamoff off, std::ios_base::seekdir way) {
    boost_in.seekg(off, way);
    return *this;
  }

  AutoGzIfstream& getline(AutoGzIfstream& in, std::string &s) {
    std::getline(in.boost_in, s);
    return in;
  }

  
  /***** AutoGzOfstream class implementation *****/

  void AutoGzOfstream::openOrExit(const std::string &file, std::ios_base::openmode mode) {
    fout.open(file.c_str(), mode);
    if (!fout) {
      cerr << "ERROR: Unable to open file: " << file << endl;
      exit(1);
    }
    if ((int) file.length() > 3 && file.substr(file.length()-3) == ".gz")
      boost_out.push(boost::iostreams::gzip_compressor());
    boost_out.push(fout);
  }

  void AutoGzOfstream::close() {
    boost_out.reset();
  }

  AutoGzOfstream& AutoGzOfstream::operator << (std::ostream&(*manip)(std::ostream&)) {
    manip(boost_out);
    return *this;
  }

  void AutoGzOfstream::unsetf(std::ios_base::fmtflags mask) {
    boost_out.unsetf(mask);
  }

  AutoGzOfstream::operator bool() const {
    return !boost_out.fail();
  }

}
