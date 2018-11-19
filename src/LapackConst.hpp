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

#ifndef LAPACKCONST_HPP
#define LAPACKCONST_HPP

#ifdef USE_MKL

#include "mkl.h"
#define DGER_MACRO dger
#define DGEMV_MACRO dgemv
#define DGEMM_MACRO dgemm
#define SGEMM_MACRO sgemm
#define DGELS_MACRO dgels
#define DGESVD_MACRO dgesvd

#else

  extern "C" int dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s,
			 double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork,
			 int *info);
  extern "C" int dgemv_(char *TRANS, int *M, int *N, double *ALPHA, double *A, int *LDA,
			const double *X, int *INCX, double *BETA, double *Y, int *INCY);
  extern "C" int dger_(int *M, int *N, double *ALPHA, double *X, int *INCX, const double *Y,
		       int *INCY, double *A, int *LDA);
  extern "C" int dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA,
			const double *A, int *LDA, const double *B, int *LDB, double *BETA,
			double *C, int *LDC);
  extern "C" int sgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA,
			const float *A, int *LDA, const float *B, int *LDB, float *BETA,
			float *C, int *LDC);
  extern "C" int dgels_(char *TRANS, int *M, int *N, int *NRHS, double *A, int *LDA, double *B,
			int *LDB, double *WORK, int *LWORK, int *INFO);

#define DGER_MACRO dger_
#define DGEMV_MACRO dgemv_
#define DGEMM_MACRO dgemm_
#define SGEMM_MACRO sgemm_
#define DGELS_MACRO dgels_
#define DGESVD_MACRO dgesvd_

#endif



/*
namespace LapackConst {

#ifndef USE_MKL
#ifdef USE_MKL
  inline CBLAS_TRANSPOSE lapackTransToMKL(char trans) {
    return (trans=='N'||trans=='n') ? CblasNoTrans : CblasTrans;
  }
#endif

  void dgemm_wrap(char TRANSA, char TRANSB, int M, int N, int K, double ALPHA,
		  const double *A, int LDA, const double *B, int LDB, double BETA,
		  double *C, int LDC);
}
*/
#endif
