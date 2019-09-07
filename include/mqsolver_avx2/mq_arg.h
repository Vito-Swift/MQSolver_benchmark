// Filename:
//       mq_arg.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       miscellaneous functions related to multivariate quadratic
//       problems based on AVX2 intrinsic set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#ifndef MQSOLVER_MQ_ARG_H
#define MQSOLVER_MQ_ARG_H

#include <stdint.h>
#include <fstream>
#include <cstring>
#include <emmintrin.h>
#include <immintrin.h>

#include "param.h"
#include "utils.h"
#include "mqsolver_avx2/avx2_intrinsic_ext.h"
#include "mqsolver_avx2/lin_func.h"

#define VEC_SIZE 8

// data structure to store a mq polynomial.
// poly is a uint64_t array, each element stands for one term.
// in each term, no more that 2 bits will be set.
struct poly {
  uint64_t *p;
  uint32_t length;
};

static uint64_t grayCodeKey = 0x80000000 | ((uint64_t) 0x1ff << RES_VARNUM);
static uint64_t binKey = 0;
static poly ResidentPoly[LIN_EQNUM][N + 1];
static poly VerifyPoly[M][N + 1];
static uint16_t partialDerivative[N + 1][RES_VARNUM] = {0};
static uint64_t ColVal[N + 1] = {0};

// add one term to a given polynomial.
void addTerm(poly &dstPoly, uint64_t term);

// read the entire MQ system from challenge file.
void file2poly(FILE *fr, poly dstPoly[M][N + 1]);

// read the first m equation from resident file.
void ffile2poly(FILE *fr, poly dstPoly[M][N + 1], int m);

// calculate the partial derivative of resident polynomial and
// store in colPD.
void loadPD(poly residentPoly[M][N + 1], uint16_t colPD[N + 1][RES_VARNUM]);

// validate 256 candidates simultaneously
__m256i validateResult(__m256i guess[N]);

// initialization of the mqsolver
void mqInit(const std::string residentFileName, const std::string verifyFileName);

// main loop of the mqsolver
void mqLoop();

#endif //MQSOLVER_MQ_ARG_H
