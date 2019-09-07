// Filename:
//       mq_arg.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       miscellaneous functions related to multivariate quadratic
//       problems based on AVX512 intrinsic set.
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
#include "mqsolver_avx512/lin_func.h"

#define VEC_SIZE 32
#define TD_PARNUM 5
#define GC_VARNUM 32
#define GC_INIT_KEY 0x80000000

// data structure to store a mq polynomial.
// poly is a uint64_t array, each element stands for one term.
// in each term, no more that 2 bits will be set.
struct poly {
    uint64_t *p;
    uint32_t length;
};

struct verifyUnit {
    uint64_t p1;
    uint64_t p2;
};

struct verifySeq {
    verifyUnit* vu;
    uint64_t length;
};

static verifySeq verifySequences[M];

static poly ResidentPoly[LIN_EQNUM][N + 1];
static poly VerifyPoly[M][N + 1];
static uint16_t ColPartialDrt[N + 1][RES_VARNUM];
static __m512i ColVal[N + 1];
const uint64_t paraseg = KEY_END / 32;
alignas(64) static uint64_t binKey[VEC_SIZE] = {0};
alignas(64) static uint64_t grayCodeKey[VEC_SIZE] = {0};

// add one term to a given polynomial.
void addTerm(poly &dstPoly, uint64_t term);

// read the entire MQ system from challenge file.
void file2poly(FILE *fr, poly dstPoly[M][N + 1]);

// read the first m equation from resident file.
void ffile2poly(FILE *fr, poly dstPoly[M][N + 1], int m);

// calculate the partial derivative of resident polynomial and
// store in colPD.
void loadPD(poly residentPoly[M][N + 1], uint32_t colPD[N + 1][RES_VARNUM]);

// validate 256 candidates simultaneously
__m512i validateResult(__m512i guess[N], poly vpoly[M][N + 1]);

// initialization of the verify sequences
void initVerifySeq();

// initialization of the mqsolver
void mqInit(const std::string residentFileName, const std::string verifyFileName);

// main loop of the mqsolver
void mqLoop();

#endif //MQSOLVER_MQ_ARG_H
