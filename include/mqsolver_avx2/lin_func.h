// Filename:
//       lin_func.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       miscellaneous functions related to linear systems over
//       GF(2) based on AVX2 intrinsic set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#ifndef MQSOLVER_LIN_FUNC_H
#define MQSOLVER_LIN_FUNC_H

#include "param.h"
#include "mqsolver_avx2/avx2_intrinsic_ext.h"

// check whether a linear system is solvable, the consistensy
// flag is stored in the mask.
// if the linear system is solvable, then after the process,
// !(clist[-1][i] & mask[i]) == 0
// this function is similar to Gaussian elimination but without
// pivoting the first set bit
void checkConsist(__m256i clist[LIN_COLNUM], __m256i& mask);
void checkConsist(uint32_t clist[LIN_COLNUM], uint32_t& mask);

// extract one solution from a linear system that being processed
// by checkConsist()
void extractSolution(const uint32_t clist[LIN_COLNUM], uint32_t sol[LIN_VARNUM]);

#endif //MQSOLVER_LIN_FUNC_H
