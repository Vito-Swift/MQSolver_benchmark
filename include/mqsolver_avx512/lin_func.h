// Filename:
//       lin_func.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       miscellaneous functions related to linear systems over
//       GF(2) based on AVX512 intrinsic set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#ifndef MQSOLVER_LIN_FUNC_H
#define MQSOLVER_LIN_FUNC_H

#include "param.h"
#include "mqsolver_avx512/avx512_intrinsic_ext.h"

// check whether a linear system is solvable, the consistensy
// flag is stored in the mask.
// if the linear system is solvable, then after the process,
// !(clist[-1][i] & mask[i]) == 0
// this function is similar to Gaussian elimination but without
// pivoting the first set bit
inline void checkConsist(__m512i clist[LIN_COLNUM], __m512i& mask) {
    for (int i = 0; i < LIN_VARNUM; i++) {
        __m512i ci = _mm512_and_si512(clist[i], mask);
        __mmask32 m = ~_mm512_cmpeq_epi16_mask(ci, mm512_const_zero);
        if (m == 0) continue;
        __m512i x = avx512_lzcnt_epi16(ci);
        __m512i xorp = _mm512_maskz_srlv_epi16(m, mm512_const_leadone_epi16, x);
        __m512i xormask = _mm512_maskz_xor_epi32(m, clist[i], xorp);
        for (int j = i; j < LIN_COLNUM; j++) {
            __m512i cj = _mm512_and_si512(clist[j], xorp);
            __mmask32 _m = _mm512_cmpeq_epi16_mask(cj, xorp);
            __m512i mxorj = _mm512_maskz_loadu_epi16(_m, &xormask);
            clist[j] = _mm512_xor_si512(clist[j], mxorj);
        }
        mask = _mm512_xor_si512(mask, xorp);
    }
}

// extract one solution from a linear system that being processed
// by checkConsist()
void extractSolution(const uint16_t clist[LIN_COLNUM], uint32_t sol[LIN_VARNUM]);
inline void extractSolution(const __m512i clist[LIN_COLNUM], __m512i sol[LIN_VARNUM], 
                            uint32_t offset) {
    for (int i = 0; i < LIN_VARNUM; i++) {
        __m512i x = avx512_lzcnt_epi16(clist[i]);
        __mmask32 m = ~_mm512_cmpeq_epi16_mask(clist[i], mm512_const_zero);
        __m512i xp = _mm512_maskz_srlv_epi16(m, mm512_const_leadone_epi16, x);
        __m512i bs = _mm512_and_si512(xp, clist[LIN_COLNUM]);
        __mmask32 _m = ~_mm512_cmpeq_epi16_mask(bs, mm512_const_zero);
        __m512i s = _mm512_maskz_loadu_epi16(_m, &mm512_const_setone_epi16);
        sol[i] = _mm512_or_si512(_mm512_slli_epi16(s, offset), sol[i]);
    }
}

#endif //MQSOLVER_LIN_FUNC_H
