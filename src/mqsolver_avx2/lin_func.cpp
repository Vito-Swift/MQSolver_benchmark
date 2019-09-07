// Filename:
//       lin_func.cpp
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       implementation of miscellaneous functions related to linear
//       systems over GF(2) based on AVX2 intrinsic set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#include "mqsolver_avx2/lin_func.h"

void checkConsist(__m256i clist[LIN_COLNUM], __m256i &mask) {
    for (int i = 0; i < LIN_VARNUM; i++) {
        __m256i ci = _mm256_and_si256(clist[i], mask);
        __m256i m = ~_mm256_cmpeq_epi32(ci, mm256_const_zero);
        __m256i x = avx2_lzcnt_epi32(ci);
        __m256i xorp = _mm256_srlv_epi32(_mm256_and_si256(m, mm256_const_leadone_epi32), x);
        __m256i cia = _mm256_and_si256(clist[i], mm256_const_fullone);
        __m256i xormask = _mm256_and_si256(m, _mm256_xor_si256(cia, xorp));
        for (int j = 0; j < LIN_COLNUM; j++) {
            __m256i cj = _mm256_and_si256(clist[j], xorp);
            __mmask8 _m = ~avx2_cmpeqzero_epi32(cj);
            clist[j] = avx2_blendv_epi32(clist[j], _mm256_xor_si256(clist[j], xormask), _m);
        }
        mask = _mm256_xor_si256(mask, xorp);
    }
}

void checkConsist(uint32_t clist[LIN_COLNUM], uint32_t &mask) {
    const uint32_t _mask = 0xffffffff;
    for (int i = 0; i < LIN_VARNUM; i++) {
        uint32_t ci = clist[i] & mask;
        if (ci == 0) continue;
        uint32_t x = __builtin_clz(ci);
        uint32_t xp = (0x80000000 >> x);
        uint32_t xormask = (clist[i] & _mask) ^xp;
        for (int j = 0; j < LIN_COLNUM; j++)
            clist[j] = (clist[j] & xp) ? (clist[j] ^ xormask) : clist[j];
        mask ^= xp;
    }
}
void extractSolution(const uint32_t clist[LIN_COLNUM], uint32_t sol[LIN_VARNUM]) {
    for (int i = 0; i < LIN_VARNUM; i++) {
        if (clist[i] == 0) continue;
        uint32_t xp = 0x80000000 >> __builtin_clz(clist[i]);
        sol[i] = (bool) (clist[LIN_COLNUM] & xp);
    }
}