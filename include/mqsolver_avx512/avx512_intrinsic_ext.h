// Filename:
//       avx512_intrinsic_ext.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       extensional utility functions for avx512 intrinsic set.
//       these functions are all inlined for performance issue.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#ifndef MQSOLVER_AVX512_INTRINSIC_EXT_H
#define MQSOLVER_AVX512_INTRINSIC_EXT_H

#include <immintrin.h>
#include <emmintrin.h>

const __m512i mm512_const_zero = _mm512_setzero_si512();
const __m512i mm512_const_fullone = _mm512_set1_epi32(0xffffffff);
const __m512i mm512_const_leadone_epi16 = _mm512_set1_epi16(0x8000);
const __m512i mm512_const_setone_epi64 = _mm512_set1_epi64(1);
const __m512i mm512_const_hut_mask = _mm512_set1_epi32(0xffff0000);
const __m512i mm512_const_lut_mask = _mm512_set1_epi32(0x0000ffff);
const __m512i mm512_const_lut_offset = _mm512_set1_epi32(16);
const __m512i mm512_asc_vec = _mm512_set_epi16(31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14,
                                               13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

const __m512i mm512_asc_lut = _mm512_set_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
const __m512i mm512_asc_hut = _mm512_set_epi32(16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31);
const __m512i mm512_const_setone_epi16 = _mm512_set1_epi16(1);

// Counts the number of leading zero bits in each packed 16-bit
// integer in a, and store the results in dst.
//      for i in [0: 32]
//          a[i*16: (i+1)*16] = _lzcnt_u32(v[i*16: (i+1)*16]);
//      return a;
inline __m512i avx512_lzcnt_epi16(const __m512i &v) {
    __m512i lzcnt_hut = _mm512_lzcnt_epi32(_mm512_and_si512(v, mm512_const_hut_mask));
    __m512i lzcnt_lut = _mm512_lzcnt_epi32(_mm512_sllv_epi32(v, mm512_const_lut_offset));
    __m512i ret_lzcnt_hut = _mm512_sllv_epi32(lzcnt_hut, mm512_const_lut_offset);
    return _mm512_or_epi32(ret_lzcnt_hut, lzcnt_lut);
}

#endif //MQSOLVER_AVX512_INTRINSIC_EXT_H
