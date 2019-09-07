// Filename:
//       avx2_intrinsic_ext.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       extensional utility functions for avx2 intrinsic set.
//       these functions are all inlined for performance issue.
//       most of functions in this file based on prototypes in
//       the AVX512 intrinsics set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#ifndef MQSOLVER_AVX2_INTRINSIC_EXT_H

#include <stdint.h>
#include <immintrin.h>
#include <emmintrin.h>

const __m256i mm256_const_zero = _mm256_setzero_si256();
const __m256i mm256_const_lut_mask = _mm256_set1_epi32(0xffff);
const __m256i mm256_const_hut_mask = _mm256_set1_epi32(0xffff0000);
const __m256i mm256_const_lut_offset = _mm256_set1_epi32(16);
const __m256i mm256_const_leadone_epi16 = _mm256_set1_epi16((short) 0x8000);
const __m256i mm256_const_leadone_epi32 = _mm256_set1_epi32(0x80000000);
const __m256i mm256_const_fullone = _mm256_set1_epi16((short) 0xffff);
const __m256i mm256_const_setone_epi64 = _mm256_set1_epi64x(1);

// Counts the number of leading zero bits in each packed 16-bit
// integer in a, and store the results in dst.
//      for i in [0: 16]
//          a[i*16: (i+1)*16] = _lzcnt_u32(v[i*16: (i+1)*16]);
//      return a;
inline __m256i avx2_lzcnt_epi16(__m256i &v) {
    const __m256i lut_lo = _mm256_set_epi8(
        4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 7, 16,
        4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 7, 16
    );
    const __m256i lut_hi = _mm256_set_epi8(
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 16,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 16
    );
    const __m256i nibble_mask = _mm256_set1_epi8(0x0F);
    const __m256i byte_offset = _mm256_set1_epi16(0x0008);
    __m256i t;

    t = _mm256_and_si256(nibble_mask, v);
    v = _mm256_and_si256(_mm256_srli_epi16(v, 4), nibble_mask);
    t = _mm256_shuffle_epi8(lut_lo, t);
    v = _mm256_shuffle_epi8(lut_hi, v);
    v = _mm256_min_epu8(v, t);

    t = _mm256_srli_epi16(v, 8);
    v = _mm256_or_si256(v, byte_offset);
    v = _mm256_min_epu8(v, t);

    return v;
}

// Counts the number of leading zero bits in each packed 32-bit
// integer in a, and store the results in dst.
inline __m256i avx2_lzcnt_epi32(__m256i& v) {
    const __m256i lzcnt_epi32 = avx2_lzcnt_epi16(v);
    return _mm256_sub_epi32(lzcnt_epi32, mm256_const_lut_offset);
}

// avx512 - _mm256_cmpeq_epi32 alternative
inline bool avx2_checkcontain_zero_epi32(const __m256i &v) {
    int arr[8];
    _mm256_maskstore_epi32(arr, mm256_const_fullone, v);
    for (auto &f: arr)
        if (f == 0)
            return true;
    return false;
}

// Compare packed 32-bit integers in a and b for
// equality, and store the results in mask vector k.
inline __mmask8 avx2_cmpeqzero_epi32(const __m256i &v) {
    __mmask8 ret = 0;

    int arr[8];
    _mm256_maskstore_epi32(arr, mm256_const_fullone, v);
    for (uint32_t i = 0; i < 8; i++)
        if (arr[i])
            ret |= 1 << i;
    return ret;
}

// Compare packed 16-bit integers in a and b for
// equality, and store the results in mask vector k.
inline __mmask16 avx2_cmpeqzero_epi16(const __m256i &v) {
    __mmask16 ret = 0;

    int arr[8];
    _mm256_maskstore_epi32(arr, mm256_const_fullone, v);
    for (uint32_t i = 0; i < 8; i++) {
        if (arr[i] & 0xffff)
            ret |= 1 << (i * 2);
        if (arr[i] & 0xffff0000)
            ret |= 1 << (i * 2 + 1);
    }
    return ret;
}

// Load 256-bits (composed of 8 packed 32-bit integers)
// from memory into dst.
// mem_addr has to be aligned to 32 bit
inline __m256i avx2_loadu_epi16(const uint16_t *mem_addr) {
    __m256i ret;
    ret = _mm256_maskload_epi32((int *) mem_addr, mm256_const_fullone);
    return ret;
}

// Load packed 16-bit integers from memory into dst using
// zeromask k (elements are zeroed out when the corresponding
// mask bit is not set). mem_addr must be aligned on a 32-byte
// boundary or a general-protection exception may be generated.
inline __m256i avx2_maskz_loadu_epi16(const __mmask16 mask, const __m256i &v) {

}

// Load packed 32-bit integers from memory into dst using
// zeromask k (elements are zeroed out when the corresponding
// mask bit is not set).
inline __m256i avx2_maskz_loadu_epi32(const __mmask8 mask, const __m256i &v) {

}

// Shift packed 16-bit integers in a right by the amount
// specified by the corresponding element in count while shifting
// in zeros, and store the results in dst using zeromask k
// (elements are zeroed out when the corresponding mask bit is
// not set).
inline __m256i avx2_maskz_srlv_epi16(const __mmask16 mask, const __m256i &v, const __m256i &count) {
    __m256i lut = _mm256_and_si256(v, mm256_const_lut_mask);
    __m256i ret_lut = _mm256_srlv_epi32(lut, _mm256_and_si256(count, mm256_const_lut_mask));
    __m256i hut = _mm256_srlv_epi32(v, _mm256_srlv_epi32(count, mm256_const_lut_offset));
    __m256i ret_hut = _mm256_and_si256(hut, mm256_const_hut_mask);

    return avx2_maskz_loadu_epi16(mask, _mm256_or_si256(ret_lut, ret_hut));
}

// Convert 8-bit control mask to vectorized control mask.
inline __m256i avx2_imm8_convert_mm256_epi32(const __mmask8& mask) {
    alignas(32) int arr[8] = {0};
    for (uint32_t i = 0; i < 8; i++)
        if ((mask >> i) & 1)
            arr[i] = (int) 0xffffffff;
    return _mm256_maskload_epi32(arr, mm256_const_fullone);
}

// Shift packed 32-bit integers in a right by the amount specified
// by the corresponding element in count while shifting in zeros,
// and store the results in dst using writemask k (elements are
// copied from src when the corresponding mask bit is not set).
inline __m256i avx2_maskz_srlv_epi32(const __mmask8 mask, const __m256i &v, const __m256i count) {
    __m256i mask_epi32 = avx2_imm8_convert_mm256_epi32(mask);
    return _mm256_srlv_epi32(_mm256_and_si256(v, mask_epi32), count);
}

// Compute the bitwise XOR of packed 32-bit integers in a and b,
// and store the results in dst using zeromask k (elements are
// zeroed out when the corresponding mask bit is not set).
inline __m256i avx2_maskz_xor_epi32(const __mmask8 mask, const __m256i &a, const __m256i &b) {
    __m256i mask_epi32 = avx2_imm8_convert_mm256_epi32(mask);
    return _mm256_and_si256(_mm256_xor_si256(a, b), mask_epi32);
}

// Blend packed 32-bit integers from a and b using mask, and store
// the results in dst.
inline __m256i avx2_blendv_epi32(const __m256i& a, const __m256i& b, const __mmask8 mask) {
    __m256i mask_epi32 = avx2_imm8_convert_mm256_epi32(mask);
    return _mm256_blendv_epi8(a, b, mask_epi32);
}

// Copy the lower 32 bit in each packed 64 bit of src and store
// in 128 bit dst.
inline void avx2_cvtsi256_si128_epi64(__m128i* dst, const __m256i& src) {
    const __m256i K_PERM = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
    __m256i permuted = _mm256_permutevar8x32_epi32(src, K_PERM);
    __m128i lo128 = _mm256_extractf128_si256(permuted, 0);
    _mm_storeu_si128((__m128i*)dst, lo128);
}

// Extract the lower 32 bit in each packed 64 bit of ma and mb
// and combine two segmentation to a 256 bit vector dst.
inline __m256i avx2_cvt_si256_2xsi128_mask(const __m256i& ma, const __m256i& mb) {
    __m128i ma_128, mb_128;
    avx2_cvtsi256_si128_epi64(&ma_128, ma);
    avx2_cvtsi256_si128_epi64(&mb_128, mb);
    __m256i ret = _mm256_castsi128_si256(ma_128);
    return _mm256_insertf128_si256(ret, mb_128, 1);
}

#define MQSOLVER_AVX2_INTRINSIC_EXT_H

#endif //MQSOLVER_AVX2_INTRINSIC_EXT_H
