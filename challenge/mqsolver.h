#ifndef MQSOLVER_MQSOLVER_H
#define MQSOLVER_MQSOLVER_H

#include <iostream>
#include <stdint.h>
#include <time.h>
#include <fstream>
#include <string>
#include <bitset>
#include <cstring>

#include <emmintrin.h>
#include <immintrin.h>

#define M 55
#define N 55
#define RES_VARNUM 46
#define TD_VECSIZE 32
#define TD_VECPOW 5
#define ITER_VARNUM 41
#define LIN_VARNUM 9
#define LIN_COLNUM 10
#define LIN_EQNUM 16

#define TEST_SPEED true

struct poly {
  uint64_t *p;
  uint32_t length;
};

struct verifyUnit {
  uint64_t p1;
  uint64_t p2;
};

struct verifySeq {
  verifyUnit *vu;
  uint64_t length;
};

verifySeq verifySequences[M];

const std::string ResidentFilename = "mq-resident4-55-2.txt";
const std::string ChallengeFilename = "mq4-55-2-f.txt";
poly ResidentPoly[M][N + 1];
poly VerifyPoly[M][N + 1];
uint16_t ColPartialDrt[N + 1][RES_VARNUM];
__m512i ColVal[N + 1];
uint64_t binKey = 0x1186279a8c3 - 0x400000; // solution = 0x1186279a8c2
uint64_t grayCodeKey = (binKey - 1) ^((binKey - 1) >> 1);

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

inline __m512i avx512_lzcnt_epi16(const __m512i &v) {
    __m512i hut = _mm512_and_si512(v, mm512_const_hut_mask);
    __m512i lut = _mm512_sllv_epi32(v, mm512_const_lut_offset);
    __m512i lzcnt_hut = _mm512_lzcnt_epi32(hut);
    __m512i lzcnt_lut = _mm512_lzcnt_epi32(lut);
    __m512i ret_lzcnt_hut = _mm512_sllv_epi32(lzcnt_hut, mm512_const_lut_offset);
    return _mm512_or_epi32(ret_lzcnt_hut, lzcnt_lut);
}

inline void ffile2poly(FILE *fr, poly dstPoly[M][N + 1], int m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < N + 1; j++) {
            fscanf(fr, "%u", &(dstPoly[i][j].length));
            dstPoly[i][j].p = (uint64_t *) malloc(dstPoly[i][j].length * sizeof(uint64_t));
            for (int k = 0; k < dstPoly[i][j].length; k++) {
                uint32_t temp[3];
                for (int l = 0; l < 3; l++)
                    fscanf(fr, "%u", &(temp[l]));
                dstPoly[i][j].p[k] = ((uint64_t) temp[1] << 32) | (uint64_t) temp[0];
            }
        }
    }
}

inline void loadPD(poly residentPoly[M][N + 1], uint16_t colPD[N + 1][RES_VARNUM]) {
    // store the partial derivative in row form
    uint64_t tmpRowPD[M][N];
    for (int eqIndex = 0; eqIndex < LIN_EQNUM; eqIndex++) {
        for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++) {
            for (int termIndex = 0; termIndex < residentPoly[eqIndex][varIndex].length; termIndex++) {
                uint64_t los = 63 - _lzcnt_u64(residentPoly[eqIndex][varIndex].p[termIndex]);
                if (los > varIndex) {
                    tmpRowPD[eqIndex][varIndex] ^= (uint64_t) 1 << los;
                    tmpRowPD[eqIndex][los] ^= (uint64_t) 1 << varIndex;
                } else
                    tmpRowPD[eqIndex][varIndex] ^= (uint64_t) 1 << N;
            }
        }
    }

    // transpose row-form into column-form
    for (int varIndex = 0; varIndex < N + 1; varIndex++)
        for (int bitIndex = 0; bitIndex < RES_VARNUM; bitIndex++)
            for (int eqIndex = 0; eqIndex < LIN_EQNUM; eqIndex++)
                colPD[varIndex][bitIndex] |= ((tmpRowPD[eqIndex][bitIndex] >> varIndex) & 1) << eqIndex;
}

inline void initVerifySeq() {
    for (int i = 0; i < M; i++) {
        uint64_t len = 0;
        for (int j = 0; j < N; j++)
            len += VerifyPoly[i][j].length;
        verifySequences[i].length = len;
        verifySequences[i].vu = (verifyUnit *) malloc(len * sizeof(verifyUnit));
        int counter = 0;
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < VerifyPoly[i][j].length; k++) {
                int los = 63 - _lzcnt_u64(VerifyPoly[i][j].p[k]);
                verifySequences[i].vu[counter].p1 = los;
                verifySequences[i].vu[counter].p2 = j;
                counter++;
            }
        }
    }
}

inline void checkConsist(__m512i clist[LIN_COLNUM], __m512i &mask) {
    for (int i = 0; i < LIN_VARNUM; i++) {
        __m512i ci = _mm512_and_si512(clist[i], mask);
        __mmask32 m = ~_mm512_cmpeq_epi16_mask(ci, mm512_const_zero);
        if (m == 0) continue;
        __m512i x = avx512_lzcnt_epi16(ci);
        __m512i xorp = _mm512_maskz_srlv_epi16(m, mm512_const_leadone_epi16, x);
        __m512i xormask = _mm512_xor_epi32(clist[i], xorp);
        xormask = _mm512_maskz_loadu_epi16(m, &xormask);
        for (int j = i; j < LIN_COLNUM; j++) {
            __m512i cj = _mm512_and_si512(clist[j], xorp);
            __mmask32 _m = ~_mm512_cmpeq_epi16_mask(cj, mm512_const_zero);
            __m512i mxorj = _mm512_maskz_loadu_epi16(_m, &xormask);
            clist[j] = _mm512_xor_si512(clist[j], mxorj);
        }
        mask = _mm512_xor_si512(mask, xorp);
    }
}

inline void extractSolution(const __m512i clist[LIN_COLNUM], __m512i sol[LIN_VARNUM],
                            uint32_t offset) {
    for (int i = 0; i < LIN_VARNUM; i++) {
        __m512i x = avx512_lzcnt_epi16(clist[i]);
        __mmask32 m = ~_mm512_cmpeq_epi16_mask(clist[i], mm512_const_zero);
        __m512i xp = _mm512_maskz_srlv_epi16(m, mm512_const_leadone_epi16, x);
        __m512i bs = _mm512_and_si512(xp, clist[LIN_VARNUM]);
        __mmask32 _m = ~_mm512_cmpeq_epi16_mask(bs, mm512_const_zero);
        __m512i s = _mm512_maskz_sllv_epi16(_m, mm512_const_setone_epi16, _mm512_set1_epi16(offset));
        sol[i] = _mm512_or_si512(s, sol[i]);
    }
}

inline __m512i validateResult(__m512i guess[N], poly vpoly[M][N + 1]) {
    __m512i ret = _mm512_setzero_epi32();
    for (int i = 0; i < M; i++) {
        __m512i v = _mm512_setzero_epi32();
        for (int j = 0; j < verifySequences[i].length; j++) {
            verifyUnit vunit = verifySequences[i].vu[j];
            __m512i vec = _mm512_and_si512(guess[vunit.p1], guess[vunit.p2]);
            v = _mm512_xor_si512(vec, v);

        }
        if (vpoly[i][N].length)
            v = _mm512_xor_si512(v, mm512_const_fullone);
        ret = _mm512_or_si512(v, ret);

        __mmask8 m = _mm512_cmpeq_epi64_mask(ret, mm512_const_fullone);
        if (m == 0xff)
            return ret;
    }
    return ret;
}

inline const std::string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}

inline void mqInit() {
    // load resident partial derivative from file
    FILE *ResidentFile = fopen(ResidentFilename.c_str(), "rb");
    ffile2poly(ResidentFile, ResidentPoly, LIN_EQNUM);
    loadPD(ResidentPoly, ColPartialDrt);
    fclose(ResidentFile);
    std::cout << "\t" << currentDateTime() << "\t[info]\tpartial derivative loaded from " << ResidentFilename
              << std::endl;

    FILE *ChallengeFile = fopen(ChallengeFilename.c_str(), "rb");
    ffile2poly(ChallengeFile, VerifyPoly, M);
    fclose(ChallengeFile);
    std::cout << "\t" << currentDateTime() << "\t[info]\tchallenge loaded from " << ChallengeFilename << std::endl;

    initVerifySeq();
    uint64_t initBinKey[TD_VECSIZE] = {0};
    uint64_t initGrayCodeKey[TD_VECSIZE] = {0};
    for (uint64_t i = 0; i < TD_VECSIZE; i++) {
        initBinKey[i] = ((uint64_t) i << ITER_VARNUM) | binKey;
        initGrayCodeKey[i] = ((uint64_t) i << ITER_VARNUM) | grayCodeKey;
    }

    alignas(64) uint16_t tempEquationVal[TD_VECSIZE][N + 1] = {0};
    alignas(64) uint16_t procEquationVal[N + 1][TD_VECSIZE] = {0};
    for (uint64_t tdIndex = 0; tdIndex < TD_VECSIZE; tdIndex++) {
        for (uint64_t eqIndex = 0; eqIndex < LIN_EQNUM; eqIndex++) {
            for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++) {
                if ((initGrayCodeKey[tdIndex] >> varIndex) & 1) {
                    for (int termIndex = 0; termIndex < ResidentPoly[eqIndex][varIndex].length; termIndex++) {
                        uint64_t los = 63 - _lzcnt_u64(ResidentPoly[eqIndex][varIndex].p[termIndex]);
                        if (los == varIndex)
                            tempEquationVal[tdIndex][N] ^= 1 << eqIndex;
                        else if (los < RES_VARNUM)
                            tempEquationVal[tdIndex][N] ^= ((initGrayCodeKey[tdIndex] >> los) & 1) << eqIndex;
                        else
                            tempEquationVal[tdIndex][los] ^= 1 << eqIndex;
                    }
                }
            }
            for (int varIndex = RES_VARNUM; varIndex < N; varIndex++)
                if (ResidentPoly[eqIndex][varIndex].length)
                    tempEquationVal[tdIndex][varIndex] ^= 1 << eqIndex;
            tempEquationVal[tdIndex][N] ^= ResidentPoly[eqIndex][N].length << eqIndex;
        }
    }
    for (int tdIndex = 0; tdIndex < TD_VECSIZE; tdIndex++)
        for (int varIndex = 0; varIndex < N + 1; varIndex++)
            procEquationVal[varIndex][tdIndex] = tempEquationVal[tdIndex][varIndex];
    for (int varIndex = 0; varIndex < N + 1; varIndex++)
        ColVal[varIndex] = _mm512_maskz_loadu_epi16(0xffffffff, procEquationVal[varIndex]);
}

void mqLoop() {
    uint64_t loopNum = 0;
    __m512i sols[N];
    for (uint i = 0; i < N; i++)
        sols[i] = mm512_const_zero;

    for (uint32_t varIndex = 0; varIndex < TD_VECPOW; varIndex++) {
        uint16_t gcKey_hut_bit[TD_VECSIZE] = {0};
        for (uint32_t i = 0; i < TD_VECSIZE; i++)
            if ((i >> varIndex) & 1) gcKey_hut_bit[i] = 0xffff;
        sols[ITER_VARNUM + varIndex] = _mm512_maskz_loadu_epi16(0xffffffff, gcKey_hut_bit);
    }
    uint32_t solOffset = 0;

    std::cout << "\t" << currentDateTime() << "\t[info]\tbruteforcing..." << std::endl;
    while (true) {
        __m512i Val[LIN_VARNUM + 1];
        std::copy(ColVal + RES_VARNUM, ColVal + N + 1, Val);

        // set matrix value
        for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++)
            Val[LIN_VARNUM] = _mm512_xor_si512(Val[LIN_VARNUM], ColVal[varIndex]);

        // check consistency
        __m512i mask = _mm512_set1_epi16((short) 0xffff);
        checkConsist(Val, mask);
        
        if (!TEST_SPEED) {
            // extract solution
            // TODO: not need to extract all the solutions
            for (uint32_t varIndex = 0; varIndex < ITER_VARNUM; varIndex++)
                if ((grayCodeKey >> varIndex) & 1)
                    sols[varIndex] = _mm512_or_si512(sols[varIndex], _mm512_slli_epi16(_mm512_set1_epi16(1), solOffset));
            extractSolution(Val, sols + RES_VARNUM, solOffset);
            solOffset++;

            // check result
            if (solOffset == 16) {
                __m512i flag = validateResult(sols, VerifyPoly);
                __mmask16 flag_mask = _mm512_cmpeq_epi32_mask(flag, mm512_const_fullone);
                if (flag_mask != 0xffff)
                    std::cout << "\t" << currentDateTime() << "\t[info]\tone solution is found" << std::endl;
                for (uint i = 0; i < ITER_VARNUM; i++)
                    sols[i] = mm512_const_zero;
                for (uint i = RES_VARNUM; i < N; i++)
                    sols[i] = mm512_const_zero;
                solOffset = 0;
            }
        }

        // update value
        uint64_t los = __builtin_ffsl(((uint64_t) binKey >> 1) ^ binKey ^ grayCodeKey) - 1;
        for (int varIndex = 0; varIndex < ITER_VARNUM; varIndex++)
            if ((grayCodeKey >> varIndex) & 1)
                ColVal[varIndex] = _mm512_xor_si512(_mm512_set1_epi16(ColPartialDrt[varIndex][los]), ColVal[varIndex]);
        for (int varIndex = 0; varIndex < TD_VECPOW; varIndex++) {
            __m512i tdVal = _mm512_srli_epi16(mm512_asc_vec, varIndex);
            __mmask32 m = ~_mm512_cmpeq_epi16_mask(_mm512_and_si512(tdVal, mm512_const_setone_epi16), mm512_const_zero);
            __m512i xormask = _mm512_maskz_set1_epi16(m, ColPartialDrt[varIndex + ITER_VARNUM][los]);
            ColVal[varIndex + ITER_VARNUM] = _mm512_xor_si512(ColVal[varIndex + ITER_VARNUM], xormask);
        }
        for (int varIndex = RES_VARNUM; varIndex < N + 1; varIndex++)
            ColVal[varIndex] = _mm512_xor_si512(_mm512_set1_epi16(ColPartialDrt[varIndex][los]), ColVal[varIndex]);

        // update gray code and bin key
        grayCodeKey = (binKey >> 1) ^ binKey;
        binKey++;

        if (++loopNum % 0x4000000 == 0)
            std::cout << "\t" << currentDateTime() << "\t[info]\tprocess: " << (uint64_t) (loopNum / 0x10000000)
                      << "/16834" << std::endl;
    }
}

#endif //MQSOLVER_MQSOLVER_H
