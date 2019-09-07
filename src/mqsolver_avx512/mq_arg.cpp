// Filename:
//       mq_arg.cpp
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       implementation of miscellaneous functions related to
//       multivariate quadratic problems based on AVX512 intrinsic
//       set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#include <bitset>
#include "mqsolver_avx512/mq_arg.h"

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")

void addTerm(poly &dstPoly, uint64_t term) {
    uint32_t len = dstPoly.length + 1;
    if ((dstPoly.p = (uint64_t *) realloc(dstPoly.p, len * sizeof(uint64_t))) == nullptr) {
        std::cerr << "\t" << currentDateTime() << "\t[error]\tOut of memory!" << std::endl;
        exit(1);
    }
    dstPoly.p[dstPoly.length] = term;
    dstPoly.length = len;
}

void file2poly(FILE *fr, poly dstPoly[M][N + 1]) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N + 1; j++) {
            dstPoly[i][j].length = 0;
            dstPoly[i][j].p = nullptr;
        }
        uint32_t temp;
        for (int j = 0; j < N; j++) {
            for (int k = 0; k <= j; k++) {
                fscanf(fr, "%u", &temp);
                if (temp) {
                    uint64_t term = (uint64_t) 1 << j;
                    addTerm(dstPoly[i][j], term);
                }
            }
        }
        for (int j = 0; j < N; j++) {
            fscanf(fr, "%u", &temp);
            if (temp) {
                uint64_t term = (uint64_t) 1 << j;
                addTerm(dstPoly[i][j], term);
            }
        }
        fscanf(fr, "%u", &temp);
        if (temp) {
            uint64_t term = (uint64_t) 1 << N;
            addTerm(dstPoly[i][N], term);
        }
    }
}

void ffile2poly(FILE *fr, poly dstPoly[M][N + 1], int m) {
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

void loadPD(poly residentPoly[M][N + 1], uint16_t colPD[N + 1][RES_VARNUM]) {
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

__m512i validateResult(__m512i guess[N], poly vpoly[M][N + 1]) {
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

void initVerifySeq() {
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

void mqInit(const std::string residentFileName, const std::string verifyFileName) {
    // load resident partial derivative from file
    FILE *ResidentFile = fopen(residentFileName.c_str(), "rb");
    ffile2poly(ResidentFile, ResidentPoly, LIN_EQNUM);
    loadPD(ResidentPoly, ColPartialDrt);
    fclose(ResidentFile);
    std::cout << "\t" << currentDateTime() << "\t[info]\tpartial derivative loaded from " << residentFileName
              << std::endl;
    // load verify polynomial from challenge file
    FILE *VerifyFile = fopen(verifyFileName.c_str(), "rb");
    ffile2poly(VerifyFile, VerifyPoly, M);
    fclose(VerifyFile);
    std::cout << "\t" << currentDateTime() << "\t[info]\tverify polynomial loaded from " << verifyFileName << std::endl;
    // enumeration initialization
    for (int tdIndex = 0; tdIndex < VEC_SIZE; tdIndex++) {
        binKey[tdIndex] = (uint64_t) tdIndex << GC_VARNUM;
        grayCodeKey[tdIndex] = (uint64_t) tdIndex << GC_VARNUM | GC_INIT_KEY;
    }

    // equation value initialization
    alignas(64) uint16_t tempEquationVal[VEC_SIZE][N + 1] = {0};
    alignas(64) uint16_t procEquationVal[N + 1][VEC_SIZE] = {0};
    for (int tdIndex = 0; tdIndex < VEC_SIZE; tdIndex++) {
        for (int eqIndex = 0; eqIndex < LIN_EQNUM; eqIndex++) {
            for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++) {
                if ((grayCodeKey[tdIndex] >> varIndex) & 1) {
                    for (int termIndex = 0; termIndex < ResidentPoly[eqIndex][varIndex].length; termIndex++) {
                        uint64_t los = 63 - _lzcnt_u64(ResidentPoly[eqIndex][varIndex].p[termIndex]);
                        if (los == varIndex)
                            tempEquationVal[tdIndex][N] ^= 1 << eqIndex;
                        else if (los < RES_VARNUM)
                            tempEquationVal[tdIndex][N] ^= ((grayCodeKey[tdIndex] >> los) & 1) << eqIndex;
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
    for (int tdIndex = 0; tdIndex < VEC_SIZE; tdIndex++)
        for (int varIndex = 0; varIndex < N + 1; varIndex++)
            procEquationVal[varIndex][tdIndex] = tempEquationVal[tdIndex][varIndex];
    for (int varIndex = 0; varIndex < N + 1; varIndex++)
        ColVal[varIndex] = _mm512_maskz_loadu_epi16(0xffffffff, procEquationVal[varIndex]);
    initVerifySeq();
}

void mqLoop() {
    uint64_t loopNum = 0;
    uint64_t gcKeyVal = grayCodeKey[0];
    uint64_t binKeyVal = 0;
    __m512i sols[N];
    for (uint32_t varIndex = 0; varIndex < TD_PARNUM; varIndex++) {
        uint16_t gcKey_hut_bit[VEC_SIZE] = {0};
        for (uint32_t i = 0; i < VEC_SIZE; i++) gcKey_hut_bit[i] = 0xffff;
        _mm512_mask_loadu_epi16(sols[GC_VARNUM + varIndex], 0xffffffff, gcKey_hut_bit);
    }
    uint32_t solOffset = 0;

    std::cout << "\t" << currentDateTime() << "\t[info]\tstart mq main loop" << std::endl;
    while (binKeyVal < KEY_START + paraseg) {
        __m512i Val[LIN_VARNUM + 1];
        std::copy(ColVal + RES_VARNUM, ColVal + N + 1, Val);

        // set matrix value
        for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++)
            Val[LIN_VARNUM] = _mm512_xor_si512(Val[LIN_VARNUM], ColVal[varIndex]);

        // check consistency
        __m512i mask = _mm512_set1_epi16((short) LIN_MASK);
        checkConsist(Val, mask);
        
        // extract solution
        for (uint32_t varIndex = 0; varIndex < GC_VARNUM; varIndex++)
            if ((gcKeyVal >> varIndex) & 1)
                sols[varIndex] = _mm512_slli_epi16(_mm512_set1_epi16(1), solOffset);
        extractSolution(Val, sols + RES_VARNUM, solOffset);
        solOffset++;
        
        // check result
        if (solOffset == 16) {
            __m512i flag = validateResult(sols, VerifyPoly);
            __mmask16 flag_mask = _mm512_cmpeq_epi32_mask(flag, mm512_const_fullone);
            if (flag_mask != 0xffff)
                std::cout << "\t" << currentDateTime() << "\t[info]\tone solution is found" << std::endl;
            solOffset = 0;
        }

        // update value
        uint64_t los = __builtin_ffsl(((uint64_t) binKeyVal >> 1) ^ binKeyVal ^ gcKeyVal) - 1;
        for (int varIndex = 0; varIndex < GC_VARNUM; varIndex++)
            if ((gcKeyVal >> varIndex) & 1)
                ColVal[varIndex] = _mm512_xor_si512(_mm512_set1_epi16(ColPartialDrt[varIndex][los]), ColVal[varIndex]);
        for (int varIndex = 0; varIndex < TD_PARNUM; varIndex++) {
            __m512i tdVal = _mm512_srli_epi16(mm512_asc_vec, varIndex);
            __mmask32 m = ~_mm512_cmpeq_epi16_mask(_mm512_and_si512(tdVal, mm512_const_fullone), mm512_const_zero);
            __m512i xormask = _mm512_maskz_set1_epi16(m, ColPartialDrt[varIndex + GC_VARNUM][los]);
            ColVal[varIndex + GC_VARNUM] = _mm512_xor_si512(ColVal[varIndex + GC_VARNUM], xormask);
        }
        for (int varIndex = RES_VARNUM; varIndex < N + 1; varIndex++)
            ColVal[varIndex] = _mm512_xor_si512(_mm512_set1_epi16(ColPartialDrt[varIndex][los]), ColVal[varIndex]);

        // update gray code and bin key
        gcKeyVal = (binKeyVal >> 1) ^ binKeyVal;
        binKeyVal++;

        if (loopNum++ % 0x2000000 == 0)
            std::cout << "\t" << currentDateTime() << "\t[info]\tbruteforcing... process: " << int(loopNum / 0x2000000) << "/128."
                      << std::endl;
    }
    std::cout << "\t" << currentDateTime() << "\t[info]\tsearch complete. exit." << std::endl;
}
