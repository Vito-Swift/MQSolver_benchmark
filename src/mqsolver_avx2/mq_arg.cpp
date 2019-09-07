// Filename:
//       mq_arg.cpp
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       implementation of miscellaneous functions related to
//       multivariate quadratic problems based on AVX2 intrinsic
//       set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#include "mqsolver_avx2/mq_arg.h"

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

__m256i validateResult(__m256i guess[N]) {
    __m256i ret = _mm256_setzero_si256();
    for (int eqIndex = 0; eqIndex < M; eqIndex++) {
        __m256i v = _mm256_setzero_si256();
        for (int varIndex = 0; varIndex < N; varIndex++) {
            for (int termIndex = 0; termIndex < VerifyPoly[eqIndex][varIndex].length; termIndex++) {
                uint64_t los = 63 - _lzcnt_u64(VerifyPoly[eqIndex][varIndex].p[termIndex]);
                v = _mm256_xor_si256(_mm256_and_si256(guess[los], guess[varIndex]), v);
            }
            if (VerifyPoly[eqIndex][N].length)
                v = _mm256_xor_si256(v, mm256_const_fullone);
            ret = _mm256_or_si256(v, ret);

            if (!avx2_checkcontain_zero_epi32(ret))
                return ret;
        }
    }
    return ret;
}

void mqInit(const std::string residentFileName, const std::string verifyFileName) {
    // load resident partial derivative from file
    FILE *ResidentFile = fopen(residentFileName.c_str(), "rb");
    ffile2poly(ResidentFile, ResidentPoly, LIN_EQNUM);
    loadPD(ResidentPoly, partialDerivative);
    fclose(ResidentFile);
    std::cout << "\t" << currentDateTime() << "\t[info]\tpartial derivative loaded from" << residentFileName
              << std::endl;

    // load verify polynomial from challenge file
    FILE *VerifyFile = fopen(verifyFileName.c_str(), "rb");
    ffile2poly(VerifyFile, VerifyPoly, M);
    fclose(VerifyFile);
    std::cout << "\t" << currentDateTime() << "\t[info]\tverify polynomial loaded from" << verifyFileName << std::endl;

    for (int eqIndex = 0; eqIndex < LIN_EQNUM; eqIndex++) {
        for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++) {
            if ((grayCodeKey >> varIndex) & 1) {
                for (int termIndex = 0; termIndex < ResidentPoly[eqIndex][varIndex].length; termIndex++) {
                    uint64_t los = 63 - _lzcnt_u64(ResidentPoly[eqIndex][varIndex].p[termIndex]);
                    if (los == varIndex)
                        ColVal[N] ^= 1 << eqIndex;
                    else if (los < RES_VARNUM)
                        ColVal[N] ^= ((grayCodeKey >> los) & 1) << eqIndex;
                    else
                        ColVal[los] ^= 1 << eqIndex;
                }
            }
        }
        for (int varIndex = RES_VARNUM; varIndex < N; varIndex++)
            if (ResidentPoly[eqIndex][varIndex].length)
                ColVal[varIndex] ^= 1 << eqIndex;
        ColVal[N] ^= ResidentPoly[eqIndex][N].length << eqIndex;
    }
}

void mqLoop() {
    uint64_t loopNum = 0;
    std::cout << "\t" << currentDateTime() << "\t[info]\tstart mq main loop" << std::endl;
    for (; binKey < KEY_END; binKey++) {
        // set matrix value
        uint32_t Val[LIN_VARNUM + 1] = {0};
        std::copy(ColVal + RES_VARNUM, ColVal + N + 1, Val);

        for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++)
            Val[LIN_VARNUM] ^= ColVal[varIndex];

        // check consistency
        uint32_t mask = 0xffff;
        checkConsist(Val, mask);

        // check result

        // update value
        uint64_t los = __builtin_ffsl(((uint64_t) binKey >> 1) ^ binKey ^ grayCodeKey) - 1;
        for (int varIndex = 0; varIndex < RES_VARNUM; varIndex++)
            if ((grayCodeKey >> varIndex) & 1)
                ColVal[varIndex] ^= partialDerivative[varIndex][los];
        for (int varIndex = RES_VARNUM; varIndex < N + 1; varIndex++)
            ColVal[varIndex] ^= partialDerivative[varIndex][los];

        // update graycode and binary key
        grayCodeKey = (binKey >> 1) ^ binKey;

        if (loopNum++ % 0x4000000 == 0)
            std::cout << "\t" << currentDateTime() << "\t[info]\tprocess: " << int(loopNum / 0x4000000) << "/4096."
                      << std::endl;
    }
}