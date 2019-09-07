// Filename:
//       lin_func.cpp
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       implementation of miscellaneous functions related to linear
//       systems over GF(2) based on AVX512 intrinsic set.
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#include "mqsolver_avx512/lin_func.h"

void extractSolution(const uint16_t clist[LIN_COLNUM], uint32_t sol[LIN_VARNUM]) {
    for (int i = 0; i < LIN_VARNUM; i++) {
        if (clist[i] == 0) continue;
        uint32_t xp = 0x80000000 >> __builtin_clz(clist[i]);
        sol[i] = (bool) (clist[LIN_COLNUM] & xp);
    }
}
