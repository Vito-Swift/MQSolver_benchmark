// Filename:
//       mqsolver.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       aggregate library header file
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#ifndef MQSOLVER_MQSOLVER_H
#define MQSOLVER_MQSOLVER_H

#ifdef MODE_AVX512
#include "mqsolver_avx512/avx512_intrinsic_ext.h"
#include "mqsolver_avx512/lin_func.h"
#include "mqsolver_avx512/mq_arg.h"
#endif // MODE_AVX512
#ifdef MODE_AVX2
#include "mqsolver_avx2/avx2_intrinsic_ext.h"
#include "mqsolver_avx2/lin_func.h"
#include "mqsolver_avx2/mq_arg.h"
#endif // MODE_AVX2

#include "param.h"
#include "utils.h"

#endif //MQSOLVER_MQSOLVER_H
