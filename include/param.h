// Filename:
//       param.h
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       header file to define basic parameters in the mq problem
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#ifndef MQSOLVER_PARAM_H
#define MQSOLVER_PARAM_H

#include <iostream>

// N = RES_VARNUM + LIN_VARNUM;
#define M 46
#define N 46
#define RES_VARNUM 37
#define LIN_VARNUM 9
#define LIN_COLNUM 10
#define LIN_EQNUM 10
#define LIN_MASK 0x3ff

#define KEY_START 0
#define KEY_END 0x2000000000

const std::string residentFileName = "mq46-resident.txt";
const std::string verifyFileName = "mq46.txt";

#endif //MQSOLVER_PARAM_H
