// Filename:
//       mqsolver.cpp
// Author:
//       Chenhao Wu (chenhaowu@link.cuhk.edu.cn)
//
// Start Date: Sep 17, 2019             Last Change: Sep 17, 2019
//
// Description:
//       Main executable for testing MQ problem on M = N = 46
//
// Remark:
//       This code is licensed under GNU license.
//       see LICENSE.txt for details.

#include "mqsolver.h"

int main() {
    mqInit("mq46_resident.txt", "mq46.txt");
    mqLoop();
    
    return EXIT_SUCCESS;
}