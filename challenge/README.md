This directory stores the implementation on mq challenge Type 4 N = 55. To see the submission, one can refer to [Submission Details](https://www.mqchallenge.org/details/details_IV_20190908_1.html).

Note: for the test program in each iteration all the solution remained in the vector will be verified once, no matter whether the linear system is solvable. If the linear system has a huge probability to be solvable, then that would not bring extra computation time significantly. But if the linear system tends to be unsolvable, it definitely slows down the computation speed. 

The linear systems in Type 4 N = 55 problems have 9 variables and 16 equations, which means it has a huge probability to be unsolvable. By checking the consistency of the linear system, one can easily boost the execution speed. The test program is to show an extreme case that each candidates *has to be* verified.

```bash
$ g++ -Ofast -march=native mqsolver.cpp mqsolver.h -o main
$ ./main
```