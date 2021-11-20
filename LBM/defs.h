/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef DEFS_H
#define DEFS_H

#include "vec.h"

void Init(void);
void Streaming(void);
void CellCalc(void);
void Stress_Measure(void);
void Heat_bath(void);
void Collision(void);
void Output(int step);


void feq_calc(_vec<double> vel, double* feq, double rho);
void Periodic(int& x, int& y);
void Boundary(void);
#endif
