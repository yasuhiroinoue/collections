//rlg.cpp
#ifndef DEFS_H
#define DEFS_H

#include "rclass.h"
void Initialize(void);
void Propagation(double& SX);
void CellCalculation(void);
void Collision(void);
void RotationMatrix(void);
//è’ìÀä÷êî
void ColorCollision(void);
void MultiColorCollision(void);
void Scalling(void);
void Collision(void);
//
void Mean(void);
void SystemMonitor(int step);
void Output(int step, double SX);
void AVR_output(int step,double& SX);
void Stress_Measure(void);
void Ptcl_Output(void);
void Read_Ptcl(void);

//bc.cpp
void Periodic(int& x, int& y);
void Periodic(double &x,double &y);
void Wall(_ptcl_base &p,double& x,double& y,double& SX);

#endif
