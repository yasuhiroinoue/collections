/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef LBM_H
#define LBM_H

#include "lbm_class.h"
#include "system.h"

_cell CELL[X][Y];
_pcell PCELL[X][Y];
_vec<int> STM[9];//<-- STreaMing
_vec<double> EVEL[9];//<-- E_i, Unit Velocity Vector
double WT[9];//<-- WeighT

_wall WALL[X][Y];
#endif
