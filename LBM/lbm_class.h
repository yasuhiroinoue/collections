/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef LBM_CLASS_H
#define LBM_CLASS_H
#include "vec.h"
#include "system.h"

struct _cell{
	double f[COL_NUM][9];
	double feq[COL_NUM][9];
	double foc[COL_NUM][9];
	double rho[COL_NUM];
	_vec<double> vel[COL_NUM];
	_vec<double> vel_total;
	double rho_total;
	int col[COL_NUM];
	int col_end;
	
	double s[3];//stress --> s
};

struct _wall{
	int act;
	double rho;
};

struct _pcell{
	double f[COL_NUM][9];
	int col[COL_NUM];
	int col_end;
};

/*
struct _cell{
	double f[3][9];
	double pf[3][9];
	double feq[3][9];
	double rho[3];
	_vec<double> vel[3];
	int col[3];
	int pcol[3];
	int col;
	int pcol;
};
*/
/*
struct _cell{
	double f[9];
	double pf[9];
	double feq[9];
	double rho;
	_vec<double> vel;
};
*/

#endif
