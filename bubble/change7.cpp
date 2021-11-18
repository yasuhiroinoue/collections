/*
// For contact information
// Yasuhiro Inoue
// inoue.yasuhiro.4n@kyoto-u.ac.jp
//
// super minor revision (ext: plt --> tec)
//
*/
//////////////////////////////////////////////
///////			Ishida Kazuki		//////////
///////	ç≈èIçXêVì˙	2008/08/05		//////////
///////			PFñ@				//////////
//////////////////////////////////////////////

////////////////////////////////////////////////////
//ëäïœâªìÒëäó¨
//collocatedäiéq
////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
//#include "mpi.h"


//#define WALLX
#define WALLY
//#define WALLZ
//#define USEMPI //MPI
#define XMAX 50.0
#define YMAX 50.0
#define XCELL_NUM 50	//cellêî
#define YCELL_NUM 50	//cellêî
#define dt	0.2	//éûä‘çèÇ›ïù
#define RHOmax 0.405	//âtëÃÇÃñßìx
#define RHOmin 0.265	//ãCëÃÇÃñßìx
#define MUL 0.2		//îSê´åWêî Re=100
#define MUG 0.2	//ãCëÃÇÃîSê´åWêî
#define gammas 0.0	//îGÇÍê´
#define vdWa 1.0	//vdW ModelÇÃíËêî
#define vdWb 1.0	//vdW ModelÇÃíËêî
#define vdWT 0.293	//vdW ModelÇÃâ∑ìx
#define Ts 0.293
#define deltaT 1.465e-2// 0.00293
#define vdWc 1.0
#define kappa1 0.01 //2.94e-2
#define kappa2 0.2
#define mack 0.1
#define kt 0.2

struct _cell_type {
	double v_x;	//xï˚å¸ó¨ë¨
	double v_y;	//yï˚å¸ó¨ë¨
	double momentum_x;
	double momentum_y;
	double momentum_xnew;
	double momentum_ynew;
	double press;		//à≥óÕ
	double pressnew;
	double v_xnew;	//âºÇÃxï˚å¸ó¨ë¨
	double v_ynew;	//âºÇÃyï˚å¸ó¨ë¨
	double x;	//xç¿ïW
	double y;	//yç¿ïW
	double rho;	//ñßìx
	double rhonew;
	double temper;	//â∑ìx
	double tempernew;
	double energy;
	double energynew;
	double mu;
};

double	dx = XMAX / (double) XCELL_NUM;
double	dy = YMAX / (double) YCELL_NUM;
double dxi = 1.0/dx;
double dyi = 1.0/dy;
int CELL_SIZE = (XCELL_NUM+6)*(YCELL_NUM+6);

int n,nout,TCELL_NUM;


void indexsetdefine(int i, int j);
void indexSORdefine(int i, int j);
void indexnseq(int i, int j);
void indexbndx(int j);
void indexbndy(int i);
void initial(void);
void hontai(int chkflag);
void maccormack(int chkflag);
void nseq(void);
void vnewbnd(void);
void pressnewbnd(void);
void momentumnewbnd(void);
void rhonewbnd(void);
void tempernewbnd(void);
void energynewbnd(void);
void vbnd(void);
void pressbnd(void);
void momentumbnd(void);
void rhobnd(void);
void temperbnd(void);
void energybnd(void);
void filewrite(void);

_cell_type* CELL;
double* ddrho;
double wetrho;
double pressstart;
double dti = 1.0/dt;

int index0;
int index_xm;
int index_xp;
int index_yp;
int index_ym;
int index_xpp;
int index_ypp;
int index_xmm;
int index_ymm;
int index_xp_yp;
int index_xm_ym;
int index_xm_yp;
int index_xp_ym;
int index_x0;
int index_x1;
int index_x2;
int index_x3;
int index_x4;
int index_x5;
int index_xcp0;
int index_xcp1;
int index_xcp2;
int index_xcp3;
int index_xcp4;
int index_xcp5;
int index_x6;
			
_cell_type* cellp_x0;
_cell_type* cellp_x1;
_cell_type* cellp_x2;
_cell_type* cellp_x3;
_cell_type* cellp_x4;
_cell_type* cellp_x5;
_cell_type* cellp_xcp0;
_cell_type* cellp_xcp1;
_cell_type* cellp_xcp2;
_cell_type* cellp_xcp3;
_cell_type* cellp_xcp4;
_cell_type* cellp_xcp5;
_cell_type* cellp_x6;
int index_y0;
int index_y1;
int index_y2;
int index_y3;
int index_y4;
int index_y5;
int index_y6;
int index_ycp0;
int index_ycp1;
int index_ycp2;
int index_ycp3;
int index_ycp4;
int index_ycp5;
			
_cell_type* cellp_y0;
_cell_type* cellp_y1;
_cell_type* cellp_y2;
_cell_type* cellp_y3;
_cell_type* cellp_y4;
_cell_type* cellp_y5;
_cell_type* cellp_y6;
_cell_type* cellp_ycp0;
_cell_type* cellp_ycp1;
_cell_type* cellp_ycp2;
_cell_type* cellp_ycp3;
_cell_type* cellp_ycp4;
_cell_type* cellp_ycp5;


_cell_type* cellp;
_cell_type* cellp_xp;
_cell_type* cellp_yp;
_cell_type* cellp_xm;
_cell_type* cellp_ym;
_cell_type* cellp_xpp;
_cell_type* cellp_ypp;
_cell_type* cellp_xmm;
_cell_type* cellp_ymm;
_cell_type* cellp_xp_yp;
_cell_type* cellp_xm_ym;
_cell_type* cellp_xm_yp;
_cell_type* cellp_xp_ym;



void indexsetdefine(int i, int j){
	index0     = j  + i*(YCELL_NUM+6);
	index_xm   = j + (i-1)*(YCELL_NUM+6);
	index_xp  =  j     + (i+1)*(YCELL_NUM+6);
	index_yp  =  (j+1) + i*(YCELL_NUM+6);
	index_ym  =  (j-1) + i*(YCELL_NUM+6);
	index_xpp =  j     + (i+2)*(YCELL_NUM+6);
	index_ypp =  (j+2) + i*(YCELL_NUM+6);
	index_xmm =  j     + (i-2)*(YCELL_NUM+6);
	index_ymm =  (j-2) + i*(YCELL_NUM+6);
	
	
	cellp= &CELL[index0];
	cellp_xp = &CELL[index_xp];
	cellp_yp = &CELL[index_yp];
	cellp_xm = &CELL[index_xm];
	cellp_ym = &CELL[index_ym];
	cellp_xpp = &CELL[index_xpp];
	cellp_ypp = &CELL[index_ypp];
	cellp_xmm = &CELL[index_xmm];
	cellp_ymm = &CELL[index_ymm];
}

void indexSORdefine(int i, int j){
	index0   = j     + i*(YCELL_NUM+6);
	index_xm = j     + (i-1)*(YCELL_NUM+6);
	index_xp = j     + (i+1)*(YCELL_NUM+6);
	index_yp = (j+1) + i*(YCELL_NUM+6);
	index_ym = (j-1) + i*(YCELL_NUM+6);
					
	cellp= &CELL[index0];
	cellp_xp = &CELL[index_xp];
	cellp_yp = &CELL[index_yp];
	cellp_xm = &CELL[index_xm];
	cellp_ym = &CELL[index_ym];
}

	
void indexnsdefine(int i, int j){
	index0      = j     + i*(YCELL_NUM+6);
	index_xp    = j     + (i+1)*(YCELL_NUM+6);
	index_xpp   = j     + (i+2)*(YCELL_NUM+6);
	index_xm    = j     + (i-1)*(YCELL_NUM+6);
	index_xmm   = j     + (i-2)*(YCELL_NUM+6);
	index_yp    = (j+1) + i*(YCELL_NUM+6);
	index_ypp   = (j+2) + i*(YCELL_NUM+6);
	index_ym    = (j-1) + i*(YCELL_NUM+6);
	index_ymm   = (j-2) + i*(YCELL_NUM+6);
	index_xp_yp = (j+1) + (i+1)*(YCELL_NUM+6);
	index_xm_ym = (j-1) + (i-1)*(YCELL_NUM+6);
	index_xm_yp = (j+1) + (i-1)*(YCELL_NUM+6);
	index_xp_ym = (j-1) + (i+1)*(YCELL_NUM+6);
				
	cellp       = &CELL[index0];
	cellp_xp    = &CELL[index_xp];
	cellp_yp    = &CELL[index_yp];
	cellp_xm    = &CELL[index_xm];
	cellp_ym    = &CELL[index_ym];
	cellp_xpp   = &CELL[index_xpp];
	cellp_xmm   = &CELL[index_xmm];
	cellp_ypp   = &CELL[index_ypp];
	cellp_ymm   = &CELL[index_ymm];
	cellp_xp_yp = &CELL[index_xp_yp];
	cellp_xm_ym = &CELL[index_xm_ym];
	cellp_xm_yp = &CELL[index_xm_yp];
	cellp_xp_ym = &CELL[index_xp_ym];
}


void indexbndx(int j){
	index_x0   = j;
	index_x1   = j +   (YCELL_NUM+6);
	index_x2   = j + 2*(YCELL_NUM+6);
	index_x3   = j + 3*(YCELL_NUM+6);
	index_x4   = j + 4*(YCELL_NUM+6);
	index_x5   = j + 5*(YCELL_NUM+6);
	index_xcp0 = j + (XCELL_NUM)*(YCELL_NUM+6);
	index_xcp1 = j + (XCELL_NUM+1)*(YCELL_NUM+6);
	index_xcp2 = j + (XCELL_NUM+2)*(YCELL_NUM+6);
	index_xcp3 = j + (XCELL_NUM+3)*(YCELL_NUM+6);
	index_xcp4 = j + (XCELL_NUM+4)*(YCELL_NUM+6);
	index_xcp5 = j + (XCELL_NUM+5)*(YCELL_NUM+6);
	index_x6   = j + 6*(YCELL_NUM+6);
				
	cellp_x0   = &CELL[index_x0];
	cellp_x1   = &CELL[index_x1];
	cellp_x2   = &CELL[index_x2];
	cellp_x3   = &CELL[index_x3];
	cellp_x4   = &CELL[index_x4];
	cellp_x5   = &CELL[index_x5];
	cellp_xcp0 = &CELL[index_xcp0];
	cellp_xcp1 = &CELL[index_xcp1];
	cellp_xcp2 = &CELL[index_xcp2];
	cellp_xcp3 = &CELL[index_xcp3];
	cellp_xcp4 = &CELL[index_xcp4];
	cellp_xcp5 = &CELL[index_xcp5];
	cellp_x6   = &CELL[index_x6];
}


void indexbndy(int i){
	index_y0   =                    i*(YCELL_NUM+6);
	index_y1   = 1 + i*(YCELL_NUM+6);
	index_y2   = 2 + i*(YCELL_NUM+6);
	index_y3   = 3 + i*(YCELL_NUM+6);
	index_y4   = 4 + i*(YCELL_NUM+6);
	index_y5   = 5 + i*(YCELL_NUM+6);
	index_y6   = 6 + i*(YCELL_NUM+6);
	index_ycp0 = YCELL_NUM + i*(YCELL_NUM+6);
	index_ycp1 = (YCELL_NUM+1) + i*(YCELL_NUM+6);
	index_ycp2 = (YCELL_NUM+2) + i*(YCELL_NUM+6);
	index_ycp3 = (YCELL_NUM+3) + i*(YCELL_NUM+6);
	index_ycp4 = (YCELL_NUM+4) + i*(YCELL_NUM+6);
	index_ycp5 = (YCELL_NUM+5) + i*(YCELL_NUM+6);
				
	cellp_y0   = &CELL[index_y0];
	cellp_y1   = &CELL[index_y1];
	cellp_y2   = &CELL[index_y2];
	cellp_y3   = &CELL[index_y3];
	cellp_y4   = &CELL[index_y4];
	cellp_y5   = &CELL[index_y5];
	cellp_y6   = &CELL[index_y6];
	cellp_ycp0 = &CELL[index_ycp0];
	cellp_ycp1 = &CELL[index_ycp1];
	cellp_ycp2 = &CELL[index_ycp2];
	cellp_ycp3 = &CELL[index_ycp3];
	cellp_ycp4 = &CELL[index_ycp4];
	cellp_ycp5 = &CELL[index_ycp5];
}


int main(void){
	int restartout;
	int chkflag;
	TCELL_NUM=100000;
	nout=1000;
	restartout=10000;
	chkflag=0;
	n=0;
	
	initial();
	filewrite();
	
	for(; n<(TCELL_NUM+1); n++) {
		hontai(chkflag);
		if(n==1) filewrite();
		if(n%nout==0) {
			filewrite();
//			masscount();
		}
		if(n%restartout==0){
//			restartwrite();
		}
	}
	
	
	delete [] CELL;
	delete [] ddrho;
	return 0;
}

			//-------------- mainä÷êîèIÇÌÇË ----------------//

			//-------------- èâä˙âª -----------------------//
void initial(void){
	CELL = new _cell_type[CELL_SIZE];
	memset(CELL,0,CELL_SIZE*sizeof(_cell_type));
	ddrho = new double[CELL_SIZE];
	memset(ddrho,0,CELL_SIZE*sizeof(double));
	
//	wetrho = gammas/kappa1;
	wetrho = 0.0;
	
	for(int i=0; i<(XCELL_NUM+6); i++){
		for(int j=0; j<(YCELL_NUM+6); j++){
			int index = j + i*(YCELL_NUM+6);
			
			_cell_type* cellp= &CELL[index];
			
			cellp->x = ((double)i-2.5)*dx;
			cellp->y = ((double)j-2.5)*dy;
		}
	}
	
	for(int i=0; i<(XCELL_NUM+6); i++) {
		for(int j=0; j<(YCELL_NUM+6); j++) {
			index0 = j + i*(YCELL_NUM+6);
			_cell_type* cellp= &CELL[index0];
			
			cellp->rho = RHOmax;
			cellp->rhonew = RHOmax;
			cellp->temper = Ts;
			cellp->mu = MUL;
		}
	}
	
	vbnd();
	temperbnd();
	rhobnd();
	
	for(int i=2; i<(XCELL_NUM+5); i++) {
		for(int j=2; j<(YCELL_NUM+5); j++) {
			indexnsdefine(i,j);
			
			cellp->momentum_x = cellp->rho*cellp->v_x;
			cellp->momentum_y = cellp->rho*cellp->v_y;
			cellp->press = cellp->rho*cellp->temper/(1.0-vdWb*cellp->rho) - vdWa*cellp->rho*cellp->rho 
								- kappa1*cellp->rho*((cellp_xp->rho-2.0*cellp->rho+cellp_xm->rho)*dxi*dxi+(cellp_yp->rho-2.0*cellp->rho+cellp_ym->rho)*dyi*dyi)
								- 0.5*kappa1*((cellp_xp->rho-cellp_xm->rho)*dxi*0.5*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5
								+(cellp_yp->rho-cellp_ym->rho)*dyi*0.5*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5);
			cellp->energy = 0.5*cellp->rho*(cellp->v_x*cellp->v_x+cellp->v_y*cellp->v_y) + cellp->rho*(vdWc*cellp->temper - vdWa*cellp->rho) 
							+ 0.5*kappa1*((cellp_xp->rho-cellp_xm->rho)*dxi*0.5*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5
							+(cellp_yp->rho-cellp_ym->rho)*dyi*0.5*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5);
		
			if((i==25)&&(j==25)) pressstart=cellp->press;
		}
	}
	
	printf("press = %f \n",pressstart);
	
	pressbnd();
	momentumbnd();
	energybnd();
}


			//-------------- èâä˙âªèIÇÌÇË -----------------//


void hontai(int chkflag){
	maccormack(chkflag);
}

void maccormack(int chkflag){
	
	nseq();
	
}


void nseq(void){
	double macFxp,macFx,macFyp,macFy,macFxm,macFym;
	double macQxp,macQx,macQyp,macQy,macQxm,macQym;
	
	
	
	for(int i=3; i<=(XCELL_NUM+3); i++){
		for(int j=3; j<=(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
			if(j!=3){
				macFxp = cellp_xp->momentum_x;
				macFx  = cellp->momentum_x;
				macFyp = cellp_yp->momentum_y;
				macFy  = cellp->momentum_y;
				macQxp = 0.125*mack*(cellp_xpp->rho-2.0*cellp_xp->rho+cellp->rho);
				macQx  = 0.125*mack*(cellp_xp->rho-2.0*cellp->rho+cellp_xm->rho);
				macQyp = 0.125*mack*(cellp_ypp->rho-2.0*cellp_yp->rho+cellp->rho);
				macQy  = 0.125*mack*(cellp_yp->rho-2.0*cellp->rho+cellp_ym->rho);
				
				cellp->rhonew = cellp->rho - dt*(dxi*(macFxp - macFx) + dyi*(macFyp - macFy) + dxi*(macQxp - macQx) + dyi*(macQyp - macQy));
			}
			else {
				cellp->rhonew = cellp->rho - 0.5*dyi*dt*(4.0*cellp_yp->momentum_y-cellp_ypp->momentum_y);
			}
//			if((cellp_xp->rho<RHOmax)&&(cellp_xp->rho>RHOmin)){
				macFxp = cellp_xp->momentum_x*cellp_xp->v_x + cellp_xp->press + kappa1*(cellp_xpp->rho-cellp->rho)*(cellp_xpp->rho-cellp->rho)*dxi*dxi*0.5*0.5
						- 4.0/3.0*cellp_xp->mu*(cellp_xpp->v_x - cellp->v_x)*dxi*0.5+2.0/3.0*cellp_xp->mu*(cellp_xp_yp->v_y - cellp_xp_ym->v_y)*dyi*0.5;
//			}
//			else{
//				macFxp = cellp_xp->momentum_x*cellp_xp->v_x + cellp_xp->press
//						- 4.0/3.0*cellp_xp->mu*(cellp_xpp->v_x - cellp->v_x)*dxi*0.5+2.0/3.0*cellp_xp->mu*(cellp_xp_yp->v_y - cellp_xp_ym->v_y)*dyi*0.5;
//			}
//			if((cellp_yp->rho<RHOmax)&&(cellp_yp->rho>RHOmin)){
				macFyp = cellp_yp->momentum_y*cellp_yp->v_x + kappa1*(cellp_xp_yp->rho-cellp_xm_yp->rho)*(cellp_ypp->rho-cellp->rho)*dxi*dyi*0.5*0.5
						- cellp_yp->mu*(cellp_ypp->v_x - cellp->v_x)*dyi*0.5-cellp_yp->mu*(cellp_xp_yp->v_y - cellp_xm_yp->v_y)*dxi*0.5;
//			}
//			else{
//				macFyp = cellp_yp->momentum_y*cellp_yp->v_x
//						- cellp_yp->mu*(cellp_ypp->v_x - cellp->v_x)*dyi*0.5-cellp_yp->mu*(cellp_xp_yp->v_y - cellp_xm_yp->v_y)*dxi*0.5;
//			}
//			if((cellp->rho<RHOmax)&&(cellp->rho>RHOmin)){
				macFx  = cellp->momentum_x*cellp->v_x + cellp->press + kappa1*(cellp_xp->rho-cellp_xm->rho)*(cellp_xp->rho-cellp_xm->rho)*dxi*dxi*0.5*0.5
						- 4.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5+2.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5;
				macFy  = cellp->momentum_y*cellp->v_x + kappa1*(cellp_xp->rho-cellp_xm->rho)*(cellp_yp->rho-cellp_ym->rho)*dxi*dyi*0.5*0.5 
						- cellp->mu*(cellp_yp->v_x - cellp_ym->v_x)*dyi*0.5-cellp->mu*(cellp_xp->v_y - cellp_xm->v_y)*dxi*0.5;
//			}
//			else{
//				macFx  = cellp->momentum_x*cellp->v_x + cellp->press
//						- 4.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5+2.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5;
//				macFy  = cellp->momentum_y*cellp->v_x
//						- cellp->mu*(cellp_yp->v_x - cellp_ym->v_x)*dyi*0.5-cellp->mu*(cellp_xp->v_y - cellp_xm->v_y)*dxi*0.5;
//			}
			macQxp = 0.125*mack*(cellp_xpp->momentum_x-2.0*cellp_xp->momentum_x+cellp->momentum_x);
			macQx  = 0.125*mack*(cellp_xp->momentum_x-2.0*cellp->momentum_x+cellp_xm->momentum_x);
			macQyp = 0.125*mack*(cellp_ypp->momentum_x-2.0*cellp_yp->momentum_x+cellp->momentum_x);
			macQy  = 0.125*mack*(cellp_yp->momentum_x-2.0*cellp->momentum_x+cellp_ym->momentum_x);
			
			cellp->momentum_xnew = cellp->momentum_x - dt*(dxi*(macFxp - macFx) + dyi*(macFyp - macFy) + dxi*(macQxp - macQx) + dyi*(macQyp - macQy));
			
//			if((cellp_xp->rho<RHOmax)&&(cellp_xp->rho>RHOmin)){
				macFxp = cellp_xp->momentum_x*cellp_xp->v_y + kappa1*(cellp_xp_yp->rho-cellp_xp_ym->rho)*(cellp_xpp->rho-cellp->rho)*dxi*dyi*0.5*0.5
						- cellp_xp->mu*(cellp_xpp->v_y - cellp->v_y)*dxi*0.5-cellp_xp->mu*(cellp_xp_yp->v_x - cellp_xp_ym->v_x)*dyi*0.5;
//			}
//			else{
//				macFxp = cellp_xp->momentum_x*cellp_xp->v_y
//						- cellp_xp->mu*(cellp_xpp->v_y - cellp->v_y)*dxi*0.5-cellp_xp->mu*(cellp_xp_yp->v_x - cellp_xp_ym->v_x)*dyi*0.5;
//			}
//			if((cellp_yp->rho<RHOmax)&&(cellp_yp->rho>RHOmin)){
				macFyp = cellp_yp->momentum_y*cellp_yp->v_y + cellp_yp->press + kappa1*(cellp_ypp->rho-cellp->rho)*(cellp_ypp->rho-cellp->rho)*dyi*dyi*0.5*0.5
						- 4.0/3.0*cellp_yp->mu*(cellp_ypp->v_y - cellp->v_y)*dyi*0.5+2.0/3.0*cellp_yp->mu*(cellp_xp_yp->v_x - cellp_xm_yp->v_x)*dxi*0.5;
//			}
//			else{
//				macFyp = cellp_yp->momentum_y*cellp_yp->v_y + cellp_yp->press
//						- 4.0/3.0*cellp_yp->mu*(cellp_ypp->v_y - cellp->v_y)*dyi*0.5+2.0/3.0*cellp_yp->mu*(cellp_xp_yp->v_x - cellp_xm_yp->v_x)*dxi*0.5;
//			}
//			if((cellp->rho<RHOmax)&&(cellp->rho>RHOmin)){
				macFx  = cellp->momentum_x*cellp->v_y + kappa1*(cellp_yp->rho-cellp_ym->rho)*(cellp_xp->rho-cellp_xm->rho)*dxi*dyi*0.5*0.5
						- cellp->mu*(cellp_xp->v_y - cellp_xm->v_y)*dxi*0.5-cellp->mu*(cellp_yp->v_x - cellp_ym->v_x)*dyi*0.5;
				macFy  = cellp->momentum_y*cellp->v_y + cellp->press + kappa1*(cellp_yp->rho-cellp_ym->rho)*(cellp_yp->rho-cellp_ym->rho)*dyi*dyi*0.5*0.5
						- 4.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5+2.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5;
//			}
//			else{
//				macFx  = cellp->momentum_x*cellp->v_y
//						- cellp->mu*(cellp_xp->v_y - cellp_xm->v_y)*dxi*0.5-cellp->mu*(cellp_yp->v_x - cellp_ym->v_x)*dyi*0.5;
//				macFy  = cellp->momentum_y*cellp->v_y + cellp->press
//						- 4.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5+2.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5;
//			}
			macQxp = 0.125*mack*(cellp_xpp->momentum_y-2.0*cellp_xp->momentum_y+cellp->momentum_y);
			macQx  = 0.125*mack*(cellp_xp->momentum_y-2.0*cellp->momentum_y+cellp_xm->momentum_y);
			macQyp = 0.125*mack*(cellp_ypp->momentum_y-2.0*cellp_yp->momentum_y+cellp->momentum_y);
			macQy  = 0.125*mack*(cellp_yp->momentum_y-2.0*cellp->momentum_y+cellp_ym->momentum_y);
			
			cellp->momentum_ynew = cellp->momentum_y - dt*(dxi*(macFxp - macFx) + dyi*(macFyp - macFy) + dxi*(macQxp - macQx) + dyi*(macQyp - macQy));
			
//			if((cellp_xp->rho<RHOmax)&&(cellp_xp->rho>RHOmin)){
				macFxp = cellp_xp->energy*cellp_xp->v_x + cellp_xp->press*cellp_xp->v_x
							+ kappa1*(cellp_xpp->rho-cellp->rho)*(cellp_xpp->rho-cellp->rho)*dxi*dxi*0.5*0.5*cellp_xp->v_x 
							+ kappa1*(cellp_xpp->rho-cellp->rho)*(cellp_xp_yp->rho-cellp_xp_ym->rho)*dxi*dyi*0.5*0.5*cellp_xp->v_y
							- ((4.0/3.0*cellp_xp->mu*(cellp_xpp->v_x - cellp->v_x)*dxi*0.5-2.0/3.0*cellp_xp->mu*(cellp_xp_yp->v_y - cellp_xp_ym->v_y)*dyi*0.5)*cellp_xp->v_x
							+(cellp_xp->mu*(cellp_xp_yp->v_x - cellp_xp_ym->v_x)*dyi*0.5+cellp_xp->mu*(cellp_xpp->v_y - cellp->v_y)*dxi*0.5)*cellp_xp->v_y)
							- kt*(cellp_xpp->temper - cellp->temper)*dxi*0.5 
							+ kappa1*cellp_xp->rho*((cellp_xpp->v_x-cellp->v_x)*dxi*0.5+(cellp_xp_yp->v_y-cellp_xp_ym->v_y)*dyi*0.5)*(cellp_xpp->rho-cellp->rho)*dxi*0.5;
//			}
//			else{
//				macFxp = cellp_xp->energy*cellp_xp->v_x + cellp_xp->press*cellp_xp->v_x
//							- (4.0/3.0*cellp_xp->mu*(cellp_xpp->v_x - cellp->v_x)*dxi*0.5-2.0/3.0*cellp_xp->mu*(cellp_xp_yp->v_y - cellp_xp_ym->v_y)*dyi*0.5)*cellp_xp->v_x
//							- kt*(cellp_xpp->temper - cellp->temper)*dxi*0.5;
//			}
//			if((cellp_yp->rho<RHOmax)&&(cellp_yp->rho>RHOmin)){
				macFyp = cellp_yp->energy*cellp_yp->v_y + cellp_yp->press*cellp_yp->v_y 
							+ kappa1*(cellp_ypp->rho-cellp->rho)*(cellp_ypp->rho-cellp->rho)*dyi*dyi*0.5*0.5*cellp_yp->v_y 
							+ kappa1*(cellp_ypp->rho-cellp->rho)*(cellp_xp_yp->rho-cellp_xm_yp->rho)*dyi*dxi*0.5*0.5*cellp_yp->v_x
							- ((4.0/3.0*cellp_yp->mu*(cellp_ypp->v_y - cellp->v_y)*dyi*0.5-2.0/3.0*cellp_yp->mu*(cellp_xp_yp->v_x - cellp_xm_yp->v_x)*dxi*0.5)*cellp_yp->v_y
							+(cellp_yp->mu*(cellp_xp_yp->v_y - cellp_xm_yp->v_y)*dxi*0.5+cellp_yp->mu*(cellp_ypp->v_x - cellp->v_x)*dyi*0.5)*cellp_yp->v_x)
							- kt*(cellp_ypp->temper - cellp->temper)*dyi*0.5 
							+ kappa1*cellp_yp->rho*((cellp_ypp->v_y-cellp->v_y)*dyi*0.5+(cellp_xp_yp->v_x-cellp_xm_yp->v_x)*dxi*0.5)*(cellp_ypp->rho-cellp->rho)*dyi*0.5;
//			}
//			else{
//				macFyp = cellp_yp->energy*cellp_yp->v_y + cellp_yp->press*cellp_yp->v_y 
//							- (4.0/3.0*cellp_yp->mu*(cellp_ypp->v_y - cellp->v_y)*dyi*0.5-2.0/3.0*cellp_yp->mu*(cellp_xp_yp->v_x - cellp_xm_yp->v_x)*dxi*0.5)*cellp_yp->v_y
//							- kt*(cellp_ypp->temper - cellp->temper)*dyi*0.5;
//			}
//			if((cellp->rho<RHOmax)&&(cellp->rho>RHOmin)){
				macFx  = cellp->energy*cellp->v_x + cellp->press*cellp->v_x
							+ kappa1*(cellp_xp->rho-cellp_xm->rho)*(cellp_xp->rho-cellp_xm->rho)*dxi*dxi*0.5*0.5*cellp->v_x 
							+ kappa1*(cellp_xp->rho-cellp_xm->rho)*(cellp_yp->rho-cellp_ym->rho)*dxi*dyi*0.5*0.5*cellp->v_y
							- ((4.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5-2.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5)*cellp->v_x
							+(cellp->mu*(cellp_yp->v_x - cellp_ym->v_x)*dyi*0.5+cellp->mu*(cellp_xp->v_y - cellp_xm->v_y)*dxi*0.5)*cellp->v_y)
							- kt*(cellp_xp->temper - cellp_xm->temper)*dxi*0.5 
							+ kappa1*cellp->rho*((cellp_xp->v_x-cellp_xm->v_x)*dxi*0.5+(cellp_yp->v_y-cellp_ym->v_y)*dyi*0.5)*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5;
				macFy  = cellp->energy*cellp->v_y + cellp->press*cellp->v_y 
							+ kappa1*(cellp_yp->rho-cellp_ym->rho)*(cellp_yp->rho-cellp_ym->rho)*dyi*dyi*0.5*0.5*cellp->v_y 
							+ kappa1*(cellp_yp->rho-cellp_ym->rho)*(cellp_xp->rho-cellp_xm->rho)*dyi*dxi*0.5*0.5*cellp->v_x
							- ((4.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5-2.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5)*cellp->v_y
							+(cellp->mu*(cellp_xp->v_y - cellp_xm->v_y)*dxi*0.5+cellp->mu*(cellp_yp->v_x - cellp_ym->v_x)*dyi*0.5)*cellp->v_x)
							- kt*(cellp_yp->temper - cellp_ym->temper)*dyi*0.5 
							+ kappa1*cellp->rho*((cellp_yp->v_y-cellp_ym->v_y)*dyi*0.5+(cellp_xp->v_x-cellp_xm->v_x)*dxi*0.5)*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5;
//			}
//			else{
//				macFx  = cellp->energy*cellp->v_x + cellp->press*cellp->v_x
//							- (4.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5-2.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5)*cellp->v_x
//							- kt*(cellp_xp->temper - cellp_xm->temper)*dxi*0.5;
//				macFy  = cellp->energy*cellp->v_y + cellp->press*cellp->v_y 
//							- (4.0/3.0*cellp->mu*(cellp_yp->v_y - cellp_ym->v_y)*dyi*0.5-2.0/3.0*cellp->mu*(cellp_xp->v_x - cellp_xm->v_x)*dxi*0.5)*cellp->v_y
//							- kt*(cellp_yp->temper - cellp_ym->temper)*dyi*0.5;
//			}
			macQxp = 0.125*mack*(cellp_xpp->energy-2.0*cellp_xp->energy+cellp->energy);
			macQx  = 0.125*mack*(cellp_xp->energy-2.0*cellp->energy+cellp_xm->energy);
			macQyp = 0.125*mack*(cellp_ypp->energy-2.0*cellp_yp->energy+cellp->energy);
			macQy  = 0.125*mack*(cellp_yp->energy-2.0*cellp->energy+cellp_ym->energy);
			
			cellp->energynew = cellp->energy - dt*(dxi*(macFxp - macFx) + dyi*(macFyp - macFy) + dxi*(macQxp - macQx) + dyi*(macQyp - macQy));
		}
	}
	
	
	rhonewbnd();
	momentumnewbnd();
	energynewbnd();
	
	
	for(int i=3; i<(XCELL_NUM+3); i++){
		for(int j=3; j<(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
			ddrho[index0] = (-cellp_xpp->rhonew+16.0*cellp_xp->rhonew-30.0*cellp->rhonew+16.0*cellp_xm->rhonew-cellp_xmm->rhonew)/12.0*dxi*dxi+(-cellp_ypp->rhonew+16.0*cellp_yp->rhonew-30.0*cellp->rhonew+16.0*cellp_ym->rhonew-cellp_ymm->rhonew)/12.0*dyi*dyi;
		}
	}
	
	for(int i=3; i<=(XCELL_NUM+3); i++){
		for(int j=3; j<=(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
			cellp->v_xnew = cellp->momentum_xnew/cellp->rhonew;
			cellp->v_ynew = cellp->momentum_ynew/cellp->rhonew;
		}
	}
	
	vnewbnd();
	
	for(int i=3; i<=(XCELL_NUM+3); i++){
		for(int j=3; j<=(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
//			if((cellp->rhonew<RHOmax)&&(cellp->rhonew>RHOmin)){
				cellp->tempernew = (cellp->energynew - 0.5*cellp->rhonew*(cellp->v_xnew*cellp->v_xnew+cellp->v_ynew*cellp->v_ynew)
								+ cellp->rhonew*vdWa*cellp->rhonew
								- 0.5*kappa1*((cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5
								+(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5*(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5))/(cellp->rhonew*vdWc);
//			}
//			else{
//				cellp->tempernew = (cellp->energynew - 0.5*cellp->rhonew*(cellp->v_xnew*cellp->v_xnew+cellp->v_ynew*cellp->v_ynew)
//								+ cellp->rhonew*vdWa*cellp->rhonew
//								)/(cellp->rhonew*vdWc);
//			}
//			if((cellp->rhonew<RHOmax)&&(cellp->rhonew>RHOmin)){
				cellp->pressnew = cellp->rhonew*cellp->tempernew/(1.0-vdWb*cellp->rhonew) - vdWa*cellp->rhonew*cellp->rhonew 
									- kappa1*cellp->rhonew*ddrho[index0]
									- 0.5*kappa1*((cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5
									+(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5*(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5);
//			}
//			else{
//				cellp->pressnew = cellp->rhonew*cellp->tempernew/(1.0-vdWb*cellp->rhonew) - vdWa*cellp->rhonew*cellp->rhonew;
//			}
			
			
		}
	}
	
	tempernewbnd();
	
	
	for(int i=3; i<=(XCELL_NUM+3); i++){
			indexnsdefine(i,3);
			
				cellp->pressnew = cellp->rhonew*cellp->tempernew/(1.0-vdWb*cellp->rhonew) - vdWa*cellp->rhonew*cellp->rhonew 
									- kappa1*cellp->rhonew*ddrho[index0]
//		- kappa1*cellp->rhonew*((cellp_xp->rhonew-2.0*cellp->rhonew+cellp_xm->rhonew)*(cellp_xp->rhonew-2.0*cellp->rhonew+cellp_xm->rhonew)+(cellp_yp->rhonew-2.0*cellp->rhonew+cellp_ym->rhonew)*(cellp_yp->rhonew-2.0*cellp->rhonew+cellp_ym->rhonew))
									- 0.5*kappa1*((cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5
									+(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5*(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5);
//			}
//			else{
//				cellp->pressnew = cellp->rhonew*cellp->tempernew/(1.0-vdWb*cellp->rhonew) - vdWa*cellp->rhonew*cellp->rhonew;
//			}
	
			cellp->energynew = 0.5*cellp->rhonew*(cellp->v_xnew*cellp->v_xnew+cellp->v_ynew*cellp->v_ynew) + cellp->rhonew*(vdWc*cellp->tempernew - vdWa*cellp->rhonew) 
							+ 0.5*kappa1*((cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5
							+(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5*(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5);
	}
	
	pressnewbnd();
	
	for(int i=3; i<=(XCELL_NUM+3); i++){
		for(int j=3; j<=(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
			if(j!=3){
				macFxm = cellp_xm->momentum_xnew;
				macFx  = cellp->momentum_xnew;
				macFym = cellp_ym->momentum_ynew;
				macFy  = cellp->momentum_ynew;
				macQxm = 0.125*mack*(cellp->rhonew-2.0*cellp_xm->rho+cellp_xmm->rhonew);
				macQx  = 0.125*mack*(cellp_xp->rhonew-2.0*cellp->rhonew+cellp_xm->rhonew);
				macQym = 0.125*mack*(cellp->rhonew-2.0*cellp_ym->rhonew+cellp_ymm->rhonew);
				macQy  = 0.125*mack*(cellp_yp->rhonew-2.0*cellp->rhonew+cellp_ym->rhonew);
				
				cellp->rho = 0.5*((cellp->rho+cellp->rhonew) - dt*(dxi*(macFx - macFxm) + dyi*(macFy - macFym) - dxi*(macQx - macQxm) - dyi*(macQy - macQym)));
			}
			else {
				cellp->rho = 0.5*(cellp->rhonew+cellp->rho - 0.5*dyi*dt*(4.0*cellp_yp->momentum_ynew-cellp_ypp->momentum_ynew));
			}
//			if((cellp_xm->rhonew<RHOmax)&&(cellp_xm->rhonew>RHOmin)){
				macFxm = cellp_xm->momentum_xnew*cellp_xm->v_xnew + cellp_xm->pressnew + kappa1*(cellp->rhonew-cellp_xmm->rhonew)*(cellp->rhonew-cellp_xmm->rhonew)*dxi*dxi*0.5*0.5
						- 4.0/3.0*cellp_xm->mu*(cellp->v_xnew - cellp_xmm->v_xnew)*dxi*0.5+2.0/3.0*cellp_xm->mu*(cellp_xm_yp->v_ynew - cellp_xm_ym->v_ynew)*dyi*0.5;
//			}
//			else{
//				macFxm = cellp_xm->momentum_xnew*cellp_xm->v_xnew + cellp_xm->pressnew
//						- 4.0/3.0*cellp_xm->mu*(cellp->v_xnew - cellp_xmm->v_xnew)*dxi*0.5+2.0/3.0*cellp_xm->mu*(cellp_xm_yp->v_ynew - cellp_xm_ym->v_ynew)*dyi*0.5;
//			}
//			if((cellp_ym->rhonew<RHOmax)&&(cellp_ym->rhonew>RHOmin)){
				macFym = cellp_ym->momentum_ynew*cellp_ym->v_xnew + kappa1*(cellp_xp_ym->rhonew-cellp_xm_ym->rhonew)*(cellp->rhonew-cellp_ymm->rhonew)*dxi*dyi*0.5*0.5
					- cellp_ym->mu*(cellp->v_xnew - cellp_ymm->v_xnew)*dyi*0.5-cellp_ym->mu*(cellp_xp_ym->v_ynew - cellp_xm_ym->v_ynew)*dxi*0.5;
//			}
//			else{
//				macFym = cellp_ym->momentum_ynew*cellp_ym->v_xnew
//					- cellp_ym->mu*(cellp->v_xnew - cellp_ymm->v_xnew)*dyi*0.5-cellp_ym->mu*(cellp_xp_ym->v_ynew - cellp_xm_ym->v_ynew)*dxi*0.5;
//			}
//			if((cellp->rhonew<RHOmax)&&(cellp->rhonew>RHOmin)){
				macFx  = cellp->momentum_xnew*cellp->v_xnew + cellp->pressnew + kappa1*(cellp_xp->rhonew-cellp_xm->rhonew)*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*dxi*0.5*0.5
						- 4.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5+2.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5;
				macFy  = cellp->momentum_ynew*cellp->v_xnew + kappa1*(cellp_xp->rhonew-cellp_xm->rhonew)*(cellp_yp->rhonew-cellp_ym->rhonew)*dxi*dyi*0.5*0.5 
						- cellp->mu*(cellp_yp->v_xnew - cellp_ym->v_xnew)*dyi*0.5-cellp->mu*(cellp_xp->v_ynew - cellp_xm->v_ynew)*dxi*0.5;
//			}
//			else{
//				macFx  = cellp->momentum_xnew*cellp->v_xnew + cellp->pressnew
//						- 4.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5+2.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5;
//				macFy  = cellp->momentum_ynew*cellp->v_xnew
//						- cellp->mu*(cellp_yp->v_xnew - cellp_ym->v_xnew)*dyi*0.5-cellp->mu*(cellp_xp->v_ynew - cellp_xm->v_ynew)*dxi*0.5;
//			}
			macQxm = 0.125*mack*(cellp->momentum_xnew-2.0*cellp_xm->momentum_xnew+cellp_xmm->momentum_xnew);
			macQx  = 0.125*mack*(cellp_xp->momentum_xnew-2.0*cellp->momentum_xnew+cellp_xm->momentum_xnew);
			macQym = 0.125*mack*(cellp->momentum_xnew-2.0*cellp_ym->momentum_xnew+cellp_ymm->momentum_xnew);
			macQy  = 0.125*mack*(cellp_yp->momentum_xnew-2.0*cellp->momentum_xnew+cellp_ym->momentum_xnew);
			
			cellp->momentum_x = 0.5*((cellp->momentum_x+cellp->momentum_xnew) - dt*(dxi*(macFx - macFxm) + dyi*(macFy - macFym) - dxi*(macQx - macQxm) - dyi*(macQy - macQym)));
			
//			if((cellp_xm->rhonew<RHOmax)&&(cellp_xm->rhonew>RHOmin)){
				macFxm = cellp_xm->momentum_xnew*cellp_xm->v_ynew + kappa1*(cellp_xm_yp->rhonew-cellp_xm_ym->rhonew)*(cellp->rhonew-cellp_xmm->rhonew)*dxi*dyi*0.5*0.5
						- cellp_xm->mu*(cellp->v_ynew - cellp_xmm->v_ynew)*dxi*0.5-cellp_xm->mu*(cellp_xm_yp->v_xnew - cellp_xm_ym->v_xnew)*dyi*0.5;
//			}
//			else{
//				macFxm = cellp_xm->momentum_xnew*cellp_xm->v_ynew
//						- cellp_xm->mu*(cellp->v_ynew - cellp_xmm->v_ynew)*dxi*0.5-cellp_xm->mu*(cellp_xm_yp->v_xnew - cellp_xm_ym->v_xnew)*dyi*0.5;
//			}
//			if((cellp_ym->rhonew<RHOmax)&&(cellp_ym->rhonew>RHOmin)){
				macFym = cellp_ym->momentum_ynew*cellp_ym->v_ynew + cellp_ym->pressnew + kappa1*(cellp->rhonew-cellp_ymm->rhonew)*(cellp->rhonew-cellp_ymm->rhonew)*dyi*dyi*0.5*0.5
					- 4.0/3.0*cellp_ym->mu*(cellp->v_ynew - cellp_ymm->v_ynew)*dyi*0.5+2.0/3.0*cellp_ym->mu*(cellp_xp_ym->v_xnew - cellp_xm_ym->v_xnew)*dxi*0.5;
//			}
//			else{
//				macFym = cellp_ym->momentum_ynew*cellp_ym->v_ynew + cellp_ym->pressnew
//					- 4.0/3.0*cellp_ym->mu*(cellp->v_ynew - cellp_ymm->v_ynew)*dyi*0.5+2.0/3.0*cellp_ym->mu*(cellp_xp_ym->v_xnew - cellp_xm_ym->v_xnew)*dxi*0.5;
//			}
//			if((cellp->rhonew<RHOmax)&&(cellp->rhonew>RHOmin)){
				macFx  = cellp->momentum_xnew*cellp->v_ynew + kappa1*(cellp_yp->rhonew-cellp_ym->rhonew)*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*dyi*0.5*0.5
						- cellp->mu*(cellp_xp->v_ynew - cellp_xm->v_ynew)*dxi*0.5-cellp->mu*(cellp_yp->v_xnew - cellp_ym->v_xnew)*dyi*0.5;
				macFy  = cellp->momentum_ynew*cellp->v_ynew + cellp->pressnew + kappa1*(cellp_yp->rhonew-cellp_ym->rhonew)*(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*dyi*0.5*0.5
						- 4.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5+2.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5;
//			}
//			else{
//				macFx  = cellp->momentum_xnew*cellp->v_ynew
//						- cellp->mu*(cellp_xp->v_ynew - cellp_xm->v_ynew)*dxi*0.5-cellp->mu*(cellp_yp->v_xnew - cellp_ym->v_xnew)*dyi*0.5;
//				macFy  = cellp->momentum_ynew*cellp->v_ynew + cellp->pressnew
//						- 4.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5+2.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5;
//			}
			macQxm = 0.125*mack*(cellp->momentum_ynew-2.0*cellp_xm->momentum_ynew+cellp_xmm->momentum_ynew);
			macQx  = 0.125*mack*(cellp_xp->momentum_ynew-2.0*cellp->momentum_ynew+cellp_xm->momentum_ynew);
			macQym = 0.125*mack*(cellp->momentum_ynew-2.0*cellp_ym->momentum_ynew+cellp_ymm->momentum_ynew);
			macQy  = 0.125*mack*(cellp_yp->momentum_ynew-2.0*cellp->momentum_ynew+cellp_ym->momentum_ynew);
			
			cellp->momentum_y = 0.5*((cellp->momentum_y+cellp->momentum_ynew) - dt*(dxi*(macFx - macFxm) + dyi*(macFy - macFym) - dxi*(macQx - macQxm) - dyi*(macQy - macQym)));
			
//			if((cellp_xm->rhonew<RHOmax)&&(cellp_xm->rhonew>RHOmin)){
				macFxm = cellp_xm->energynew*cellp_xm->v_xnew + cellp_xm->press*cellp_xm->v_xnew
							+ kappa1*(cellp->rhonew-cellp_xmm->rhonew)*(cellp->rhonew-cellp_xmm->rhonew)*dxi*dxi*0.5*0.5*cellp_xm->v_xnew 
							+ kappa1*(cellp->rhonew-cellp_xmm->rhonew)*(cellp_xm_yp->rhonew-cellp_xm_ym->rhonew)*dxi*dyi*0.5*0.5*cellp_xm->v_ynew
							- ((4.0/3.0*cellp_xm->mu*(cellp->v_xnew - cellp_xmm->v_xnew)*dxi*0.5-2.0/3.0*cellp_xm->mu*(cellp_xm_yp->v_ynew - cellp_xm_ym->v_ynew)*dyi*0.5)*cellp_xm->v_xnew
							+(cellp_xm->mu*(cellp_xm_yp->v_xnew - cellp_xm_ym->v_xnew)*dyi*0.5+cellp_xm->mu*(cellp->v_ynew - cellp_xmm->v_ynew)*dxi*0.5)*cellp_xm->v_ynew)
							- kt*(cellp->tempernew - cellp_xmm->tempernew)*dxi*0.5 
							+ kappa1*cellp_xm->rhonew*((cellp->v_xnew-cellp_xmm->v_xnew)*dxi*0.5+(cellp_xm_yp->v_ynew-cellp_xm_ym->v_ynew)*dyi*0.5)*(cellp->rhonew-cellp_xmm->rhonew)*dxi*0.5;
//			}
//			else{
//				macFxm = cellp_xm->energynew*cellp_xm->v_xnew + cellp_xm->press*cellp_xm->v_xnew
//							- (4.0/3.0*cellp_xm->mu*(cellp->v_xnew - cellp_xmm->v_xnew)*dxi*0.5-2.0/3.0*cellp_xm->mu*(cellp_xm_yp->v_ynew - cellp_xm_ym->v_ynew)*dyi*0.5)*cellp_xm->v_xnew
//							- kt*(cellp->tempernew - cellp_xmm->tempernew)*dxi*0.5;
//			}
//			if((cellp_ym->rhonew<RHOmax)&&(cellp_ym->rhonew>RHOmin)){
				macFym = cellp_ym->energynew*cellp_ym->v_ynew + cellp_ym->press*cellp_ym->v_ynew 
						+ kappa1*(cellp->rhonew-cellp_ymm->rhonew)*(cellp->rhonew-cellp_ymm->rhonew)*dyi*dyi*0.5*0.5*cellp_ym->v_ynew 
						+ kappa1*(cellp->rhonew-cellp_ymm->rhonew)*(cellp_xp_ym->rhonew-cellp_xm_ym->rhonew)*dyi*dxi*0.5*0.5*cellp_ym->v_xnew
						- ((4.0/3.0*cellp_ym->mu*(cellp->v_ynew - cellp_ymm->v_ynew)*dyi*0.5-2.0/3.0*cellp_ym->mu*(cellp_xp_ym->v_xnew - cellp_xm_ym->v_xnew)*dxi*0.5)*cellp_ym->v_ynew
						+(cellp_ym->mu*(cellp_xp_ym->v_ynew - cellp_xm_ym->v_ynew)*dxi*0.5+cellp_ym->mu*(cellp->v_xnew - cellp_ymm->v_xnew)*dyi*0.5)*cellp_ym->v_xnew)
						- kt*(cellp->tempernew - cellp_ymm->tempernew)*dyi*0.5 
						+ kappa1*cellp_ym->rhonew*((cellp->v_ynew-cellp_ymm->v_ynew)*dyi*0.5+(cellp_xp_ym->v_xnew-cellp_xm_ym->v_xnew)*dxi*0.5)*(cellp->rhonew-cellp_ymm->rhonew)*dyi*0.5;
//			}
//			else{
//				macFym = cellp_ym->energynew*cellp_ym->v_ynew + cellp_ym->press*cellp_ym->v_ynew 
//						- (4.0/3.0*cellp_ym->mu*(cellp->v_ynew - cellp_ymm->v_ynew)*dyi*0.5-2.0/3.0*cellp_ym->mu*(cellp_xp_ym->v_xnew - cellp_xm_ym->v_xnew)*dxi*0.5)*cellp_ym->v_ynew
//						- kt*(cellp->tempernew - cellp_ymm->tempernew)*dyi*0.5;
//			}
//			if((cellp->rhonew<RHOmax)&&(cellp->rhonew>RHOmin)){
				macFx  = cellp->energynew*cellp->v_xnew + cellp->press*cellp->v_xnew
							+ kappa1*(cellp_xp->rhonew-cellp_xm->rhonew)*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*dxi*0.5*0.5*cellp->v_xnew 
							+ kappa1*(cellp_xp->rhonew-cellp_xm->rhonew)*(cellp_yp->rhonew-cellp_ym->rhonew)*dxi*dyi*0.5*0.5*cellp->v_ynew
							- ((4.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5-2.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5)*cellp->v_xnew
							+(cellp->mu*(cellp_yp->v_xnew - cellp_ym->v_xnew)*dyi*0.5+cellp->mu*(cellp_xp->v_ynew - cellp_xm->v_ynew)*dxi*0.5)*cellp->v_ynew)
							- kt*(cellp_xp->tempernew - cellp_xm->tempernew)*dxi*0.5 
							+ kappa1*cellp->rhonew*((cellp_xp->v_xnew-cellp_xm->v_xnew)*dxi*0.5+(cellp_yp->v_ynew-cellp_ym->v_ynew)*dyi*0.5)*(cellp_xp->rhonew-cellp_xm->rhonew)*dxi*0.5;
				macFy  = cellp->energynew*cellp->v_ynew + cellp->press*cellp->v_ynew 
							+ kappa1*(cellp_yp->rhonew-cellp_ym->rhonew)*(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*dyi*0.5*0.5*cellp->v_ynew 
							+ kappa1*(cellp_yp->rhonew-cellp_ym->rhonew)*(cellp_xp->rhonew-cellp_xm->rhonew)*dyi*dxi*0.5*0.5*cellp->v_xnew
							- ((4.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5-2.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5)*cellp->v_ynew
							+(cellp->mu*(cellp_xp->v_ynew - cellp_xm->v_ynew)*dxi*0.5+cellp->mu*(cellp_yp->v_xnew - cellp_ym->v_xnew)*dyi*0.5)*cellp->v_xnew)
							- kt*(cellp_yp->tempernew - cellp_ym->tempernew)*dyi*0.5 
							+ kappa1*cellp->rhonew*((cellp_yp->v_ynew-cellp_ym->v_ynew)*dyi*0.5+(cellp_xp->v_xnew-cellp_xm->v_xnew)*dxi*0.5)*(cellp_yp->rhonew-cellp_ym->rhonew)*dyi*0.5;
//			}
//			else{
//				macFx  = cellp->energynew*cellp->v_xnew + cellp->press*cellp->v_xnew
//							- (4.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5-2.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5)*cellp->v_xnew
//							- kt*(cellp_xp->tempernew - cellp_xm->tempernew)*dxi*0.5;
//				macFy  = cellp->energynew*cellp->v_ynew + cellp->press*cellp->v_ynew 
//							- (4.0/3.0*cellp->mu*(cellp_yp->v_ynew - cellp_ym->v_ynew)*dyi*0.5-2.0/3.0*cellp->mu*(cellp_xp->v_xnew - cellp_xm->v_xnew)*dxi*0.5)*cellp->v_ynew
//							- kt*(cellp_yp->tempernew - cellp_ym->tempernew)*dyi*0.5;
//			}
			macQxm = 0.125*mack*(cellp->energynew-2.0*cellp_xm->energynew+cellp_xmm->energynew);
			macQx  = 0.125*mack*(cellp_xp->energynew-2.0*cellp->energynew+cellp_xm->energynew);
			macQym = 0.125*mack*(cellp->energynew-2.0*cellp_ym->energynew+cellp_ymm->energynew);
			macQy  = 0.125*mack*(cellp_yp->energynew-2.0*cellp->energynew+cellp_ym->energynew);
			
			cellp->energy = 0.5*((cellp->energy+cellp->energynew) - dt*(dxi*(macFx - macFxm) + dyi*(macFy - macFym) - dxi*(macQx - macQxm) - dyi*(macQy - macQym)));
		}
	}
	
	
	rhobnd();
	momentumbnd();
	energybnd();
	
	
	for(int i=3; i<(XCELL_NUM+3); i++){
		for(int j=3; j<(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
			ddrho[index0] = (-cellp_xpp->rho+16.0*cellp_xp->rho-30.0*cellp->rho+16.0*cellp_xm->rho-cellp_xmm->rho)/12.0*dxi*dxi+(-cellp_ypp->rho+16.0*cellp_yp->rho-30.0*cellp->rho+16.0*cellp_ym->rho-cellp_ymm->rho)/12.0*dyi*dyi;
		}
	}
		
	for(int i=3; i<=(XCELL_NUM+3); i++){
		for(int j=3; j<=(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
			cellp->v_x = cellp->momentum_x/cellp->rho;
			cellp->v_y = cellp->momentum_y/cellp->rho;
		}
	}
	
	vbnd();
	
	for(int i=3; i<=(XCELL_NUM+3); i++){
		for(int j=3; j<=(YCELL_NUM+3); j++){
			indexnsdefine(i,j);
			
				cellp->temper = (cellp->energy - 0.5*cellp->rho*(cellp->v_x*cellp->v_x+cellp->v_y*cellp->v_y) 
								+ cellp->rho*vdWa*cellp->rho-0.5*kappa1*((cellp_xp->rho-cellp_xm->rho)*dxi*0.5*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5
								+(cellp_yp->rho-cellp_ym->rho)*dyi*0.5*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5))/(cellp->rho*vdWc);
//			}
//			else{
//				cellp->temper = (cellp->energy - 0.5*cellp->rho*(cellp->v_x*cellp->v_x+cellp->v_y*cellp->v_y) 
//								+ cellp->rho*vdWa*cellp->rho)/(cellp->rhonew*vdWc);
//			}
//			if((cellp->rho<RHOmax)&&(cellp->rho>RHOmin)){
				cellp->press = cellp->rho*cellp->temper/(1.0-vdWb*cellp->rho) - vdWa*cellp->rho*cellp->rho 
								- kappa1*cellp->rho*ddrho[index0]
//	- kappa1*cellp->rho*((cellp_xp->rho-2.0*cellp->rho+cellp_xm->rho)*(cellp_xp->rho-2.0*cellp->rho+cellp_xm->rho)+(cellp_yp->rho-2.0*cellp->rho+cellp_ym->rho)*(cellp_yp->rho-2.0*cellp->rho+cellp_ym->rho))
								- 0.5*kappa1*((cellp_xp->rho-cellp_xm->rho)*dxi*0.5*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5
								+(cellp_yp->rho-cellp_ym->rho)*dyi*0.5*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5);
//			}
//			else{
//				cellp->press = cellp->rho*cellp->temper/(1.0-vdWb*cellp->rho) - vdWa*cellp->rho*cellp->rho;
//			}
		}
	}
	
	temperbnd();
	
	
	for(int i=3; i<=(XCELL_NUM+3); i++){
			indexnsdefine(i,3);
			
				cellp->press = cellp->rho*cellp->temper/(1.0-vdWb*cellp->rho) - vdWa*cellp->rho*cellp->rho 
								- kappa1*cellp->rho*ddrho[index0]
								- 0.5*kappa1*((cellp_xp->rho-cellp_xm->rho)*dxi*0.5*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5
								+(cellp_yp->rho-cellp_ym->rho)*dyi*0.5*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5);
//			}
//			else{
//				cellp->press = cellp->rho*cellp->temper/(1.0-vdWb*cellp->rho) - vdWa*cellp->rho*cellp->rho;
//			}
			cellp->energy = 0.5*cellp->rho*(cellp->v_x*cellp->v_x+cellp->v_y*cellp->v_y) + cellp->rho*(vdWc*cellp->temper - vdWa*cellp->rho) 
							+ 0.5*kappa1*((cellp_xp->rho-cellp_xm->rho)*dxi*0.5*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5
							+(cellp_yp->rho-cellp_ym->rho)*dyi*0.5*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5);
	}
	
	pressbnd();
}


void momentumnewbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->momentum_xnew   = cellp_ycp0->momentum_xnew;
		cellp_y1->momentum_xnew   = cellp_ycp1->momentum_xnew;
		cellp_y2->momentum_xnew   = cellp_ycp2->momentum_xnew;
		cellp_ycp3->momentum_xnew = cellp_y3->momentum_xnew;
		cellp_ycp4->momentum_xnew = cellp_y4->momentum_xnew;
		cellp_ycp5->momentum_xnew = cellp_y5->momentum_xnew;
		cellp_y0->momentum_ynew   = cellp_ycp0->momentum_ynew;
		cellp_y1->momentum_ynew   = cellp_ycp1->momentum_ynew;
		cellp_y2->momentum_ynew   = cellp_ycp2->momentum_ynew;
		cellp_ycp3->momentum_ynew = cellp_y3->momentum_ynew;
		cellp_ycp4->momentum_ynew = cellp_y4->momentum_ynew;
		cellp_ycp5->momentum_ynew = cellp_y5->momentum_ynew;
	}
#endif
#ifndef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
	
		cellp_x0->momentum_xnew   = cellp_xcp0->momentum_xnew;
		cellp_x1->momentum_xnew   = cellp_xcp1->momentum_xnew;
		cellp_x2->momentum_xnew   = cellp_xcp2->momentum_xnew;
		cellp_xcp3->momentum_xnew = cellp_x3->momentum_xnew;
		cellp_xcp4->momentum_xnew = cellp_x4->momentum_xnew;
		cellp_xcp5->momentum_xnew = cellp_x5->momentum_xnew;
		cellp_x0->momentum_ynew   = cellp_xcp0->momentum_ynew;
		cellp_x1->momentum_ynew   = cellp_xcp1->momentum_ynew;
		cellp_x2->momentum_ynew   = cellp_xcp2->momentum_ynew;
		cellp_xcp3->momentum_ynew = cellp_x3->momentum_ynew;
		cellp_xcp4->momentum_ynew = cellp_x4->momentum_ynew;
		cellp_xcp5->momentum_ynew = cellp_x5->momentum_ynew;
	}
#endif
/*
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->momentum_xnew   = -cellp_y5->momentum_xnew;
		cellp_y1->momentum_xnew   = -cellp_y4->momentum_xnew;
		cellp_y2->momentum_xnew   = -cellp_y3->momentum_xnew;
		cellp_ycp3->momentum_xnew = -cellp_ycp2->momentum_xnew;
		cellp_ycp4->momentum_xnew = -cellp_ycp1->momentum_xnew;
		cellp_ycp5->momentum_xnew = -cellp_ycp0->momentum_xnew;
		cellp_y0->momentum_ynew   = -cellp_y5->momentum_ynew;
		cellp_y1->momentum_ynew   = -cellp_y4->momentum_ynew;
		cellp_y2->momentum_ynew   = -cellp_y3->momentum_ynew;
		cellp_ycp3->momentum_ynew = -cellp_ycp2->momentum_ynew;
		cellp_ycp4->momentum_ynew = -cellp_ycp1->momentum_ynew;
		cellp_ycp5->momentum_ynew = -cellp_ycp0->momentum_ynew;
	}
#endif*/
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->momentum_xnew   = -cellp_y6->momentum_xnew;
		cellp_y1->momentum_xnew   = -cellp_y5->momentum_xnew;
		cellp_y2->momentum_xnew   = -cellp_y4->momentum_xnew;
		cellp_y3->momentum_xnew   = 0.0;
		cellp_ycp3->momentum_xnew = cellp_ycp2->momentum_xnew;
		cellp_ycp4->momentum_xnew = cellp_ycp1->momentum_xnew;
		cellp_ycp5->momentum_xnew = cellp_ycp0->momentum_xnew;
		cellp_y0->momentum_ynew   = -cellp_y6->momentum_ynew;
		cellp_y1->momentum_ynew   = -cellp_y5->momentum_ynew;
		cellp_y2->momentum_ynew   = -cellp_y4->momentum_ynew;
		cellp_y3->momentum_ynew   = 0.0;
		cellp_ycp3->momentum_ynew = cellp_ycp2->momentum_ynew;
		cellp_ycp4->momentum_ynew = cellp_ycp1->momentum_ynew;
		cellp_ycp5->momentum_ynew = cellp_ycp0->momentum_ynew;
	}
#endif
#ifdef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->momentum_xnew   = -cellp_x5->momentum_xnew;
		cellp_x1->momentum_xnew   = -cellp_x4->momentum_xnew;
		cellp_x2->momentum_xnew   = -cellp_x3->momentum_xnew;
		cellp_xcp3->momentum_xnew = -cellp_xcp2->momentum_xnew;
		cellp_xcp4->momentum_xnew = -cellp_xcp1->momentum_xnew;
		cellp_xcp5->momentum_xnew = -cellp_xcp0->momentum_xnew;
		cellp_x0->momentum_ynew   = -cellp_x5->momentum_ynew;
		cellp_x1->momentum_ynew   = -cellp_x4->momentum_ynew;
		cellp_x2->momentum_ynew   = -cellp_x3->momentum_ynew;
		cellp_xcp3->momentum_ynew = -cellp_xcp2->momentum_ynew;
		cellp_xcp4->momentum_ynew = -cellp_xcp1->momentum_ynew;
		cellp_xcp5->momentum_ynew = -cellp_xcp0->momentum_ynew;
	}
#endif
}


void vnewbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->v_xnew   = cellp_ycp0->v_xnew;
		cellp_y1->v_xnew   = cellp_ycp1->v_xnew;
		cellp_y2->v_xnew   = cellp_ycp2->v_xnew;
		cellp_ycp3->v_xnew = cellp_y3->v_xnew;
		cellp_ycp4->v_xnew = cellp_y4->v_xnew;
		cellp_ycp5->v_xnew = cellp_y5->v_xnew;
		cellp_y0->v_ynew   = cellp_ycp0->v_ynew;
		cellp_y1->v_ynew   = cellp_ycp1->v_ynew;
		cellp_y2->v_ynew   = cellp_ycp2->v_ynew;
		cellp_ycp3->v_ynew = cellp_y3->v_ynew;
		cellp_ycp4->v_ynew = cellp_y4->v_ynew;
		cellp_ycp5->v_ynew = cellp_y5->v_ynew;
	}
#endif
#ifndef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
	
		cellp_x0->v_xnew   = cellp_xcp0->v_xnew;
		cellp_x1->v_xnew   = cellp_xcp1->v_xnew;
		cellp_x2->v_xnew   = cellp_xcp2->v_xnew;
		cellp_xcp3->v_xnew = cellp_x3->v_xnew;
		cellp_xcp4->v_xnew = cellp_x4->v_xnew;
		cellp_xcp5->v_xnew = cellp_x5->v_xnew;
		cellp_x0->v_ynew   = cellp_xcp0->v_ynew;
		cellp_x1->v_ynew   = cellp_xcp1->v_ynew;
		cellp_x2->v_ynew   = cellp_xcp2->v_ynew;
		cellp_xcp3->v_ynew = cellp_x3->v_ynew;
		cellp_xcp4->v_ynew = cellp_x4->v_ynew;
		cellp_xcp5->v_ynew = cellp_x5->v_ynew;
	}
#endif
/*
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->v_xnew   = -cellp_y5->v_xnew;
		cellp_y1->v_xnew   = -cellp_y4->v_xnew;
		cellp_y2->v_xnew   = -cellp_y3->v_xnew;
		cellp_ycp3->v_xnew = -cellp_ycp2->v_xnew;
		cellp_ycp4->v_xnew = -cellp_ycp1->v_xnew;
		cellp_ycp5->v_xnew = -cellp_ycp0->v_xnew;
		cellp_y0->v_ynew   = -cellp_y5->v_ynew;
		cellp_y1->v_ynew   = -cellp_y4->v_ynew;
		cellp_y2->v_ynew   = -cellp_y3->v_ynew;
		cellp_ycp3->v_ynew = -cellp_ycp2->v_ynew;
		cellp_ycp4->v_ynew = -cellp_ycp1->v_ynew;
		cellp_ycp5->v_ynew = -cellp_ycp0->v_ynew;
	}
#endif*/
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->v_xnew   = -cellp_y6->v_xnew;
		cellp_y1->v_xnew   = -cellp_y5->v_xnew;
		cellp_y2->v_xnew   = -cellp_y4->v_xnew;
		cellp_y3->v_xnew   = 0.0;
		cellp_ycp3->v_xnew = cellp_ycp2->v_xnew;
		cellp_ycp4->v_xnew = cellp_ycp1->v_xnew;
		cellp_ycp5->v_xnew = cellp_ycp0->v_xnew;
		cellp_y0->v_ynew   = -cellp_y6->v_ynew;
		cellp_y1->v_ynew   = -cellp_y5->v_ynew;
		cellp_y2->v_ynew   = -cellp_y4->v_ynew;
		cellp_y3->v_ynew   = 0.0;
		cellp_ycp3->v_ynew = cellp_ycp2->v_ynew;
		cellp_ycp4->v_ynew = cellp_ycp1->v_ynew;
		cellp_ycp5->v_ynew = cellp_ycp0->v_ynew;
	}
#endif
#ifdef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->v_xnew   = -cellp_x5->v_xnew;
		cellp_x1->v_xnew   = -cellp_x4->v_xnew;
		cellp_x2->v_xnew   = -cellp_x3->v_xnew;
		cellp_xcp3->v_xnew = -cellp_xcp2->v_xnew;
		cellp_xcp4->v_xnew = -cellp_xcp1->v_xnew;
		cellp_xcp5->v_xnew = -cellp_xcp0->v_xnew;
		cellp_x0->v_ynew   = -cellp_x5->v_ynew;
		cellp_x1->v_ynew   = -cellp_x4->v_ynew;
		cellp_x2->v_ynew   = -cellp_x3->v_ynew;
		cellp_xcp3->v_ynew = -cellp_xcp2->v_ynew;
		cellp_xcp4->v_ynew = -cellp_xcp1->v_ynew;
		cellp_xcp5->v_ynew = -cellp_xcp0->v_ynew;
	}
#endif
}
			//-------------- âºÇÃó¨ë¨ã´äEèåèèIÇÌÇË ---------------//


//-------------- à≥óÕã´äEèåè -----------------//
void pressnewbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->pressnew   = cellp_ycp0->pressnew;
		cellp_y1->pressnew   = cellp_ycp1->pressnew;
		cellp_y2->pressnew   = cellp_ycp2->pressnew;
		cellp_ycp3->pressnew = cellp_y3->pressnew;
		cellp_ycp4->pressnew = cellp_y4->pressnew;
		cellp_ycp5->pressnew = cellp_y5->pressnew;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->pressnew   = cellp_xcp0->pressnew;
		cellp_x1->pressnew   = cellp_xcp1->pressnew;
		cellp_x2->pressnew   = cellp_xcp2->pressnew;
		cellp_xcp3->pressnew = cellp_x3->pressnew;
		cellp_xcp4->pressnew = cellp_x4->pressnew;
		cellp_xcp5->pressnew = cellp_x5->pressnew;
	}
#endif
/*
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->pressnew   = cellp_y5->pressnew;
		cellp_y1->pressnew   = cellp_y4->pressnew;
		cellp_y2->pressnew   = cellp_y3->pressnew;
		cellp_ycp3->pressnew = cellp_ycp2->pressnew;
		cellp_ycp4->pressnew = cellp_ycp1->pressnew;
		cellp_ycp5->pressnew = cellp_ycp0->pressnew;
	}
#endif
*/
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->pressnew   = cellp_y3->pressnew;
		cellp_y1->pressnew   = cellp_y3->pressnew;
		cellp_y2->pressnew   = cellp_y3->pressnew;
		cellp_ycp3->pressnew = pressstart;
		cellp_ycp4->pressnew = pressstart;
		cellp_ycp5->pressnew = pressstart;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->pressnew   = cellp_x5->pressnew;
		cellp_x1->pressnew   = cellp_x4->pressnew;
		cellp_x2->pressnew   = cellp_x3->pressnew;
		cellp_xcp3->pressnew = cellp_xcp2->pressnew;
		cellp_xcp4->pressnew = cellp_xcp1->pressnew;
		cellp_xcp5->pressnew = cellp_xcp0->pressnew;
	}
#endif
}


//-------------- â∑ìxã´äEèåè -----------------//
void tempernewbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->tempernew   = cellp_ycp0->tempernew;
		cellp_y1->tempernew   = cellp_ycp1->tempernew;
		cellp_y2->tempernew   = cellp_ycp2->tempernew;
		cellp_ycp3->tempernew = cellp_y3->tempernew;
		cellp_ycp4->tempernew = cellp_y4->tempernew;
		cellp_ycp5->tempernew = cellp_y5->tempernew;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->tempernew   = cellp_xcp0->tempernew;
		cellp_x1->tempernew   = cellp_xcp1->tempernew;
		cellp_x2->tempernew   = cellp_xcp2->tempernew;
		cellp_xcp3->tempernew = cellp_x3->tempernew;
		cellp_xcp4->tempernew = cellp_x4->tempernew;
		cellp_xcp5->tempernew = cellp_x5->tempernew;
	}
#endif
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
//		cellp_y0->tempernew   = cellp_y5->tempernew;
//		cellp_y1->tempernew   = cellp_y4->tempernew;
//		cellp_y2->tempernew   = cellp_y3->tempernew;
//		cellp_ycp3->tempernew = cellp_ycp2->tempernew;
//		cellp_ycp4->tempernew = cellp_ycp1->tempernew;
//		cellp_ycp5->tempernew = cellp_ycp0->tempernew;
		cellp_y0->tempernew   = Ts;
		cellp_y1->tempernew   = Ts;
		cellp_y2->tempernew   = Ts;
		cellp_y3->tempernew   = Ts;
		cellp_ycp3->tempernew = Ts;
		cellp_ycp4->tempernew = Ts;
		cellp_ycp5->tempernew = Ts;
	}
	for(int i=23; i<33; i++){
//	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->tempernew   = Ts+deltaT;
		cellp_y1->tempernew   = Ts+deltaT;
		cellp_y2->tempernew   = Ts+deltaT;
		cellp_y3->tempernew   = Ts+deltaT;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
//		cellp_x0->tempernew   = cellp_x5->tempernew;
//		cellp_x1->tempernew   = cellp_x4->tempernew;
//		cellp_x2->tempernew   = cellp_x3->tempernew;
//		cellp_xcp3->tempernew = cellp_xcp2->tempernew;
//		cellp_xcp4->tempernew = cellp_xcp1->tempernew;
//		cellp_xcp5->tempernew = cellp_xcp0->tempernew;
		
		cellp_x0->tempernew   = Ts;
		cellp_x1->tempernew   = Ts;
		cellp_x2->tempernew   = Ts;
		cellp_xcp3->tempernew = Ts;
		cellp_xcp4->tempernew = Ts;
		cellp_xcp5->tempernew = Ts;
	}
	
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->tempernew   = Ts+deltaT;
		cellp_x1->tempernew   = Ts+deltaT;
		cellp_x2->tempernew   = Ts+deltaT;
	}
	
#endif
}


void rhonewbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->rhonew   = cellp_ycp0->rhonew;
		cellp_y1->rhonew   = cellp_ycp1->rhonew;
		cellp_y2->rhonew   = cellp_ycp2->rhonew;
		cellp_ycp3->rhonew = cellp_y3->rhonew;
		cellp_ycp4->rhonew = cellp_y4->rhonew;
		cellp_ycp5->rhonew = cellp_y5->rhonew;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->rhonew   = cellp_xcp0->rhonew;
		cellp_x1->rhonew   = cellp_xcp1->rhonew;
		cellp_x2->rhonew   = cellp_xcp2->rhonew;
		cellp_xcp3->rhonew = cellp_x3->rhonew;
		cellp_xcp4->rhonew = cellp_x4->rhonew;
		cellp_xcp5->rhonew = cellp_x5->rhonew;
	}
#endif
	
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
			
		if((cellp_y3->rhonew>RHOmin)&&(cellp_y3->rhonew<RHOmax)){
			cellp_y0->rhonew   = 2.0*cellp_y3->rhonew-cellp_y6->rhonew;
			cellp_y1->rhonew   = 2.0*cellp_y3->rhonew-cellp_y5->rhonew;
			cellp_y2->rhonew   = 2.0*cellp_y3->rhonew-cellp_y4->rhonew;
//			cellp_y3->rhonew   = cellp_y4->rhonew;
		}
		else {
			cellp_y0->rhonew   = 2.0*cellp_y3->rhonew-cellp_y6->rhonew;
			cellp_y1->rhonew   = 2.0*cellp_y3->rhonew-cellp_y5->rhonew;
			cellp_y2->rhonew   = 2.0*cellp_y3->rhonew-cellp_y4->rhonew;
//			cellp_y3->rhonew   = cellp_y4->rhonew;
		}
		if((cellp_ycp2->rhonew>RHOmin)&&(cellp_ycp2->rhonew<RHOmax)){
			cellp_ycp3->rhonew = cellp_ycp2->rhonew-wetrho*dy*1.0;
			cellp_ycp4->rhonew = cellp_ycp1->rhonew-wetrho*dy*3.0;
			cellp_ycp5->rhonew = cellp_ycp0->rhonew-wetrho*dy*5.0;
		}
		else {
			cellp_ycp3->rhonew = cellp_ycp2->rhonew;
			cellp_ycp4->rhonew = cellp_ycp1->rhonew;
			cellp_ycp5->rhonew = cellp_ycp0->rhonew;
		}
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		if((cellp_x3->rhonew>RHOmin)&&(cellp_x3->rhonew<RHOmax)){
			cellp_x0->rhonew   = cellp_x5->rhonew-wetrho*dx*5.0;
			cellp_x1->rhonew   = cellp_x4->rhonew-wetrho*dx*3.0;
			cellp_x2->rhonew   = cellp_x3->rhonew-wetrho*dx*1.0;
		}
		else{
			cellp_x0->rhonew   = cellp_x5->rhonew;
			cellp_x1->rhonew   = cellp_x4->rhonew;
			cellp_x2->rhonew   = cellp_x3->rhonew;
			}
		if((cellp_xcp2->rhonew>RHOmin)&&(cellp_xcp2->rhonew<RHOmax)){
			cellp_xcp3->rhonew = cellp_xcp2->rhonew-wetrho*dx*1.0;
			cellp_xcp4->rhonew = cellp_xcp1->rhonew-wetrho*dx*3.0;
			cellp_xcp5->rhonew = cellp_xcp0->rhonew-wetrho*dx*5.0;
		}
		else{
			cellp_xcp3->rhonew = cellp_xcp2->rhonew;
			cellp_xcp4->rhonew = cellp_xcp1->rhonew;
			cellp_xcp5->rhonew = cellp_xcp0->rhonew;
		}
	}
#endif
}


void energynewbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->energynew   = cellp_ycp0->energynew;
		cellp_y1->energynew   = cellp_ycp1->energynew;
		cellp_y2->energynew   = cellp_ycp2->energynew;
		cellp_ycp3->energynew = cellp_y3->energynew;
		cellp_ycp4->energynew = cellp_y4->energynew;
		cellp_ycp5->energynew = cellp_y5->energynew;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->energynew   = cellp_xcp0->energynew;
		cellp_x1->energynew   = cellp_xcp1->energynew;
		cellp_x2->energynew   = cellp_xcp2->energynew;
		cellp_xcp3->energynew = cellp_x3->energynew;
		cellp_xcp4->energynew = cellp_x4->energynew;
		cellp_xcp5->energynew = cellp_x5->energynew;
	}
#endif
	
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
			
		cellp_y0->energynew   = cellp_y4->energynew;
		cellp_y1->energynew   = cellp_y4->energynew;
		cellp_y2->energynew   = 2.0*cellp_y3->energynew-cellp_y4->energynew;
//		cellp_y3->energynew   = cellp_y4->energynew;
		cellp_ycp3->energynew = cellp_ycp2->energynew;
		cellp_ycp4->energynew = cellp_ycp1->energynew;
		cellp_ycp5->energynew = cellp_ycp0->energynew;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->energynew   = cellp_x5->energynew;
		cellp_x1->energynew   = cellp_x4->energynew;
		cellp_x2->energynew   = cellp_x3->energynew;
		cellp_xcp3->energynew = cellp_xcp2->energynew;
		cellp_xcp4->energynew = cellp_xcp1->energynew;
		cellp_xcp5->energynew = cellp_xcp0->energynew;
	}
#endif
}


void momentumbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->momentum_x   = cellp_ycp0->momentum_x;
		cellp_y1->momentum_x   = cellp_ycp1->momentum_x;
		cellp_y2->momentum_x   = cellp_ycp2->momentum_x;
		cellp_ycp3->momentum_x = cellp_y3->momentum_x;
		cellp_ycp4->momentum_x = cellp_y4->momentum_x;
		cellp_ycp5->momentum_x = cellp_y5->momentum_x;
		cellp_y0->momentum_y   = cellp_ycp0->momentum_y;
		cellp_y1->momentum_y   = cellp_ycp1->momentum_y;
		cellp_y2->momentum_y   = cellp_ycp2->momentum_y;
		cellp_ycp3->momentum_y = cellp_y3->momentum_y;
		cellp_ycp4->momentum_y = cellp_y4->momentum_y;
		cellp_ycp5->momentum_y = cellp_y5->momentum_y;
	}
#endif
#ifndef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
	
		cellp_x0->momentum_x   = cellp_xcp0->momentum_x;
		cellp_x1->momentum_x   = cellp_xcp1->momentum_x;
		cellp_x2->momentum_x   = cellp_xcp2->momentum_x;
		cellp_xcp3->momentum_x = cellp_x3->momentum_x;
		cellp_xcp4->momentum_x = cellp_x4->momentum_x;
		cellp_xcp5->momentum_x = cellp_x5->momentum_x;
		cellp_x0->momentum_y   = cellp_xcp0->momentum_y;
		cellp_x1->momentum_y   = cellp_xcp1->momentum_y;
		cellp_x2->momentum_y   = cellp_xcp2->momentum_y;
		cellp_xcp3->momentum_y = cellp_x3->momentum_y;
		cellp_xcp4->momentum_y = cellp_x4->momentum_y;
		cellp_xcp5->momentum_y = cellp_x5->momentum_y;
	}
#endif
/*
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->momentum_x   = -cellp_y5->momentum_x;
		cellp_y1->momentum_x   = -cellp_y4->momentum_x;
		cellp_y2->momentum_x   = -cellp_y3->momentum_x;
		cellp_ycp3->momentum_x = -cellp_ycp2->momentum_x;
		cellp_ycp4->momentum_x = -cellp_ycp1->momentum_x;
		cellp_ycp5->momentum_x = -cellp_ycp0->momentum_x;
		cellp_y0->momentum_y   = -cellp_y5->momentum_y;
		cellp_y1->momentum_y   = -cellp_y4->momentum_y;
		cellp_y2->momentum_y   = -cellp_y3->momentum_y;
		cellp_ycp3->momentum_y = -cellp_ycp2->momentum_y;
		cellp_ycp4->momentum_y = -cellp_ycp1->momentum_y;
		cellp_ycp5->momentum_y = -cellp_ycp0->momentum_y;
	}
#endif*/
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->momentum_x   = -cellp_y6->momentum_x;
		cellp_y1->momentum_x   = -cellp_y5->momentum_x;
		cellp_y2->momentum_x   = -cellp_y4->momentum_x;
		cellp_y3->momentum_x   = 0.0;
		cellp_ycp3->momentum_x = cellp_ycp2->momentum_x;
		cellp_ycp4->momentum_x = cellp_ycp1->momentum_x;
		cellp_ycp5->momentum_x = cellp_ycp0->momentum_x;
		cellp_y0->momentum_y   = -cellp_y6->momentum_y;
		cellp_y1->momentum_y   = -cellp_y5->momentum_y;
		cellp_y2->momentum_y   = -cellp_y4->momentum_y;
		cellp_y3->momentum_y   = 0.0;
		cellp_ycp3->momentum_y = cellp_ycp2->momentum_y;
		cellp_ycp4->momentum_y = cellp_ycp1->momentum_y;
		cellp_ycp5->momentum_y = cellp_ycp0->momentum_y;
	}
#endif
#ifdef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->momentum_x   = -cellp_x5->momentum_x;
		cellp_x1->momentum_x   = -cellp_x4->momentum_x;
		cellp_x2->momentum_x   = -cellp_x3->momentum_x;
		cellp_xcp3->momentum_x = -cellp_xcp2->momentum_x;
		cellp_xcp4->momentum_x = -cellp_xcp1->momentum_x;
		cellp_xcp5->momentum_x = -cellp_xcp0->momentum_x;
		cellp_x0->momentum_y   = -cellp_x5->momentum_y;
		cellp_x1->momentum_y   = -cellp_x4->momentum_y;
		cellp_x2->momentum_y   = -cellp_x3->momentum_y;
		cellp_xcp3->momentum_y = -cellp_xcp2->momentum_y;
		cellp_xcp4->momentum_y = -cellp_xcp1->momentum_y;
		cellp_xcp5->momentum_y = -cellp_xcp0->momentum_y;
	}
#endif
}



void vbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->v_x   = cellp_ycp0->v_x;
		cellp_y1->v_x   = cellp_ycp1->v_x;
		cellp_y2->v_x   = cellp_ycp2->v_x;
		cellp_ycp3->v_x = cellp_y3->v_x;
		cellp_ycp4->v_x = cellp_y4->v_x;
		cellp_ycp5->v_x = cellp_y5->v_x;
		cellp_y0->v_y   = cellp_ycp0->v_y;
		cellp_y1->v_y   = cellp_ycp1->v_y;
		cellp_y2->v_y   = cellp_ycp2->v_y;
		cellp_ycp3->v_y = cellp_y3->v_y;
		cellp_ycp4->v_y = cellp_y4->v_y;
		cellp_ycp5->v_y = cellp_y5->v_y;
	}
#endif
#ifndef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
	
		cellp_x0->v_x   = cellp_xcp0->v_x;
		cellp_x1->v_x   = cellp_xcp1->v_x;
		cellp_x2->v_x   = cellp_xcp2->v_x;
		cellp_xcp3->v_x = cellp_x3->v_x;
		cellp_xcp4->v_x = cellp_x4->v_x;
		cellp_xcp5->v_x = cellp_x5->v_x;
		cellp_x0->v_y   = cellp_xcp0->v_y;
		cellp_x1->v_y   = cellp_xcp1->v_y;
		cellp_x2->v_y   = cellp_xcp2->v_y;
		cellp_xcp3->v_y = cellp_x3->v_y;
		cellp_xcp4->v_y = cellp_x4->v_y;
		cellp_xcp5->v_y = cellp_x5->v_y;
	}
#endif
/*
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->v_x   = -cellp_y5->v_x;
		cellp_y1->v_x   = -cellp_y4->v_x;
		cellp_y2->v_x   = -cellp_y3->v_x;
		cellp_ycp3->v_x = -cellp_ycp2->v_x;
		cellp_ycp4->v_x = -cellp_ycp1->v_x;
		cellp_ycp5->v_x = -cellp_ycp0->v_x;
		cellp_y0->v_y   = -cellp_y5->v_y;
		cellp_y1->v_y   = -cellp_y4->v_y;
		cellp_y2->v_y   = -cellp_y3->v_y;
		cellp_ycp3->v_y = -cellp_ycp2->v_y;
		cellp_ycp4->v_y = -cellp_ycp1->v_y;
		cellp_ycp5->v_y = -cellp_ycp0->v_y;
	}
#endif*/
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){										//ï«ñ 
		indexbndy(i);
			
		cellp_y0->v_x   = -cellp_y6->v_x;
		cellp_y1->v_x   = -cellp_y5->v_x;
		cellp_y2->v_x   = -cellp_y4->v_x;
		cellp_y3->v_x   = 0.0;
		cellp_ycp3->v_x = cellp_ycp2->v_x;
		cellp_ycp4->v_x = cellp_ycp1->v_x;
		cellp_ycp5->v_x = cellp_ycp0->v_x;
		cellp_y0->v_y   = -cellp_y6->v_y;
		cellp_y1->v_y   = -cellp_y5->v_y;
		cellp_y2->v_y   = -cellp_y4->v_y;
		cellp_y3->v_y   = 0.0;
		cellp_ycp3->v_y = cellp_ycp2->v_y;
		cellp_ycp4->v_y = cellp_ycp1->v_y;
		cellp_ycp5->v_y = cellp_ycp0->v_y;
	}
#endif
#ifdef WALLX
	for (int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->v_x   = -cellp_x5->v_x;
		cellp_x1->v_x   = -cellp_x4->v_x;
		cellp_x2->v_x   = -cellp_x3->v_x;
		cellp_xcp3->v_x = -cellp_xcp2->v_x;
		cellp_xcp4->v_x = -cellp_xcp1->v_x;
		cellp_xcp5->v_x = -cellp_xcp0->v_x;
		cellp_x0->v_y   = -cellp_x5->v_y;
		cellp_x1->v_y   = -cellp_x4->v_y;
		cellp_x2->v_y   = -cellp_x3->v_y;
		cellp_xcp3->v_y = -cellp_xcp2->v_y;
		cellp_xcp4->v_y = -cellp_xcp1->v_y;
		cellp_xcp5->v_y = -cellp_xcp0->v_y;
	}
#endif
}



//-------------- à≥óÕã´äEèåè -----------------//
void pressbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->press   = cellp_ycp0->press;
		cellp_y1->press   = cellp_ycp1->press;
		cellp_y2->press   = cellp_ycp2->press;
		cellp_ycp3->press = cellp_y3->press;
		cellp_ycp4->press = cellp_y4->press;
		cellp_ycp5->press = cellp_y5->press;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->press   = cellp_xcp0->press;
		cellp_x1->press   = cellp_xcp1->press;
		cellp_x2->press   = cellp_xcp2->press;
		cellp_xcp3->press = cellp_x3->press;
		cellp_xcp4->press = cellp_x4->press;
		cellp_xcp5->press = cellp_x5->press;
	}
#endif
/*
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->press   = cellp_y5->press;
		cellp_y1->press   = cellp_y4->press;
		cellp_y2->press   = cellp_y3->press;
		cellp_ycp3->press = cellp_ycp2->press;
		cellp_ycp4->press = cellp_ycp1->press;
		cellp_ycp5->press = cellp_ycp0->press;
	}
#endif*/
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->press   = cellp_y3->press;
		cellp_y1->press   = cellp_y3->press;
		cellp_y2->press   = cellp_y3->press;
		cellp_ycp3->press = pressstart;
		cellp_ycp4->press = pressstart;
		cellp_ycp5->press = pressstart;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->press   = cellp_x5->press;
		cellp_x1->press   = cellp_x4->press;
		cellp_x2->press   = cellp_x3->press;
		cellp_xcp3->press = cellp_xcp2->press;
		cellp_xcp4->press = cellp_xcp1->press;
		cellp_xcp5->press = cellp_xcp0->press;
	}
#endif
}


//-------------- â∑ìxã´äEèåè -----------------//
void temperbnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){											//é¸ä˙ã´äE
		indexbndy(i);
			
		cellp_y0->temper   = cellp_ycp0->temper;
		cellp_y1->temper   = cellp_ycp1->temper;
		cellp_y2->temper   = cellp_ycp2->temper;
		cellp_ycp3->temper = cellp_y3->temper;
		cellp_ycp4->temper = cellp_y4->temper;
		cellp_ycp5->temper = cellp_y5->temper;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->temper   = cellp_xcp0->temper;
		cellp_x1->temper   = cellp_xcp1->temper;
		cellp_x2->temper   = cellp_xcp2->temper;
		cellp_xcp3->temper = cellp_x3->temper;
		cellp_xcp4->temper = cellp_x4->temper;
		cellp_xcp5->temper = cellp_x5->temper;
	}
#endif

#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
//		cellp_y0->temper   = cellp_y5->temper;
//		cellp_y1->temper   = cellp_y4->temper;
//		cellp_y2->temper   = cellp_y3->temper;
//		cellp_ycp3->temper = cellp_ycp2->temper;
//		cellp_ycp4->temper = cellp_ycp1->temper;
//		cellp_ycp5->temper = cellp_ycp0->temper;
		cellp_y0->temper   = Ts;
		cellp_y1->temper   = Ts;
		cellp_y2->temper   = Ts;
		cellp_y3->temper   = Ts;
		cellp_ycp3->temper = Ts;
		cellp_ycp4->temper = Ts;
		cellp_ycp5->temper = Ts;
	}
	for(int i=23; i<33; i++){
//	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->temper   = Ts+deltaT;
		cellp_y1->temper   = Ts+deltaT;
		cellp_y2->temper   = Ts+deltaT;
		cellp_y3->temper   = Ts+deltaT;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
//		cellp_x0->temper   = cellp_x5->temper;
//		cellp_x1->temper   = cellp_x4->temper;
//		cellp_x2->temper   = cellp_x3->temper;
//		cellp_xcp3->temper = cellp_xcp2->temper;
//		cellp_xcp4->temper = cellp_xcp1->temper;
//		cellp_xcp5->temper = cellp_xcp0->temper;
		
		cellp_x0->temper   = Ts;
		cellp_x1->temper   = Ts;
		cellp_x2->temper   = Ts;
		cellp_xcp3->temper = Ts;
		cellp_xcp4->temper = Ts;
		cellp_xcp5->temper = Ts;
	}
	
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
		
		cellp_x0->temper   = Ts+deltaT;
		cellp_x1->temper   = Ts+deltaT;
		cellp_x2->temper   = Ts+deltaT;
	}
#endif
}


void rhobnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->rho   = cellp_ycp0->rho;
		cellp_y1->rho   = cellp_ycp1->rho;
		cellp_y2->rho   = cellp_ycp2->rho;
		cellp_ycp3->rho = cellp_y3->rho;
		cellp_ycp4->rho = cellp_y4->rho;
		cellp_ycp5->rho = cellp_y5->rho;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->rho   = cellp_xcp0->rho;
		cellp_x1->rho   = cellp_xcp1->rho;
		cellp_x2->rho   = cellp_xcp2->rho;
		cellp_xcp3->rho = cellp_x3->rho;
		cellp_xcp4->rho = cellp_x4->rho;
		cellp_xcp5->rho = cellp_x5->rho;
	}
#endif
	
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
			
		if((cellp_y3->rho>RHOmin)&&(cellp_y3->rho<RHOmax)){
			cellp_y0->rho   = 2.0*cellp_y3->rho-cellp_y6->rho;
			cellp_y1->rho   = 2.0*cellp_y3->rho-cellp_y5->rho;
			cellp_y2->rho   = 2.0*cellp_y3->rho-cellp_y4->rho;
//			cellp_y3->rho   = cellp_y4->rho;
		}
		else {
			cellp_y0->rho   = 2.0*cellp_y3->rho-cellp_y6->rho;
			cellp_y1->rho   = 2.0*cellp_y3->rho-cellp_y5->rho;
			cellp_y2->rho   = 2.0*cellp_y3->rho-cellp_y4->rho;
//			cellp_y3->rho   = cellp_y4->rho;
		}
		if((cellp_ycp2->rho>RHOmin)&&(cellp_ycp2->rho<RHOmax)){
			cellp_ycp3->rho = cellp_ycp2->rho-wetrho*dy*1.0;
			cellp_ycp4->rho = cellp_ycp1->rho-wetrho*dy*3.0;
				cellp_ycp5->rho = cellp_ycp0->rho-wetrho*dy*5.0;
		}
		else {
			cellp_ycp3->rho = cellp_ycp2->rho;
			cellp_ycp4->rho = cellp_ycp1->rho;
			cellp_ycp5->rho = cellp_ycp0->rho;
		}
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		if((cellp_x3->rho>RHOmin)&&(cellp_x3->rho<RHOmax)){
			cellp_x0->rho   = cellp_x5->rho-wetrho*dx*5.0;
			cellp_x1->rho   = cellp_x4->rho-wetrho*dx*3.0;
			cellp_x2->rho   = cellp_x3->rho-wetrho*dx*1.0;
		}
		else{
			cellp_x0->rho   = cellp_x5->rho;
			cellp_x1->rho   = cellp_x4->rho;
			cellp_x2->rho   = cellp_x3->rho;
			}
		if((cellp_xcp2->rho>RHOmin)&&(cellp_xcp2->rho<RHOmax)){
			cellp_xcp3->rho = cellp_xcp2->rho-wetrho*dx*1.0;
			cellp_xcp4->rho = cellp_xcp1->rho-wetrho*dx*3.0;
			cellp_xcp5->rho = cellp_xcp0->rho-wetrho*dx*5.0;
		}
		else{
			cellp_xcp3->rho = cellp_xcp2->rho;
			cellp_xcp4->rho = cellp_xcp1->rho;
			cellp_xcp5->rho = cellp_xcp0->rho;
		}
	}
#endif
}


void energybnd(void){
#ifndef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
		
		cellp_y0->energy   = cellp_ycp0->energy;
		cellp_y1->energy   = cellp_ycp1->energy;
		cellp_y2->energy   = cellp_ycp2->energy;
		cellp_ycp3->energy = cellp_y3->energy;
		cellp_ycp4->energy = cellp_y4->energy;
		cellp_ycp5->energy = cellp_y5->energy;
	}
#endif
#ifndef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->energy   = cellp_xcp0->energy;
		cellp_x1->energy   = cellp_xcp1->energy;
		cellp_x2->energy   = cellp_xcp2->energy;
		cellp_xcp3->energy = cellp_x3->energy;
		cellp_xcp4->energy = cellp_x4->energy;
		cellp_xcp5->energy = cellp_x5->energy;
	}
#endif
	
#ifdef WALLY
	for(int i=0; i<(XCELL_NUM+6); i++){
		indexbndy(i);
			
		cellp_y0->energy   = cellp_y4->energy;
		cellp_y1->energy   = cellp_y4->energy;
		cellp_y2->energy   = 2.0*cellp_y3->energy-cellp_y4->energy;
//		cellp_y3->energy   = cellp_y4->energy;
		cellp_ycp3->energy = cellp_ycp2->energy;
		cellp_ycp4->energy = cellp_ycp1->energy;
		cellp_ycp5->energy = cellp_ycp0->energy;
	}
#endif
#ifdef WALLX
	for(int j=0; j<(YCELL_NUM+6); j++){
		indexbndx(j);
			
		cellp_x0->energy   = cellp_x5->energy;
		cellp_x1->energy   = cellp_x4->energy;
		cellp_x2->energy   = cellp_x3->energy;
		cellp_xcp3->energy = cellp_xcp2->energy;
		cellp_xcp4->energy = cellp_xcp1->energy;
		cellp_xcp5->energy = cellp_xcp0->energy;
	}
#endif
}

void filewrite(void){
	double tecx,tecy,tecu,tecv,tecp,tect,tecrho,tece;
	char filename_tec[30];
	char filename_dat[30];
	FILE *tecout;
	FILE *energyout;
	double energyall = 0.0;

	/* Tecplotpo */
	printf("Storing data for Tecplot\n");
	sprintf(filename_tec,"d-%d.tec",n);
	tecout = fopen(filename_tec,"wb");
	fprintf(tecout,"variables = \"x\",\"y\",\"u\",\"v\",\"pres\",\"rho\",\"temper\",\"energy\"\n");
//	fprintf(tecout,"variables = \"x\",\"y\",\"u\",\"v\",\"pres\",\"divv\",\"rho\"\n");
//	fprintf(tecout,"zone i=%d j=%d f=point\n",XCELL_NUM+1,YCELL_NUM+3);
	fprintf(tecout,"zone i=%d j=%d f=point\n",XCELL_NUM+1,YCELL_NUM+1);
//	fe = freeeg();
	for(int j=3;j<=YCELL_NUM+3;j++){
		for(int i=3;i<=XCELL_NUM+3;i++){
			index_xm_ym = (j-1) + (i-1)*(YCELL_NUM+6);
			
			
			cellp_xm_ym = &CELL[index_xm_ym];
			
		
			index0      = j     + i*(YCELL_NUM+6);
			index_xp    = j     + (i+1)*(YCELL_NUM+6);
			index_yp    = (j+1) + i*(YCELL_NUM+6);
			index_xm    = j     + (i-1)*(YCELL_NUM+6);
			index_ym    = (j-1) + i*(YCELL_NUM+6);
			
			cellp= &CELL[index0];
			cellp_xp = &CELL[index_xp];
			cellp_yp = &CELL[index_yp];
			cellp_xm = &CELL[index_xm];
			cellp_ym = &CELL[index_ym];
			
			tecx = cellp->x-0.5*dx;
			tecy = cellp->y-0.5*dy;
			
/*			tecu = 0.25*(cellp->v_x + cellp_xm->v_x + cellp_ym->v_x + cellp_xm_ym->v_x);
			tecv = 0.25*(cellp->v_y + cellp_xm->v_y + cellp_ym->v_y + cellp_xm_ym->v_y);
			tecp = 0.25*(cellp->press + cellp_xm->press + cellp_ym->press + cellp_xm_ym->press);
//			tecp = cellp->press;
			tecrho = 0.25*(cellp->rho + cellp_xm->rho + cellp_ym->rho + cellp_xm_ym->rho);
			tect = 0.25*(cellp->temper + cellp_xm->temper + cellp_ym->temper + cellp_xm_ym->temper);
			tece = 0.25*(cellp->energy + cellp_xm->energy + cellp_ym->energy + cellp_xm_ym->energy);
*/			
			tecu = cellp->v_x;
			tecv = cellp->v_y;
			tecp = cellp->press;
			tecrho = cellp->rho;
//			tecrho = 0.5*kappa1*((cellp_xp->rho-cellp_xm->rho)*dxi*0.5*(cellp_xp->rho-cellp_xm->rho)*dxi*0.5
//								+(cellp_yp->rho-cellp_ym->rho)*dyi*0.5*(cellp_yp->rho-cellp_ym->rho)*dyi*0.5);
			tect = cellp->temper;
			tece = cellp->energy;
			
			
			energyall+=cellp->temper;
			
			fprintf(tecout,"%6.4lf %6.4lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf %20.16lf \n",tecx,tecy,tecu,tecv,tecp,tecrho,tect,tece);
		}
	}
	
	sprintf(filename_dat,"energy.dat");
	energyout = fopen(filename_dat,"a+");
	fprintf(energyout,"%d %f \n",n,energyall);
	
	fclose(tecout);
	fclose(energyout);
}
