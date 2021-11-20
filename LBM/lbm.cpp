/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "system.h"
#include "lbm_class.h"
#include "defs.h"
#include "lbm.h"
#include "mymath.h"

/* カラーの番号として０（ゼロ）は使わない */

void Init(){

	memset(CELL,0,X*Y*sizeof(_cell));
	memset(PCELL,0,X*Y*sizeof(_pcell));
	memset(WALL,0,X*Y*sizeof(_wall));
	
	//
	for( int x = 0; x < X; x++){
		WALL[x][0].act = 1;
		WALL[x][0].rho = 1.0;

		WALL[x][Y-1].act = 1;
		WALL[x][Y-1].rho = 1.0;
	}
	
	for( int x = X/2 - 5; x < (X/2 + 6); x++){
		for( int y = Y/2 - 5; y < (Y/2 + 6); y++){
			WALL[x][y].act = 1;
			WALL[x][y].rho = 1.0;
		}
	}
	
	
	//
	
	//STM --> STreaMing
	
	STM[0].x = 0;
	STM[0].y = 0;
	
	STM[1].x = 1;
	STM[1].y = 0;
	
	STM[2].x = 1;
	STM[2].y = 1;
	
	STM[3].x = 0;
	STM[3].y = 1;
	
	STM[4].x = -1;
	STM[4].y = 1;
	
	STM[5].x = -1;
	STM[5].y = 0;
	
	STM[6].x = -1;
	STM[6].y = -1;
	
	STM[7].x = 0;
	STM[7].y = -1;
	
	STM[8].x = 1;
	STM[8].y = -1;
	
	EVEL[0].x = 0.00F;
	EVEL[0].y = 0.00F;
	
	EVEL[1].x = 1.00F;
	EVEL[1].y = 0.00F;
	
	EVEL[2].x = 1.00F;
	EVEL[2].y = 1.00F;
	
	EVEL[3].x = 0.00F;
	EVEL[3].y = 1.00F;
	
	EVEL[4].x = -1.00F;
	EVEL[4].y =  1.00F;
	
	EVEL[5].x = -1.00F;
	EVEL[5].y = 0.00F;
	
	EVEL[6].x = -1.00F;
	EVEL[6].y = -1.00F;
	
	EVEL[7].x = 0.00F;
	EVEL[7].y = -1.00F;
	
	EVEL[8].x = 1.00F;
	EVEL[8].y = -1.00F;
	
	WT[0] = 4.0/9.0;
	for( int i = 1; i < 9; i += 2){
		WT[i] = 1.0/9.0;
	}
	for( int i = 2; i < 9; i += 2){
		WT[i] = 1.0/36.0;
	}
	
	
	/*
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
			for( int col = 0; col < COL_NUM; col++){
				for( int i = 0; i < 9; i++){
					CELL[x][y].f[col][i] = WT[i];
				}
				CELL[x][y].col[col] = (col + 1 );
				CELL[x][y].col_end = col;
			}
		}
	}
	*/
	
	/*
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
			for( int i = 0; i < 9; i++){
				CELL[x][y].f[0][i] = WT[i];
			}
		
			if( x < (int)(X/2.0) && y < (int)(Y/2.0) ){
				CELL[x][y].col[0] = 1;
				CELL[x][y].col_end = 0;
			}else if( x < (int)(X/2.0) ){
				CELL[x][y].col[0] = 2;
				CELL[x][y].col_end = 0;

			}else if( y < (int)(Y/2.0)){
				CELL[x][y].col[0] = 2;
				CELL[x][y].col_end = 0;
			}else{
				CELL[x][y].col[0] = 1;
				CELL[x][y].col_end = 0;
			}
			
		}
	}
	*/
	
	_vec<double> vel;
	double rho = 1.0;
	double feq[9];
	int cx0 = (int)(X/4.0); int cy0 = (int)(Y/4.0);
	int cx1 = (int)((3.0*X)/2.0); int cy1 = (int)(Y/2.0);
	double R = 12.0;
	
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
		if( WALL[x][y].act != 1 ){
			double R20 = (double)( SQR((x-cx0)) + SQR((y-cy0)) );
			double R21 = (double)( SQR((x-cx1)) + SQR((y-cy1)) );
			
			if( R20 < (R*R) ){
				vel.IN(0.0, 0.0);
				(void)feq_calc(vel,&feq[0],rho);
				for( int i = 0; i < 9; i++){
					CELL[x][y].f[0][i] = feq[i];
				}
				CELL[x][y].col[0] = 2;
				CELL[x][y].col_end = 0;
			}else if( R21 < /*(R*R)*/ -1.0 ){
				vel.IN(0.0, 0.0);
				(void)feq_calc(vel,&feq[0],rho);
				for( int i = 0; i < 9; i++){
					CELL[x][y].f[0][i] = feq[i];
				}
				CELL[x][y].col[0] = 3;
				CELL[x][y].col_end = 0;
			}else{
				vel.IN(0.0,0.0);
				(void)feq_calc(vel,&feq[0],rho);
				for( int i = 0; i < 9; i++){
					CELL[x][y].f[0][i] = feq[i];
				}
				CELL[x][y].col[0] = 1;
				CELL[x][y].col_end = 0;
			}
			
			/*
			vel.IN(0.0,0.0);
			(void)feq_calc(vel,&feq[0],rho);
			for( int i = 0; i < 9; i++){
				CELL[x][y].f[0][i] = feq[i];
			}
			
			if( y < (int)(Y/2.0) ){
				CELL[x][y].col[0] = 1;
				CELL[x][y].col_end = 0;
			}else{
				CELL[x][y].col[0] = 2;
				CELL[x][y].col_end = 0;
			}
			*/
		}
		}
	}

}

void feq_calc(_vec<double> vel, double* feq, double rho){

	double u2_cs2_2 = (vel*vel)/(2.0*Cs2);
	feq[0] = rho*WT[0]*(1.0 - u2_cs2_2 );
	for( int i = 1; i < 9; i++){
		double uiu_cs2= (EVEL[i]*vel)/Cs2;
		feq[i] = rho*WT[i]*(1.0+ uiu_cs2 + 0.50*SQR(uiu_cs2) - u2_cs2_2 );
	}
}

void Heat_Bath(){
	/*
	for( int x = 0; x < X; x++){
		_cell* cellp = &CELL[x][0];
		double rho = 1.0;

		double feq[9];
		(void)feq_calc(_vec<double>(0.0,0.0),&feq[0],rho);
		
		for( int i = 0; i < 9; i++){
			cellp->f[0][i] = feq[i];
		}
		cellp->col[0] = 1;
		cellp->col_end = 0;
	}
	
	
	for( int x = 0; x < X; x++){
		_cell* cellp = &CELL[x][Y-1];
		double rho = 1.0;

		double feq[9];
		(void)feq_calc(_vec<double>(0.0,0.0),&feq[0],rho);
		
		for( int i = 0; i < 9; i++){
			cellp->f[0][i] = feq[i];
		}
		cellp->col[0] = 1;
		cellp->col_end = 0;
	}
	
	for( int y = 0; y < Y; y++){
		_cell* cellp = &CELL[0][y];
		double rho = 1.0;

		double feq[9];
		(void)feq_calc(_vec<double>(0.0,0.0),&feq[0],rho);
		
		for( int i = 0; i < 9; i++){
			cellp->f[0][i] = feq[i];
		}
		cellp->col[0] = 1;
		cellp->col_end = 0;
	}
	
	for( int y = 0; y < Y; y++){
		_cell* cellp = &CELL[X-1][y];
		double rho = 1.0;

		double feq[9];
		(void)feq_calc(_vec<double>(0.0,0.0),&feq[0],rho);
		
		for( int i = 0; i < 9; i++){
			cellp->f[0][i] = feq[i];
		}
		cellp->col[0] = 1;
		cellp->col_end = 0;
	}
	*/
	for( int y = 1; y < Y-1; y++){
		_cell* cellp = &CELL[0][y];
		double rho = 1.0;

		double feq[9];
		(void)feq_calc(_vec<double>(0.0,0.0),&feq[0],rho);
		
		for( int i = 0; i < 9; i++){
			cellp->f[0][i] = feq[i];
		}
		cellp->col[0] = 1;
		cellp->col_end = 0;
	}
	
	for( int y = 1; y < Y-1; y++){
		_cell* cellp = &CELL[X-1][y];
		double rho = 0.90;

		double feq[9];
		(void)feq_calc(_vec<double>(0.0,0.0),&feq[0],rho);
		
		for( int i = 0; i < 9; i++){
			cellp->f[0][i] = feq[i];
		}
		cellp->col[0] = 1;
		cellp->col_end = 0;
	}
}

void Streaming(){
	
	memset(CELL,0,X*Y*sizeof(_cell));
	
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){

			for( int i = 0; i < 9; i++){
				int px = x + STM[i].x;
				int py = y + STM[i].y;
				(void)Periodic(px,py);
				_cell* cellp = &CELL[px][py];
				_pcell* pcell = &PCELL[x][y];
				
				if( cellp->col[0] == 0 ){
					for( int col = 0; col <= pcell->col_end; col++){
						cellp->f[col][i] = pcell->f[col][i];
						cellp->col[col] = pcell->col[col];
					}
					cellp->col_end = pcell->col_end;
				//	std::cout << cellp->col_end << std::endl;
				}else{
					for( int col = 0; col <= pcell->col_end; col++){
						int flag = 0;
						for( int c = 0; c <= cellp->col_end; c++){
							
							if( cellp->col[c] == pcell->col[col] ){
								cellp->f[c][i] = pcell->f[col][i];
								flag = 1;
							}
							
						}
						
						if( flag == 0 ){
							if( pcell->col[col] != 0 ){
								cellp->col_end++;
								if( cellp->col_end >= COL_NUM ){
									std::cout << "memory violation" << std::endl;
									std::cout << cellp->col[0] << std::endl;
									std::cout << pcell->col[0] << std::endl;
									std::cout << x << " " << y << std::endl;
									std::cout << i << std::endl;
									exit(0);
								}
							
								int c = cellp->col_end;
								cellp->f[c][i] = pcell->f[col][i];
								cellp->col[c] = pcell->col[col];
							}
						}
					}
				}
				
				
			}
			
			
			
		}
	}
	
	(void)Boundary();
	
	/*
	for( int x = 0; x < X; x++){
		for( int y = 1; y < (Y-1); y++){
			for( int i = 0; i < 9; i++){
				int px = x + STM[i].x;
				int py = y + STM[i].y;
				int dummy = 1;
				(void)Periodic(px,dummy);
				CELL[px][py].f[i] = CELL[x][y].pf[i];
			}
		}
	}
	*/
	(void)Heat_Bath();
}

void Boundary(){
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){

			if( WALL[x][y].act == 1 ){
				_cell* cellp = &CELL[x][y];
				
				
				for( int c = 0; c <= cellp->col_end; c++){
						double tmp;
						for( int i = 1; i <= 4; i++){
							tmp = cellp->f[c][i+4];
							cellp->f[c][i+4] = cellp->f[c][i];
							cellp->f[c][i] = tmp;
						}
				}
			}
		}
	}
}



void CellCalc(){
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
			/*
			int col = COL_NUM - 1;
			while( col > 0 && CELL[x][y].col[col] == 0 ){
				col--;
			}
			CELL[x][y].col_end = col;*/
			
			for( int c = 0; c <= CELL[x][y].col_end; c++){
				for( int i = 0; i < 9; i++){
					CELL[x][y].rho[c] += CELL[x][y].f[c][i];
					CELL[x][y].rho_total += CELL[x][y].f[c][i];
					
				//	CELL[x][y].vel[c] += ( CELL[x][y].f[c][i]*EVEL[i] );
					CELL[x][y].vel_total += ( CELL[x][y].f[c][i]*EVEL[i] );
				}
				/*
				if( CELL[x][y].rho[c] > 0.0 ){
					CELL[x][y].vel[c] /= CELL[x][y].rho[c];
				}
				*/
			}
			
			if( CELL[x][y].rho_total > 0 ){
				CELL[x][y].vel_total /= CELL[x][y].rho_total;
			}
		}
	}
	
/* BUG CHECK	
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
			double col = 0.0;
			for( int c = 0; c <= CELL[x][y].col[c]; c++){
				col += CELL[x][y].rho[c] * CELL[x][y].col[c];
			}
			if( CELL[x][y].rho_total > 0.0){
				col /= CELL[x][y].rho_total;
			}
			std::cout << x << " "<< y << " " << col << std::endl;
		}
	}
*/
	
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
			if( WALL[x][y].act != 1 ){
			/* カラー計算始め */
			double DL = 4.000;
			double GAA[4];
			GAA[0] = G11; GAA[1] = G12;
			GAA[2] = G22; GAA[3] = G23;
			
			int px = x + 1; int nx = x - 1;
			int py = y + 1; int ny = y - 1;
			(void)Periodic(px,py); (void)Periodic(nx,ny);
			
			_vec<double> FLD[(CELL[x][y].col_end + 1)];
			
			
			for( int ci = 0; ci <= CELL[x][y].col_end; ci++){
				double rho = CELL[x][y].rho[ci];
				
				if( CELL[x][y].col[ci] == 1 ){
					
				}else{
					FLD[ci].x += WL * rho * WALL[px][y].rho;
					FLD[ci].y += WL * rho * WALL[x][py].rho;
					FLD[ci].x -= WL * rho * WALL[nx][y].rho;
					FLD[ci].y -= WL * rho * WALL[x][ny].rho;
					
					FLD[ci].x += WL * rho * WALL[px][py].rho/DL;
					FLD[ci].y += WL * rho * WALL[px][py].rho/DL;
					
					FLD[ci].x -= WL * rho * WALL[nx][py].rho/DL;
					FLD[ci].y += WL * rho * WALL[nx][py].rho/DL;
					
					FLD[ci].x -= WL * rho * WALL[nx][ny].rho/DL;
					FLD[ci].y -= WL * rho * WALL[nx][ny].rho/DL;
					
					FLD[ci].x += WL * rho * WALL[px][ny].rho/DL;
					FLD[ci].y -= WL * rho * WALL[px][ny].rho/DL;
				}
			}
			
			for( int ci = 0; ci <= CELL[x][y].col_end; ci++){
				if( CELL[x][y].col[ci] == 1 ){
				/**/
					for( int cj = 0; cj <= CELL[px][y].col_end; cj++){
						if( CELL[px][y].col[cj] == 1 ){
							FLD[ci].x += GAA[0] * CELL[x][y].rho[ci]*CELL[px][y].rho[cj];
						}else{
							FLD[ci].x += GAA[1] * CELL[x][y].rho[ci]*CELL[px][y].rho[cj];
						}
					}
					
					for( int cj = 0; cj <= CELL[nx][y].col_end; cj++){
						if( CELL[nx][y].col[cj] == 1 ){
							FLD[ci].x -= GAA[0] * CELL[x][y].rho[ci]*CELL[nx][y].rho[cj];
						}else{
							FLD[ci].x -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][y].rho[cj];
						}
					}
					
					for( int cj = 0; cj <= CELL[x][py].col_end; cj++){
						if( CELL[x][py].col[cj] == 1 ){
							FLD[ci].y += GAA[0] * CELL[x][y].rho[ci]*CELL[x][py].rho[cj];
						}else{
							FLD[ci].y += GAA[1] * CELL[x][y].rho[ci]*CELL[x][py].rho[cj];
						}
					}
					
					for( int cj = 0; cj <= CELL[x][ny].col_end; cj++){
						if( CELL[x][ny].col[cj] == 1 ){
							FLD[ci].y -= GAA[0] * CELL[x][y].rho[ci]*CELL[x][ny].rho[cj];
						}else{
							FLD[ci].y -= GAA[1] * CELL[x][y].rho[ci]*CELL[x][ny].rho[cj];
						}
					}
					
					
					
					//
					for( int cj = 0; cj <= CELL[px][py].col_end; cj++){
						if( CELL[px][py].col[cj] == 1 ){
							FLD[ci].x += GAA[0] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
							FLD[ci].y += GAA[0] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
						}else{
							FLD[ci].x += GAA[1] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
							FLD[ci].y += GAA[1] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
						}
					}
					
					for( int cj = 0; cj <= CELL[nx][py].col_end; cj++){
						if( CELL[nx][py].col[cj] == 1 ){
							FLD[ci].x -= GAA[0] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
							FLD[ci].y += GAA[0] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
						}else{
							FLD[ci].x -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
							FLD[ci].y += GAA[1] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
						}
					}

					for( int cj = 0; cj <= CELL[nx][ny].col_end; cj++){
						if( CELL[nx][ny].col[cj] == 1 ){
							FLD[ci].x -= GAA[0] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[0] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
						}else{
							FLD[ci].x -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
						}
					}

					for( int cj = 0; cj <= CELL[px][ny].col_end; cj++){
						if( CELL[px][ny].col[cj] == 1 ){
							FLD[ci].x += GAA[0] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[0] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
						}else{
							FLD[ci].x += GAA[1] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[1] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
						}
					}

					
				/**/
				}else{
				/**/
					for( int cj = 0; cj <= CELL[px][y].col_end; cj++){
						if( CELL[px][y].col[cj] == 1 ){
							FLD[ci].x += GAA[1] * CELL[x][y].rho[ci]*CELL[px][y].rho[cj];
						}else if( CELL[x][y].col[ci] == CELL[px][y].col[cj] ){
							FLD[ci].x += GAA[2] * CELL[x][y].rho[ci]*CELL[px][y].rho[cj];
						}else{
							FLD[ci].x += GAA[3] * CELL[x][y].rho[ci]*CELL[px][y].rho[cj];
						}
					}
					
					
					for( int cj = 0; cj <= CELL[nx][y].col_end; cj++){
						if( CELL[nx][y].col[cj] == 1 ){
							FLD[ci].x -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][y].rho[cj];
						}else if( CELL[x][y].col[ci] == CELL[nx][y].col[cj] ){
							FLD[ci].x -= GAA[2] * CELL[x][y].rho[ci]*CELL[nx][y].rho[cj];
						}else{
							FLD[ci].x -= GAA[3] * CELL[x][y].rho[ci]*CELL[nx][y].rho[cj];
						}
					}
					
					
					for( int cj = 0; cj <= CELL[x][py].col_end; cj++){
						if( CELL[x][py].col[cj] == 1 ){
							FLD[ci].y += GAA[1] * CELL[x][y].rho[ci]*CELL[x][py].rho[cj];
						}else if( CELL[x][y].col[ci] == CELL[x][py].col[cj] ){
							FLD[ci].y += GAA[2] * CELL[x][y].rho[ci]*CELL[x][py].rho[cj];
						}else{
							FLD[ci].y += GAA[3] * CELL[x][y].rho[ci]*CELL[x][py].rho[cj];
						}
					}
					
					for( int cj = 0; cj <= CELL[x][ny].col_end; cj++){
						if( CELL[x][ny].col[cj] == 1 ){
							FLD[ci].y -= GAA[1] * CELL[x][y].rho[ci]*CELL[x][ny].rho[cj];
						}else if( CELL[x][y].col[ci] == CELL[x][ny].col[cj] ){
							FLD[ci].y -= GAA[2] * CELL[x][y].rho[ci]*CELL[x][ny].rho[cj];
						}else{
							FLD[ci].y -= GAA[3] * CELL[x][y].rho[ci]*CELL[x][ny].rho[cj];
						}
					}
					

					//
					for( int cj = 0; cj <= CELL[px][py].col_end; cj++){
						if( CELL[px][py].col[cj] == 1 ){
							FLD[ci].x += GAA[1] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
							FLD[ci].y += GAA[1] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
						}else if( CELL[x][y].col[ci] == CELL[px][py].col[cj] ){
							FLD[ci].x += GAA[2] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
							FLD[ci].y += GAA[2] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
						}else{
							FLD[ci].x += GAA[3] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
							FLD[ci].y += GAA[3] * CELL[x][y].rho[ci]*CELL[px][py].rho[cj]/DL;
						}
					}

					for( int cj = 0; cj <= CELL[nx][py].col_end; cj++){
						if( CELL[nx][py].col[cj] == 1 ){
							FLD[ci].x -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
							FLD[ci].y += GAA[1] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
						}else if( CELL[x][y].col[ci] == CELL[nx][py].col[cj] ){
							FLD[ci].x -= GAA[2] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
							FLD[ci].y += GAA[2] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
						}else{
							FLD[ci].x -= GAA[3] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
							FLD[ci].y += GAA[3] * CELL[x][y].rho[ci]*CELL[nx][py].rho[cj]/DL;
						}
					}

					for( int cj = 0; cj <= CELL[nx][ny].col_end; cj++){
						if( CELL[nx][ny].col[cj] == 1 ){
							FLD[ci].x -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[1] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
						}else if( CELL[x][y].col[ci] == CELL[nx][ny].col[cj] ){
							FLD[ci].x -= GAA[2] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[2] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
						}else{
							FLD[ci].x -= GAA[3] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[3] * CELL[x][y].rho[ci]*CELL[nx][ny].rho[cj]/DL;
						}
					}

					for( int cj = 0; cj <= CELL[px][ny].col_end; cj++){
						if( CELL[px][ny].col[cj] == 1 ){
							FLD[ci].x += GAA[1] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[1] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
						}else if( CELL[x][y].col[ci] == CELL[px][ny].col[cj] ){
							FLD[ci].x += GAA[2] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[2] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
						}else{
							FLD[ci].x += GAA[3] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
							FLD[ci].y -= GAA[3] * CELL[x][y].rho[ci]*CELL[px][ny].rho[cj]/DL;
						}
					}

				/**/
				}
			}
			
			for( int ci = 0; ci <= CELL[x][y].col_end; ci++){
				if( CELL[x][y].rho[ci] > 0.0){
					FLD[ci] /= CELL[x][y].rho[ci];
					CELL[x][y].vel[ci] = CELL[x][y].vel_total + (TAU*FLD[ci]);
/* Luo's accelaration */
//					CELL[x][y].vel[ci] = CELL[x][y].vel_total;
				}else{
					CELL[x][y].vel[ci] = _vec<double>(0.0,0.0);
					FLD[ci].IN(0.0,0.0);
				}
			}
			
			/* カラー計算終わり */
			
			
			for( int c = 0; c <= CELL[x][y].col_end; c++){
				double u2_cs2_2 = (CELL[x][y].vel[c] * CELL[x][y].vel[c])/(2.0*Cs2);
				CELL[x][y].feq[c][0] = CELL[x][y].rho[c]*WT[0]*(1.0 - u2_cs2_2 );
				
				/* correction */
				double FF = (FLD[c]*FLD[c]);
				CELL[x][y].feq[c][0] += (1.50*WT[0]*CELL[x][y].rho[c]*FF*TAU*TAU);
				
				
				for( int i = 1; i < 9; i++){
					double uiu_cs2= (EVEL[i]*CELL[x][y].vel[c])/Cs2;
					CELL[x][y].feq[c][i] = CELL[x][y].rho[c]*WT[i]*(1.0+ uiu_cs2 + 0.50*SQR(uiu_cs2) - u2_cs2_2 );
					
					/* correction */
					double FE2 = (FLD[c]*EVEL[i]); FE2 *= FE2;
					CELL[x][y].feq[c][i] -= (0.50*WT[i]*CELL[x][y].rho[c]*TAU*TAU*(9.0*FE2 - 3.0*FF));
					
/* Luo's accelaration */
//					CELL[x][y].foc[c][i] = (WT[i]*CELL[x][y].rho[c])*(FLD[c]*EVEL[i])/(Cs2);
/*					double Wa = (3.0 * WT[i]*CELL[x][y].rho[c]);					
					CELL[x][y].foc[c][i] = (EVEL[i]-CELL[x][y].vel_total)*FLD[c];
					CELL[x][y].foc[c][i] += (3.0 * ( EVEL[i]*CELL[x][y].vel_total) ) *(EVEL[i]*FLD[c]);
					CELL[x][y].foc[c][i] *= Wa;
					*/
				}
			}
			
		  }//if( WALL[x][y].act !=1 ) に対応する括弧
		}
	}

}

void Stress_Measure(){
	
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
			
			_cell* cellp = &CELL[x][y];
			
			for( int c = 0; c <= cellp->col_end;c++){
				for( int i = 1; i < 9; i++){
					cellp->s[0] += cellp->f[c][i];
					cellp->s[2] += cellp->f[c][i];
				}
				cellp->s[0] -= ( cellp->f[c][3] + cellp->f[c][7]);
				cellp->s[2] -= ( cellp->f[c][1] + cellp->f[c][5]);
				
				cellp->s[1] += cellp->f[c][2];
				cellp->s[1] += cellp->f[c][6];
				cellp->s[1] -= cellp->f[c][4];
				cellp->s[1] -= cellp->f[c][8];
			}
		}
	}
}

void Collision(){
	memset(PCELL,0,X*Y*sizeof(_pcell));

	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){

			_cell* cellp = &CELL[x][y];
			_pcell* pcell = &PCELL[x][y];
			pcell->col_end = cellp->col_end;
			
			if( WALL[x][y].act != 1 ){
				for( int c = 0; c <= cellp->col_end; c++){
					pcell->col[c] = cellp->col[c];
				
					for( int i = 0; i < 9; i++){
						pcell->f[c][i] = cellp->feq[c][i]/TAU + (1.0 - 1.0/TAU)*cellp->f[c][i];
/* Luo's accelaration */
//					pcell->f[c][i] = cellp->feq[c][i]/TAU + (1.0 - 1.0/TAU)*cellp->f[c][i] + cellp->foc[c][i];
					}
				}
			}else{
				for( int c = 0; c <= cellp->col_end; c++){
					pcell->col[c] = cellp->col[c];
				
					for( int i = 0; i < 9; i++){
						pcell->f[c][i] = cellp->f[c][i];
					}
				}
			}
		}
	}
}

void Output(int step){
	char fna[30];
	sprintf(fna,"d%07d.tec",step);
	FILE* fptr = fopen(fna,"w");
	
	fprintf(fptr,"variables = \"x\",\"y\",\"rho\",\"c\",\"u\",\"v\",\"s0\",\"s1\",\"s2\"\n");
	fprintf(fptr,"zone i=%d j=%d f=point\n",X,Y);
	
	for( int y = 0; y < Y; y++){
		for( int x = 0; x < X; x++){
			double col = 0.0;
			
			/*
			for( int c = 0; c <= CELL[x][y].col_end; c++){
				int color = CELL[x][y].col[c];
				col +=  (double)color * CELL[x][y].rho[c];
			}
			*/
			
			col = CELL[x][y].rho_total;
			for( int c = 0; c <= CELL[x][y].col_end; c++){
				int color = CELL[x][y].col[c];
				if( color == 1 ){
					col -=  CELL[x][y].rho[c];
				}
			}
			
			if( CELL[x][y].rho_total > 0.0){
				col /= CELL[x][y].rho_total;
			}
			
			_vec<double> vel = CELL[x][y].vel_total;
			if( WALL[x][y].act == 1 ){ vel.IN(0.0,0.0); col = -1.0; }
			
			fprintf(fptr,"%d %d %lf %lf %lf %lf",x,y,CELL[x][y].rho_total,col,vel.x, vel.y);
			for( int j = 0; j < 3; j++)
			fprintf(fptr," %lf",CELL[x][y].s[j]);
/*			for( int i = 0; i < 9; i++){
				fprintf(fptr," %lf",CELL[x][y].f[i]);
			}
*/
			fprintf(fptr,"\n");
		}
	}
	
	fclose(fptr);
}
