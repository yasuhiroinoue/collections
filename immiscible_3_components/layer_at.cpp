/*
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "vec.h"
#include "mymath.h"

using namespace std;


class _avr{
	public:
	_vec<double> vel;
	double col;
	_avr :: _avr(void){};
};


main(){
	int X=96;
	int Y=64;
	_avr AVR[96][64];
	int num[20];
	
	num[0] = 800;
	num[1] = 1050;
	num[2] = 1300;
	num[3] = 1549;
	num[4] = 1799;
	num[5] = 2049;
	num[6] = 2299;
	num[7] = 2550;
	num[8] = 2800;
	num[9] = 3050;
	num[10] = 3300;
	num[11] = 3550;
	num[12] = 3800;
	num[13] = 4050;
	num[14] = 4300;
	num[15] = 4550;
	num[16] = 4800;
	num[17] = 5050;
	
	
	
	for( int inum = 1; inum < 18; inum++){
	memset(AVR,0,96*64*sizeof(_avr));
	char fna[30];
	memset(fna,0,30*sizeof(char));
	cout << "file name:" << endl;
	sprintf(fna,"avr_U0%d_0030000.dat",num[inum]);
	
	FILE *fptr;
	
	fptr = fopen(fna,"r");
	_vec<double> vel;
	_vec<double> loc;
	double col;
	
	cout << "now reading:" << fna << endl;
	
	char buf[1024];
	int n_r = 0;
	while( fgets( buf, 1024, fptr ) != NULL && n_r < 2){
		n_r++;
		cout << buf << endl;
	}
	float tx,ty,tz;
	float m;
	
	while(fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&loc.x, &loc.y, &vel.x, &vel.y, &m, &col, &tx, &ty, &tz) !=EOF){
		int x,y;
		
		x = (int)loc.x;
		y = (int)loc.y;
		
		AVR[x][y].vel += vel;
		AVR[x][y].col += col;
		
	}
	fclose(fptr);
	cout << "reading done" << endl;
	
	for( int x= 1; x < X; x++){
		for( int y =0; y < Y; y++){
			AVR[0][y].vel += AVR[x][y].vel;
			AVR[0][y].col += AVR[x][y].col;
		}
	}
	
	double dvd = 0.0;

	for( int y = 0; y < Y; y++){
		dvd += AVR[0][y].col;
	}
	
	char wfn[30];
	FILE *wfp;
	sprintf(wfn,"c%s",fna);
	wfp = fopen(wfn,"w");
	for( int y = 0; y < Y; y++ ){
		AVR[0][y].vel /= (double)X;
		AVR[0][y].col /= (dvd*(double)X);
		fprintf(wfp,"%d %lf %lf %lf\n",y,AVR[0][y].vel.x,AVR[0][y].vel.y,AVR[0][y].col);
	}
	fclose(wfp);
	
}
}


