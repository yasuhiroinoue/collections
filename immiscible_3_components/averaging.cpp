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
	int count;
	double col;
	_vec<double> vel;
	_avr :: _avr(void){};
};


main(){
	int X=96;
	int Y=64;
	_avr AVR[96][64];
	int num[20];
	
	
	//cout << "START, END STEP?" << endl;
	
	int step_end,step_start;
	step_end = 20000;
	step_start = 10000;
	//scanf("%d %d",&step_start,&step_end);
	
	
	num[0] = 3;
	num[1] = 5;
	num[2] = 8;
	num[3] = 10;
	num[4] = 13;
	num[5] = 15;
	num[6] = 18;
	num[7] = 20;
	num[8] = 23;
	num[9] = 25;
	num[10] = 28;
	num[11] = 30;
	num[12] = 33;
	num[13] = 35;
	num[14] = 38;
	num[15] = 40;
	num[16] = 43;
	num[17] = 45;
	num[18] = 48;
	num[19] = 50;
	
	
	for(int i = 2; i < 3; i++){
	memset(AVR,0,96*64*sizeof(_avr));
	char fn[30];
	//cout << "file name:" << endl;
	sprintf(fn,"d_U%05d",num[i]);
	
	
	
	for( int n = step_start; n <= step_end; n += 500){
			
			FILE *fptr;
			char fna[30];
			
			sprintf(fna,"%s_%07d.dat",fn,n);
			fptr = fopen(fna,"r");
			_vec<double> vel;
			_vec<int> loc;
			double col;
			
			cout << "now reading:" << fna << endl;
			
			char buf[1024];
			int n_r = 0;
			while( fgets( buf, 1024, fptr ) != NULL && n_r < 2){
				n_r++;
				cout << buf << endl;
			}
			
			while(fscanf(fptr,"%d %d %lf %lf %lf\n",&loc.x,&loc.y,&vel.x,&vel.y,&col) !=EOF){
				
				AVR[loc.x][loc.y].vel += vel;
				AVR[loc.x][loc.y].col += col;
				
			}
			fclose(fptr);
			cout << "reading done" << endl;
		
	}
	
	
	
	FILE *fp;
	char f[30];
	sprintf(f,"%s_AVC.dat",fn);
	
	fp = fopen(f,"w");
	double sample;
	
	
	sample = (double)((step_end - step_start)/500 + 1);
	cout << "sample=" << sample << endl;
	
	cout << "ok" << endl;

	for( int y = 0; y < Y; y++)
	for( int x = 0; x < X; x++){
		
		AVR[x][y].vel /= sample;
		AVR[x][y].col /= sample;
		
		if( x != 0 ){
			AVR[0][y].vel += AVR[x][y].vel;
			AVR[0][y].col += AVR[x][y].col;
		}
	
	}
	
	int y_0;
	int y_Y;
	
	y_0 = 0;//(int)( 18.0/sqrt(3.14159) + 1.0 );
	y_Y = Y-1;//(int)( 62 - 18.0/sqrt(3.14159) );
	double S=0.0;
	
	for(int y=y_0; y <= y_Y; y++){
		AVR[0][y].vel /= (double)X;
		AVR[0][y].col /= (double)X;
		S += AVR[0][y].col;
	}


	for(int y=y_0; y <= y_Y; y++){
		fprintf(fp,"%d %f %f %f\n",y,AVR[0][y].vel.x,AVR[0][y].vel.y,AVR[0][y].col/S);
	}
	
	
	fclose(fp);
	FILE *fp_plot;
	char fna_plot[30];
	sprintf(fna_plot,"plot.plt");
	fp_plot = fopen(fna_plot,"a");
	fprintf(fp_plot,"rep \"%s\" u 1:4 title \"%d\"\n",f,num[i]);
	fclose(fp_plot);
	}
}
