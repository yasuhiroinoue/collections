#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "defs.h"
#include "system.h"

int main(void){
	
	printf("read ptcl.dat? Y-->1, N-->0\n");
	int i;
	scanf("%d",&i);
	
	if( i == 1 ){
		printf("Reading data\n");
		(void)Read_Ptcl();
	}else{
		printf("Initializing\n");
		(void)Initialize();
	}
	printf("A12 = %f\n",A12);
	(void)CellCalculation();
	(void)SystemMonitor(0);
	(void)Output(0,0);
	
	double SX=0.0;
	printf("start step: ");
	int step_start=1;
	scanf("%d",&step_start);
	for( int step = step_start; step < STEP_END; step++){
		(void)Propagation(SX);
		(void)CellCalculation();
		
		if( step > START_STEP ){
			 (void)Mean();
			 (void)Stress_Measure();
		
		
			if( step % MTIME == 0 ){
				(void)AVR_output(step,SX);
			}
		
		}
		
		(void)MultiColorCollision();
		(void)Collision();
		(void)SystemMonitor(step);

		if( step % STEP_OUT == 0 )
			(void)Output(step,SX);
	
		if( step % PTCL_OUT == 0 )
			(void)Ptcl_Output();
		
		cout << step << " step\n";
	}
	(void)Ptcl_Output();
}
