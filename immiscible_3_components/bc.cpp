/*
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#include <iostream>
#include <math.h>
#include "mymath.h"
#include "system.h"
#include "rclass.h"


void Periodic(int& x, int& y){
	while( x < 0 ) x += (int)X;
	while( x >= (int)X ) x -= (int)X;
	while( y < 0 ) y += (int)Y;
	while( y >= (int)Y ) y -= (int)Y;
}

void Periodic(double& x, double& y){
	while( x < 0.0 ) x += X;
	while( x >= X ) x -= X;
	while( y < 0.0 ) y += Y;
	while( y >= Y ) y -= Y;
}

void Wall(_ptcl_base& p,double& x,double& y,double& SX){
	if( y < 1.0 ){
		double tau = (1.0 - p.loc.y) / p.vel.y;
		x -= p.vel.x * tau;
		p.vel.x = (double)gauss(T/p.mas,VX);
		p.vel.y = (double)gamma_noflow(T/p.mas);
		
		tau = RAND;
		x += p.vel.x * tau;
		y =  1.0 + p.vel.y * tau;
	}else if( y >= (Y - 1) ){
		double tau = (Y - 1.0 - p.loc.y) / p.vel.y;
		x -= p.vel.x * tau;
		p.vel.x =  (double)gauss(T/p.mas,(VX+SX));
		p.vel.y = -(double)gamma_noflow(T/p.mas);
		
		tau = RAND;
		x += p.vel.x * tau;
		y =  Y - 1.0 - DELTA + p.vel.y * tau;
	}
	(void)Periodic(x,y);
}
