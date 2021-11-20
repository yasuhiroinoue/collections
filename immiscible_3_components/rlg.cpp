/*
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
//rlg.cpp
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "mymath.h"
#include "system.h"
#include "rlg.h"
#include "defs.h"
#include <map>
using namespace std;


void Read_Ptcl(){
	ptcl_size = PN;
	cell_size = X * Y;
	
	
	ptcl = new _ptcl_base[ptcl_size];
	memset( ptcl, 0 , ptcl_size * sizeof(_ptcl_base));
	
	memset( cell, 0, cell_size * sizeof(_cell_base));
	memset( node, 0, cell_size * sizeof(_node));
	
	
	FILE *fptr;
	char fna[30];
	
	sprintf(fna,"ptcl.dat");
	fptr = fopen(fna,"r");
	
	_vec<double> loc;
	_vec<double> vel;
	double mas;
	int col;
	int i =0;
	
	while(fscanf(fptr,"%d %lf %lf %lf %lf %lf",&col, &mas, &loc.x, &loc.y, &vel.x, &vel.y  ) != EOF ){
		if( i >= ptcl_size ){
		cout << "over ptcl_size\n" << endl;
			exit(0);
		}
		ptcl[i].col = col;
		ptcl[i].mas = mas;
		ptcl[i].loc.IN(loc);
		ptcl[i].vel.IN(vel);
		_vec<int> iloc = ptcl[i].loc.icast();
		
		if( iloc.x < 0 || iloc.x >= X || iloc.y < 0 || iloc.y >= Y ){
			cout << "reading error\n" << endl;
			exit(0);
		}
		
		list[iloc.x][iloc.y].pmap.insert(std::make_pair(ptcl[i].col,&ptcl[i]));
		
		i++;
	}
	fclose(fptr);
	
	for( int x = 0; x < X; x++)
	for( int y = 0; y < Y; y++){
		int xp = x + 1; int yp = y + 1; (void)Periodic(xp,yp);
		int xn = x - 1; int yn = y - 1; (void)Periodic(xn,yn);
		list[x][y].lptr[0] = &list[xp][y];
		list[x][y].lptr[1] = &list[x][yp];
		list[x][y].lptr[2] = &list[xn][y];
		list[x][y].lptr[3] = &list[x][yn];
		//         1
		//        ��
		//   2 ��   �� 0
		//        ��
		//         3
		list[x][y].lptr[4] = &list[xp][yp];
		list[x][y].lptr[5] = &list[xn][yp];
		list[x][y].lptr[6] = &list[xn][yn];
		list[x][y].lptr[7] = &list[xp][yn];
		//       5 1  4
		//       ������
		//   2 ��    �� 0
		//       ������
		//      6  3  7
		
		
		
	}
	cout << "ok\n" << endl;
	
}



void Initialize(void){
	ptcl_size = PN;
	cell_size = X * Y;
	
	
	ptcl = new _ptcl_base[ptcl_size];
	memset( ptcl, 0 , ptcl_size * sizeof(_ptcl_base));
	

	memset( cell, 0, cell_size * sizeof(_cell_base));

	memset( node, 0, cell_size * sizeof(_node));
	
	_vec<int> loc;
	for( int i = 0; i < ptcl_size; i++){
		ptcl[i].mas = Mass;
	//	ptcl[i].loc.IN(RAND*X,RAND*(Y - 2.0) + 1.0);
		ptcl[i].loc.IN(RAND*X,RAND*Y);
		ptcl[i].vel.IN((double)gauss(T/Mass,VY),(double)gauss(T/Mass,VY));
	
	
	
	
		ptcl[i].col = (int)(RAND*SPCS);
		loc = ptcl[i].loc.icast();
		list[loc.x][loc.y].pmap.insert(std::make_pair(ptcl[i].col,&ptcl[i]));
	}
	
	
	for( int x = 0; x < X; x++)
	for( int y = 0; y < Y; y++){
		int xp = x + 1; int yp = y + 1; (void)Periodic(xp,yp);
		int xn = x - 1; int yn = y - 1; (void)Periodic(xn,yn);
		list[x][y].lptr[0] = &list[xp][y];
		list[x][y].lptr[1] = &list[x][yp];
		list[x][y].lptr[2] = &list[xn][y];
		list[x][y].lptr[3] = &list[x][yn];
		//         1
		//        ��
		//   2 ��   �� 0
		//        ��
		//         3
		list[x][y].lptr[4] = &list[xp][yp];
		list[x][y].lptr[5] = &list[xn][yp];
		list[x][y].lptr[6] = &list[xn][yn];
		list[x][y].lptr[7] = &list[xp][yn];
		//       5 1  4
		//       ������
		//   2 ��    �� 0
		//       ������
		//      6  3  7
		
		
		
	}
	
	/*
	int col=1;
	
	for( int x = 5; x < X - 4; x+= 10){
		for( int y = 7; y < Y - 6; y+=11){
			_list* lst = &list[x][y];
			multimap<int, _ptcl_base*>::iterator pos;
			
			for( pos = lst->pmap.begin(); pos != lst->pmap.end(); ++pos){
				pos->second->mas = Mass_H;
				pos->second->col = col;
				pos->second->vel.IN((double)gauss(T/pos->second->mas,VX),(double)gauss(T/pos->second->mas,VY));
			}
			
			for( int i = 0; i < 8; i++){
				for( pos = lst->lptr[i]->pmap.begin(); pos != lst->lptr[i]->pmap.end(); ++pos ){
				pos->second->mas = Mass_H;
				pos->second->col = col;
				pos->second->vel.IN((double)gauss(T/pos->second->mas,VX),(double)gauss(T/pos->second->mas,VY));
				}
			}
			
			for( int i = 0; i < 8; i++){
				for( int j = 0; j < 8; j++){
					for( pos = lst->lptr[i]->lptr[j]->pmap.begin(); pos != lst->lptr[i]->lptr[j]->pmap.end(); ++pos ){
						pos->second->mas = Mass_H;
						pos->second->col = col;
						pos->second->vel.IN((double)gauss(T/pos->second->mas,VX),(double)gauss(T/pos->second->mas,VY));
					}
				}
			}
			
			for( int i = 0; i < 8; i++){
				for( int j = 0; j < 8; j++){
					for( int k = 0; k < 8; k++){
						for( pos = lst->lptr[i]->lptr[j]->lptr[k]->pmap.begin(); pos != lst->lptr[i]->lptr[j]->lptr[k]->pmap.end(); ++pos ){
						pos->second->mas = Mass_H;
						pos->second->col = col;
						pos->second->vel.IN((double)gauss(T/pos->second->mas,VX),(double)gauss(T/pos->second->mas,VY));
						}
					}
				}
			}
			
			
			col++;
		}
	}
	
	*/
	
	double mas=0.0;
	_vec<double> vel;
	for( int i = 0; i < ptcl_size; i++){
		mas += ptcl[i].mas;
		vel += ptcl[i].mas * ptcl[i].vel;
	}
	vel /= mas;
	
	for( int i = 0; i < ptcl_size; i++){
		ptcl[i].vel -= vel;
	}
	
//	cout << col << endl; exit(0);
//	cout << "LINE113-115" << "color=0 as a flow field" << endl;
//	cout << "total number of colors = " << col << endl;
	
}


void Propagation(double& SX)
{
	_vec<double> loc;
	
	for( int i = 0; i < ptcl_size; i++){
		loc =( ptcl[i].loc + ptcl[i].vel);
		
		
		(void)Periodic(loc.x,loc.y);
	//	(void)Wall(*(ptcl + i),loc.x,loc.y,SX);
		
		ptcl[i].loc = loc;
		
	}
}

void CellCalculation()
{
	int i;
	
	memset( cell, 0 ,cell_size * sizeof(_cell_base));
	
	for( int x = 0; x < X; x++)
	for( int y = 0; y < Y; y++){
		list[x][y].pmap.clear();
	}
	
	_vec<int> loc;
	for(i = 0; i < ptcl_size; i++){
		loc.IN(ptcl[i].loc.icast());
		if( loc.x < 0.0 || loc.x >= X || loc.y < 0.0 || loc.y >= Y ){
			cout << i << endl;
			cout << loc.x << " " << loc.y <<  endl;
			exit(0);
		}
		_cell_base *cellp;
		
		cellp = &cell[loc.x][loc.y];
		
		cellp->mas += ptcl[i].mas;
		cellp->vel += (ptcl[i].mas * ptcl[i].vel);
		cellp->ene += ptcl[i].mas *SQR( ptcl[i].vel);
		cellp->count++;
		list[loc.x][loc.y].pmap.insert(std::make_pair(ptcl[i].col,&ptcl[i]));
		//}
	}
	
	
	for(int x = 0; x < X; x++){
		for(int y = 0; y < Y; y++){
			if(cell[x][y].mas > 0.0 ){
				cell[x][y].vel /= cell[x][y].mas;
				cell[x][y].ene -= cell[x][y].mas * SQR( cell[x][y].vel );
			}
		}
	}
}


void RotationMatrix()
{
	int x,y;
	double th;
	
	for(x = 0; x < X; x++){
		for(y = 0; y < Y; y++){
			if( RAND < 0.50 ){
				cell[x][y].rot[0] = 0.0;
				cell[x][y].rot[1] = -1.0;
			}else{
				cell[x][y].rot[0] = 0.0;
				cell[x][y].rot[1] = 1.0;
			}
		}
	}
}



void MultiColorCollision(){
	
	
	for( int x = 0; x < X; x++)
	for( int y = 0; y < Y; y++){
		double inner = 0.0;
		double outer = 0.0;
		
		_cell_base* cellp = &cell[x][y];
		
		_list* lst = &list[x][y];
		
		//
		if( (int)( cellp->mas ) != 0 ){
		// �������p�O����
		int count=0;
		int TMCOL=0;
		
		while( count == 0 && TMCOL < SPCS ){
			count = (int)( lst->number_color(TMCOL) );
			TMCOL++;
		}
		
		if( cellp->count != count ){
		
		
		std::map< int,_vec<double> > CFLD;
		_vec<double> MFLD(0.0,0.0);
		_vec<double> WBC(0.0,0.0);
		
		/*
		if( y == 1 ){
			WBC.IN(0.0,-2.0*WALL_COLOR);
		}else if( y == (Y - 2) ){
			WBC.IN(0.0, 2.0*WALL_COLOR);
		}
		*/
		
		int pcol = -1;
		int col;
		std::multimap<int, _ptcl_base*>::iterator pos;
		for( pos = lst->lptr[0]->pmap.begin(); pos != lst->lptr[0]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].x +=  lst->lptr[0]->number_color(col);
				pcol = col;
			}
		}
		
		pcol = -1;
		for( pos = lst->lptr[2]->pmap.begin(); pos != lst->lptr[2]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].x -=  lst->lptr[2]->number_color(col);
				pcol = col;
			}
		}
		
		pcol = -1;
		for( pos = lst->lptr[1]->pmap.begin(); pos != lst->lptr[1]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].y +=  lst->lptr[1]->number_color(col);
				pcol = col;
			}
		}
		
		pcol = -1;
		for( pos = lst->lptr[3]->pmap.begin(); pos != lst->lptr[3]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].y -=  lst->lptr[3]->number_color(col);
				pcol = col;
			}
		}
		
		pcol = -1;
		for( pos = lst->lptr[4]->pmap.begin(); pos != lst->lptr[4]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].x += (lst->lptr[4]->number_color(col)*0.707106781186548);
				CFLD[col].y += (lst->lptr[4]->number_color(col)*0.707106781186548);
				pcol = col;
			}
		}
		
		pcol = -1;
		for( pos = lst->lptr[7]->pmap.begin(); pos != lst->lptr[7]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].x += (lst->lptr[7]->number_color(col)*0.707106781186548);
				CFLD[col].y -= (lst->lptr[7]->number_color(col)*0.707106781186548);
				pcol = col;
			}
		}
		
		pcol = -1;
		for( pos = lst->lptr[5]->pmap.begin(); pos != lst->lptr[5]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].x -= (lst->lptr[5]->number_color(col)*0.707106781186548);
				CFLD[col].y += (lst->lptr[5]->number_color(col)*0.707106781186548);
				pcol = col;
			}
		}
		
		pcol = -1;
		for( pos = lst->lptr[6]->pmap.begin(); pos != lst->lptr[6]->pmap.end(); ++pos){
			col = pos->second->col;
			if( col != pcol ){
				CFLD[col].x -= (lst->lptr[6]->number_color(col)*0.707106781186548);
				CFLD[col].y -= (lst->lptr[6]->number_color(col)*0.707106781186548);
				pcol = col;
			}
		}
		
		
		//
		for( int col = 0; col < SPCS; col++){
			MFLD += CFLD[col];
		}
		///
		
		
		_vec<double> SFLD(0.0,0.0);
		SFLD = MFLD - CFLD[0];
		
		for( int col = 0; col < SPCS; col++){
			if( lst->pmap.count(col) != 0 ){
				_vec<double> fld;
				_vec<double> flx;
				
				
				if( col == 0 ){
					fld =  A33 * CFLD[0] - A13 * SFLD - B3 * WBC;
				}else{
					fld =  A11 * CFLD[col] - A13 * CFLD[0] - A12 * ( SFLD - CFLD[col] ) - B1 * WBC;
				}
				
				
				multimap<int, _ptcl_base*>::iterator pos;
				for( pos = lst->pmap.lower_bound(col); pos != lst->pmap.upper_bound(col); ++pos){
					flx +=  ( pos->second->vel - cellp->vel );
					
					//�F�ꂩ��󂯂�͂́A�F�݂̂Ō��肳���ׂ��B
					//���ʂ́A�����^���̕ω��̂��ɂ����Ȃ̂ŁA�͂��󂯂���̉^����Ԃ����肷�邱�Ƃ͂����Ă�
					//�F��ƐF���q�̑��ݍ�p�̑傫�������߂Ă͂Ȃ�Ȃ��B
					//flx += qc * ( pos->second->vel - cellp->vel );
					
					
				}
				
				inner += flx * fld;
				outer += flx % fld;
				
			}
		}

		CFLD.clear();
		//�����˾��Ͳ�ž���٤η׻������롣
		}
		
		double th;
		int jiku_chk = 0;
		
		if( RAND < 0.50 ){ jiku_chk = 1; outer *= -1.0; }

		
		if( (inner * inner) > 0.00010 ){
			th = atan2(outer,inner);
			
			double alpha = cos(th);
			if( inner < 0.0 ){
				if( alpha > 0.0 ){ th += PI;}
			}else{
				if( alpha < 0.0 ) th += PI;
			}
			
			cellp->rot[0] = cos(th);
			cellp->rot[1] = sin(th);

		}else{
			cellp->rot[0] = 0.00;
			cellp->rot[1] = 1.00;
		}
		
		if( jiku_chk == 1 ) cellp->rot[1] *= -1.0;
	}
	
	}
}



void Stress_Measure(){
	
	for(int i = 0; i < ptcl_size; i++){
		_vec<int> loc;
		loc = ptcl[i].loc.icast();
		_node* nodep = &node[loc.x][loc.y];
		_cell_base* cellp = &cell[loc.x][loc.y];
		_ptcl_base* p = &ptcl[i];
		
		nodep->s[0] += p->mas * p->vel.x * p->vel.x;
		nodep->s[1] += p->mas * p->vel.x * p->vel.y;
		nodep->s[2] += p->mas * p->vel.y * p->vel.y;
		
		
	}
	
	/*
	
	for( int x = 0; x < X; x++){
		for( int y = 0; y < Y; y++){
			_cell_base* cellp = &cell[x][y];
			_node* nodep = &node[x][y];
			
			nodep->s[0] += cellp->mas * cellp->vel.x * cellp->vel.x;
			nodep->s[1] += cellp->mas * cellp->vel.x * cellp->vel.y;
			nodep->s[2] += cellp->mas * cellp->vel.y * cellp->vel.y;
			
			nodep->mas += cellp->mas;
			nodep->mom += (nodep->mas * cellp->vel);
		}
	}
	*/
	
}

void Mean(){
	
	for( int i=0; i < ptcl_size; i++){
		_vec<int> loc = ptcl[i].loc.icast();
		_node* nodep = &node[loc.x][loc.y];
		
		nodep->mas += ptcl[i].mas;
		nodep->mom += (ptcl[i].mas * ptcl[i].vel);
		if( ptcl[i].col != 0 ) nodep->col += 1.0;
	}
	
}

void AVR_output(int step, double& SX){
	
	char fna[30];
	sprintf(fna,"avr_U%05d_%07d.dat",(int)(10000.0*SX),step);
	
	FILE *fptr = fopen(fna,"w");
	fprintf(fptr,"variables = \"x\",\"y\",\"u\",\"v\",\"m\", \"c\", \"sxx\", \"sxy\", \"syy\"\n");
	fprintf(fptr,"zone i=%d j=%d f=point\n",X,Y);
	
	
	for(int y = 0; y < Y; y++)
	for(int x = 0; x < X; x++){
		_node* nodep = &node[x][y];
		if( nodep->mas > 0 ){ 
			nodep->mom /= nodep->mas;
			nodep->mas /= (double)MTIME;
			nodep->col /= (double)MTIME;
			
			nodep->s[0] /= (double)MTIME;
			nodep->s[1] /= (double)MTIME;
			nodep->s[2] /= (double)MTIME;
			
			
			nodep->s[0] -= nodep->mas * SQR( nodep->mom.x );
			nodep->s[1] -= nodep->mas * nodep->mom.x * nodep->mom.y;
			nodep->s[2] -= nodep->mas * SQR( nodep->mom.y );
			
		}
		fprintf(fptr,"%f %f %f %f %f %f %f %f %f\n",x+0.50,y+0.50,nodep->mom.x,nodep->mom.y, nodep->mas, nodep->col,nodep->s[0], nodep->s[1], nodep->s[2]);
	}
	fclose(fptr);
	
	
	double sxx = 0.0;
	double sxy = 0.0;
	double syy = 0.0;
	
	for(int x = 0; x < X; x++)
	for(int y = 6; y < (Y-6); y++){
		_node* nodep = &node[x][y];
		sxx += nodep->s[0];
		sxy += nodep->s[1];
		syy += nodep->s[2];
	}
	
	double SAMPLE = (double)(X * (Y - 12));
	sxx /= SAMPLE;
	sxy /= SAMPLE;
	syy /= SAMPLE;
	
	char fnaa[30];
	FILE *fptra;
	sprintf(fnaa,"shear_strain.dat");
	fptra = fopen(fnaa,"a");
	fprintf(fptra,"%f %f %f %f\n",SX, sxx, sxy, syy);
	fclose(fptra);
	memset( node, 0, cell_size * sizeof(_node));
}


void Collision()
{
//	(void)RotationMatrix();
	
	_vec<int> loc;
	_vec<double> vel;
	_vec<double> rot;
	
	for(int i = 0; i < ptcl_size; i++){
		loc = ptcl[i].loc.icast();
		vel = ptcl[i].vel - cell[loc.x][loc.y].vel;
		rot.x = (vel.x * cell[loc.x][loc.y].rot[0] - vel.y * cell[loc.x][loc.y].rot[1]);
		rot.y = (vel.x * cell[loc.x][loc.y].rot[1] + vel.y * cell[loc.x][loc.y].rot[0]);
		
		ptcl[i].vel = rot + cell[loc.x][loc.y].vel;
	//	if( ptcl[i].col == 0 ) ptcl[i].vel.x += G;
  }
}


void SystemMonitor(int step)
{
  int i;
  _total total;
  FILE *fptr;
  char fna[30];
  
  sprintf(fna,"system.dat");
  fptr = fopen(fna,"a");
  memset(&total, 0, sizeof(_total));

	
	for( i = 0; i < ptcl_size; i++){
		total.mas += ptcl[i].mas;
		total.mom.x += ptcl[i].mas * ptcl[i].vel.x;
		total.mom.y += ptcl[i].mas * ptcl[i].vel.y;
		total.ene += ptcl[i].mas * (SQR( ptcl[i].vel.x ) + SQR( ptcl[i].vel.y ));
	}
  total.mom.x /= total.mas;
  total.mom.y /= total.mas;
  
  total.ene -= total.mas * (SQR( total.mom.x ) +SQR( total.mom.y ));
	
	
  fprintf(fptr, "%d %f %f %f %f\n",
	  /* ����     */ step,
	  /* ������   */ total.mas,
	  /* ʿ�Ѳ��� */ total.ene/(2.0*(double)(ptcl_size)),
	  /* ʿ��ή® */ total.mom.x, total.mom.y);
  fclose(fptr);

}


void Output(int step,double SX){
	
	char fna[30];
	sprintf(fna,"d%07d.tec",step);
	FILE *fptr = fopen(fna,"w");
	fprintf(fptr,"variables = \"x\",\"y\",\"u\",\"v\",\"c\"\n");
	fprintf(fptr,"zone i=%d j=%d f=point\n",X,Y);
	
	float col=0.0;
	
	for(int y = 0; y < Y; y++)
	for(int x = 0; x < X; x++){
		/*col = 1.0;
		if( (int)cell[x][y].mas > 0 ){
			float nc_0 = list[x][y].number_color(0);
			nc_0 /= (float)cell[x][y].count;
			col -= nc_0;
		}else col = 0.0;
		*/
		
		//start
		
		col = 0.0;
		if( (int)cell[x][y].mas > 0 ){ 
			for(int i = 0; i < SPCS; i++)
					col += (float)( i * list[x][y].number_color(i) );
			col /= (float)cell[x][y].count;
		}else col = 0.0;
		
		//
		//end
		//
		//if( col > 100 ){cerr<< col <<" "<<list[x][y].pmap.count(0) << " " << list[x][y].pmap.count(1) <<" "<< cell[x][y].mas << endl;exit(0);}
		fprintf(fptr,"%d %d %f %f %f\n",x,y,cell[x][y].vel.x,cell[x][y].vel.y,col);
	}
	fclose(fptr);
}

void Ptcl_Output(){
	FILE *fptr;
	char fna[30];
	
	sprintf(fna,"ptcl.dat");
	
	fptr = fopen(fna,"w");
	
	for( int i = 0; i < ptcl_size; i++){
		_ptcl_base* p = &ptcl[i];
		fprintf(fptr,"%d %3.13f %3.13f %3.13f %3.13f %3.13f\n",p->col, p->mas, p->loc.x, p->loc.y, p->vel.x, p->vel.y);
	}
	fclose(fptr);
}

