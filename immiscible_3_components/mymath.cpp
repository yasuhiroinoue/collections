//mylib.cpp
#include <math.h>
#include "mymath.h"

//ÇQéüå≥çsóÒÇ…ä÷Ç∑ÇÈíËã`
//MATRIX matrix_new(int col, int row){
//	MATRIX a = new double *[col];
//	for( int i = 0; i < row; i++)
//			a[i] = new double[row];
//
//	return a;
//}
//
//void matrix_delete(MATRIX a){
//	MATRIX b = a;
//	
//	while( *b != 0 ){
//		delete [] *b++;
//	}
//	delete [] a;
//}

//êîäwä÷êî
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

double gauss(double t, double u){
  double ran, gatai;
  
	do ran = RAND; while(ran == 0.0);
	gatai = (double)(sqrt(-2.0 * t * log(ran)) * sin(PI2 * RAND)) + u;
	return(gatai);
}

double gamma_noflow(double t)
{
  double ran, gamma;
  
  do ran = RAND; while(ran == 0.0);
  gamma = sqrt(-2.0 * t * log(ran));
  return(gamma);
}

void sgenrand(unsigned long seed)
{
    int i;

    for (i=0;i<MTN;i++) {
         mt[i] = seed & 0xffff0000;
         seed = 69069 * seed + 1;
         mt[i] |= (seed & 0xffff0000) >> 16;
         seed = 69069 * seed + 1;
    }
    mti = MTN;
}

/* Initialization by "sgenrand()" is an example. Theoretically,      */
/* there are 2^19937-1 possible states as an intial state.           */
/* This function allows to choose any of 2^19937-1 ones.             */
/* Essential bits in "seed_array[]" is following 19937 bits:         */
/*  (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]. */
/* (seed_array[0]&LOWER_MASK) is discarded.                          */ 
/* Theoretically,                                                    */
/*  (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]  */
/* can take any values except all zeros.                             */

double genrand(){
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= MTN) { /* generate N words at one time */
        int kk;

        if (mti == MTN+1)   /* if sgenrand() has not been called, */
	  sgenrand(4357); /* a default initial seed is used   */

        for (kk=0;kk<MTN-MTM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MTM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<MTN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MTM-MTN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[MTN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MTN-1] = mt[MTM-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y * 2.3283064365386963e-10 ); /* reals: [0,1)-interval */
    /* return y; */ /* for integer generation */
}
