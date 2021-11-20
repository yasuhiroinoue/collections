//mylib.h
/*
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#ifndef MYLIB_H 
#define MYLIB_H 

//�Q�����s��Ɋւ����`
//typedef double** MATRIX;
//MATRIX matrix_new(int, int);
//void matrix_delete(MATRIX);

//���w�֐�
static double HPI = 1.57079632679490;
static double PI =  3.14159265358979;
static double PI2 = 6.28318530717959;
#define TANE 10

#define SQR(x) (x * x)
#define ABS(x) sqrt(SQR(x))

/* ran2()�Ɠ����̐��x������ran2()�̂Q�{���� Mersenne Twister */
/* http://www.math.keio.ac.jp/~matumoto/mt.html */
/* A C-program for MT19937: Real number version([0,1)-interval) */
/* (1999/10/28)                                                 */
/*   genrand() generates one pseudorandom real number (double)  */
/* which is uniformly distributed on [0,1)-interval, for each   */
/* call. sgenrand(seed) sets initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be       */
/* called once. (seed is any 32-bit integer.)                   */
/* Integer generator is obtained by modifying two lines.        */
/*   Coded by Takuji Nishimura, considering the suggestions by  */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.            */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura. */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

/* Period parameters */  
#define MTN 624
#define MTM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[MTN]; /* the array for the state vector  */
static int mti=MTN+1; /* mti==N+1 means mt[N] is not initialized */

double genrand(void);
void sgenrand( unsigned long seed );
double gauss( double t, double u);
double gamma_noflow( double t);
#define RAND genrand()


#endif
