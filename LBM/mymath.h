/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
//mylib.h



#ifndef MYMATH_H 

#define MYMATH_H 



//‚QŽŸŒ³s—ñ‚ÉŠÖ‚·‚é’è‹`

//typedef double** MATRIX;

//MATRIX matrix_new(int, int);

//void matrix_delete(MATRIX);



//”ŠwŠÖ”

static double HPI = 1.57079632679490;

static double PI =  3.14159265358979;

static double PI2 = 6.28318530717959;

#define TANE 10



#define SQR(x) (x * x)

#define ABS(x) sqrt(SQR(x))



/* Period parameters */  

#define N_MT 624

#define M_MT 397

#define MATRIX_A 0x9908b0dfUL   /* constant vector a */

#define UMASK 0x80000000UL /* most significant w-r bits */

#define LMASK 0x7fffffffUL /* least significant r bits */

#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )

#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))



static unsigned long state[N_MT]; /* the array for the state vector  */

static int left = 1;

static int initf = 0;

static unsigned long *next;





double genrand_real3(void);

void init_genrand(unsigned long s);

double gauss( double t, double u);

double gamma_noflow( double t);

#define RAND (double)genrand_real3()





#endif

