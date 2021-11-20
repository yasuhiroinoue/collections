/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#include <math.h>
#include "system.h"

void Periodic(int &x, int &y){
	if( x < 0 ) x += X; if( x >= X) x -= X;
	if( y < 0 ) y += Y; if( y >= Y) y -= Y;
}
