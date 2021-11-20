/* 
Yasuhiro Inoue
inoue.yasuhiro.4n@kyoto-u.ac.jp
*/
#include <stdio.h>
#include "defs.h"
#include "system.h"

int main(){
	(void)Init();
	for( int step = 0; step <= STEP_END; step++){
		(void)CellCalc();
		(void)Stress_Measure();
		(void)Collision();

		if( step % STEP_OUT == 0){
			(void)Output(step);
		}

		(void)Streaming();
	}
}
