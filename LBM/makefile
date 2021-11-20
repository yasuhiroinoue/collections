all :	a.out
a.out : main.o lbm.o bc.o mymath.o
	g++ -lm main.o lbm.o bc.o mymath.o 
lbm.o:	lbm.cpp  lbm.h vec.h  lbm_class.h system.h defs.h mymath.h
	g++ -c -O3 lbm.cpp 
mymath.o : mymath.cpp mymath.h
	g++ -c  -O3 mymath.cpp
bc.o : bc.cpp vec.h  lbm_class.h system.h
	g++ -c  -O3 bc.cpp 
main.o : main.cpp defs.h system.h lbm_class.h vec.h
	g++ -c  -O3 main.cpp

clean : 
	rm *.o
dclean :
	rm *.dat
