all :	a.out
a.out : main.o rlg.o bc.o mymath.o
	g++	 main.o rlg.o bc.o mymath.o 
rlg.o:	rlg.cpp  vec.h rlg.h rclass.h system.h defs.h mymath.h
	g++ -c rlg.cpp  -ffast-math -O3
mymath.o : mymath.cpp mymath.h
	g++ -c mymath.cpp  -ffast-math -O3
bc.o : bc.cpp vec.h  rclass.h system.h
	g++ -c bc.cpp -ffast-math -O3
main.o : main.cpp defs.h system.h rclass.h vec.h
	g++ -c main.cpp -ffast-math -O3

clean : 
	rm *.o
dclean :
	rm *.dat
