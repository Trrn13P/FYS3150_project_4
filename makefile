CPPflags=  clang++ -Xpreprocessor -fopenmp -std=c++14
LIB = -larmadillo -llapack -lblas
#CPPflags= g++-10

all: compile execute

compile:
	${CPPflags} main.cpp MC_solver.cpp -o ./main.out ${LIB}

execute:
	./main.out
