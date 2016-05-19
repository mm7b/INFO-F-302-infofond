CXX = g++-4.9 
CXXFLAGS = -Wall -g
LDFLAGS = -lm
TARGET = solver

${TARGET}: Solver.o OrthogonalPackingSolver.o Main.o
	${CXX} -o $@ $^ ${LDFLAGS}

clean:
	rm -f *.o ${TARGET}


