G++ = g++-4.9 -Wall  -g

sudosat: Main.o Solver.o
	$(G++) *.o -o path

Solver.o: Solver.cpp Solver.hpp SolverTypes.hpp Global.hpp Heap.hpp
	$(G++) -c Solver.cpp 

Main.o: Main.cpp  Solver.hpp 
	$(G++) -c Main.cpp

clean:
	$(RM) *.o path
