G++ = g++-4.9 -Wall  -g

sudosat: Main.o Solver.o
	$(G++) *.o -o orthogonal_packing

Solver.o: Solver.cpp Solver.hpp SolverTypes.hpp Global.hpp Heap.hpp
	$(G++) -c Solver.cpp 

Main.o: Main.cpp  Solver.hpp OrthogonalPackingSolver.hpp
	$(G++) -c Main.cpp

clean:
	$(RM) *.o orthogonal_packing
