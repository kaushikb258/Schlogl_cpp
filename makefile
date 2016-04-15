all:
	g++ main.cpp func.cpp vtk_cpp.cpp -o schlogl.exe
clean:
	rm *~ *.exe *.o *.vtk
