CC = g++
CPPFLAGS = -g -Wall -fopenmp -lpthread -O2
LDFLAGS =

OBJ = main.o TopOpt.o nsProblem.o assemble.o adjointNS.o assembleAdjoint.o  OPTIMIZATION/Optimizer.o OPTIMIZATION/Mean_Diffusion_Filter.o OPTIMIZATION/Mean_Diffusion_Filter_Node.o MATH_C++/parallel.o

%.o: %.cpp
		$(CC) $(CPPFLAGS) -c $^ -o $@

main: $(OBJ)
		$(CC) $(CPPFLAGS) $(OBJ) -L/usr/lib/gcc/x86_64-linux-gnu/11/  PARDISO/libpardiso600-GNU800-X86-64.so -lgfortran -llapack -lblas -o main  
		
clean:
		rm -rf $(OBJ) main

.PHONY: clean