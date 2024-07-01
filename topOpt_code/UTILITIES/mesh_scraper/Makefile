CC = g++
CPPFLAGS = -g -Wall -fopenmp -lpthread -O2
LDFLAGS =

OBJ = scrape_meshes.o scrape_main.o

%.o: %.cpp
		$(CC) $(CPPFLAGS) -c $^ -o $@

main: $(OBJ)
		$(CC) $(CPPFLAGS) $(OBJ) -L/usr/lib/gcc/x86_64-linux-gnu/11/ -lgfortran -llapack -lblas -o main  
		
clean:
		rm -rf $(OBJ) main

.PHONY: clean