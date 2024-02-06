CC = g++
CPPFLAGS = -g -fopenmp -lpthread -Wall -O2
LDFLAGS = 

OBJ = mainLoadFromComsol.o

%.o: %.cpp
		$(CC) $(CPPFLAGS) -c $^ -o $@

main: $(OBJ)
		$(CC) $(CPPFLAGS) $(OBJ) -o load
		
clean:
		rm -rf $(OBJ) load

.PHONY: clean