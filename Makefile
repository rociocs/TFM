COMP=gcc
LIBS=-lm -llapack -lblas -lpthread
OMP=-fopenmp -liomp5

all: Sequential_LU OpenMP_LU

Sequential_LU: Sequential_LU.c 
	$(COMP) $(LIBS) Sequential_LU.c -o Sequential_LU.o

OpenMP_LU: OpenMP_LU.c
	$(COMP) $(LIBS) $(OMP) OpenMP_LU.c -o OpenMP_LU.o

clean:
	rm -f *.o

