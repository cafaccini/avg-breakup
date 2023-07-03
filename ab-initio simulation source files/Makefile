all: stokes_test

#include ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables
#include ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscrules

# On Quartz, the lib directory is not under the ${PETSC_ARCH} directory
include ${PETSC_DIR}/lib/petsc/conf/petscvariables
include ${PETSC_DIR}/lib/petsc/conf/petscrules

CC = cc
CFLAGS = -g 

#CLEANFILES = stokes
SOURCE = main_test.c $(wildcard ./src/*.c)
OBJS = $(SOURCE:.c=.o)

stokes_test: $(OBJS)
	$(CC) -o $@ $^ ${PETSC_LIB}

clean:
	rm -f stokes_test *.o src/*.o
