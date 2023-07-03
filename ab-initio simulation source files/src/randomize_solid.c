#include <stdio.h>
#include <petscksp.h>
#include <petscdmda.h>
#include <math.h>
#include "stokes.h"
#include <getopt.h>

void randomize_solid(levelset_vec *G, parameter *para) {

    PetscInt llr, llz, lsizer, lsizez, rank;
    DMDAGetCorners(G->da, &llr, &llz, 0, &lsizer, &lsizez, 0);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscScalar **array;
    DMDAVecGetArray(G->da, G->data, &array);

    PetscInt i, j;
    PetscScalar z, rnum;
    for (j=llz; j<llz+lsizez; j++) {
        z = j * para->dz;
        if (z > para->maxz-para->solidwidth/2) {
            for (i=llr; i<llr+lsizer; i++) {
                rnum = ((double) random() / (RAND_MAX)) - 0.5;
                array[j][i] += rnum * para->pertb;
            }

        }
    }

    DMDAVecRestoreArray(G->da, G->data, &array);

    return;
    
}
