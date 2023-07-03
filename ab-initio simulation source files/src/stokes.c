#include <stdio.h>
#include <petscksp.h>
#include <petscdmda.h>
#include <math.h>
#include "stokes.h"
#include <getopt.h>

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif


levelset_vec create_levelset(parameter *para) {
  levelset_vec X;
  PetscInt nr, nz;
  DM da;
  
  nr = round( para->maxr/para->dr );
  nz = round( para->maxz/para->dz );
  DMDACreate2d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_BOX,
               nr, nz,
               1, PETSC_DECIDE, 1,
               para->nghostlayer, NULL, NULL, &da);
  
  X.nr = nr;
  X.nz = nz;
  X.da = da;
  X.t = 0;

  DMSetFromOptions(X.da);
  DMSetUp(X.da);
  DMCreateGlobalVector(X.da,&X.data);
  VecSet(X.data,0.0);
  
  return X;
}

void destroy_levelset(levelset_vec X) {
  free(X.data);
}


stokes_matrix create_stokes_matrix(parameter *para) {
  
  stokes_matrix B;
  DM da;
  int dof = 3;
  
  DMDACreate2d(PETSC_COMM_WORLD, 
	       DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
	       DMDA_STENCIL_BOX,
	       para->nr, para->nz,
	       1, PETSC_DECIDE, dof,
	       1, NULL, NULL, &da);
  
  B.nz = para->nz;
  B.nr = para->nr;
  B.da = da;

  return B;
}

void assembly_stokes_matrix(stokes_matrix *B) {

  DMSetMatrixPreallocateOnly(B->da, PETSC_TRUE);
  DMSetUp(B->da);
  //  DMCreateMatrix(B->da, MATAIJ, &B->data);
  DMCreateMatrix(B->da, &B->data);
  
  //MatZeroEntries(B->data);

}


stokes_force create_stokes_force(parameter *para, DM da) {
  stokes_force X;

  X.nr = para->nr;
  X.nz = para->nz;
  X.da = da;

  DMCreateGlobalVector(X.da, &X.data);
  VecSet(X.data, 0.1);
  
  return X;
  
}


/*******************************************************/

void initial_levelset(levelset_vec *G, parameter *para) {
  
  PetscInt llr, llz, lsizer, lsizez, rank;

  DMDAGetCorners(G->da, &llr, &llz, 0, &lsizer, &lsizez, 0);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscScalar interfacei, interfaceR;
  PetscScalar z, lambdaz;
  
  interfacei = para->r0/para->dr;
  
  PetscInt i, j, j0;
  PetscScalar number;
  PetscScalar **array;
  DMDAVecGetArray(G->da, G->data, &array);
  
  if( strcmp(para->mode,"periodic") == 0 ) {
    lambdaz = para->maxz - para->dz/2;
    srandom(rank);
    
    for(j=llz; j<llz+lsizez; j++) {
      z = j * para->dz;
      interfaceR = interfacei + (para->pertb/para->dr) * cos(z/lambdaz*2*M_PI);
      for(i=llr; i<llr+lsizer; i++) {
	array[j][i] = (i-interfaceR) * para->dr;
      }
    }
    DMDAVecRestoreArray(G->da, G->data, &array);
        
    return;
  }
  else if( strcmp(para->mode,"end") == 0 ) {
    j0 = G->nz - (int) (para->init_thread/para->dz);
    srandom(rank);
    //    printf("j0 = %d\n",j0);
    for( j=llz; j<llz+lsizez; j++ ) {
      if(j<j0) {
	for( i=llr; i<llr+lsizer; i++) {
	  array[j][i] = sqrt( (i)*(i)*para->dr*para->dr + (j-j0)*(j-j0)*para->dz*para->dz ) - para->r0;
	}
      }
      else {
	number = ((double) random() / (RAND_MAX)) - 0.5;
	for( i=llr; i<llr+lsizer; i++) {
	  interfaceR = interfacei + number * para->pertb;
	  array[j][i] = (i - interfaceR) * para->dr;
	}
      }
    }
    DMDAVecRestoreArray(G->da, G->data, &array);
    
    return;
  }
  else {
    if(rank==0) printf("Please check the initial condition configuration: periodic or end \n");
    DMDAVecRestoreArray(G->da, G->data, &array);
    
    return;
  }
  
  DMDAVecRestoreArray(G->da, G->data, &array);
  
  
  return;
}

/******* viscosity_mu2 function *******/
PetscScalar viscosity_mu2(parameter *para, PetscInt i) {
    PetscScalar zposi;
    zposi = i*para->dz;
    //    return func_1(zposi, para);
}

PetscScalar func_1(PetscScalar zposi, parameter *para) {
    PetscScalar muo1 = 1000;
    PetscScalar muo2 = 10;
    PetscScalar mu2_val;
    mu2_val = 10 + zposi/(para->maxz-para->solidwidth) * (1000-10);
    return mu2_val;
}

/****************************************************************************/
/********************** get the temperature *******************************/
void get_temperature(PetscScalar *t_z, PetscScalar zposi, parameter *para) {
  PetscScalar dT;
  PetscScalar tsi;
  PetscScalar troom;
  /* PetscScalar shift = para->shift;*/
  PetscScalar maxz = para->maxz;
  PetscScalar tsolid = para->tsolid;
  PetscScalar thigh = para->thigh;
  PetscScalar tlow = para->tlow;

  PetscScalar solidwidth = para->solidwidth;
  PetscScalar twidth1 = para->twidth1;
  PetscScalar twidth2 = para->twidth2;
  PetscScalar thigh_width = para->thigh_width;
  PetscScalar z2 = para->z2;
  PetscScalar E = para->E;
  
  tsi = 1414;
  troom = 20;
  dT =(thigh - tsi);
  /*  970 means 2mm */
  /* shift = arctanh((2*tsi - thigh - tlow)/(thigh-tlow));*/
  
  if (para->tanh_profile == 1) {
    if (zposi > (maxz - solidwidth)) *t_z = tsolid;
    else
      /* *t_z = tlow + dT * tanh((para->maxz-solidwidth-zposi+para->startx)/twidth1);*/
     *t_z = thigh - dT*exp(-(para->maxz - zposi + para->startx)/twidth1);
      /* *t_z = thigh - dT*exp(-(para->maxz - zposi + para->startx + z2)/twidth1) - E*exp(-(para->maxz - zposi + para->startx + z2)/twidth2); */
    return;
  }
 
  if (zposi > (maxz - solidwidth))
    *t_z = tsolid;
  else if (zposi > (maxz - solidwidth - twidth1))
    *t_z = tlow + dT/(-twidth1) * (zposi - (maxz-solidwidth));
  else if (zposi > (maxz - solidwidth - twidth1 - thigh_width))
    *t_z = thigh;
  else if (zposi > (maxz - solidwidth - twidth1 - thigh_width - twidth2))
    *t_z = thigh + dT/(twidth2) * (zposi - (maxz - solidwidth - twidth1 - thigh_width));
  else
    *t_z = tlow;

  return;

}

/****************************************************************************/

void get_input(int argc, char **args, viscosity *mu, parameter *para) {
  
  int opt;
  
  struct option opts[] = {
    { "maxr",          1, NULL, 1 },
    { "maxz",          1, NULL, 2 },
    { "dr",            1, NULL, 3 },
    { "dz",            1, NULL, 4 },
    { "r0",            1, NULL, 5 },
    { "tension",       1, NULL, 6 },
    { "pertb",         1, NULL, 7 },
    { "period_or_end", 1, NULL, 8 },
    { "mui",           1, NULL, 9 },
    { "muo",           1, NULL, 10 },
    { "vf",            1, NULL, 11 },
    { "temp_profile",  1, NULL, 12 },
    { "tlow",          1, NULL, 13 },
    { "thigh",         1, NULL, 14 },
    { "twidth",        1, NULL, 15 },
    { "lowtwidth",     1, NULL, 16 },
    { "restart",       1, NULL, 17 },
    { "trestart",      1, NULL, 18 },
    { "outputdt",      1, NULL, 19 },
    { "nghostlayer",   1, NULL, 20 },
    { "epsilon",       1, NULL, 21 },
    { "reinitstep",    1, NULL, 22 },
    { "tsolid",        1, NULL, 23 },
    { "solidwidth",    1, NULL, 24 },
    { "twidth1",       1, NULL, 25 },
    { "twidth2",       1, NULL, 26 },
    { "thigh_width",   1, NULL, 27 },
    { "init_thread",   1, NULL, 28 },
    { "tanh_profile",  1, NULL, 29 },
    { "startx",        1, NULL, 30 },
    { NULL,   0, NULL, 0 }
  };

  while ((opt = getopt_long(argc, args, "h", opts, NULL)) != -1) {
    switch (opt) {
    case 1:
      para->maxr = atof(optarg);
      break;
    case 2:
      para->maxz = atof(optarg);
      break;
    case 3:
      para->dr = atof(optarg);
      break;
    case 4:
      para->dz = atof(optarg);
      break;
    case 5:
      para->r0 = atof(optarg);
      break;
    case 6:
      para->tension = atof(optarg);
      break;
    case 7:
      para->pertb = atof(optarg);
      break;
    case 8:
      para->period_or_end = atof(optarg);
      break;
    case 9:
      para->mui = atof(optarg);
      break;
    case 10:
      para->muo = atof(optarg);
      break;
    case 11:
      para->vf = atof(optarg);
      break;
    case 12:
      para->temp_profile = atof(optarg);
      break;
    case 13:
      para->tlow = atof(optarg);
      break;
    case 14:
      para->thigh = atof(optarg);
      break;
    case 15:
      para->twidth = atof(optarg);
      break;
    case 16:
      para->lowtwidth  = atof(optarg);
      break;
    case 17:
      para->restart = atof(optarg);
      break;
    case 18:
      para->trestart = atof(optarg);
      break;
    case 19:
      para->outputdt = atof(optarg);
      break;
    case 20:
      para->nghostlayer = atof(optarg);
      break;
    case 21:
      para->epsilon = atof(optarg);
      break;
    case 22:
      para->reinitstep = atof(optarg);
      break;
    case 23:
      para->tsolid = atof(optarg);
      break;
    case 24:
      para->solidwidth = atof(optarg);
      break;
    case 25:
      para->twidth1 = atof(optarg);
      break;
    case 26:
      para->twidth2 = atof(optarg);
      break;
    case 27:
      para->thigh_width = atof(optarg);
      break;
    case 28:
      para->init_thread = atof(optarg);
      break;
    case 29:
      para->tanh_profile = atof(optarg);
      break;
    case 30:
      para->startx = atof(optarg);
      break;
    }
  }
  
  para->nr = (int)(para->maxr/para->dr);
  para->nz = (int)(para->maxz/para->dz);

  if(para->period_or_end) strcpy(para->mode, "periodic");
  else strcpy(para->mode, "end");

  int i, nz;
  PetscScalar t_z, z_width;
  nz = para->nz;
  z_width = para->twidth / para->dz;
  mu->mu1 = malloc( nz * sizeof(PetscScalar) );
  mu->mu2 = malloc( nz * sizeof(PetscScalar) );

  PetscScalar zposi;

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  for(i=0; i<nz; i++) {
    zposi = i * para->dz;
    if(para->temp_profile == 1) {
      get_temperature(&t_z, zposi, para);
      //t_z = (para->thigh - para->tlow) * tanh( (nz + 1.0*para->lowtwidth/para->dz - i) / z_width ) + para->tlow;
      mu->mu2[i] = pow(10, 26909.0/(t_z+273)-7.2348);
      mu->mu2[i] = mu->mu2[i]/1.0e3; // outer
      mu->mu1[i] = pow(10, 819.0/(t_z+273)-3.727);
      mu->mu1[i] /= 1.0e3; // inner
    }
    else {
        if (i >= para->nz - para->solidwidth/para->dz) mu->mu2[i] = 1.0e20;
        else {
            mu->mu2[i] = para->muo;
            //mu->mu2[i] = func_1(zposi, para);
            mu->mu1[i] = para->mui;
        }
        //       
        //mu->mu2[i] = para->muo;
    }
    
    if(rank==0 && i%(nz/100)==0) printf("z = %f, T = %f, mu2 = %f\n", zposi, t_z, mu->mu2[i]);
  }
  
  if(para->temp_profile == 1) para->tension = para->tension * 1000;
   
  if (rank==0)
    printf("\
running simulation with \n \
maxr = %g, maxz = %g, dr=%g, dz = %g, nr = %d, nz = %d, \n \
r0 = %g, \n \
tension = %g, \n \
pertb = %g, \n \
mode = %s, \n \
mu1 = %g, mu2 = %g, \n \
vf = %g, \n\
temp_profile = %d, \n\
tlow = %g, thigh = %g, tsolid = %g \n\
solidwidth = %g, twidth1 = %g, thigh_width = %g, twidth2 = %g, \n\
restart = %g, trestart = %g, \n\
outputdt = %g, reinitstep = %d, \n\
tanh_profile = %d, startx = %g \n\
\n",
           para->maxr, para->maxz, para->dr, para->dz, para->nr, para->nz,
           para->r0,
           para->tension,
           para->pertb,
           para->mode,
           mu->mu1[0], mu->mu2[0],
           para->vf,
           para->temp_profile,
           para->tlow, para->thigh, para->tsolid,
           para->solidwidth, para->twidth1, para->thigh_width, para->twidth2,
           para->restart, para->trestart,
           para->outputdt, para->reinitstep,
	   para->tanh_profile, para->startx);
  
  return;
}

