#include <petscksp.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <math.h>
#include <time.h>
#include "./src/stokes.h"
#include <stdio.h>
#include <petscviewerhdf5.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

int main(int argc, char **args) {

  levelset_vec G;

  PetscInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  viscosity *mu;
  parameter *para;
  mu = malloc(sizeof(viscosity));
  para = malloc(sizeof(parameter));

  get_input(argc, args, mu, para);

  G = create_levelset(para);
  PetscObjectSetName((PetscObject)G.data, "levelset");

  int output = 1;
  if(output) mk_dir(para);

  if(para->restart == 0) {
    initial_levelset(&G, para);

    if(output) output_levelset(&G, para);

    reinit(&G, para);

    para->trestart = 0;
  }
  else {
    load_levelset(&G, para);
    G.t = para->trestart;
  }
  
  stokes_matrix B;
  B = create_stokes_matrix(para);  
  assembly_stokes_matrix(&B);
 
  stokes_force force, uwp;
  force = create_stokes_force(para, B.da);  
  uwp = create_stokes_force(para, B.da);

  MYKSP solver;
  initialize_MYKSP(&solver);

  PetscScalar dt, t1, t2, t_end;
  PetscInt output_count = 1, evolve_count = 0;
  t1 = para->trestart;
  t2 = t1;
  t_end = para->maxz / para->vf * 4;
  // t_end = 0.2;
  
  double solve_t1, solve_t2;

  PetscViewer viewer;
  
  while (t2 < t_end) {
    MPI_Barrier( MPI_COMM_WORLD );
    
    if(rank==0)
      solve_t1 = MPI_Wtime();

    evolve_count++;

    get_B_from_G(&B, &G, para, mu);
    get_force_from_G(&force, &G, para);

    stokes_solve(&uwp, 
		 &solver, 
		 &B, &force, 
		 &G, 
		 para, mu);
    
    lab_frame_shift_w(&uwp, para->vf);
    
    dt = advection_timestep(&uwp, para);
    //    dt = 0.001;
    
    advection_evolve(&G, &uwp, para, dt);
    
    if( evolve_count % 5 == 0) 
      reinit(&G, para);

    t2 = t1 + dt;
    G.t = t2;
    uwp.t = t2;
    force.t = t2;

    // output_uwp(&uwp, para);
    if( t1 < output_count * para->outputdt + para->trestart
	&& t2 >= output_count * para->outputdt + para->trestart ) {
      G.t = output_count * para->outputdt + para->trestart;
      uwp.t = G.t;
      force.t = G.t;
      
      if(output) {
          output_levelset(&G, para);
      }

          
      // output_force(&force, para);
      MPI_Barrier( MPI_COMM_WORLD );
      /*
      if(rank==0) {
        printf("t1 = %g, t2 = %g, dt = %g\n", t1, t2, dt);
        solve_t1 = MPI_Wtime() - solve_t1;
        printf("solving time = %g\n", solve_t1);
      }
      */
      output_count++;
    }

    if(rank==0) {
      printf("t1 = %g, t2 = %g, dt = %g\n", t1, t2, dt);
      solve_t1 = MPI_Wtime() - solve_t1;
      printf("solving time = %g\n", solve_t1);
    }
      
    t1 = t2;
  }
  
  /** free memory **/
  free(mu->mu1);
  free(mu->mu2);
  free(mu);
  free(para);

  VecDestroy(&G.data);
  VecDestroy(&force.data);
  VecDestroy(&uwp.data);
  MatDestroy(&B.data);
  KSPDestroy(&solver.ksp);

  PetscFinalize();
  return 0;
}
