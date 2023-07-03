#include <petscksp.h>
#include <sys/stat.h>
#include <string.h>
#include "stokes.h"

void mk_dir(parameter *para){
  char buffer[1000];
  
  sprintf(buffer, "./data");
  mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  if (para->temp_profile==1)
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth1_%.1f_thigh_width_%.1f_twidth2_%.1f_reinit_%d", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->tlow, para->thigh, para->twidth1, para->thigh_width, para->twidth2, para->reinitstep);
  else
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_mui_%.7f_muo_%.7f_reinit_%d", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->mui, para->muo, para->reinitstep);

  mkdir(buffer, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  return;
}


void output_levelset(levelset_vec *G, parameter *para) {
  
  PetscViewer viewer;
  char buffer[1000];
  
  if (para->temp_profile==1)
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth1_%.1f_thigh_width_%.1f_twidth2_%.1f_reinit_%d/outputG_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->tlow, para->thigh, para->twidth1, para->thigh_width, para->twidth2, para->reinitstep, G->t);
  else
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_mui_%.7f_muo_%.7f_reinit_%d/outputG_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->mui, para->muo, para->reinitstep, G->t);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(G->data, viewer);
  PetscViewerDestroy(&viewer);
  
  return;
}

void output_uwp(stokes_force *uwp, parameter *para){
  
  PetscViewer viewer;
  char buffer[1000];
  if (para->temp_profile==1)
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth1_%.1f_thigh_width_%.1f_twidth2_%.1f_reinit_%d/outputUWP_t_%.6f", para->maxz, para->maxr, para->dz, para->vf,para->tension, para->r0, para->tlow, para->thigh, para->twidth1, para->thigh_width, para->twidth2, para->reinitstep, uwp->t);
  else
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_mui_%.7f_muo_%.7f_reinit_%d/outputUWP_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->mui, para->muo, para->reinitstep, uwp->t);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(uwp->data, viewer);
  PetscViewerDestroy(&viewer);
  printf("output uwp \n");  
  return;
}

void output_force(stokes_force *force, parameter *para){
  
  PetscViewer viewer;
  char buffer[1000];

    if (para->temp_profile==1)
        sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth1_%.1f_thigh_width_%.1f_twidth2_%.1f_reinit_%d/outputforce_t_%.6f", para->maxz, para->maxr, para->dz, para->vf,para->tension, para->r0, para->tlow, para->thigh, para->twidth1, para->thigh_width, para->twidth2, para->reinitstep, force->t);
    else
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_mui_%.7f_muo_%.7f_reinit_%d/outputforce_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->mui, para->muo, para->reinitstep, force->t);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_WRITE, &viewer);
  VecView(force->data, viewer);
  PetscViewerDestroy(&viewer);
  //  printf("output uwp \n");  
  return;
}


void load_levelset(levelset_vec *G, parameter *para) {
  char buffer[1000];
  PetscViewer viewer;

  if (para->temp_profile==1)
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_tlow_%.1f_thigh_%.1f_twidth1_%.1f_thigh_width_%.1f_twidth2_%.1f_reinit_%d/outputG_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->tlow, para->thigh, para->twidth1, para->thigh_width, para->twidth2, para->reinitstep, para->trestart);
    else
      sprintf(buffer, "./data/maxz_%.1f_maxr_%.1f_dz_%.2f_vf_%.5f_tension_%.2f_r0_%.2f_mui_%.7f_muo_%.7f_reinit_%d/outputG_t_%.6f", para->maxz, para->maxr, para->dz, para->vf, para->tension, para->r0, para->mui, para->muo, para->reinitstep, para->trestart);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, buffer, FILE_MODE_READ, &viewer);
  VecLoad(G->data, viewer);
  PetscViewerDestroy(&viewer);
}
