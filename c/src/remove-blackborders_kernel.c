#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* internal modules */
int remove_blackborders(double**,int*,int*,double*,double*,double*,int,int,int,int,int,double,char); 

/* remove_blackborders: crop a high resolution image to remove borders
                        caused by apodization
  
   Module Parameters
   =================
   
   + u : graylevels of the input high-resolution image 
 
   + Nx : width of the high-resolution image domain
 
   + Ny : height of the high-resolution image domain 
 
   + nx : width of the low-resolution image domain
  
   + ny : height of the low-resolution image domain 
 
   + dx : (double array) horizontal displacement sequence (must be the
          same sequence as that used for apodization of the
          low-resolution sequence)
 
   + dy : (double array) vertical displacement sequence (must be the
          same sequence as that used for apodization of the
          low-resolution sequence)
  
   + ntranslations : number of displacement vectors (number of
                     elements in dx and dy arrays)
  
   + sigma : sharpness parameter of the apodization profiles (must be
             the same as that used for apodization of the low-resolution
             sequence)
  
   + vflag : verbose flag (set vflag = 0 to disable verbose mode, or
             vflag != 0 to enable verbose mode)
  
   + ucrop : on entry -> non allocated double array, on exit ->
             allocated double array containing the graylevels of the
             output cropped image
 
   + Nxcrop : on exit -> width of the output cropped image
  
   + Nycrop : on exit -> height of the output cropped image
 
   Module output
   =============

   This module returns EXIT_SUCCESS if no error occured during its
   execution, it returns EXIT_FAILURE otherwise.
 
*/
int remove_blackborders(double **ucrop,int *Nxcrop, int *Nycrop, double *u, double *dx, double *dy, int Nx, int Ny, int nx, int ny, int ntranslations, double sigma, char vflag)
{
  int adr,x,y,delta_x,delta_y; 
  double dxmax=dx[0],dymax=dy[0];

  // retrieve maximal displacements along both directions // 
  for(adr=0;adr<ntranslations;adr++) {
    dxmax = fmax(dxmax,dx[adr]); 
    dymax = fmax(dymax,dy[adr]); 
  }

  // compute thickness of the black borders along both directions //  
  delta_x = (int)ceil(((double)Nx/(double)nx)*dxmax + 10.*sigma);
  delta_y = (int)ceil(((double)Ny/(double)ny)*dymax + 10.*sigma);

  // compute dimensions of the cropped image //
  *Nxcrop = Nx-2*delta_x; 
  *Nycrop = Ny-2*delta_y;
  if((*Nxcrop <= 0) || (*Nycrop <= 0)) {
    printf("could not remove black-borders from the input image (no pixel would remain left after cropping)\n");
    return EXIT_FAILURE; 
  }

  // deal with verbose mode // 
  if(vflag) {
    printf("Remove black-borders caused by apodization\n");
    printf("==========================================\n");
    printf("remove %d pixels from the left side\n",delta_x);
    printf("remove %d pixels from the right side\n",delta_x);
    printf("remove %d pixels from the top side\n",delta_y);
    printf("remove %d pixels from the bottom side\n",delta_y);
    printf("\n");
  }
  
  // memory allocation //
  ASSERT_ALLOC(*ucrop = (double*) malloc ((*Nxcrop)*(*Nycrop)*sizeof(double)));

  // image crop //
  for(adr=0,y=delta_y;y<=Ny-delta_y-1;y++)    
    for(x=delta_x;x<=Nx-delta_x-1;x++,adr++) 
      (*ucrop)[adr] = u[y*Nx+x]; 
  
  return EXIT_SUCCESS;
  
}
