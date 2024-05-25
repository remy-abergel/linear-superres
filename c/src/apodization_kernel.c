#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* apodization: compute low and/or high resolution multiplicative
                apodization filters, and/or compute the apodized
                low-resolution sequence.

   The construction of those apodization filters is discussed in
   detail in Section 2.3 of the IPOL paper. We briefly recall here
   that the goal of those filters is to avoid artifcats due to the non
   periodicity (in the sense that the borders of those images are not
   alike) of the real life low-resolution sequences.

   Module parameters
   =================

   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              denoted M and Ny is denoted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              denoted m and ny is denoted n)

   + nimages : number of images contained in the input stack of
               low-resolution images (in the IPOL paper, nimages is
               denoted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + vflag : a flag to enable/disable the verbose mode (set vflag=0 to
             disable the verbose mode, or vflag!=0 to enable the
             verbose mode)

   + sigma : sharpness parameter for the 1D apodization profiles (see
             Equation (14) of the IPOL paper). This parameter is such
             that 10*sigma is the length of the transition intervals
             where the profiles increase from 0 to 1 or decrease from
             1 to 0

   + u0 : on entry -> the stack of low-resolution input images, on
          exit -> the apodized sequence of low-resolution images
          (correspond to the sequence "u0 tilde" in Equation (16) of
          the IPOL paper)

   + apod_lr : on entry -> a preallocated double array with same size
               of u0, on exit -> the stack of low-resolution filters
               involved in the right-hand term of Equation (16) of the
               IPOL paper

   + apod_hr : on entry -> a preallocated double array with size
               Nx*Ny, on exit -> the high-resolution filter (width =
               Nx, height = Ny) involved in the right-hand term of
               Equation (15)

   NOTE: calling this module with apod_lr = NULL or/and apod_hr = NULL
   or/and u0 = NULL is possible. When one of those pointers is NULL,
   its corresponding values are never computed (it remains NULL on
   exit).

   Module output
   =============

   This module should return EXIT_SUCCESS.

 */
int apodization(u0,apod_lr,apod_hr,dx,dy,sigma,nx,ny,Nx,Ny,nimages,vflag)
     double **u0,**apod_lr,*apod_hr,*dx,*dy,sigma;
     int nx,ny,Nx,Ny,nimages;
     char vflag;
{
  int x,y,adr,id;
  double xx,yy,zx,zy,tx,ty,Dx,Dy,cof,half_Nx,half_Ny,NNx,NNy,vx,vy;

  // initialization //
  zx = (double)Nx/(double)nx;
  zy = (double)Ny/(double)ny;
  Dx = fabs(dx[0]);
  Dy = fabs(dy[0]);
  for(id=1;id<nimages;id++) { // compute highest elements of dx and dy
    if(Dx < fabs(dx[id])) Dx = fabs(dx[id]);
    if(Dy < fabs(dy[id])) Dy = fabs(dy[id]);
  }
  Dx *= zx; // highest displacement along the horizontal direction
  Dy *= zy; // highest displacement along the vertical direction
  cof = 1./(sigma*sqrt(2.));
  half_Nx = 0.5*((double)Nx-1.);
  half_Ny = 0.5*((double)Ny-1.);
  NNx = half_Nx-Dx-5.*sigma;
  NNy = half_Ny-Dy-5.*sigma;

  // deal with the low-resolution apodization filters //
  if(u0 || apod_lr){
    if(vflag) {
      printf("Compute low-resolution filters\n");
      printf("==============================\n");
    }
    for(id=0;id<nimages;id++) {
      tx = dx[id];
      ty = dy[id];
      for(adr=y=0;y<ny;y++) {
	yy = zy*((double)y+ty);
	vy = 0.5*erfc(cof*(fabs(half_Ny-yy)-NNy));
	for(x=0;x<nx;x++,adr++){
	  xx = zx*((double)x+tx);
	  vx = 0.5*erfc(cof*(fabs(half_Nx-xx)-NNx));
	  if(u0) u0[id][adr] *= vx*vy;
	  if(apod_lr) apod_lr[id][adr] = vx*vy;
	}
      }
    }
    if(vflag && apod_lr) printf("compute low-resolution filters: done\n");
    if(vflag && u0) printf("apodization of the stack u0: done\n");
    if(vflag) printf("\n");
  }

  // deal with the high-resolution apodization filter //
  if(apod_hr) {
    if(vflag) {
      printf("Compute high-resolution filter\n");
      printf("==============================\n");
    }
    for(adr=y=0;y<Ny;y++) {
      vy = 0.5*erfc(cof*(fabs(half_Ny-(double)y)-NNy));
      for(x=0;x<Nx;x++,adr++) {
	vx = 0.5*erfc(cof*(fabs(half_Nx-(double)x)-NNx));
	apod_hr[adr] = vx*vy;
      }
    }
    if(vflag) printf("done\n\n");
  }

  return EXIT_SUCCESS;

}
