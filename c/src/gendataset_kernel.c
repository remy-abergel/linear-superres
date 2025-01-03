#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }
#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

/* external modules */
extern int direct_stack_operator(fftw_complex**,fftw_complex*,double*,double*,int,int,int,int,int); // see source file 'operators_kernel.c'

/* internal modules */
double getval(fftw_complex*,double,double,int,int);
int perdecomp(fftw_complex*,fftw_complex*,fftw_complex*,int,int);
int gendataset(fftw_complex*, fftw_complex***, fftw_complex**, double*, double*, int*, int*, int*, int*, int, double*, double*, int, int, int, int, char, char);

/* getval: get the subpixellic value of the input image (u) using
   bilinear interpolation.

   Note: the input image is complex but only its real part is being
   processed.

   Module parameters
   =================

   + u : (fftw_complex array) contains the graylevels of the input
         image stored in ROW-MAJOR order
   
   + nx, ny : horizontal and vertical dimensions (width and height) of
              the input image (u)
   
   + x, y : subpixellic horizontal (x) and vertical (y) coordinates

   Module output
   =============

   + out : bilinear interpolation of the real part of u sampled at the
           subpixellic location (x,y)
	   
*/
double getval(u,x,y,nx,ny)
  fftw_complex *u;
  double x,y;
  int nx,ny;
{
  int k,l,x0,y0,adr;
  double u00,u01,u10,u11,dx,dy,out;

  // decompose (x,y) = (x0,y0) + (dx,dy)
  x0 = (int) floor(x);
  y0 = (int) floor(y);
  dx = x-(double)x0; 
  dy = y-(double)y0; 

  // get top-left value
  k = max(0,min(x0,nx-1));
  l = max(0,min(y0,ny-1));
  adr = l*nx+k;
  u00 = u[adr][0];

  // get top-right value
  k = max(0,min(x0+1,nx-1));
  l = max(0,min(y0,ny-1));
  adr = l*nx+k;
  u10 = u[adr][0];
  
  // get bottom-left value
  k = max(0,min(x0,nx-1));
  l = max(0,min(y0+1,ny-1));
  adr = l*nx+k;
  u01 = u[adr][0];

  // get bottom-right value
  k = max(0,min(x0+1,nx-1));
  l = max(0,min(y0+1,ny-1));
  adr = l*nx+k;
  u11 = u[adr][0];

  // compute bilinear interpolation of real(u) at location (x,y)
  out = (u00*(1.-dx) + u10*dx)*(1.-dy) + (u01*(1.-dx) + u11*dx)*(dy);

  return out;
}

/* perdecomp: periodic plus smooth decomposition (see [1])

   [1] L. Moisan. "Periodic plus smooth image decomposition." Journal
       of Mathematical Imaging and Vision 39 (2011): 161-179.
       
   Note: the input image is complex but only its real part is being
   processed.

   Module parameters
   =================

   + u : (fftw_complex array) contains the graylevels of the input
         image stored in ROW-MAJOR order
	 
   + nx, ny : horizontal and vertical dimensions (width and height) of
              the input image (u)
	      
   + p : (fftw_complex array) on entry, a preallocated array with size
         nx*ny, on exit, the real part contains the graylevels of the
         periodic component of the real part of u and the imaginary
         part is filled with zeros

   + s : (fftw_complex array) on entry, a preallocated array with size
         nx*ny, on exit, the real part contains the graylevels of the
         smooth component of the real part of u and the imaginary part
         is filled with zeros
   
   Module output
   =============
   
   This module should return EXIT_SUCCESS if no error occured during
   its execution.
   
 */
int perdecomp(u,p,s,nx,ny)
     fftw_complex *u,*p,*s;
     int nx,ny;
{
  fftw_complex *dft_s;
  fftw_plan plan_s,iplan_s;
  int adr,x,y;
  double a,cx,cy,b,cof;

  /* memory allocations */
  ASSERT_ALLOC(dft_s = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(plan_s = fftw_plan_dft_2d(ny,nx,s,dft_s,FFTW_FORWARD,FFTW_ESTIMATE));
  ASSERT_ALLOC(iplan_s = fftw_plan_dft_2d(ny,nx,dft_s,s,FFTW_BACKWARD,FFTW_ESTIMATE));

  /* fill "boundary image" */
  memset(s,0.,nx*ny*sizeof(fftw_complex));
  for (x=0;x<nx;x++) {
    b = u[x][0]-u[(ny-1)*nx+x][0];
    s[x][0] += b;
    s[(ny-1)*nx+x][0] += -b;
  }
  for (y=0;y<ny;y++) {
    b = u[y*nx][0]-u[y*nx+nx-1][0];
    s[y*nx][0] += b;
    s[y*nx+nx-1][0] += -b;
  }

  /* Fourier transform */
  fftw_execute(plan_s); // set dft_s = DFT(s);
  
  /* (0,0) frequency */
  dft_s[0][0] = dft_s[0][1] = 0.;

  cx = 2.*M_PI/(double)nx;
  cy = 2.*M_PI/(double)ny;

  /* other frequencies */
  y = 1; // to avoid (0,0) frequency
  cof = 1/(double)(nx*ny);
  for (x=0;x<nx;x++) {
    a = cos(cx*(double)x);
    for (;y<ny;y++) {
      adr = y*nx+x;
      b = 0.5/(double)(2.-a-cos(cy*(double)y));
      dft_s[adr][0] *= b*cof;
      dft_s[adr][1] *= b*cof;
    }
    y=0;
  }

  /* inverse Fourier transform */
  fftw_execute(iplan_s); // now s = IDFT(dft_u)
  
  /* get p=u-s */
  for (adr=nx*ny;adr--;) {
    p[adr][0] = u[adr][0]-s[adr][0];
    p[adr][1] = 0.;
  }
  
  /* free memory */
  fftw_free(dft_s);
  fftw_destroy_plan(plan_s);
  fftw_destroy_plan(iplan_s);
  
  return EXIT_SUCCESS;
}

/* gendataset: generate a synthetic low-resolution stack from the
               input high-resolution image (in), avoiding
               periodic-like boundaries

   Module parameters
   =================

   + in : (fftw_complex array) contains the input (high-resolution)
          image graylevels stored in ROW-MAJOR order.

   + Nx, Ny : on entry, *Nx and *Ny correspond to the width and height
              of the input high-resolution image in, on exit, *Nx and
              *Ny are equal to the width and height of the cropped
              high-resolution reference image

   + nx, ny : on exit, *nx and *ny correspond to the width and height
              of the generated low-resolution images stored in u0

   + nimages : number of low-resolution images to generate

   + dx : double array containing the (nimages) horizontal shifts
   
   + dy : double array containing the (nimages) vertical shifts
   
   + zx, zy : on entry, *zx and *zy correspond to the requested
              horizontal and vertical subsampling factor, on exit *zx
              and *zy correspond to the actually applied horizontal
              and vertical subsampling factors

   + cx : (integer) allow discarding up to cx columns to the input
          high-resolution image in order to minimize the distance
          between the requested and actually used value of the
          horizontal subsampling factor (*zx)

   + cy : (integer) allow discarding up to cy rows to the input
          high-resolution image in order to minimize the distance
          between the requested and actually used value of the
          vertical subsampling factor (*zy)
	  
   + p, q : (integers) the final high-resolution ROI size will be
            equal to (q/p) * the dimensions of the cropped HR image
            (i.e., after the cropping step of up to cy rows and cx
            columns)

   + ref_needed : (char) if ref_needed != NULL, ref will be allocated
                  and filled with the graylevel values of the
                  high-resolution image associated to the generated
                  low-resolution sequence (obtained by cropping again
                  the input image, the cropping ROI is integer by
                  construction so no interpolation is needed).

   + vflag : (char) enable verbose mode when vflag != NULL

   + u0 : on entry, u0 is a non allocated fftw_complex***, on exit,
          (*u0) is sequence of low-resolution images synthesized from
          in (--> (*u0)[k] corresponds the the k-th low-resolution
          image of the sequence)

   + ref : on entry, ref is a non allocated fftw_complex**, on exit,
           ref is either unchanged (when ref_needed is NULL) or an
           allocated array with size (*Nx)*(*Ny) containing the
           graylevel values of the reference high-resolution image
           associated to the sequence u0
   
   Module output
   =============
   
   This module should return EXIT_SUCCESS if no error occured during
   its execution.
       
*/
int gendataset(in,u0,ref,dx,dy,Nx,Ny,nx,ny,nimages,zx,zy,cx,cy,p,q,ref_needed,vflag)
     fftw_complex *in,**ref,***u0; 
     double *dx,*dy,*zx,*zy;
     int *Nx,*Ny,*nx,*ny,cx,cy,p,q,nimages;
     char vflag,ref_needed;
{
  fftw_complex *u, *up, *us, **p0, **s0; 
  double zx_tmp, zy_tmp, smin, s, xx0, yy0,tx,ty;
  int k, kopt, Ax, Ay, ax, ay, Nx2, Ny2, x0, y0, X0, Y0, x, y, adr, adr0, err;
  
  // first cropping step along the X-axis (we will remove up to cx
  // columns from the input HR image in order fulfill as best as
  // possible the zx requirements)
  zx_tmp = *zx;
  smin = INFINITY;
  for(k=0;k<=cx;k++) {
    Ax = (int)floor((*Nx-k)/p);
    ax = (int)round(Ax/zx_tmp);
    s = (double)(*zx) - (double)Ax/(double)ax;
    s *= s;
    if(s < smin) { smin = s; kopt = k; }
  }
  Ax = (int)floor((*Nx-kopt)/p);
  ax = (int)round(Ax/zx_tmp);
  *zx = (double)Ax/(double)ax;
  Nx2 = Ax*p; 
  if(zx_tmp != *zx) printf("WARNING: zx has been changed into %.17g (instead of zx=%g) to allow integer width for the output image\n",*zx,zx_tmp);
  
  // first cropping step along the Y-axis (remove up to cy rows from
  // the input HR image in order fulfill as best as possible the zx
  // requirements)
  zy_tmp = *zy;
  smin = INFINITY;
  for(k=0;k<=cy;k++) {
    Ay = (int)floor((*Ny-k)/p);
    ay = (int)round(Ay/zy_tmp);
    s = (double)(*zy) - (double)Ay/(double)ay;
    s *= s;
    if(s < smin) {smin = s; kopt = k; }
  }
  Ay = (int)floor((*Ny-kopt)/p);
  ay = (int)round(Ay/zy_tmp);
  *zy = (double)Ay/(double)ay;
  Ny2 = Ay*p; 
  if(zy_tmp != *zy) printf("WARNING: zy has been changed into %.17g (instead of zy=%g) to allow integer height for the output image\n",*zy,zy_tmp);
  
  // perform first cropping step
  x0 = (*Nx-Nx2)/2; // integer division here
  y0 = (*Ny-Ny2)/2; // integer division here
  ASSERT_ALLOC(u = (fftw_complex*) fftw_malloc (Nx2*Ny2*sizeof(fftw_complex)));
  for(y=y0;y<y0+Ny2;y++)
    for(x=x0;x<x0+Nx2;x++) {
      adr0 = (y-y0)*Nx2+(x-x0);
      u[adr0][0] = in[y*(*Nx)+x][0];
      u[adr0][1] = 0.;
    }
  
  // perform periodic + smooth decomposition of the HR image (u)
  ASSERT_ALLOC(up = (fftw_complex*) fftw_malloc (Nx2*Ny2*sizeof(fftw_complex)));
  ASSERT_ALLOC(us = (fftw_complex*) fftw_malloc (Nx2*Ny2*sizeof(fftw_complex)));
  perdecomp(u,up,us,Nx2,Ny2);
  
  // compute recropping parameters & shift input translations to allow
  // integer positions of the top-left corner of the cropping area in
  // the upcoming LR stacks (p0 & s0)
  *Nx = q*Ax;
  *nx = q*ax;
  *Ny = q*Ay;
  *ny = q*ay;
  X0 = (Nx2 - *Nx)/2;
  Y0 = (Ny2 - *Ny)/2;
  xx0 = X0/(*zx); 
  yy0 = Y0/(*zy);
  x0 = (int)floor(xx0);
  y0 = (int)floor(yy0);
  tx = xx0 - (double)x0;
  ty = yy0 - (double)y0;
  for(k=0;k<nimages;k++) {
    dx[k] += tx;
    dy[k] += ty;
  }

  // deal with verbose mode
  if(vflag) {
    printf("Retrieve image size parameters\n");
    printf("==============================\n");
    printf("width of the low-resolution domain: m = %d\n",*nx);
    printf("height of the low-resolution domain: n = %d\n",*ny);
    printf("width of the high-resolution domain: M = %d\n",*Nx);
    printf("height of the high-resolution domain: N = %d\n",*Ny);
    printf("horizontal subsampling factor: zx = M/m = %g\n",*zx);
    printf("vertical subsampling factor: zy = N/n = %g\n",*zy);
    printf("\n");
  }
  
  // memory allocation for p0, s0 & u0
  ASSERT_ALLOC((*u0) = (fftw_complex**) malloc (nimages*sizeof(fftw_complex*)));
  ASSERT_ALLOC(p0 = (fftw_complex**) malloc (nimages*sizeof(fftw_complex*)));
  ASSERT_ALLOC(s0 = (fftw_complex**) malloc (nimages*sizeof(fftw_complex*)));
  for(k=0;k<nimages;k++) {
    ASSERT_ALLOC((*u0)[k] = (fftw_complex*) fftw_malloc ((*nx)*(*ny)*sizeof(fftw_complex)));
    ASSERT_ALLOC(p0[k] = (fftw_complex*) fftw_malloc (ax*p*ay*p*sizeof(fftw_complex)));
    ASSERT_ALLOC(s0[k] = (fftw_complex*) fftw_malloc (ax*p*ay*p*sizeof(fftw_complex)));
  }
  
  // compute low-resolution stack p0
  err = direct_stack_operator(p0,up,dx,dy,ax*p,ay*p,Ax*p,Ay*p,nimages);
  
  if(err != EXIT_SUCCESS) {
    for(k=0;k<nimages;k++) { fftw_free((*u0)[k]); fftw_free(p0[k]); fftw_free(s0[k]); }
    free(*u0); free(p0); free(s0); fftw_free(u); fftw_free(up); fftw_free(us); 
    return EXIT_FAILURE;
  }
  
  fftw_free(up);
  
  // compute low-resolution stack s0 (shift & subsample us using bilinear interpolation)
  for(k=0;k<nimages;k++)
    for(adr=y=0;y<ay*p;y++)
      for(x=0;x<ax*p;x++,adr++){
	s0[k][adr][0] = getval(us,(*zx)*((double)(x)+dx[k]),(*zy)*((double)(y)+dy[k]),Ax*p,Ay*p);
	s0[k][adr][1] = 0.; 
      }
  
  fftw_free(us);
  
  // compute (if needed) the corresponding HR image
  if(ref_needed) {
    ASSERT_ALLOC((*ref) = (fftw_complex*) fftw_malloc((*Nx)*(*Ny)*sizeof(fftw_complex)));
    for(y=Y0;y<Y0+(*Ny);y++)
      for(x=X0;x<X0+(*Nx);x++) {
	adr0 = (y-Y0)*(*Nx)+(x-X0);
	(*ref)[adr0][0] = u[y*(Ax*p)+x][0]; 
	(*ref)[adr0][1] = 0; 
      }
  }
  
  // compute & crop the low-resolution stack p0+s0
  for(k=0;k<nimages;k++)
    for(y=y0;y<y0+(*ny);y++)
      for(x=x0;x<x0+(*nx);x++) {
	adr0 = (y-y0)*(*nx)+(x-x0);
	adr = y*(ax*p)+x;
	(*u0)[k][adr0][0] = p0[k][adr][0] + s0[k][adr][0];
	(*u0)[k][adr0][1] = 0.;
      }
  
  // free memory (everything except u0 and ref)
  for(k=0;k<nimages;k++) { fftw_free(p0[k]); fftw_free(s0[k]); }
  free(p0); free(s0);
  fftw_free(u);

  return EXIT_SUCCESS;
}
