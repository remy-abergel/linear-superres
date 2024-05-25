#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

// internal modules
int complex_shannon_translation(fftw_complex*,fftw_complex*,double,double,int,int);
int complex_shannon_zooming(fftw_complex*,fftw_complex*,int,int,int,int);
int fuse_stack(double*,double*,fftw_complex**,int,int,int);
int metrics(double*,double*,double*,double*,double*,int,int,double*);
static int compar(const void*,const void*);

/* complex_shannon_translation: subpixellic shift of a 2D image (using
                                the complex Shannon interpolate,
                                referred as Uc in the IPOL companion
                                paper).

   Module parameters
   =================

   + width : width of the input/output images

   + height : height of the input/output images

   + dx : horizontal subpixel shift (double precision)

   + dy : vertical subpixel shift (double precision)

   + in : graylevel values of the input image (fftw_complex format,
          row-major order)

   + 'out': on entry -> a preallocated fftw_complex array with size
            width*height, on exit -> graylevel values of the output
            shifted image (row-major order)

   Module output 
   =============

   This module should return EXIT_SUCCESS.
   
*/
int complex_shannon_translation(out,in,dx,dy,width,height)
     fftw_complex *out,*in;
     double dx,dy;
     int width,height;
{
  fftw_complex *dft_out,*dft_in;
  fftw_plan plan_in,iplan_out;
  int a,b,x,y,adr;
  double nrm,theta,c,s,re,im;

  // memory allocation
  ASSERT_ALLOC(dft_in = (fftw_complex*) fftw_malloc (width*height*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_out = (fftw_complex*) fftw_malloc (width*height*sizeof(fftw_complex)));
  ASSERT_ALLOC(plan_in = fftw_plan_dft_2d(height,width,in,dft_in,FFTW_FORWARD,FFTW_ESTIMATE));
  ASSERT_ALLOC(iplan_out = fftw_plan_dft_2d(height,width,dft_out,out,FFTW_BACKWARD,FFTW_ESTIMATE));

  // kernel
  fftw_execute(plan_in); // set dft_in = DFT(in);
  nrm = 1./((double)width*(double)height);
  for(b=-(height/2),y=0;y<height;y++,b++) {
    for(a=-(width/2),x=0;x<width;x++,a++) {
      theta = -2.*M_PI*((double)a*dx/(double)width + (double)b*dy/(double)height);
      c = cos(theta);
      s = sin(theta);
      adr = ((b+height)%height)*width+(a+width)%width; // = location of freq. (a,b) in dft_in & dft_out
      re = dft_in[adr][0];
      im = dft_in[adr][1];
      dft_out[adr][0] = nrm * (re*c - im*s);
      dft_out[adr][1] = nrm * (im*c + re*s);
    }
  }
  fftw_execute(iplan_out);

  // free memory
  fftw_free(dft_in);
  fftw_free(dft_out);
  fftw_destroy_plan(plan_in);
  fftw_destroy_plan(iplan_out);

  return EXIT_SUCCESS;
}


/* complex_shannon_zooming: zooming of a 2D image (using the complex
                            Shannon interpolate, referred as Uc in the
                            IPOL companion manuscript).

   Module parameters
   =================

   + nx: width of the input image

   + ny: height of the input image

   + Nx: width of the output image (the zooming factor along the
         horizontal direction is Nx/nx)

   + Ny: height of the output image (the zooming factor along the
         vertical direction is Ny/ny)

   + in: graylevel values of the input image (fftw_complex format,
         row-major order)

   + 'out': graylevel values of the output image (fftw_complex format,
            row-major order)
   
   Module output 
   =============

   This module should return EXIT_SUCCESS.
   
*/
int complex_shannon_zooming(out,in,nx,ny,Nx,Ny)
     fftw_complex *out,*in;
     int nx,ny,Nx,Ny;
{
  fftw_complex *dft_in,*dft_out;
  fftw_plan plan_in,iplan_out;
  double nrm;
  int a,b,x,y,adr1,adr2;

  // memory allocation
  ASSERT_ALLOC(dft_in = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_out = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(plan_in = fftw_plan_dft_2d(ny,nx,in,dft_in,FFTW_FORWARD,FFTW_ESTIMATE));
  ASSERT_ALLOC(iplan_out = fftw_plan_dft_2d(Ny,Nx,dft_out,out,FFTW_BACKWARD,FFTW_ESTIMATE));

  // kernel
  fftw_execute(plan_in);
  memset(dft_out,0.,Nx*Ny*sizeof(fftw_complex));
  nrm = 1./((double)nx*(double)ny);
  for(b=-(ny/2),y=0;y<ny;y++,b++) {
    for(a=-(nx/2),x=0;x<nx;x++,a++) {
      adr1 = ((b+ny)%ny)*nx+(a+nx)%nx; // = location of freq. (a,b) in dft_in
      adr2 = ((b+Ny)%Ny)*Nx+(a+Nx)%Nx; // = location of freq. (a,b) in dft_out
      dft_out[adr2][0] = nrm * dft_in[adr1][0];
      dft_out[adr2][1] = nrm * dft_in[adr1][1];
    }
  }
  fftw_execute(iplan_out);

  // free memory
  fftw_free(dft_in);
  fftw_free(dft_out);
  fftw_destroy_plan(plan_in);
  fftw_destroy_plan(iplan_out);

  return EXIT_SUCCESS;
}


/* metrics: evaluate MSE/SNR/PSNR metrics as described below.

   Given two images im1, im2, containing N pixels, compute: 

   MSE(im1,im2) = ||im1-im2|| / N
   SNR(im1,im2) = 10*log10(VAR(im1)/MSE(im1,im2))
   PSNR(im1,im2) = 10*log10(DYN(im1)^2/MSE(im1,im2))

   denoting by DYN(im1) the dynamic of im1 (can be specified by the
   user, of set by default using DYN(im1) = MAX(im1)-MIN(im1))

   Module parameters
   =================

   + width: width of the input images im1 & im2

   + height: height of the input images im1 & im2

   + dyn: dynamic value to use in PSNR evaluation (when dyn == NULL,
          we use dyn = MAX(im1)-MIN(im1))

   + im1: graylevel values of the input image im1 (row-major order)

   + im2: graylevel values of the input image im2 (row-major order)

   + psnr: on exit, *psnr = computed value of PSNR(im1,im2)

   + snr: on exit, *snr = computed value of SNR(im1,im2)

   + mse: on exit, *mse = computed value of MSE(im1,im2)

   Module output 
   =============

   This module should return EXIT_SUCCESS.
   
*/
int metrics(psnr,snr,mse,im1,im2,width,height,dyn)
     double *psnr,*snr,*mse,*im1,*im2,*dyn;
     int width,height;
{
  double var1,mean1,diff,size,max1,min1;
  int adr;

  size = (double)width*(double)height;

  // compute the mean, max and min of im1 //
  for(mean1=0.,max1=min1=im1[adr],adr=0;adr<width*height;adr++) {
    mean1 += im1[adr];
    if(max1<im1[adr]) max1 = im1[adr];
    if(min1>im1[adr]) min1 = im1[adr];
  }
  mean1 /= size;

  // compute the variance of im1 //
  for(var1=0.,adr=0;adr<width*height;adr++) {
    diff = (im1[adr]-mean1);
    var1 += diff*diff;
  }
  var1 /= size - 1.;

  // compute the metrics //
  for(*mse=0.,adr=0;adr<width*height;adr++) {
    diff = im1[adr]-im2[adr];
    (*mse) += diff*diff;
  }
  (*mse) /= size;
  (*snr) = 10.*log10(var1/(*mse));
  (*psnr) = 10.*log10(((dyn)?((*dyn)*(*dyn)):((max1-min1)*(max1-min1)))/(*mse));

  return EXIT_SUCCESS;
}



/* compar: comparison function for sorting double values in ascending order using qsort() */
static int compar(const void *a, const void *b)
{
  double v = *((double const *)a)-*((double const *)b);
  return (v>0 ? 1 : ((v<0) ? -1 : 0));
}


/* fuse_stack: fuse a stack of 2D images by pixelwise averaging or
               median along the temporal (third) dimension (useful for
               shift-and-add like algorithms)

   Module parameters
   =================

   + width: width of the input stack

   + height: height of the input stack

   + nimages: number of images contained in the input stack

   + in: input stack, such as in[i] corresponds to the graylevel
         values of the i-th image (0 <= i < nimages) of the stack,
         stored in row-major order

   + av: on enty, a preallocated double array with size width*height,
         on exit, graylevel values of the pixelwise temporal averaging
         of the input stack

   + med: on entry, a preallocated double array with size
          width*height, on exit, graylevel values of the pixelwise
          temporal median of the input stack

   Module output 
   =============

   This module should return EXIT_SUCCESS.

*/
int fuse_stack(med,av,in,width,height,nimages)
     double *med,*av;
     fftw_complex **in;
     int width,height,nimages;
{
  double *tmp,m,val;
  int adr,id;

  // memory allocation
  if(med) { ASSERT_ALLOC(tmp = (double*) malloc (nimages*sizeof(double))); }

  // main loop
  for(adr=0;adr<width*height;adr++) {
    m = 0.;
    for(id=0;id<nimages;id++) {
      val = in[id][adr][0];
      if(med) tmp[id] = val;
      if(av) m += val;
    }
    if(med) {
      qsort((void *)tmp,nimages,sizeof(double),compar);
      med[adr] = ((nimages%2) ? tmp[nimages/2] : 0.5*(tmp[nimages/2] + tmp[-1+nimages/2]));
    }
    if(av) av[adr] = m/(double)nimages;
  }

  // free memory
  if(med) free(tmp);

  return EXIT_SUCCESS;
}
