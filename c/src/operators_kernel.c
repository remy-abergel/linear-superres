#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* internal modules */
int direct_operator_dft(fftw_complex*,fftw_complex*,double,double,int,int,int,int,char);
int adjoint_operator_dft(fftw_complex*,fftw_complex*,double,double,double,int,int,int,int,char);
int direct_stack_operator(fftw_complex**,fftw_complex*,double*,double*,int,int,int,int,int);
int adjoint_stack_operator_dft(fftw_complex*,fftw_complex**,double*,double*,int,int,int,int,int,char);
int adjoint_weighted_stack_operator_dft(fftw_complex*,fftw_complex**,double*,double*,double*,int,int,int,int,int,char);

/*
   General Note about the 'normalize' parameter in the *_dft modules:
   ==================================================================

   The FFTW library provides routines for forward and backward Discrete Fourier
   Transforms (DFTs).

   However, given a signal 's' containing N elements, applying successively the
   forward and backward FFTW routine to 's', we do not retrieve 's' but the
   signal 'N*s'.

   This mean that when computing directly in 'dft_s' the DFT coefficients of the
   signal 's', we can be interested in applying the normalization dft_s /= N
   before returning 'dft_s'. When such normalization is done, 'dft_s' does not
   contain the DFT coefficients of 's' anymore but the normalized DFT
   coefficients, so that, if we further apply the Backward routine in FFTW to
   'dft_s' we get the signal 's' (instead of 'N*s');

   In all modules dedicated to computation of signals directly in the Frequency
   domain (see the modules *_dft below), we use the 'normalize' parameter to
   specify wether the output signal should be normalized or not.

*/

/*******************************************************************************/
/*                                                                             */
/*     MODULES DEDICATED TO THE COMPUTATION OF THE DIRECT OPERATORS (A,Aj)     */
/*                                                                             */
/*******************************************************************************/

/* direct_operator_dft: shift and subsample an image (the operation is
                        done in the Fourier domain)

   Apply sub-pixellic translation and subsampling to the input image
   (see mathematical description below). This operator corresponds to
   the "Aj" operator described the IPOL paper, and its computation is
   based on Equation (27).

   In this module, the input (high-resolution) image is given in the
   frequency domain and the output (low-resolution) image is computed
   in the frequency domain as well (concretly, this means that we
   compute the DFT coefficents of the output image 'out' from the DFT
   coefficients of the input image 'in').

   Module parameters
   =================

   + dx, dy : the horizontal and vertical components of the
              translation vector ('dx' and 'dy' can be any real
              number).

   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + dft_in : a ffw_complex array containing the (unnormalized) Nx*Ny
              DFT coefficients of the input image (ROW-MAJOR order).

   + normalize : a flag that specifies wether dft_out should be
                 normalized (i.e. divided by (nx*ny)) or not on
                 exit. Set normalize=0 to disable normalization, or
                 normalized!=0 to enable normalization.

   + dft_out : on entry -> a preallocated a fftw_complex array
               containing nx*ny ellements, on exit -> contains the
               (normalized or unnormalized) DFT coefficients of the
               shifted and subsampled image 'out' described below.

   Module output
   =============

   This module should return EXIT_SUCCESS.

   Mathematical description:
   =========================

   Denoting by zx = Nx/nx and zy = Ny/ny, in the spatial domain, we
   have, for any pixel location (x,y) in {0,...,nx-1} x {0,...,ny-1},

   out[y*nx+x][0] = real part of Uc(zx*(x+dx),zy*(y+dy)),
   out[y*nx+x][1] = imaginary part of Uc(zx*(x+dx),zy*(y+dy)),

   where Uc denotes the "complex variant" of the Shannon interpolate of the
   input image 'in', defined in Equation (10) of the IPOL paper.

*/
int direct_operator_dft(dft_out,dft_in,dx,dy,nx,ny,Nx,Ny,normalize)
     fftw_complex *dft_in,*dft_out;
     double dx,dy;
     int nx,ny,Nx,Ny;
     char normalize;
{
  int alf,bet,a,b,p,q,p0,p1,q0,q1,adr_alf,adr_bet,adr_a,adr_b,adr;
  double zx,zy,cof,cx,cy,nrm,s_re,s_im,phi_re,phi_im,dft_re,dft_im,theta;

  // initialize some local variables //
  zx = (double)Nx/(double)nx;
  zy = (double)Ny/(double)ny;
  cof = 1./(zx*zy);
  cx = 2.*M_PI*dx/(double)nx;
  cy = 2.*M_PI*dy/(double)ny;
  nrm = (normalize) ? 1./((double)nx*(double)ny) : 1.;

  // main loop: range the low-frequency domain in order to fill 'dft_out' frequency by frequency //
  for(b=-ny/2; b<=(ny-1)/2; b++) {
    adr_b = (b+ny)%ny;                               // vertical position of frequency b in 'dft_out'
    q0 = (int)ceil(-zy/2.-(double)b/(double)ny);     // first element of the index set I_{Ny,ny}(b)
    q1 = (int)ceil( zy/2.-(double)b/(double)ny)-1;   // last element of the index set I_{Ny,ny}(b)
    for(a=-nx/2; a<=(nx-1)/2; a++) {
      adr_a = (a+nx)%nx;                             // horizontal position of frequency a in 'dft_out'
      p0 = (int)ceil(-zx/2.-(double)a/(double)nx);   // first element of the index set I_{Nx,nx}(a)
      p1 = (int)ceil( zx/2.-(double)a/(double)nx)-1; // last element of the index set I_{Nx,nx}(a)
      s_re = s_im = 0.;
      for(q=q0;q<=q1;q++) {
	bet = b + q*ny;
	adr_bet = (Ny+bet)%Ny;                       // vertical position of frequency bet in 'dft_in'
	for(p=p0;p<=p1;p++) {
	  alf = a + p*nx;
	  adr_alf = (Nx+alf)%Nx;                     // horizontal position of frequency alf in 'dft_in'
	  adr = adr_bet*Nx+adr_alf;                  // position of the frequency (alf,bet) in 'dft_in'
	  dft_re = dft_in[adr][0];                   // real part of the DFT coefficient stored in 'dft_in' at location (alf,bet)
	  dft_im = dft_in[adr][1];                   // imaginary part of the DFT coefficient stored in 'dft_in' at location (alf,bet)
	  theta = cx*alf+cy*bet;
	  phi_re = cos(theta);                       // real part of the ramp-phase at location (alf,bet)
	  phi_im = sin(theta);                       // imaginary part of the ramp-phase at location (alf,bet)
	  s_re += (dft_re*phi_re - dft_im*phi_im);
	  s_im += (dft_re*phi_im + dft_im*phi_re);
	}
	// fill 'dft_out' at the frequency location (a,b) //
	adr = adr_b*nx+adr_a;
	dft_out[adr][0] = nrm * cof * s_re;
	dft_out[adr][1] = nrm * cof * s_im;
      }
    }
  }

  return EXIT_SUCCESS;
}


/* direct_stack_operator: (companion module of 'direct_operator_dft')

   In this module, we compute a sequence (or stack) of translated and
   subsampled images. Each image of the ouptut (low-resolution) stack
   is computed by applying the Aj operator to the input
   (high-resolution) image. This operator corresponds to the "A"
   operator described the IPOL paper (see Equation (12) for the formal
   definition).

   In this module, the input (high-resolution) image is given in the
   spatial domain and the output (low-resolution) stack is computed in
   the spatial domain.

   Module parameters
   =================

   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + nimages : number of images contained in the input stack of
               low-resolution images (in the IPOL paper, nimages is
               noted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + in : (fftw_complex array) contains the input (high-resolution)
          image graylevels stored in ROW-MAJOR order.

   + 'out' : contains the graylevels of the output stack of
             low-resolution images which is described below.

   Module output
   =============

   This module should return EXIT_SUCCESS.

   Mathematical description:
   =========================

   Denoting by zx = Nx/nx and zy = Ny/ny, in the spatial domain, we
   have,

   for any pixel location (x,y) in {0,...,nx-1} x {0,...,ny-1},
   for any j in {0,...,nimages-1},

   out[j][y*nx+x][0] = real part of Uc(zx*(x+dx[j]),zy*(y+dy[j])),
   out[j][y*nx+x][1] = imaginary part of Uc(zx*(x+dx[j]),zy*(y+dy[j])),

   where Uc denotes the "complex variant" of the Shannon interpolate of the
   input image 'in'.

*/
int direct_stack_operator(out,in,dx,dy,nx,ny,Nx,Ny,nimages)
     fftw_complex **out,*in;
     double *dx,*dy;
     int nx,ny,Nx,Ny,nimages;
{
  fftw_complex *dft_in,*dft_tmp;
  fftw_plan plan_in,iplan_tmp;
  int j,err=EXIT_SUCCESS;
  char normalize=1;

  /* memory allocation */
  ASSERT_ALLOC(dft_in = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_tmp = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(plan_in = fftw_plan_dft_2d(Ny,Nx,in,dft_in,FFTW_FORWARD,FFTW_ESTIMATE));

  /* CORE: loop over j and to apply each Aj operator */
  fftw_execute(plan_in); // now dft_in contains the DFT coefficients of 'in'
  for(j=0;j<nimages;j++) {

    // compute in dft_tmp the (normalized) DFT coefficients of Aj(in[j])
    err = err || direct_operator_dft(dft_tmp,dft_in,dx[j],dy[j],nx,ny,Nx,Ny,normalize);

    // apply inverse DFT to dft_tmp and store the resulting signal in out[j]
    ASSERT_ALLOC(iplan_tmp = fftw_plan_dft_2d(ny,nx,dft_tmp,out[j],FFTW_BACKWARD,FFTW_ESTIMATE));
    fftw_execute(iplan_tmp);
    fftw_destroy_plan(iplan_tmp);

  }

  /* free memory */
  fftw_free(dft_in);
  fftw_free(dft_tmp);
  fftw_destroy_plan(plan_in);

  return EXIT_SUCCESS;
}


/******************************************************************************/
/*                                                                            */
/*   MODULES DEDICATED TO THE COMPUTATION OF THE ADJOINT OPERATORS (A*,Aj*)   */
/*                                                                            */
/******************************************************************************/

/* adjoint_operator_dft:

   This module computes dft_out = weight * DFT(Aj*(in)), where "Aj*"
   corresponds to the adjoint of the "Aj" operator described in the
   companion IPOL paper. The computation of DFT(Aj*(in)) is described
   in Equation (28)

   In this module, the input image is given in the frequency domain
   and the output image is computed in the frequency domain as well.

   Notice that, the weight parameter is helpful for computing weighted
   least squares and the iteratively reweighted least-squares (IRLS)
   algorithm.

   Inputs description:
   ===================

   + dx, dy : the horizontal and vertical components of the
              translation vector ('dx' and 'dy' can be any real
              number).

   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + weight : weight parameter

   + dft_in : a ffw_complex array containing the (unnormalized) nx*ny
              DFT coefficients of the input image (ROW-MAJOR order).

   + normalize : a flag that specifies wether dft_out should be
                 normalized (i.e. divided by (Nx*Ny)) or not on
                 exit. Set normalize=0 to disable normalization, or
                 normalized!=0 to enable normalization.

   + dft_out : on entry -> a preallocated a fftw_complex array
               containing Nx*Ny ellements, on exit -> contains the
               (normalized or unnormalized) DFT coefficients of weight
               * Aj*(in).

   Module output
   =============

   This module should return EXIT_SUCCESS.

*/
int adjoint_operator_dft(dft_out,dft_in,dx,dy,weight,nx,ny,Nx,Ny,normalize)
     fftw_complex *dft_out,*dft_in;
     double dx,dy,weight;
     int nx,ny,Nx,Ny;
     char normalize;
{

  int alf,bet,a,b,adr_a,adr_b,adr_alf,adr_bet,adr,adr2;
  double zx,zy,cx,cy,nrm,dft_re,dft_im,phi_re,phi_im,theta,half_nx,half_ny;

  /* initialize some local variables */
  zx = (double)Nx/(double)nx;
  zy = (double)Ny/(double)ny;
  cx = 2.*M_PI*dx/(double)nx;
  cy = 2.*M_PI*dy/(double)ny;
  nrm = (normalize) ? 1./((double)Nx*(double)Ny) : 1.;
  nrm *= weight;
  half_nx = (double)nx/2.;
  half_ny = (double)ny/2.;

  /* main loop: range the high-frequency domain in order to fill 'dft_out' frequency by frequency */
  for(bet=-Ny/2; bet<=(Ny-1)/2; bet++) {
    b = (bet%ny); if(b>=half_ny) b -= ny;   // aliased position of bet in the low-frequency domain
    adr_b = (b+ny)%ny;                      // vertical position of frequency b in 'dft_in'
    adr_bet = (bet+Ny)%Ny;                  // vertical position of frequency bet in 'dft_out'
    for(alf=-Nx/2; alf<=(Nx-1)/2; alf++) {
      a = (alf%nx); if(a>=half_nx) a -= nx; // aliased position of alf in the low-frequency domain
      adr_a = (a+nx)%nx;                    // horizontal position of frequency a in 'dft_in'
      adr_alf = (alf+Nx)%Nx;                // horizontal position of frequency alf in 'dft_out'
      adr = adr_b*nx+adr_a;                 // position of the frequency (a,b) in 'dft_in'
      adr2 = adr_bet*Nx+adr_alf;            // position of the frequency (alf,bet) in 'dft_out'
      dft_re = dft_in[adr][0];              // real part of the DFT coefficient stored in 'dft_in' at location (a,b)
      dft_im = dft_in[adr][1];              // imaginary part of the DFT coefficient stored in 'dft_in' at location (a,b)
      theta = -(cx*alf+cy*bet);
      phi_re = cos(theta);                  // real part of the ramp-phase at location (alf,bet)
      phi_im = sin(theta);                  // imaginary part of the ramp-phase at location (alf,bet)
      dft_out[adr2][0] = nrm * (dft_re*phi_re - dft_im*phi_im);
      dft_out[adr2][1] = nrm * (dft_re*phi_im + dft_im*phi_re);
    }
  }

  return EXIT_SUCCESS;

}


/* adjoint_stack_operator_dft:

   This module applies the A* operator (denoting by A* the adjoint of
   the "A" operator) to a sequence of low-resolution images. We refer
   to Equation (12) of the companion IPOL paper for the formal
   definition of "A")

   In this module, the input stack of (low-resolution) images is given
   in the spatial domain, and the output (high-resolution) image is
   computed in the frequency domain.

   Inputs description:
   ===================

   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + nimages : number of images contained in the input stack of
               low-resolution images (in the IPOL paper, nimages is
               noted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + in : contains the graylevels of the stack of low-resolution input
          images.

   + normalize : a flag that specifies wether dft_out should be
                 normalized (i.e. divided by (Nx*Ny)) or not on
                 exit. Set normalize=0 to disable normalization, or
                 normalized!=0 to enable normalization.

   + dft_out : on entry -> a preallocated a fftw_complex array
               containing Nx*Ny ellements, on exit -> contains the
               (normalized or unnormalized) DFT coefficients A*(in)

   Module output
   =============
   
   This module returns EXIT_SUCCESS if no error occured and
   EXIT_FAILURE otherwise.

   Practical computation:
   ======================

   Using the 'adjoint_operator_dft' module, the present module computes

   DFT(out) = sum_{j=0}^{nimages-1} DFT(Aj*(in[j])),

   which corresponds to the left-side of Equation (28) in the IPOL
   paper.

*/
int adjoint_stack_operator_dft(dft_out,in,dx,dy,nx,ny,Nx,Ny,nimages,normalize)
     fftw_complex *dft_out,**in;
     double *dx,*dy;
     int nx,ny,Nx,Ny,nimages;
     char normalize;
{
  fftw_complex *dft_in,*dft_tmp;
  fftw_plan plan_in;
  int j,adr,err=EXIT_SUCCESS,e;

  /* memory allocation */
  ASSERT_ALLOC(dft_in = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_tmp = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));

  /* CORE : compute sum_{j} DFT(Aj*(in[j])) */
  memset(dft_out,0.,Nx*Ny*sizeof(fftw_complex));
  for(j=0;j<nimages;j++) {

    // compute the unnormalized DFT coefficients of in[j]
    ASSERT_ALLOC(plan_in = fftw_plan_dft_2d(ny,nx,in[j],dft_in,FFTW_FORWARD,FFTW_ESTIMATE));
    fftw_execute(plan_in); // now 'dft_in' contains the DFT coefficients of in[j]
    fftw_destroy_plan(plan_in);

    // compute in dft_tmp the (normalized or unnormalized) DFT coefficients of Aj*(in[j])
    e = adjoint_operator_dft(dft_tmp,dft_in,dx[j],dy[j],1.,nx,ny,Nx,Ny,normalize);
    err = err || e;

    // add dft_tmp to dft_out
    for(adr=0;adr<Nx*Ny;adr++) {
      dft_out[adr][0] += dft_tmp[adr][0];
      dft_out[adr][1] += dft_tmp[adr][1];
    }

  }

  /* free memory */
  fftw_free(dft_in);
  fftw_free(dft_tmp);

  return err;
}


/* adjoint_weighted_stack_operator_dft:

   This module does a similar operation as that provided in the module
   adjoint_stack_operator_dft, but replacing each Aj operator into Aj
   * sqrt(weight[j]) and each low-resolution image in[j] into in[j] *
   sqrt(weight[j]).

   In this module, the input stack of (low-resolution) images is given
   in the frequency domain, and the output (high-resolution) image is
   computed in the frequency domain.

   Inputs description:
   ===================


   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + nimages : number of images contained in the input stack of
               low-resolution images (in the IPOL paper, nimages is
               noted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + weights : (double array containing nimages elements) weighting
               terms

   + dft_in : a fftw_complex array containing the DFT coefficients of
              the low-resolution images (dft_in[j] = DFT(in[j]),
              denoting by in[j] the j-th low resolution image).

   + normalize : a flag that specifies wether dft_out should be
                 normalized (i.e. divided by (Nx*Ny)) or not on
                 exit. Set normalize=0 to disable normalization, or
                 normalized!=0 to enable normalization.

   + dft_out : on entry -> a preallocated a fftw_complex array
               containing Nx*Ny ellements, on exit -> contains the
               (normalized or unnormalized) DFT coefficients of the
               output images (see below)

   Module output
   =============
   
   This module returns EXIT_SUCCESS if no error occured and
   EXIT_FAILURE otherwise.

   Practical computation:
   ======================

   Using the 'adjoint_operator_dft' module, the present module computes

   DFT(out) = sum_{j=0}^{nimages-1} weights[j] * DFT(Aj*(in[j])),

   which is useful to compute weighted least-squares (see Section 6 of
   the companion IPOL paper for more details).

*/
int adjoint_weighted_stack_operator_dft(dft_out,dft_in,dx,dy,weights,nx,ny,Nx,Ny,nimages,normalize)
     fftw_complex *dft_out,**dft_in;
     double *dx,*dy,*weights;
     int nx,ny,Nx,Ny,nimages;
     char normalize;
{
  fftw_complex *dft_tmp;
  int j,adr,err=EXIT_SUCCESS,e;

  /* memory allocation */
  ASSERT_ALLOC(dft_tmp = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));

  /* CORE: sum weights[j] * DFT(Aj*(in)) for j = 1..nimages */
  memset(dft_out,0.,Nx*Ny*sizeof(fftw_complex));
  for(j=0;j<nimages;j++) {
    if(weights[j] != 0.) {

      // compute in dft_tmp the (normalized or unnormalized) DFT coefficients of weights[j] * Aj*(in[j])
      e = adjoint_operator_dft(dft_tmp,dft_in[j],dx[j],dy[j],weights[j],nx,ny,Nx,Ny,normalize);
      err = err || e;

      // add dft_tmp to dft_out
      for(adr=0;adr<Nx*Ny;adr++) {
	dft_out[adr][0] += dft_tmp[adr][0];
	dft_out[adr][1] += dft_tmp[adr][1];
      }

    }
  }

  /* free memory */
  fftw_free(dft_tmp);

  return err;
}
