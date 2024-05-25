#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <lapacke.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }
#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

/* internal structures */
typedef struct valueandindex {
  double value;   // value
  int j;          // index
} *Valueandindex;

/* internal modules */
static int weightedleastsquares_dft(fftw_complex*,fftw_complex*,double*,double*,double*,int,int,int,int,int,double,char,char); 
int leastsquares(fftw_complex*,fftw_complex**,double*,double*,int,int,int,int,int,double,char); 
int irls(fftw_complex*,double*,fftw_complex**,double*,double*,int,int,int,int,int,double,double,char,int,double); 
int luckyimaging(fftw_complex*,fftw_complex**,double*,fftw_complex**,double*,double*,int,int,int,int,int,int,int,int,int,int*,double,char);

/* external modules */
extern int adjoint_stack_operator_dft(fftw_complex*,fftw_complex**,double*,double*,int,int,int,int,int,char); // see source file 'operators_kernel.c'
extern int adjoint_weighted_stack_operator_dft(fftw_complex*,fftw_complex**,double*,double*,double*,int,int,int,int,int,char); // see source file 'operators_kernel.c'
extern int direct_operator_dft(fftw_complex*,fftw_complex*,double,double,int,int,int,int,char); // see source file 'operators_kernel.c'
extern int pseudoinverse_rowmajor(lapack_complex_double*,lapack_int,lapack_int,double,double*); // see source file 'blockmatrix_kernel.c'
extern int compute_blockmatrix(lapack_complex_double*,double*,double*,double*,double,int,int,int); // see source file 'blockmatrix_kernel.c'
extern int cmvprod(int,int,lapack_complex_double*,lapack_complex_double*,lapack_complex_double*); // see source file 'blockmatrix_kernel.c'

/* weightedleastsquares_dft: super-resolution using the least-squares
                             (Fourier domain).

   This module computes DFT(u_ls) where u_ls denotes the weighted
   least-squares reconstruction from the stack of low-resolution image
   u0, i.e.,

   u_ls = argmin_{u} sum_{j=0..nimages-1} weights[j]*||Aj(u)-u0^{(j)}||^2,
   
   denoting Aj the operator defined in Equation (11) of the companion
   IPOL paper.

   In this module, the computation of u_ls is done from the quantity

   dft_v = sum_{j=0..nimages-1} weights[j] * DFT(Aj*(u0^{(j)})),
   denoting by Aj* the adjoint of the Aj operator.

   Notice that dft_v can be computed from the low-resolution stack u0
   using the 'adjoint_weighted_stack_operator_dft' module (see source
   file 'operator_kernel.c').

   This module is closely related to the pseudocode Algorithm 2 of the
   companion IPOL paper. The main differences with the pseudocode
   Algorithm 2 and this module are the following:

   + this module returns the DFT coeffcients of the u_ls image
   (instead of u_ls itslef)

   + this module is a bit more general due to the use of the weighting
   terms (the pseudocode Algorithm 2 of the IPOL paper corresponds to
   the particular case where weights[j] = 1 for j=0..nimages-1).

   Module parameters
   =================

   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + nimages : number of low resolution images (in the IPOL paper,
               nimages is noted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + weights : a pointer containing 'nimages' elements in double
               precision, setting weights=NULL is equivalent to set
               weights[j] = 1 for j=0..nimages-1

   + eps : a double number used as a threshold for the singular values
           of A in the pseudo inversion process (only the singular
           values > eps will be inverted)
   
   + normalize : specifies wether output dft_out should be normalized
                 (i.e. divided by (Nx*Ny)) or not. Set normalize=0 to
                 disable normalization, or normalized!=0 to enable
                 normalization.

   + vflag : set vflag!=0 or vflag=0 to enable or disable the verbose mode

   + dft_v : a fftw_complex array containing Nx*Ny DFT coefficient
             (see description above)

   + dft_out: the (normalized or unnormalized) DFT coefficients of the
              output high-resolution image u_ls.

   Module Output 
   =============

   This module should return EXIT_SUCCESS if no error occured during
   its execution, or EXIT_FAILURE otherwise.

*/
static int weightedleastsquares_dft(fftw_complex *dft_out, fftw_complex *dft_v, double *weights, double *dx, double *dy, int nx, int ny, int Nx, int Ny, int nimages, double eps, char normalize, char vflag)
{
  lapack_complex_double *M1pinv=NULL,*M2pinv=NULL,*M3pinv=NULL,*M4pinv=NULL,*Mpinv=NULL,*left,*right;
  lapack_int Z1,Z2,Z3,Z4,Z;
  double zx,zy,re,im,nrm,cof,cond1,cond2,cond3,cond4;
  int a,b,alf,bet,p,q,adr_alf,adr_bet,adr,adr2,pmin,qmin,pmax,qmax,size_p,size_q,nalloc;
  char zxisinteger,zyisinteger; 

  // initialize some local variables //
  zx = (double)Nx/(double)nx; 
  zy = (double)Ny/(double)ny; 
  cof = 1./(zx*zy);
  nrm = (normalize) ? 1./((double)Nx*(double)Ny) : 1.;
  Z1 = (lapack_int) (floor(zx)*floor(zy)); // Z1 = number of lines = number of columns of the square block matrix M1 (and its pseudoinverse)
  Z2 = (lapack_int) (ceil(zx)*floor(zy));  // Z2 = number of lines = number of columns of the square block matrix M2 (and its pseudoinverse)
  Z3 = (lapack_int) (floor(zx)*ceil(zy));  // Z3 = number of lines = number of columns of the square block matrix M3 (and its pseudoinverse)
  Z4 = (lapack_int) (ceil(zx)*ceil(zy));   // Z4 = number of lines = number of columns of the square block matrix M4 (and its pseudoinverse)
  zxisinteger = (zx == floor(zx)); 
  zyisinteger = (zy == floor(zy));
  
  // memory allocation //
  ASSERT_ALLOC(left = (lapack_complex_double*) LAPACKE_malloc (Z4*sizeof(lapack_complex_double))); // remark that Z4 = max_{1 <= i <= 4} Zi 
  ASSERT_ALLOC(right = (lapack_complex_double*) LAPACKE_malloc (Z4*sizeof(lapack_complex_double))); // remark that Z4 = max_{1 <= i <= 4} Zi

  // compute the pseudo-inverse of the block matrix M1 //
  ASSERT_ALLOC(M1pinv = (lapack_complex_double*) LAPACKE_malloc (Z1*Z1*sizeof(lapack_complex_double)));
  if(EXIT_SUCCESS != compute_blockmatrix(M1pinv,dx,dy,weights,cof,(int)floor(zx),(int)floor(zy),nimages)) { // set M1pinv = M1 (block matrix is not pseudo-inverted yet)
    printf("An error occurred during the computation of the block matrix M1.");
    LAPACKE_free(M1pinv); LAPACKE_free(left); LAPACKE_free(right);
    return EXIT_FAILURE;
  }
  if(EXIT_SUCCESS != pseudoinverse_rowmajor(M1pinv,Z1,Z1,eps,&cond1)) { // set M1pinv = pseudo inverse of the block matrix M1
    printf("An error occurred during the pseudo-inversion of the block matrix M1.");
    LAPACKE_free(M1pinv); LAPACKE_free(left); LAPACKE_free(right);
  }
  if(vflag) printf("compute and pseudo-invert the block matrix M1: size = %d x %d (condition number in l2-norm = %g)\n",Z1,Z1,cond1);

  // compute the pseudo-inverse of the other block matrices (M2, M3 and M4) //
  if (zxisinteger && zyisinteger) { // zx and zy are both integers (M4 = M3 = M2 = M1)
    M2pinv = M3pinv = M4pinv = M1pinv; 
    if(vflag) printf("other block matrices: we have M4 = M3 = M2 = M1\n\n");
    nalloc = 1; 
  }
  else { // at least one super-resolution factor is not an integer

    // compute the pseudo-inverse of M4
    ASSERT_ALLOC(M4pinv = (lapack_complex_double*) LAPACKE_malloc (Z4*Z4*sizeof(lapack_complex_double)));
    if(EXIT_SUCCESS != compute_blockmatrix(M4pinv,dx,dy,weights,cof,(int)ceil(zx),(int)ceil(zy),nimages)) { // set M4pinv = M4 (block matrix is not pseudo-inverted yet)
      printf("An error occurred during the computation of the block matrix M4.");
      LAPACKE_free(left); LAPACKE_free(right);
      if(NULL != M1pinv) LAPACKE_free(M1pinv); 
      if(NULL != M2pinv) LAPACKE_free(M2pinv); 
      if(NULL != M3pinv) LAPACKE_free(M3pinv); 
      if(NULL != M4pinv) LAPACKE_free(M4pinv); 
      return EXIT_FAILURE;
    }
    if(EXIT_SUCCESS != pseudoinverse_rowmajor(M4pinv,Z4,Z4,eps,&cond4)) { // set M4pinv = pseudo inverse of the block matrix M4
      printf("An error occurred during the pseudo-inversion of the block matrix M4.");
      LAPACKE_free(left); LAPACKE_free(right);
      if(NULL != M1pinv) LAPACKE_free(M1pinv); 
      if(NULL != M2pinv) LAPACKE_free(M2pinv); 
      if(NULL != M3pinv) LAPACKE_free(M3pinv); 
      if(NULL != M4pinv) LAPACKE_free(M4pinv); 
    }
    
    if(zxisinteger) { // zx is an integer but zy is not an integer (M2 = M1 and M3 = M4)
      M2pinv = M1pinv;
      M3pinv = M4pinv;
      if(vflag) {
	printf("compute and pseudo-invert the block matrix M4: size = %d x %d (condition number in l2-norm = %g)\n",Z4,Z4,cond4);
	printf("other block matrices: we have M2 = M1 and M3 = M4\n\n");
      }
      nalloc = 2; 
    }    
    else if(zyisinteger) { // zx is not an integer but zy is an integer (M2 = M4 and M3 = M1)
      M2pinv = M4pinv;
      M3pinv = M1pinv;
      if(vflag) {
	printf("compute and pseudo-invert the block matrix M4: size = %d x %d (condition number in l2-norm = %g)\n",Z4,Z4,cond4);
	printf("other block matrices: we have M2 = M4 and M3 = M1\n\n");
      }
      nalloc = 2; 
    }
    else { // zx and zy are both not integer, the computation of the pseudo-inverses of M2 and M3 is needed

      // compute the pseudo-inverse of M2
      ASSERT_ALLOC(M2pinv = (lapack_complex_double*) LAPACKE_malloc (Z2*Z2*sizeof(lapack_complex_double)));
      if(EXIT_SUCCESS != compute_blockmatrix(M2pinv,dx,dy,weights,cof,(int)ceil(zx),(int)floor(zy),nimages)) { // set M2pinv = M2 (block matrix is not pseudo-inverted yet)
	printf("An error occurred during the computation of the block matrix M2.");
	LAPACKE_free(left); LAPACKE_free(right);
	if(NULL != M1pinv) LAPACKE_free(M1pinv); 
	if(NULL != M2pinv) LAPACKE_free(M2pinv); 
	if(NULL != M3pinv) LAPACKE_free(M3pinv); 
	if(NULL != M4pinv) LAPACKE_free(M4pinv); 
	return EXIT_FAILURE;
      }
      if(EXIT_SUCCESS != pseudoinverse_rowmajor(M2pinv,Z2,Z2,eps,&cond2)) { // set M2pinv = pseudo inverse of the block matrix M2
	printf("An error occurred during the pseudo-inversion of the block matrix M2.");
	LAPACKE_free(left); LAPACKE_free(right);
	if(NULL != M1pinv) LAPACKE_free(M1pinv); 
	if(NULL != M2pinv) LAPACKE_free(M2pinv); 
	if(NULL != M3pinv) LAPACKE_free(M3pinv); 
	if(NULL != M4pinv) LAPACKE_free(M4pinv); 
      }
      
      // compute the pseudo-inverse of M3
      ASSERT_ALLOC(M3pinv = (lapack_complex_double*) LAPACKE_malloc (Z3*Z3*sizeof(lapack_complex_double)));
      if(EXIT_SUCCESS != compute_blockmatrix(M3pinv,dx,dy,weights,cof,(int)floor(zx),(int)ceil(zy),nimages)) { // set M3pinv = M3 (block matrix is not pseudo-inverted yet)
	printf("An error occurred during the computation of the block matrix M3.");
	LAPACKE_free(left); LAPACKE_free(right);
	if(NULL != M1pinv) LAPACKE_free(M1pinv); 
	if(NULL != M2pinv) LAPACKE_free(M2pinv); 
	if(NULL != M3pinv) LAPACKE_free(M3pinv); 
	if(NULL != M4pinv) LAPACKE_free(M4pinv); 
	return EXIT_FAILURE;
      }
      if(EXIT_SUCCESS != pseudoinverse_rowmajor(M3pinv,Z3,Z3,eps,&cond3)) { // set M3pinv = pseudo inverse of the block matrix M3
	printf("An error occurred during the pseudo-inversion of the block matrix M3.");
	LAPACKE_free(left); LAPACKE_free(right);
	if(NULL != M1pinv) LAPACKE_free(M1pinv); 
	if(NULL != M2pinv) LAPACKE_free(M2pinv); 
	if(NULL != M3pinv) LAPACKE_free(M3pinv); 
	if(NULL != M4pinv) LAPACKE_free(M4pinv); 
      }
      
      if(vflag) {
	printf("compute and pseudo-invert the block matrix M2: size = %d x %d (condition number in l2-norm = %g)\n",Z2,Z2,cond2);
	printf("compute and pseudo-invert the block matrix M3: size = %d x %d (condition number in l2-norm = %g)\n",Z3,Z3,cond3);
	printf("compute and pseudo-invert the block matrix M4: size = %d x %d (condition number in l2-norm = %g)\n\n",Z4,Z4,cond4);
      }
      nalloc = 4; 
      
    }
  }
  
  // main loop //
  for(b=-ny/2;b<=(ny-1)/2;b++) {
    
    qmin = (int)ceil(-(Ny+2.*(double)b)/(2.*(double)ny));  // first element of the index set Q(b)
    qmax = (int)ceil((Ny-2.*(double)b)/(2.*(double)ny))-1; // last element of the index set Q(b)
    size_q = qmax-qmin+1; // number of elements in Q(b)

    for(a=-nx/2;a<=(nx-1)/2;a++) {

      pmin = (int)ceil(-(Nx+2.*(double)a)/(2.*(double)nx));  // first element of the index set P(a)
      pmax = (int)ceil((Nx-2.*(double)a)/(2.*(double)nx))-1; // last element of the index set P(a)
      size_p = pmax-pmin+1; // number of elements in P(a)
      Z = (lapack_int) (size_p * size_q); // order of the block matrix M(a,b)

      // range the set of frequencies (alf,bet) that are aliased at
      // position (a,b) in the low-frequency domain and compute the
      // right-hand term of the linear system (notice that we range
      // P(a) x Q(b) in lexicographical order)
      for(adr=0,p=pmin;p<=pmax;p++) {
	alf = a + p*nx;
	adr_alf = (Nx+alf)%Nx; // horizontal position of frequency alf in 'dft_out' and 'dft_v'
	for(q=qmin;q<=qmax;q++,adr++) {
	  bet = b + q*ny;
	  adr_bet = (Ny+bet)%Ny; // vertical position of frequency bet in 'dft_out' and 'dft_v'
	  adr2 = adr_bet*Nx + adr_alf; // position of frequency (alf,bet) in 'dft_out' and 'dft_v'
	  re = dft_v[adr2][0];
	  im = dft_v[adr2][1];
	  right[adr] = lapack_make_complex_double(re,im);
	}
      }

      // set Mpinv = pseudo inverse of M(a,b)
      if ((size_p == (int)floor(zx)) && (size_q == (int)floor(zy))) Mpinv = M1pinv; 
      else if ((size_p == (int)ceil(zx)) && (size_q == (int)floor(zy))) Mpinv = M2pinv;
      else if ((size_p == (int)floor(zx)) && (size_q == (int)ceil(zy))) Mpinv = M3pinv;
      else Mpinv = M4pinv;
      
      // Apply the pseudoinverse to the right-hand term of the linear
      // system (complex matrix vector multiplication left =
      // Mpinv*right)
      cmvprod(Z,Z,Mpinv,right,left); // the same operation can be done using CBLAS_zgemv module using the three next commented lines below.
      //lapack_complex_double one = lapack_make_complex_double(1.,0.);
      //lapack_complex_double zero = lapack_make_complex_double(0.,0.);
      //cblas_zgemv(CblasRowMajor,CblasNoTrans,Z,Z,&one,Mpinv,Z,right,1,&zero,left,1);
      
      // replace frequencies coefficients at their correct positions
      // in dft_out (notice that we range P(a) x Q(b) in
      // lexicographical order again)
      for(adr=0,p=pmin;p<=pmax;p++) {
	alf = a + p*nx;
	adr_alf = (Nx+alf)%Nx; // horizontal position of frequency alf in 'dft_out' and 'dft_v'
	for(q=qmin;q<=qmax;q++,adr++) {
	  bet = b + q*ny;
	  adr_bet = (Ny+bet)%Ny; // vertical position of frequency bet in 'dft_out' and 'dft_v'
	  adr2 = adr_bet*Nx + adr_alf; // position of frequency (alf,bet) in 'dft_out' and 'dft_v'
	  dft_out[adr2][0] = nrm * lapack_complex_double_real(left[adr]);
	  dft_out[adr2][1] = nrm * lapack_complex_double_imag(left[adr]);
	}
      }

    }
  }

  // free memory //
  LAPACKE_free(left);
  LAPACKE_free(right);
  LAPACKE_free(M1pinv);
  if(nalloc >= 2) LAPACKE_free(M4pinv);
  if(nalloc == 4) {
    LAPACKE_free(M2pinv);
    LAPACKE_free(M3pinv);
  }

  return EXIT_SUCCESS;
}


/* leastsquares: super-resolution using the least-squares.

   This module implements the pseudocode Algorithm 2 of the companion
   IPOL paper.

   Module parameters
   =================
   
   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + nimages : number of low resolution images (in the IPOL paper,
               nimages is noted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + eps : a double number used as a threshold for the singular values
           of A in the pseudo inversion process (only the singular
           values > eps will be inverted)
   
   + vflag : set vflag!=0 or vflag=0 to enable or disable the verbose mode

   + u0: input stack of low-resolution images, such as u0[i]
         corresponds to the graylevel values of the i-th
         low-resolution image of the stack (0 <= i < nimages), stored
         in row-major order

   + out : on entry, a preallocated fftw_complex array with size
           Nx*Ny, on exit, graylevels of the output high_resolution
           image

   Module Output 
   =============

   This module should return EXIT_SUCCESS if no error occured during
   its execution, or EXIT_FAILURE otherwise.

*/
int leastsquares(fftw_complex *out, fftw_complex **u0, double *dx, double *dy, int nx, int ny, int Nx, int Ny, int nimages, double eps, char vflag)
{
  fftw_complex *dft_out,*dft_Asu0;
  fftw_plan iplan_out;
  char normalize = 1;
  int err; 
  
  /* memory allocation */
  ASSERT_ALLOC(dft_Asu0 = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_out = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(iplan_out = fftw_plan_dft_2d(Ny,Nx,dft_out,out,FFTW_BACKWARD,FFTW_ESTIMATE));
  
  /* fill dft_Asu0 with the unnormalized DFT coefficients of A*(u0) */
  if(EXIT_SUCCESS != adjoint_stack_operator_dft(dft_Asu0,u0,dx,dy,nx,ny,Nx,Ny,nimages,(char)0)) {
    printf("An error occurred when executing the internal 'adjoint_stack_operator_dft' module.\n");
    fftw_free(dft_Asu0);
    fftw_free(dft_out); 
    fftw_destroy_plan(iplan_out);
    return EXIT_FAILURE;
  }
  
  /* fill dft_out with the normalized DFT coefficients of the output
     (high-resolution) image and compute inverse DFT to get the output
     image in the spatial domain */
  if(EXIT_SUCCESS != weightedleastsquares_dft(dft_out,dft_Asu0,NULL,dx,dy,nx,ny,Nx,Ny,nimages,eps,normalize,vflag)) { // calling module with weights == NULL to compute the standard (non-weighted) least-squares estimator
    fftw_free(dft_Asu0);
    fftw_free(dft_out); 
    fftw_destroy_plan(iplan_out);
    return EXIT_FAILURE;
  }
  fftw_execute(iplan_out);

  /* free memory */
  fftw_free(dft_Asu0);
  fftw_free(dft_out); 
  fftw_destroy_plan(iplan_out);

  return EXIT_SUCCESS;
}


/* irls: iteratively reweighted least-squares.
   
   General description
   ===================

   This modules iterates the IRLS scheme corresponding to Equation
   (56) of the IPOL companion paper, i.e., for k >= 0, compute 

   u^{k+1} = argmin_{u} H(u,eta^{k})
   eta[j]^{k+1} = max(eps,||Aj(u) - u0^{(j)}||) for j=0..nimages-1

   denoting 

   H(u,eta) = sum_{j=0..nimages-1} ||Aj(u) - u0^{(j)}||^2 / (2*eta[j]) + eta[j]/2

   Notice that in this module we have set weights[j] = 1/eta[j] for
   j=0..nimages-1.

   The IRLS sheme described above intents to compute a minimizer
   u_l1l2 of the l1-l2 energy defined below,

   u_l1l2 = argmin_{u} E_{l1l2}(u) := sum_{j=0..nimages-1} ||Aj(u)-u0^{(j)}||^2.

   The iterations of this scheme are stopped when the maximal number
   of iteration (niter) is attained, or when 

   |E_l1l2(u^{k}) - E_l1l2(u^{k-1})| <= r*|E_l1l2(u^{k-1})|

   See Section 6 of the companion IPOL paper for more details.
   

   Module parameters
   =================
   
   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + nimages : number of low resolution images (in the IPOL paper,
               nimages is noted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + eps_pinv : a double number used as a threshold for the singular
                values of A in the pseudo inversion process (only the
                singular values > eps will be inverted)
   
   + vflag : set vflag!=0 or vflag=0 to enable or disable the verbose mode

   + u0: input stack of low-resolution images, such as u0[i]
         corresponds to the graylevel values of the i-th
         low-resolution image of the stack (0 <= i < nimages), stored
         in row-major order

   + eps : minimal value for the entries eta[j]

   + niter : maximal number of iterations

   + r : stop the iterations when |E_l1l2(u^{k}) - E_l1l2(u^{k-1})| <= 
         r*|E_l1l2(u^{k-1})|

   + weights : on entry, a preallocated double array containing
               nimages elements, on exit, output weights (=1/eta^{k})
               at the end of the iterations

   + out : on entry, a preallocated fftw_complex array with size
           Nx*Ny, on exit, graylevels of the output high_resolution
           image (noted u_l1l2 above)

   Module Output 
   =============

   This module should return EXIT_SUCCESS if no error occured during
   its execution, or EXIT_FAILURE otherwise.

*/
int irls(fftw_complex *out, double *weights, fftw_complex **u0, double *dx, double *dy, int nx, int ny, int Nx, int Ny, int nimages, double eps, double eps_pinv, char vflag, int niter, double r)
{
  fftw_complex **dft_u0=NULL,*dft_u=NULL,*dft_Aju=NULL,*dft_v=NULL;
  fftw_plan plan_u0,iplan_out;
  double z,E,E_prev,re,im;
  int iter,err,j,adr;
  char stop,normalize;

  /* memory allocation */
  ASSERT_ALLOC(dft_u = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_v = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_Aju = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_u0 = (fftw_complex**) malloc (nimages*sizeof(fftw_complex*)));
  for(j=0;j<nimages;j++) {
    ASSERT_ALLOC(dft_u0[j] = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  }

  /* initialization */
  for(j=0;j<nimages;j++) {
    ASSERT_ALLOC(plan_u0 = fftw_plan_dft_2d(ny,nx,u0[j],dft_u0[j],FFTW_FORWARD,FFTW_ESTIMATE));
    fftw_execute(plan_u0); // now dft_u0[j] contains the DFT coefficients of u0[j]
    fftw_destroy_plan(plan_u0);
    weights[j] = 1.;
  }
  iter = 1; stop = 0; normalize = 0; 

  /*** main loop ***/
  while(!stop) {

    // update u in the frequency domain (i.e., update dft_u = DFT(u)) //
    adjoint_weighted_stack_operator_dft(dft_v,dft_u0,dx,dy,weights,nx,ny,Nx,Ny,nimages,(char)0); // compute adjoint of the weighted stack operator
    weightedleastsquares_dft(dft_u,dft_v,weights,dx,dy,nx,ny,Nx,Ny,nimages,eps_pinv,normalize,(char)0); 
    
    // update z and compute the l1l2 energy E //
    E_prev = E; E = 0.;
    for(j=0;j<nimages;j++){

      // set dft_Aju = DFT(Aj(u)) //
      err = direct_operator_dft(dft_Aju,dft_u,dx[j],dy[j],nx,ny,Nx,Ny,normalize);

      // compute weights[j] = 1/max(eps,||Aj(u)-u0^{(j)}||) denoting by u0^{(j)} the j-th image of the input stack u0 //
      for(z=0.,adr=0;adr<nx*ny;adr++) {
	re = dft_Aju[adr][0] - dft_u0[j][adr][0];
	im = dft_Aju[adr][1] - dft_u0[j][adr][1];
	z += re*re+im*im;
      }
      z = sqrt(z/(nx*ny));
      weights[j] = 1. / max(eps,z); 
      E += z;

    }

    // deal with verbose mode //
    if(vflag) printf("iteration %-2d: l1-l2 energy = %-.17e\n",iter,E);

    // check stopping criterion //
    stop = (iter >= niter) || ((iter > 1) && (fabs(E-E_prev) <= r*fabs(E_prev)));
    iter++;

  }
  //if(vflag) printf("regularization parameter: eps = %-.5e\n",*eps);

  /* if needed, compute output image in the spatial domain (manual renormalization needed) */
  if(NULL != out) {
    ASSERT_ALLOC(iplan_out = fftw_plan_dft_2d(Ny,Nx,dft_u,out,FFTW_BACKWARD,FFTW_ESTIMATE));
    fftw_execute(iplan_out);
    fftw_destroy_plan(iplan_out);
    for(adr=0;adr<Nx*Ny;adr++) {
      out[adr][0] /= (double)Nx*(double)Ny;
      out[adr][1] /= (double)Nx*(double)Ny;
    }
  }

  /* free memory */
  for(j=0;j<nimages;j++) fftw_free(dft_u0[j]);
  free(dft_u0);
  fftw_free(dft_u);
  fftw_free(dft_Aju);
  fftw_free(dft_v); 

  return EXIT_SUCCESS;
}


/* comparison function for qsort (sort the structures in descending order of their 'value' attribute) */
static int compare_valueandindex(const void *a,const void *b)
{
  struct valueandindex const *pa = a;
  struct valueandindex const *pb = b;
  double diff = pa->value - pb->value;
  return ((diff < 0) ? 1 : ((diff > 0) ? -1 : 0));
}



/* luckyimaging: 

   This modules implements the luckyimaging procedure proposed in
   Section 6 of the companion IPOL paper. This procedure can be
   summarized as follows.

   Given a sequence 'u0' containing 'nimages' low-resolution images,
   and a sequence 'weights' containing 'nimages' real numbers (in
   practice this weights sequence is obtained by iterating the IRLS
   scheme implemented in the 'irls' module), do:

   + keep from 'u0' the 'n' images associated to the 'n' highest
     values in the 'weights' sequence

   + compute a super-resolved image 'u_{lucky}^{n}' using the
     least-square estimator over this restricted sequence of
     low-resolution images.

   In practice, we do not explicitly remove any image from 'u0' but we
   ignore them in the leastsquares reconstruction procedure: we use
   the 'weigthedleastsquares_dft' module wich computes the DFT of the
   weighted least squares estimator, setting a weight equal to 1 for
   images to be kept and setting a weight to 0 for images to be
   ignored, in the weighted least-squares procedure.

   See section 6 of the companion IPOL paper for more details about
   the luckyimaging procedure.

   This module is able to compute a sequence containing l images
   u_lucky^{n_i} (for l differents values of n_i between n0 and n1) or
   a single image u_lucky^{N} (with N given as input).

   In case of the computation of a sequence of images u_lucky^{n_i}, we use 

   n_i = round(a*i+b) with a = (n0-n1)/(l-1) and b = n0-a*l

   so that (n_i) is a sequence of l quasi-regularly spaced integers
   bewteen n0 and n1 in decreasing order.

   Module parameters 
   =================

   + Nx, Ny : horizontal and vertical dimensions (width and height) of
              the high-resolution image (in the IPOL paper, Nx is
              noted M and Ny is noted N)

   + nx, ny : horizontal and vertical dimensions (width and height) of
              the low-resolution images (in the IPOL paper, nx is
              noted m and ny is noted n)

   + nimages : number of low resolution images (in the IPOL paper,
               nimages is noted L)

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + eps_pinv : a double number used as a threshold for the singular
                values of A in the pseudo inversion process (only the
                singular values > eps will be inverted)
   
   + vflag : set vflag!=0 or vflag=0 to enable or disable the verbose mode

   + u0: input stack of low-resolution images, such as u0[i]
         corresponds to the graylevel values of the i-th
         low-resolution image of the stack (0 <= i < nimages), stored
         in row-major order

   + N: number of images to keep from the stack u0 (will be used only
        when mov == NULL)

   + n0: minimal value for the n_i (will be used only when mov != NULL) 

   + n1: maximal value for the n_i (will be used only when mov != NULL) 

   + l : number of n_i (will be used only when mov != NULL) 
   
   + weights: input sequence of weights containing 'nimages' real
              numbers (for j=0..nimages-1, weights[j] is associated to
              j-th low resolution image of the stack 'u0')

   + out: if mov == NULL, then, on entry, out should be a preallocated
          fftw_complex array with size Nx*Ny, and on exit out =
          u_{lucky}^{N} (if mov != NULL, out will be unchanged)
   
   + mov: if mov != NULL, then, on entry, mov should be a preallocated
          sequence containing l fftw_complex arrays with size Nx*Ny
          each, and on exit, mov[i-1] contains u_lucky^{n_i} for (1 <=
          i <= l)

   + nlist : if mov != NULL, then on entry nlist must be a
             preallocated integer array with size l, and on exit
             n_list[i-1] = n_i for 1 <= i <= l. Otherwise, when mov ==
             NULL, nlist is unchanged (setting nlist = NULL is
             allowed)
 
   WARNING
   =======
 
   notice that the values stored in 'weights' will be changed during the execution
   of the module (this input is not preserved in this implementation). 
   
   Module Output 
   =============

   This module should return EXIT_SUCCESS if no error occured during
   its execution, or EXIT_FAILURE otherwise.

*/
int luckyimaging(fftw_complex *out, fftw_complex **mov, double *weights, fftw_complex **u0, double *dx, double *dy, int nx, int ny, int Nx, int Ny, int nimages, int N, int n0, int n1, int l, int *nlist, double eps_pinv, char vflag)
{
  fftw_complex **dft_u0=NULL,*dft_out=NULL,*dft_v=NULL;
  fftw_plan plan_u0,iplan_out;
  int i,j,err,n,nold;
  double a,b;
  struct valueandindex *w;

  /* memory allocation */
  ASSERT_ALLOC(w = (struct valueandindex *) malloc (nimages*sizeof(struct valueandindex)));
  ASSERT_ALLOC(dft_out = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_v = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_u0 = (fftw_complex**) malloc (nimages*sizeof(fftw_complex*)));
  for(j=0;j<nimages;j++) {
    ASSERT_ALLOC(dft_u0[j] = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  }

  /* initialization */
  for(j=0;j<nimages;j++) {
    ASSERT_ALLOC(plan_u0 = fftw_plan_dft_2d(ny,nx,u0[j],dft_u0[j],FFTW_FORWARD,FFTW_ESTIMATE));
    fftw_execute(plan_u0); // now dft_u0[j] contains the DFT coefficients of u0[j]
    fftw_destroy_plan(plan_u0);
  }

  /* recopy the weights into w and sort them in ascending order */
  for(j=0;j<nimages;j++) {
    w[j].value = weights[j];
    w[j].j = j;
  }
  qsort((void*)w,nimages,sizeof(struct valueandindex),compare_valueandindex); // sort w in descending order

  /* lucky imaging procedure */
  if(!mov) {

    // set to 1 the N heighest weights, set to 0 the others
    for(j=0;j<nimages;j++) weights[w[j].j] = (j<N) ? 1. : 0.;

    // compute u_lucky^{(N)} using weighted least-squares with binary
    // weights (amounts to compute the unweighted least-squares after
    // removing from u0 the images associated to a zero-valued weight)
    adjoint_weighted_stack_operator_dft(dft_v,dft_u0,dx,dy,weights,nx,ny,Nx,Ny,nimages,(char)0); // compute adjoint of the weighted stack operator
    weightedleastsquares_dft(dft_out,dft_v,weights,dx,dy,nx,ny,Nx,Ny,nimages,eps_pinv,(char)1,(char)0); 
    ASSERT_ALLOC(iplan_out = fftw_plan_dft_2d(Ny,Nx,dft_out,out,FFTW_BACKWARD,FFTW_ESTIMATE));
    fftw_execute(iplan_out);
    fftw_destroy_plan(iplan_out);

    // deal with verbose mode
    if(vflag) printf("computed u_lucky^{(%d)} by removing %d images from the input low-resolution stack\n",N,nimages-N);

  }
  else {    

    // set to 1 all weights
    for(j=0;j<nimages;j++) weights[j] = 1.; 

    // compute (a,b) used to define the nj below
    a = (double)(n0-n1)/(double)(l-1);
    b = n0-a*l;    
    nold = nimages;
    
    // compute u_lucky^{ni} for i = 1 .. l
    for(i=1;i<=l;i++) {

      // compute ni = number of image to keep in the lucky-imaging procedure
      nlist[i-1] = n = (int)round(a*i+b);
      
      // set to 0 the weights with sorted index between n and nold
      for(j=n;j<nold;j++) weights[w[j].j] = 0.; 

      // compute u_lucky^{(ni)} using weighted least-squares with
      // binary weights (amounts to compute the unweighted
      // least-squares after removing from u0 the images associated to
      // a zero-valued weight)
      adjoint_weighted_stack_operator_dft(dft_v,dft_u0,dx,dy,weights,nx,ny,Nx,Ny,nimages,(char)0); // compute adjoint of the weighted stack operator
      weightedleastsquares_dft(dft_out,dft_v,weights,dx,dy,nx,ny,Nx,Ny,nimages,eps_pinv,(char)1,(char)0); 
      ASSERT_ALLOC(iplan_out = fftw_plan_dft_2d(Ny,Nx,dft_out,mov[i-1],FFTW_BACKWARD,FFTW_ESTIMATE));
      fftw_execute(iplan_out);
      fftw_destroy_plan(iplan_out);

      // update nold
      nold = n; 

      // deal with verbose mode
      if(vflag) printf("frame #%d/%d : computed u_lucky^{(%d)} by removing %d images from the input low-resolution stack\n",i,l,n,nimages-n);

    }

  }

  /* free memory */
  for(j=0;j<nimages;j++) fftw_free(dft_u0[j]);
  fftw_free(dft_out);
  free(dft_u0);
  free(w);

  return EXIT_SUCCESS;
}


