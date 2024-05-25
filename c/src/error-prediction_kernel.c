#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* internal modules */
int error_prediction(double*,double*,double*,double*,double*,double*,int,int,int,int,int,double,double,char,char); 

/* external modules */
extern int pseudoinverse_rowmajor(lapack_complex_double*,lapack_int,lapack_int,double,double*); // see source file 'blockmatrix_kernel.c'
extern int compute_blockmatrix(lapack_complex_double*,double*,double*,double*,double,int,int,int); // see source file 'blockmatrix_kernel.c'
extern int cmvprod(int,int,lapack_complex_double*,lapack_complex_double*,lapack_complex_double*); // see source file 'blockmatrix_kernel.c'

/* error_prediction: computes the error amplification factor (Fourier
                     domain) defined in Equation (41) of the IPOL
                     paper

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

   + sigma : double pointer, *sigma = noise level (standard deviation)
             of the low-resolution sequence (set sigma = NULL if
             unknown)

   + peakval : peak-value for the PSNR (this parameter will be ignored
               when psnr = NULL on entry)

   + eps : threshold for the singular values (will be passed as input
           for the pseudoinverse_rowmajor module only)

   + vflag : a flag to enable/disable the verbose mode (set vflag=0 to
             disable the verbose mode, or vflag!=0 to enable the
             verbose mode)

   + noshift : set noshift = 1 to not shift output A of half the image
               size (put the (0,0) frequency at the top-left corner of
               A instead of its center)

   + A : on entry -> a preallocated double array containing Nx*Ny
         elements *initialized with zero*, on exit -> the error
         amplification map defined in Equation (41) of the IPOL paper

   + mse : (double pointer) predicted MSE (this value is computed only
           when input sigma is given, i.e., when sigma != NULL), on
           exit -> *mse = sigma^2 * ||A||^2 / (Nx*Ny) (see Equation
           (45) of the IPOL paper)
	   
   + psnr : (double pointer) predicted PSNR (this value is computed
            only when input sigma is given, i.e., when sigma != NULL),
            on exit -> *psnr = predicted PSNR value =
            10*log10(peakval^2/*mse) (see Equation (45) of the IPOL
            paper)

   *** IMPORTANT ***

   on entry, all entries of A must be initialized with 0

   Module output
   =============

   This module returns EXIT_SUCCESS if no error occured and
   EXIT_FAILURE otherwise.
                                                                       */
int error_prediction(double *A, double *psnr, double *mse, double *dx, double *dy, double *sigma, int nx, int ny, int Nx, int Ny, int nimages, double peakval, double eps, char vflag, char noshift)
{
  lapack_complex_double *M1pinv=NULL,*M2pinv=NULL,*M3pinv=NULL,*M4pinv=NULL,*Mpinv=NULL,*left,*right;
  lapack_int Z1,Z2,Z3,Z4,Z;
  double cof,zx,zy,cond1,cond2,cond3,cond4,fx,fy,theta,nrm2,re,im;
  int nalloc,a,b,pmin,pmax,qmin,qmax,adr_alf,adr_bet,adr,j,size_p,size_q,p,q,alf,bet,adr2; 
  char zxisinteger,zyisinteger; 
  
  // initialize some local variables //
  zx = (double)Nx/(double)nx; 
  zy = (double)Ny/(double)ny; 
  cof = 1./(zx*zy);
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
  if(vflag) {
    printf("Precompute block matrices (M1,M2,M3,M4)\n");
    printf("=======================================\n");
  }
  ASSERT_ALLOC(M1pinv = (lapack_complex_double*) LAPACKE_malloc (Z1*Z1*sizeof(lapack_complex_double)));
  if(EXIT_SUCCESS != compute_blockmatrix(M1pinv,dx,dy,NULL,cof,(int)floor(zx),(int)floor(zy),nimages)) { // set M1pinv = M1 (block matrix is not pseudo-inverted yet)
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
    if(EXIT_SUCCESS != compute_blockmatrix(M4pinv,dx,dy,NULL,cof,(int)ceil(zx),(int)ceil(zy),nimages)) { // set M4pinv = M4 (block matrix is not pseudo-inverted yet)
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
      if(EXIT_SUCCESS != compute_blockmatrix(M2pinv,dx,dy,NULL,cof,(int)ceil(zx),(int)floor(zy),nimages)) { // set M2pinv = M2 (block matrix is not pseudo-inverted yet)
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
      if(EXIT_SUCCESS != compute_blockmatrix(M3pinv,dx,dy,NULL,cof,(int)floor(zx),(int)ceil(zy),nimages)) { // set M3pinv = M3 (block matrix is not pseudo-inverted yet)
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

      // set Mpinv = pseudo inverse of M(a,b)
      if ((size_p == (int)floor(zx)) && (size_q == (int)floor(zy))) Mpinv = M1pinv; 
      else if ((size_p == (int)ceil(zx)) && (size_q == (int)floor(zy))) Mpinv = M2pinv;
      else if ((size_p == (int)floor(zx)) && (size_q == (int)ceil(zy))) Mpinv = M3pinv;
      else Mpinv = M4pinv;
      
      // compute c^{j}_k coefficients (matrix-vector product between
      // Mpinv and complex exponential vector)
      for(j=0;j<nimages;j++) {
	
	// range the set of frequencies (alf,bet) that are aliased at
	// position (a,b) in the low-frequency domain and compute the
	// entries of the comple-exponential vector (notice that we
	// range P(a) x Q(b) in lexicographical order)
	for(adr=0,p=pmin;p<=pmax;p++) {
	  alf = a + p*nx;
	  fx = (double)alf*dx[j]/(double)nx;
	  for(q=qmin;q<=qmax;q++,adr++) {
	    bet = b + q*ny;
	    fy = (double)bet*dy[j]/(double)ny;
	    theta = -2.*M_PI*(fx+fy);
	    right[adr] = lapack_make_complex_double(cos(theta),sin(theta));
	  }
	}
	
	// Multiply the pseudoinverse to the complex-exponential term
	// to retrieve the c^{j} vector entries (complex matrix vector
	// multiplication left = Mpinv*right)
	cmvprod(Z,Z,Mpinv,right,left); // the same operation can be done using CBLAS_zgemv module using the three next commented lines below.
	//lapack_complex_double one = lapack_make_complex_double(1.,0.);
	//lapack_complex_double zero = lapack_make_complex_double(0.,0.);
	//cblas_zgemv(CblasRowMajor,CblasNoTrans,Z,Z,&one,Mpinv,Z,right,1,&zero,left,1);
	
	// update the error amplification map (add |c^{j}_k(a,b)|^2 to
	// each (alf_k,bet_k) position      	
	for(adr=0,p=pmin;p<=pmax;p++) {
	  alf = a + p*nx;
	  adr_alf = noshift ? (Nx+alf)%Nx : (alf+Nx/2); // horizontal position of frequency alf in 'A'
	  for(q=qmin;q<=qmax;q++,adr++) {
	    bet = b + q*ny;
	    adr_bet = noshift ? (Ny+bet)%Ny : (bet+Ny/2); // vertical position of frequency bet in 'A'
	    adr2 = adr_bet*Nx + adr_alf; // position of frequency (alf,bet) in 'A'
	    re = lapack_complex_double_real(left[adr]);
	    im = lapack_complex_double_imag(left[adr]);
	    A[adr2] += (re*re + im*im); 
	  }
	}
      } // end for j
    } // end for a
  } // end for b
  
  // multiply A by 1/(zx*zy) and takes the square root
  for(adr=0;adr<Nx*Ny;adr++) A[adr] = sqrt(cof*A[adr]); 
      
  // if necessary compute MSE & PSNR predictions
  if(sigma) {
    for(nrm2=0.,adr=0;adr<Nx*Ny;adr++) nrm2 += pow(A[adr],2);
    *mse = pow(*sigma,2)*nrm2/((double)Nx*(double)Ny);
    *psnr = 10.*log10(pow(peakval,2)/(*mse)); 
  }

  return EXIT_SUCCESS;
}

