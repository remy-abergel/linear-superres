#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }
#define min(a,b) ((a)>(b)?(b):(a))
#define max(a,b) ((a)<(b)?(b):(a))

/* internal modules */
int pseudoinverse_rowmajor(lapack_complex_double*,lapack_int,lapack_int,double,double*);
int compute_blockmatrix(lapack_complex_double*,double*,double*,double*,double,int,int,int);
int cmvprod(int,int,lapack_complex_double*,lapack_complex_double*,lapack_complex_double*);
static int zdscal(int,double,lapack_complex_double*,int); 
static int ctctprod(int,int,int,lapack_complex_double*,lapack_complex_double*,lapack_complex_double*); 

/* pseudoinverse_rowmajor: change matrix A into its pseudoinverse.

   Module parameters
   =================

   + nrow, ncol : the number of row and column of the input matrix A

   + A : (lapack_complex_double array containing nrow*ncol elements)
         on entry -> the entries of A stored in ROW-MAJOR order, which
         means that the entry of A at line k (1<=k<=nrow) and column l
         (1<=l<=ncol) is A[(k-1)*ncol+(l-1)], on exit -> the pseudo
         inverted matrix entries stored in ROW-MAJOR order

   + eps : a double number used as a threshold for the singular values
           of A in the pseudo inversion process (only the singular
           values > eps will be inverted)

   + cond : on entry -> double pointer, on exit -> *cond is equal to
            the condition number (in l2 norm) of the input matrix A
            (also equal to the condition number in l2 norm of the
            output matrix A if eps = 0)

   Module output
   =============
   
   This module returns EXIT_SUCCESS if no error occured and
   EXIT_FAILURE otherwise.

   Description of the pseudo-inversion procedure
   =============================================

   Let us denote m = min(nrow,ncol)

   STEP 1.

   Use LAPACKE_ZGESVD to compute the singular value decomposition
   (SVD) of the nrow-by-ncol matrix A, leading to

   A = U * SIGMA * V^H

   where SIGMA is a m-by-m diagonal matrix, U is an nrow-by-m matrix
   and V^H (V conjugate transposed) is an m-by-ncol matrix. The
   diagonal elements of SIGMA are the first m singular values of A
   sorted in descending order; they are real and non-negative.

   Notice that U and VH are not square matrices because we set the
   'jobu' and 'jobvt' input arguments of the LAPACKE_ZGESVD routine
   equal to 'S' (the routine returns the first m column of U and the
   first m rows of VH, instead of the traditional nrow-by-nrow unitary
   matrix U and ncol-by-ncol unitary matrix V).

   Note also that the LAPACK_ZGESVD routine returns V^H, not V.

   STEP 2.

   Compute SIGMA^+, the pseudoinverse of SIGMA, by inverting all the
   elements of SIGMA that are strictly higher than eps and keeping the
   zero-values in place

   STEP 3.

   Compute the pseudoinverse of A (noted A^+ below) using the relation

   A^+ = V * (SIGMA^+) * U^H.

   More precisely, the computation of A^+ is done in the following
   order

   A^+ = (V^H)^H * (U * SIGMA^+)^H

   where

   + (U * SIGMA^+) is computed from SIGMA^+ and U by means of complex
     scalar-vector multiplications using the zdscal routine (the l-th
     column of U is multiplied by the l-th diagonal element of
     SIGMA^+)

   + (V^H)^H * (U * SIGMA^+)^H is computed from V^H and (U * SIGMA^+)
     using the ctctprod routine.

*/
int pseudoinverse_rowmajor(lapack_complex_double *A, lapack_int nrow, lapack_int ncol, double eps, double *cond)
{
  lapack_int lda,ldu,ldvh,m,info;
  lapack_complex_double *U,*VH;
  double *SIGMA,*superb,scal;
  int k,l;

  // initalize dimensions & local variables //
  m = min(nrow,ncol);
  lda = ncol;  // leading dimension for a
  ldu = m;     // leading dimension for u
  ldvh = ncol; // leading dimension for vh

  /* memory allocation */
  ASSERT_ALLOC(U = (lapack_complex_double*) LAPACKE_malloc (ldu*nrow*sizeof(lapack_complex_double)));
  ASSERT_ALLOC(VH = (lapack_complex_double*) LAPACKE_malloc (ldvh*m*sizeof(lapack_complex_double)));
  ASSERT_ALLOC(SIGMA = (double*) malloc (m*sizeof(double)));
  ASSERT_ALLOC(superb = (double*) malloc ((m-1)*sizeof(double)));

  //-------------------------------//
  //      CORE OF THE MODULE       //
  //-------------------------------//

  // STEP 1. Compute [U,SIGMA,VH] the SVD of A (A = U*SIGMA*VH) and check for convergence (also compute cond) //
  info = LAPACKE_zgesvd(LAPACK_ROW_MAJOR,'S','S',nrow,ncol,A,lda,SIGMA,U,ldu,VH,ldvh,superb);
  if(info > 0) {
    printf("The LAPACKE SVD algorithm failed to converge.\n");
    LAPACKE_free(U); LAPACKE_free(VH); free(SIGMA); free(superb);
    exit(EXIT_FAILURE);
  }
  *cond = (SIGMA[m-1] == 0) ? INFINITY : SIGMA[0]/SIGMA[m-1];
  // STEP 2. Set U <- U * (SIGMA^+) by means of complex scalar-vector multiplications using the zdscal routine //
  for(l=0;l<m;l++) {
    scal = (SIGMA[l] > eps) ? 1./SIGMA[l] : SIGMA[l];
    zdscal(nrow, scal, &U[l], ldu); // same as cblas_zdscal(nrow, scal, &U[l], ldu); (we prefered use our own implementation of zdscal to avoid libopenblas dependency)
  }
  // STEP 3. Set A <- (VH)^H * U^h using the complex matrix-matrix multiplication routine cblas_zgemm //
  ctctprod(ncol,nrow,m,VH,U,A); // the same operation can be done using CBLAS_zgemm module using the three next commented lines below.
  //lapack_complex_double one = lapack_make_complex_double(1.,0.);
  //lapack_complex_double zero = lapack_make_complex_double(0.,0.);
  //cblas_zgemm(CblasRowMajor,CblasConjTrans,CblasConjTrans,ncol,nrow,m,&one,VH,ldvh,U,ldu,&zero,A,nrow);

  //-------------------------------//
  // END OF THE CORE OF THE MODULE //
  //-------------------------------//

  /* free memory */
  LAPACKE_free(U);
  LAPACKE_free(VH);
  free(SIGMA);
  free(superb);

  return EXIT_SUCCESS;
}

/* compute_blockmatrix: this module corresponds to the pseudocode
   Algorithm 1 of the IPOL paper (with a slight generalization due to
   the presence of weighting terms).

   This module computes the block matrix M with size Z x Z (Z = size_p
   x size_q) such that, denoting by M(k,l) the entry of M at row k and
   column l (for k=1..Z, l=1..Z), we have
   
   M(k,l) = sum_{j=0..nimages-1} weights[j] * Mj(k,l)

   Mj(k,l) = 

     (1/cof) *
     exp(2*i*pi*(floor((l-1)/size_q) - floor((k-1)/size_q))*dx[j]) *
     exp(2*i*pi*((l-1)%size_q - (k-1)%size_q)*dy[j]).
   
   Notice that:
   
   + when used with weights=NULL, or, equivalently, with weights[j]=1
     (for j = 0..nimages-1), this module implements the pseudocode
     Algorithm 1 of the companion IPOL paper.

   + when used with weight[j] = 1/eta[j] (for j = 0..nimages-1), this
     module computes the block matrix M_eta defined in Section 6 of
     the IPOL paper.

   Module parameters
   =================

   + nimages : a positive integer

   + dx, dy : (double arrays containing nimages elements each)
              horizontal and vertical components of the sequence of
              translation vectors (the translation associated to the
              j-th image is (dx[j],dy[j]))

   + weights : (double array containing nimages elements) weighting
               parameters (note that when weights = NULL this module
               assumes that weights[j] = 1 for j=0..nimages-1)

   + cof : a double number (noted gamma in the pseudocode Algorithm 1
           of the IPOL paper)

   + size_p : a positive integer

   + size_q : a positive integer

   + M : on entry -> a preallocated lapack_complex_double array
         containing ZxZ elements (with Z=size_p*size_q), on exit ->
         the entries of M stored in ROW-MAJOR order, which means that
         the entry of M(k,l) described above is M[(k-1)*ncol+(l-1)].

   Module output
   =============

   This module should return EXIT_SUCCESS.
   
*/
int compute_blockmatrix(lapack_complex_double *M, double *dx, double *dy, double *weights, double cof, int size_p, int size_q, int nimages)
{
  double cx,cy,re,im,theta;
  int j,k,l,adr;
  lapack_int z=(lapack_int)(size_p*size_q);

  memset(M,0.,z*z*sizeof(lapack_complex_double));
  for(j=0;j<nimages;j++) {
    cx = 2.*M_PI*dx[j];
    cy = 2.*M_PI*dy[j];
    for(adr=k=0;k<z;k++)
      for(l=0;l<z;l++,adr++) {
	theta = cx*(double)((l/size_q)-(k/size_q)) + cy*(double)((l%size_q)-(k%size_q));
	re = cof * ((NULL==weights) ? 1. : weights[j]) * cos(theta);
	im = cof * ((NULL==weights) ? 1. : weights[j]) * sin(theta);
	M[adr] = lapack_make_complex_double(lapack_complex_double_real(M[adr])+re,lapack_complex_double_imag(M[adr])+im);
      }
  }
  
  return EXIT_SUCCESS;
}

/* zdscal : multiplies each element of a vector by a constant
            (double-precision complex), similar to the cblas_zdscal
            function from the libopenblas library.

   Set X <- alpha * X

   Module parameters
   =================

   + N     : number of elements in the vector X

   + alpha : scalar factor

   + X     : lapack_complex_double array 

   + incX  : leading dimension for X (increment between two consecutive
             elements in vector X)
 
   Module output
   =============
   
   This module should return EXIT_SUCCESS. 

*/
static int zdscal(int N, double alpha, lapack_complex_double *X, int incX)
{
  int adr1,adr2;
  
  for(adr1=adr2=0; adr1<N; adr1++,adr2+=incX)
    X[adr2] = lapack_make_complex_double(alpha*lapack_complex_double_real(X[adr2]),alpha*lapack_complex_double_imag(X[adr2]));
  
  return EXIT_SUCCESS;
}

/* chhprod : conjugate transpose - conjugate transpose complex matrix
             product.

   This module computes

   C = A^H * B^H 

   where A is a K by M complex matrix, B is a N by K complex matrix
   and A^H and B^H denote their complex transpose.

   C is a M by N matrix equal to the product of A^H with B^H. 

   Note : the same matrix product can also be handled using the (more
   powerful) module cblas_zgemm from the libopenblas library. The
   present module is an alternative proposed to avoid the unnecessary
   dependency.

   Module parameters
   =================
   
   + M : number of rows of the matrix C, also equal to the number of
         columns of the matrix A

   + N : number of columns of the matrix C, also equal to the number
         of rows of the matrix B

   + K : number of rows of the matrix A, also equal to the number of
         columns of the matrix B

   + A : input left-hand term of the product (lapack_complex_double
         array containing the K*N matrix elements stored in ROW MAJOR
         order)

   + B : input righ-hand term of the product (lapack_complex_double
         array containing the N*K matrix elements stored in ROW MAJOR
         order) 

   + C : output matrix, on entry -> preallocated lapack_complex_double
         array, on exit -> the M*N matrix elements stored in ROW MAJOR
         order

   Module output
   =============

   This module should return EXIT_SUCCESS.

*/
static int ctctprod(int M, int N, int K, lapack_complex_double *A, lapack_complex_double *B, lapack_complex_double *C)
{
  int i,j,k,adr;
  double are,aim,bre,bim,cre,cim;

  for(i=0; i<M; i++)
    for(j=0; j<N; j++) {
      cre = cim = 0.;
      for(k=0; k<K; k++) {
	// retrieve (are,aim), the real and imaginary parts of the element of A at position (k,i)
	adr = k*M+i; 
	are = lapack_complex_double_real(A[adr]); 
	aim = lapack_complex_double_imag(A[adr]);
	// retrieve (bre,bim), the real and imaginary parts of the element of B at position (j,k)
	adr = j*K+k;
	bre = lapack_complex_double_real(B[adr]); 
	bim = lapack_complex_double_imag(B[adr]);
	// increment (cre,cim), the real and imaginary parts of the element of C at position (i,j)
	cre += are*bre - aim*bim;
	cim -= are*bim + aim*bre; 
      }
      C[i*N+j] = lapack_make_complex_double(cre,cim); 
    }

  return EXIT_SUCCESS; 
}


/* cmvprod : complex matrix vector product 

   This module computes

   Y = A * X 

   where A is a M by N complex matrix, X is a complex column vector
   containing N elements and Y is the output complex column vector
   containing M elements.

   Note : the same matrix - vector operation can also be handled using
   the (more powerful) module cblas_zgemv from the libopenblas
   library. The present module is an alternative proposed to avoid the
   unnecessary dependency.

   Module parameters
   =================

   + M : number of rows of the matrix A, also equal to the number of
         elements of the column vector Y

   + N : number of columns of the matrix A, also equal to the number
         of elements of the column vector X

   + A : lapack_complex_double array containing the M*N matrix
         entries stored in ROW MAJOR order

   + X : lapack_complex_double array containing the N elements of N
         entries of the vector X

   + Y : on entry -> preallocated lapack_complex_double array
         containing M elements, on exit -> contains the M entries of
         the matrix-vector product A*X
   
   Module output
   =============

   This module should return EXIT_SUCCESS.

 */
int cmvprod(int M, int N, lapack_complex_double *A, lapack_complex_double *X, lapack_complex_double *Y)
{ 
  int i,k,adr;
  double are,aim,xre,xim,yre,yim;
  for(adr=i=0; i<M; i++) {
    yre = yim = 0.;
    for(k=0; k<N; k++,adr++) {
      are = lapack_complex_double_real(A[adr]);
      aim = lapack_complex_double_imag(A[adr]);
      xre = lapack_complex_double_real(X[k]); 
      xim = lapack_complex_double_imag(X[k]); 
      yre += are*xre - aim*xim;
      yim += are*xim + aim*xre;
    }
    Y[i] = lapack_make_complex_double(yre,yim);
  }
}
