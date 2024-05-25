#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_fftw_complex(fftw_complex***,int*,int*,int*,char*,char); // see source file 'tiffreadcomplex.c'
extern int save_monopage_striped_tiff_from_fftw_complex(char*,fftw_complex*,int,int,char,int,char*,char*,char*); // see source file 'tiffwritecomplex.c'
extern int getasciitranslations(double**,double**,int*,char*,char); // see source file 'ascii.c'
extern int getasciiweights(double**,int*,char*,char); // see source file 'ascii.c'
extern int luckyimaging(fftw_complex*,fftw_complex**,double*,fftw_complex**,double*,double*,int,int,int,int,int,int,int,int,int,int*,double,char); // see source file 'leastsquares_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nLucky imaging procedure for least-squares image super-resolution.\n\n");
  printf("Usage: luckyimaging [-ftype type] [-M M] [-zx zx] [-N N] [-zy zy] [-n0 n0] [-n1 n1] [-l l] [-n n] [-v] in T weights out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -M M        : (default 2 * width of in) set width of the ouptut image...\n");
  printf(" -zx zx      : ... or set X-axis super-resolution factor (double >= 1.)\n");
  printf(" -N N        : (default 2 * height of in) set height of the output image...\n");
  printf(" -zy zy      : ... or set Y-axis super-resolution factor (double >= 1.)\n");
  printf(" -n0 n0      : (default min(L,ceil(zx)*ceil(zy))) minimal value for n (see below)\n"); 
  printf(" -n1 n1      : (default L) maximal value for n (see below)\n"); 
  printf(" -l l        : (default min(30,n1-n0+1)) number of n values (see below)\n");
  printf(" -n n        : force out = u_lucky^{n} (instead of computing a movie)\n"); 
  printf(" -v          : enable verbose mode\n");
  printf(" in          : input sequence (TIFF-stack) containing L low-resolution images\n");
  printf(" T           : input sequence of displacements (ASCII format)\n");
  printf(" weights     : IRLS weights in ASCII format (see the irls module)\n");
  printf(" out         : output sequence (TIFF-stack) of high-resolution images\n\n");
  printf("Details about the output sequence of high-resolution images : \n\n");
  printf("We note u_lucky^{n} the high-resolution image computed from the n low-\n");
  printf("resolution images associated to the n highest input weights.\n\n");
  printf("Unless option -n is used, this module returns a movie containing l images\n");
  printf("indexed from #1 to #l and such as the #i-th image of the movie corresponds\n");
  printf("to the image u_lucky^{n_i}, denoting by (n_i)_{1 <= i <= l} a decreasing \n");
  printf("sequence of l quasi-regularly spaced integers in [n0,n1] defined by\n");
  printf("n_i = round(a*i+b) where a = (n0-n1)/(l-1) and b = n0-a*l.\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_T=NULL,*fname_out=NULL,*fname_weights=NULL,*M_value=NULL,*zx_value=NULL,*N_value=NULL,*zy_value=NULL,*n_value=NULL,*n0_value=NULL,*n1_value=NULL,*l_value=NULL,vflag=0;
  double *dx=NULL,*dy=NULL,*weights=NULL,zx,zy,zx_tmp,zy_tmp;
  fftw_complex **u0=NULL,**mov=NULL,*u=NULL;
  int err,k,type,nx,ny,nimages,ntranslations,nweights,Nx,Ny,n,id,n0,n1,l,*nlist=NULL;
  char pname[1000],desc[1000]; 

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-v") == 0) vflag = 1;
    else if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-M") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -M requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	M_value = argv[k+1]; err = sscanf(M_value,"%d",&Nx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -M\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-zx") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -zx requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	zx_value = argv[k+1]; err = sscanf(zx_value,"%lf",&zx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -zx\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-N") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -N requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	N_value = argv[k+1]; err = sscanf(N_value,"%d",&Ny); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -N\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-zy") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -zy requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	zy_value = argv[k+1]; err = sscanf(zy_value,"%lf",&zy); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -zy\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-n") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -n requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	n_value = argv[k+1]; err = sscanf(n_value,"%d",&n); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -n\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-n0") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -n0 requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	n0_value = argv[k+1]; err = sscanf(n0_value,"%d",&n0); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -n0\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-n1") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -n1 requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	n1_value = argv[k+1]; err = sscanf(n1_value,"%d",&n1); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -n1\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-l") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -l requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	l_value = argv[k+1]; err = sscanf(l_value,"%d",&l); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -l\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_T) fname_T = argv[k];
    else if (NULL == fname_weights) fname_weights = argv[k];
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

  /********************************************/
  /* check consistency and set default values */
  /********************************************/
  if(NULL != M_value && NULL != zx_value) { display_usage(); printf("Error: options -M and -zx cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != N_value && NULL != zy_value) { display_usage(); printf("Error: options -N and -zy cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != zx_value && zx < 1.) { display_usage(); printf("Error: input 'zx' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(NULL != zy_value && zy < 1.) { display_usage(); printf("Error: input 'zy' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }

  if ((NULL == fname_in) || (NULL == fname_T) || (NULL == fname_weights) || (NULL == fname_out)) {
    display_usage();
    if (NULL == fname_in) printf("Error: input 'in' is missing\n\n");
    else if (NULL == fname_T) printf("Error: input 'T' is missing\n\n");
    else if (NULL == fname_weights) printf("Error: input 'weights' is missing\n\n");
    else printf("Error: output 'out' is missing\n\n");
    return EXIT_FAILURE;
  }

  if((NULL == datatype)||(strcmp(datatype,"float")==0)) type = 0;
  else if (strcmp(datatype,"uint8")==0) type = 8;
  else if (strcmp(datatype,"int8")==0) type = -8;
  else if (strcmp(datatype,"uint16")==0) type = 16;
  else if (strcmp(datatype,"int16")==0) type = -16;
  else if (strcmp(datatype,"uint32")==0) type = 32;
  else if (strcmp(datatype,"int32")==0) type = -32;
  else { display_usage(); printf("Error: unrecognized value for option -ftype.\n\n"); return EXIT_FAILURE; }

  /*******************************************************************/
  /* load input TIFF-Stack image (sequence of low-resolution images) */
  /*******************************************************************/
  if(EXIT_FAILURE == load_tiff_into_fftw_complex(&u0,&nx,&ny,&nimages,fname_in,0)){
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }

  if(M_value && Nx < nx) {
    display_usage();
    printf("Error: input 'M' must be greater than or equal to the width of the input sequence (M >= %d).\n\n",nx);
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    free(u0); 
    return EXIT_FAILURE;
  }
  if(N_value && Ny < ny) {
    display_usage();
    printf("Error: input 'N' must be greater than or equal to the height of the input sequence (N >= %d).\n\n",ny);
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    free(u0); 
    return EXIT_FAILURE;
  }

  /*****************************************************************************/
  /* retrieve & control I/O image dimensions (set default values if necessary) */
  /*****************************************************************************/
  if(NULL == zx_value && NULL == M_value) { zx = 2.; Nx = 2*nx; }
  else if(NULL == zx_value) zx = (double)Nx/(double)nx;
  else {
    Nx = (int)round(nx*zx);
    zx_tmp = zx;
    zx = (double)Nx/(double)nx;
    if(zx_tmp != zx) printf("WARNING: zx has been changed into %.17g (instead of zx=%g) to allow integer width of the output image\n",zx,zx_tmp);
  }

  if(NULL == zy_value && NULL == N_value) { zy = 2.; Ny = 2*ny; }
  else if(NULL == zy_value) zy = (double)Ny/(double)ny;
  else {
    Ny = (int)round(ny*zy);
    zy_tmp = zy;
    zy = (double)Ny/(double)ny;
    if(zy_tmp != zy) printf("WARNING: zy has been changed into %.17g (instead of zy=%g) to allow integer height of the output image\n",zy,zy_tmp);
  }

  // deal with -n0 option (set default value or check consistency)
  if(NULL == n0_value) n0 = (int)fmin(nimages,ceil(zx)*ceil(zy));
  else if(n0 < 1 || n0 > nimages) {
    display_usage();
    printf("Error: inconsistent value for option -n0 (must have 1 <= n0 <= L = %d).\n\n",nimages);
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    return EXIT_FAILURE;
  }

  // deal with -n1 option (set default value or check consistency)
  if(NULL == n1_value) n1 = nimages; 
  else if(n1 < n0 || n1 > nimages) {
    display_usage();
    printf("Error: inconsistent value for option -n1 (must have %d = n0 <= n1 <= L = %d).\n\n",n0,nimages);
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    return EXIT_FAILURE;
  }
    
  // deal with -l option (set default value or check consistency)
  if(NULL == l_value) l = (int)fmin(30,n1-n0+1); 
  else if(l < 1) {
    display_usage();
    printf("Error: inconsistent value for option -l (must have l >= 1).\n\n");
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    return EXIT_FAILURE;
  }
  else if(l > n1-n0+1) {
    printf("Warning: the value of l has been set to %d instead of %d (we must have l <= n1-n0+1 = %d).\n",n1-n0+1,l,n1-n0+1);
    l = n1-n0+1;
  }

  // deal with -n option (check consistency)
  if(n_value && (l_value || n0_value || n1_value)) {
    display_usage();
    printf("Error: option -n cannot be used jointly with options -n0, -n1 and -l.\n\n");
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    return EXIT_FAILURE;
  }
  if(n_value && (n <= 0 || n > nimages)) {
    display_usage();
    printf("Error: option -n must be used with 1 <= n <= L (here, the input stack contains L = %d images).\n\n",nimages);
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    return EXIT_FAILURE;
  }
  
  if(vflag) {
    printf("\nRetrieve image size parameters\n");
    printf("==============================\n");
    printf("width of the low-resolution domain: m = %d\n",nx);
    printf("height of the low-resolution domain: n = %d\n",ny);
    printf("width of the high-resolution domain: M = %d\n",Nx);
    printf("height of the high-resolution domain: N = %d\n",Ny);
    printf("horizontal super-resolution factor: zx = M/m = %g\n",zx);
    printf("vertical super-resolution factor: zy = N/n = %g\n",zy);
    printf("\n");
  }


  /**************************************/
  /* load the sequence of displacements */
  /**************************************/
  getasciitranslations(&dx,&dy,&ntranslations,fname_T,vflag);

  if(nimages != ntranslations) {
    printf("Error: number of images in the input multipage TIFF image '%s' should be the same\n",fname_in);
    printf("as the number of translations found in the input translation file '%s' (found %d images and %d translations).\n\n",fname_T,nimages,ntranslations);
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    free(u0); free(dx); free(dy);
    return EXIT_FAILURE;
  }

  /**************************/
  /* load the input weights */
  /**************************/
  getasciiweights(&weights,&nweights,fname_weights,vflag);

  if(nimages != nweights) {
    printf("Error: number of images in the input multipage TIFF image '%s' should be the same\n",fname_in);
    printf("as the number of weights found in the input weights file '%s' (found %d images and %d weights).\n\n",fname_weights,nimages,nweights);
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    free(u0); free(dx); free(dy); free(weights);
    return EXIT_FAILURE;
  }

  /*********************/
  /* memory allocation */
  /*********************/
  if(n_value) {
    ASSERT_ALLOC(u = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex)));
  }
  else {
    ASSERT_ALLOC(mov = (fftw_complex**) malloc (l*sizeof(fftw_complex*)));
    for(k=0;k<l;k++) { ASSERT_ALLOC(mov[k] = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex))); }
    ASSERT_ALLOC(nlist = (int*) malloc (l*sizeof(int)));
  }

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  if(vflag) {
    printf("Lucky imaging procedure:\n");
    printf("========================\n");
    fflush(stdout);
  }
  if(EXIT_SUCCESS != luckyimaging(u,mov,weights,u0,dx,dy,nx,ny,Nx,Ny,nimages,n,n0,n1,l,nlist,1e-9,vflag)) {
    for(k=0;k<nimages;k++) { fftw_free(u0[k]); }
    if(mov) {
      for(k=0;k<l;k++) { fftw_free(mov[k]); }
      free(mov);
    }
    if(u) fftw_free(u);
    if(nlist) free(nlist); 
    free(u0); free(dx); free(dy); free(weights);
    return EXIT_FAILURE;
  }
  if(vflag) {
    printf("\n");
  }

  /********************************************/
  /* Save output (high-resolution) TIFF image */
  /********************************************/
  if(!mov) { // output is monopage
    if(vflag) {
      printf("Save output high-resolution TIFF image and free memory\n");
      printf("======================================================\n");
    }
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_fftw_complex(fname_out,u,Nx,Ny,0,type,"w",NULL,NULL)) {
      printf("Error: failed to save output image '%s'\n",fname_out);
      for(k=0;k<nimages;k++) fftw_free(u0[k]);
      fftw_free(u); free(u0); free(dx); free(dy); free(weights);
      return EXIT_FAILURE;
    }
    if(vflag) {
      printf("save '%s' (precision %s): done",fname_out,((NULL == datatype)?"float":datatype));
      fflush(stdout);
    }
  }
  else { // output is multipage
    if(vflag) {
      printf("Save outputs high-resolution TIFF images and free memory\n");
      printf("========================================================\n");
    }
    for(k=0;k<l;k++) {
      sprintf(pname,"u_lucky^{(%d)}",nlist[k]);
      sprintf(desc,"super-resolved image computed after removing %d - %d = %d images from the input low-resolution stack",nimages,nlist[k],nimages-nlist[k]);
      if(EXIT_FAILURE == save_monopage_striped_tiff_from_fftw_complex(fname_out,mov[k],Nx,Ny,0,type,(k==0)?"w":"a",pname,desc)) {
	printf("Error: failed to save output image '%s'\n",fname_out);
	for(id=0;id<nimages;id++) { fftw_free(u0[id]); }
	for(id=0;id<l;id++) { fftw_free(mov[id]); }
	free(mov); free(u0); free(dx); free(dy); free(weights); free(nlist);
      }
      if(vflag) {
	printf("\33[2K\rsave '%s' page %d/%d (precision %s): done",fname_out,k+1,l,((NULL == datatype)?"float":datatype));
	fflush(stdout);
      }
    }
  }
  if(vflag) printf("\n\n");

  /***************/
  /* free memory */
  /***************/
  for(k=0;k<nimages;k++) fftw_free(u0[k]);
  if(mov) {
    for(k=0;k<l;k++) fftw_free(mov[k]);
    free(mov);
    free(nlist);
  }
  else { fftw_free(u); }
  free(u0);
  free(dx);
  free(dy);
  free(weights);

  return EXIT_SUCCESS;
}
