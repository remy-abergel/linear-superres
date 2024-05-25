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
extern int irls(fftw_complex*,double*,fftw_complex**,double*,double*,int,int,int,int,int,double,double,char,int,double); // see source file 'leastsquares_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nIteratively Reweighted Least-Squares.\n\n");
  printf("Usage: irls [-ftype type] [-M M] [-zx zx] [-N N] [-zy zy] [-n niter] [-e e] [-w weights] [-v] in T out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -M M        : (default 2 * width of in) set width of the ouptut image...\n");
  printf(" -zx zx      : ... or set X-axis super-resolution factor (double >= 1.)\n");
  printf(" -N N        : (default 2 * height of in) set height of the output image...\n");
  printf(" -zy zy      : ... or set Y-axis super-resolution factor (double >= 1.)\n");
  printf(" -n niter    : (default 50) maximal number of IRLS iterations\n");
  printf(" -r r        : (default 1e-5) stop the iterations when |E(n)-E(n-1)| <= r*|E(n)|,\n");
  printf("               denoting by E(n) the l1-l2 energy of the residual at iteration n\n");
  printf(" -e e        : (default 0) impose eta >= e in the IRLS scheme\n");
  printf(" -w weights  : output file containing the final weights in ASCII format (the j-th\n");
  printf("               weight is displayed at the j-th line of the file)\n");
  printf(" -v          : enable verbose mode\n");
  printf(" in          : input sequence of low-resolution images (TIFF-stack)\n");
  printf(" T           : input sequence of displacements (ASCII format)\n");
  printf(" out         : output high-resolution TIFF image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_T=NULL,*fname_out=NULL,*fname_weights=NULL,*M_value=NULL,*zx_value=NULL,*N_value=NULL,*zy_value=NULL,*n_value=NULL,*r_value=NULL,*e_value=NULL,vflag=0;
  double *dx=NULL,*dy=NULL,*weights=NULL,zx,zy,zx_tmp,zy_tmp,r,e;
  fftw_complex **u0=NULL,*u=NULL;
  int niter,err,k,type,nx,ny,nimages,ntranslations,Nx,Ny;
  FILE *fweights=NULL;

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
    else if (strcmp(argv[k],"-r") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -r requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	r_value = argv[k+1]; err = sscanf(r_value,"%lf",&r); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -r\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-n") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -n requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	n_value = argv[k+1]; err = sscanf(n_value,"%d",&niter); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -n\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-e") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -e requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	e_value = argv[k+1]; err = sscanf(e_value,"%lf",&e); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -e\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-w") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -w requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_weights = argv[k+1]; k++; }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_T) fname_T = argv[k];
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
  if(NULL == e_value) e = 1e-2;
  if(NULL == n_value) niter = 50;
  if(NULL == r_value) r = 1e-5;

  if ((NULL == fname_in) || (NULL == fname_T) || (NULL == fname_out)) {
    display_usage();
    if (NULL == fname_in) printf("Error: input 'in' is missing\n\n");
    else if (NULL == fname_T) printf("Error: input 'T' is missing\n\n");
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

  /*******************************************/
  /* if necessary, open output fweights file */
  /*******************************************/
  if(fname_weights) {
    if(NULL == (fweights = fopen(fname_weights,"w"))) {
      printf("Error: could not open file '%s'\n",fname_weights);
      return EXIT_FAILURE;
    }
  }

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

  /*********************/
  /* memory allocation */
  /*********************/
  ASSERT_ALLOC(u = (fftw_complex*) malloc (Nx*Ny*sizeof(fftw_complex)));
  ASSERT_ALLOC(weights = (double *) malloc (nimages*sizeof(double)));

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  if(vflag) {
    printf("Iteratively Reweighted Least-Squares iterations:\n");
    printf("================================================\n");
    fflush(stdout);
  }
  // run the IRLS Algorithm //
  if(EXIT_SUCCESS != irls(u,weights,u0,dx,dy,nx,ny,Nx,Ny,nimages,e,1e-9,vflag,niter,r)) {
    for(k=0;k<nimages;k++) fftw_free(u0[k]);
    fftw_free(u); free(u0); free(dx); free(dy); free(weights);
    return EXIT_FAILURE;
  }
  if(vflag) {
    printf("\n");
  }

  /*************************************/
  /* if necessary, save output weights */
  /*************************************/
  if(fweights) {
    for(k=0;k<nimages;k++) fprintf(fweights,"%.17e\n",weights[k]);
    fclose(fweights);
    fweights = NULL;
  }

  /********************************************/
  /* Save output (high-resolution) TIFF image */
  /********************************************/
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
    printf("save '%s' (precision %s): done\n",fname_out,((NULL == datatype)?"float":datatype));
    fflush(stdout);
  }
  if(vflag) printf("\n");

  /***************/
  /* free memory */
  /***************/
  for(k=0;k<nimages;k++) fftw_free(u0[k]);
  fftw_free(u);
  free(u0);
  free(dx);
  free(dy);
  free(weights);

  return EXIT_SUCCESS;
}
