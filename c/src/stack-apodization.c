#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'
extern int getasciitranslations(double**,double**,int*,char*,char); // see source file 'ascii.c'
extern int apodization(double**,double**,double*,double*,double*,double,int,int,int,int,int,char); // see source file 'apodization_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nCompute low/high resolution multiplicative apodization filters.\n\n");
  printf("Usage: stack-apodization [-ftype type] [-s sigma] [-M M] [-zx zx] [-N N] [-zy zy] [-g apod_lr] [-G apod_hr] [-a apod_u0] [-v] u0 T\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -s sigma    : (default 1.) sharpness parameter for the profile (positive double)\n");
  printf("               such as 10*sigma = length of the transition intervals\n");
  printf(" -M M        : (default 2 * width of u0) set width of the high-resolution domain...\n");
  printf(" -zx zx      : ... or set X-axis super-resolution factor (double >= 1.)\n");
  printf(" -N N        : (default 2 * height of u0) set height of the high-resolution domain...\n");
  printf(" -zy zx      : ... or set Y-axis super-resolution factor (double >= 1.)\n");
  printf(" -G apod_hr  : output high-resolution apodization filter (monopage TIFF-float image)\n");
  printf(" -g apod_lr  : output stack of low-resolution apodization filters (multipage TIFF-float image)\n");
  printf(" -a apod_u0  : output apodized stack of low-resolution images (multipage TIFF image)\n");
  printf(" -v          : enable verbose mode\n");
  printf(" u0          : input low-resolution stack (multipage TIFF image)\n");
  printf(" T           : input sequence of displacements (ASCII format)\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_apod_hr=NULL,*fname_apod_lr=NULL,*fname_apod_u0=NULL,*fname_u0=NULL,*fname_T=NULL,*M_value=NULL,*N_value=NULL,*zx_value=NULL,*zy_value=NULL,*sigma_value=NULL,vflag=0;
  double **u0=NULL,**apod_lr=NULL,*apod_hr=NULL,*dx=NULL,*dy=NULL,sigma,zx,zy,zx_tmp,zy_tmp;
  int err,k,type,nx,ny,nimages,ntranslations,Nx,Ny;

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
    else if (strcmp(argv[k],"-s") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -s requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	sigma_value = argv[k+1]; err = sscanf(sigma_value,"%lf",&sigma); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -s\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-G") == 0) { fname_apod_hr = argv[k+1]; k++; }
    else if (strcmp(argv[k],"-g") == 0) { fname_apod_lr = argv[k+1]; k++; }
    else if (strcmp(argv[k],"-a") == 0) { fname_apod_u0 = argv[k+1]; k++; }
    else if (NULL == fname_u0) fname_u0 = argv[k];
    else if (NULL == fname_T) fname_T = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

  /********************************************/
  /* check consistency and set default values */
  /********************************************/
  if(NULL != M_value && NULL != zx_value) { display_usage(); printf("Error: options -M and -zx cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != N_value && NULL != zy_value) { display_usage(); printf("Error: options -N and -zy cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != zx_value && zx < 1.) { display_usage(); printf("Error: input 'zx' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(NULL != zy_value && zy < 1.) { display_usage(); printf("Error: input 'zy' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(NULL == sigma_value) sigma = 1.;

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == fname_u0) || (NULL == fname_T)) {
    display_usage();
    if (NULL == fname_T) printf("Error: input 'T' is missing\n\n");
    else printf("Error: input 'u0' is missing\n\n");
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

  if(sigma <= 0.){ display_usage(); printf("Error: input 'sigma' must be positive.\n\n"); return EXIT_FAILURE; }

  if((NULL == fname_apod_lr) && (NULL == fname_apod_hr) && (NULL == fname_apod_u0)) {
    display_usage(); printf("Error: you must select at least one option in {-a,-g,-G}.\n\n"); return EXIT_FAILURE;
  }

  /*************************/
  /* load input TIFF image */
  /*************************/
  if(EXIT_FAILURE == load_tiff_into_double(&u0,&nx,&ny,&nimages,fname_u0,0)) {
    printf("Error: failed to open input image '%s'\n",fname_u0);
    return EXIT_FAILURE;
  }

  if(M_value && Nx < nx) {
    display_usage();
    printf("Error: input 'M' must be greater than or equal to the width of the input sequence (M >= %d).\n\n",nx);
    for(k=0;k<nimages;k++) free(u0[k]);
    free(u0); 
    return EXIT_FAILURE;
  }
  if(N_value && Ny < ny) {
    display_usage();
    printf("Error: input 'N' must be greater than or equal to the height of the input sequence (N >= %d).\n\n",ny);
    for(k=0;k<nimages;k++) free(u0[k]);
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
    if(zx_tmp != zx) printf("WARNING: zx has been changed into %g (instead of zx=%g) to allow integer width of the output image\n",zx,zx_tmp);
  }

  if(NULL == zy_value && NULL == N_value) { zy = 2.; Ny = 2*ny; }
  else if(NULL == zy_value) zy = (double)Ny/(double)ny;
  else {
    Ny = (int)round(ny*zy);
    zy_tmp = zy;
    zy = (double)Ny/(double)ny;
    if(zy_tmp != zy) printf("WARNING: zy has been changed into %g (instead of zy=%g) to allow integer height of the output image\n",zy,zy_tmp);
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
    printf("Error: number of images in the input multipage TIFF image '%s' should be the same\n",fname_u0);
    printf("as the number of translations found in the input translation file '%s' (found %d images and %d translations).\n\n",fname_T,nimages,ntranslations);
    for(k=0;k<nimages;k++) free(u0[k]); free(u0); free(dx); free(dy);
    return EXIT_FAILURE;
  }

  /*********************/
  /* memory allocation */
  /*********************/
  if(fname_apod_hr) { ASSERT_ALLOC(apod_hr = (double*) malloc (Nx*Ny*sizeof(double))); }
  if(fname_apod_lr) {
    ASSERT_ALLOC(apod_lr = (double**) malloc (nimages*sizeof(double*)));
    for(k=0;k<nimages;k++) ASSERT_ALLOC(apod_lr[k] = (double*) malloc (nx*ny*sizeof(double)));
  }

  /***********************/
  /* CORE OF THE ROUTINE */
  /***********************/
  apodization((fname_apod_u0?u0:NULL),apod_lr,apod_hr,dx,dy,sigma,nx,ny,Nx,Ny,nimages,vflag);

  /******************************/
  /* write output TIFF image(s) */
  /******************************/
  if(vflag) {
    printf("Save output TIFF images and free memory\n");
    printf("=======================================\n");
  }

  if(apod_hr) { // high-resolution apodization filter
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_apod_hr,apod_hr,Nx,Ny,0,0,"w")) {
      printf("Error: failed to save output image '%s'\n",fname_apod_hr);
      free(dx);
      free(dy);
      for(k=0;k<nimages;k++) {
	free(u0[k]);
	if(apod_lr) free(apod_lr[k]);
      }
      free(u0);
      if(apod_lr) free(apod_lr);
      if(apod_hr) free(apod_hr);
      return EXIT_FAILURE;
    }
    if(vflag) printf("save '%s' (precision float): done \n",fname_apod_hr);
  }

  if(apod_lr) { // low-resolution apodization filters
    for(k=0;k<nimages;k++) {
      if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_apod_lr,apod_lr[k],nx,ny,0,0,((k==0)?"w":"a"))) {
	printf("Error: failed to save output image '%s'\n",fname_apod_lr);
	free(dx);
	free(dy);
	for(k=0;k<nimages;k++) {
	  free(u0[k]);
	  if(apod_lr) free(apod_lr[k]);
	}
	free(u0);
	if(apod_lr) free(apod_lr);
	if(apod_hr) free(apod_hr);
	return EXIT_FAILURE;
      }
      if(vflag) printf("\33[2K\rsave '%s' page %d/%d (precision float): done ",fname_apod_lr,k+1,nimages);
    }
    if(vflag) printf("\33[2K\rsave '%s' page %d/%d (precision float): done \n",fname_apod_lr,nimages,nimages);
  }

  if(fname_apod_u0) { // apodized low-resolution stack
    for(k=0;k<nimages;k++) {
      if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_apod_u0,u0[k],nx,ny,0,type,((k==0)?"w":"a"))) {
	printf("Error: failed to save output image '%s'\n",fname_apod_u0);
	free(dx);
	free(dy);
	for(k=0;k<nimages;k++) {
	  free(u0[k]);
	  if(apod_lr) free(apod_lr[k]);
	}
	free(u0);
	if(apod_lr) free(apod_lr);
	if(apod_hr) free(apod_hr);
	return EXIT_FAILURE;
      }
      if(vflag) printf("\33[2K\rsave '%s' page %d/%d (precision %s): done ",fname_apod_u0,k+1,nimages,((NULL == datatype)?"float":datatype));
    }
    if(vflag) printf("\33[2K\rsave '%s' page %d/%d (precision %s): done \n",fname_apod_u0,nimages,nimages,((NULL == datatype)?"float":datatype));
  }
  if(vflag) printf("\n");

  /***************/
  /* free memory */
  /***************/
  for(k=0;k<nimages;k++) {
    free(u0[k]);
    if(apod_lr) free(apod_lr[k]);
  }
  free(u0);
  free(dx);
  free(dy);
  if(apod_lr) free(apod_lr);
  if(apod_hr) free(apod_hr);

  return EXIT_SUCCESS;
}
