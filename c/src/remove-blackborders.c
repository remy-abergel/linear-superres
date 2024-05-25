#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_monopage_tiff_into_double(double**,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'
extern int getasciitranslations(double**,double**,int*,char*,char); // see source file 'ascii.c'
extern int remove_blackborders(double**,int*,int*,double*,double*,double*,int,int,int,int,int,double,char); // see source file 'remove-blackborders_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nSupress borders caused by apodization in a high-resolution image.\n\n");
  printf("Usage: remove-blackborders [-ftype type] [-m m] [-zx zx] [-n n] [-zy zy] [-v] T in out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -m m        : (default round(0.5 * width of in)) set width of the low-resolution domain...\n");
  printf(" -zx zx      : ... or set X-axis super-resolution factor (double >= 1.)\n");
  printf(" -n n        : (default round(0.5 * height of in)) set height of the low-resolution domain...\n");
  printf(" -zy zy      : ... or set Y-axis super-resolution factor (double >= 1.)\n");
  printf(" -s sigma    : (default 1.) sharpness parameter for the apodization profiles\n"); 
  printf(" -v          : enable verbose mode\n");
  printf(" T           : input sequence of displacements (ASCII format)\n");
  printf(" in          : input high-resolution TIFF image\n");
  printf(" out         : output (cropped) image TIFF image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_T=NULL,*fname_out=NULL,*m_value=NULL,*n_value=NULL,*zx_value=NULL,*zy_value=NULL,*s_value=NULL,vflag=0;
  double *dx=NULL,*dy=NULL,zx,zy,zx_tmp,zy_tmp,sigma,*u=NULL,*ucrop=NULL;
  int err,k,type,nx,ny,ntranslations,Nx,Ny,Nxcrop,Nycrop;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-v") == 0) vflag = 1;
    else if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-m") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -m requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	m_value = argv[k+1]; err = sscanf(m_value,"%d",&nx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -m\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-zx") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -zx requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	zx_value = argv[k+1]; err = sscanf(zx_value,"%lf",&zx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -zx\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-n") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -n requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	n_value = argv[k+1]; err = sscanf(n_value,"%d",&ny); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -n\n\n"); return EXIT_FAILURE; }
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
	s_value = argv[k+1]; err = sscanf(s_value,"%lf",&sigma); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -s\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_T) fname_T = argv[k];
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

  /*********************/
  /* check consistency */
  /*********************/
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

  if(NULL != m_value && NULL != zx_value) { display_usage(); printf("Error: options -m and -zx cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != n_value && NULL != zy_value) { display_usage(); printf("Error: options -n and -zy cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != zx_value && zx < 1.) { display_usage(); printf("Error: input 'zx' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(NULL != zy_value && zy < 1.) { display_usage(); printf("Error: input 'zy' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(NULL == s_value) sigma = 1.; 
  
  /*****************************************************************/
  /* load input TIFF image & perform dimensions consistency checks */
  /*****************************************************************/
  if(EXIT_FAILURE == load_monopage_tiff_into_double(&u,&Nx,&Ny,fname_in,0)) {
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }
  
  if(Nx <= 1) {
    printf("Error: the width of the input image 'in' must be at least 2\n\n");
    free(u); 
    return EXIT_FAILURE; 
  } 
  if(Ny <= 1) {
    printf("Error: the height of the input image 'in' must be at least 2\n\n");
    free(u); 
    return EXIT_FAILURE; 
  } 

  if(m_value && Nx < nx) {
    display_usage();
    printf("Error: input 'm' must be smaller than or equal to the width of the input sequence (m <= %d).\n\n",Nx);
    free(u); 
    return EXIT_FAILURE;
  }
  if(n_value && Ny < ny) {
    display_usage();
    printf("Error: input 'n' must be smaller than or equal to the height of the input sequence (n <= %d).\n\n",Ny);
    free(u); 
    return EXIT_FAILURE;
  }

  if(m_value && nx <= 0) {
    display_usage();
    printf("Error: input 'm' must be positive.\n\n");
    free(u); 
    return EXIT_FAILURE;
  }
  if(n_value && ny <= 0) {
    display_usage();
    printf("Error: input 'n' must be positive.\n\n");
    free(u); 
    return EXIT_FAILURE;
  }

  /*****************************************************************************/
  /* retrieve & control I/O image dimensions (set default values if necessary) */
  /*****************************************************************************/
  if(NULL == zx_value && NULL == m_value) { nx = (int)round(0.5*(double)Nx); zx = (double)Nx/(double)nx; }
  else if(NULL == zx_value) zx = (double)Nx/(double)nx;
  else {
    nx = (int)round((double)Nx/zx);
    zx_tmp = zx;
    zx = (double)Nx/(double)nx;
    if(zx_tmp != zx) printf("WARNING: zx has been changed into %.17g (instead of zx=%g) to allow integer width of the output image\n",zx,zx_tmp);
  }

  if(NULL == zy_value && NULL == n_value) { ny = (int)round(0.5*(double)Ny); zy = (double)Ny/(double)ny; }
  else if(NULL == zy_value) zy = (double)Ny/(double)ny;
  else {
    ny = (int)round((double)Ny/zy);
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
    printf("horizontal subsampling factor: zx = M/m = %g\n",zx);
    printf("vertical subsampling factor: zy = N/n = %g\n",zy);
    printf("\n");
  }

  /**************************************/
  /* load the sequence of displacements */
  /**************************************/
  getasciitranslations(&dx,&dy,&ntranslations,fname_T,vflag);

  /***********************/
  /* CORE OF THE ROUTINE */
  /***********************/
  if(EXIT_FAILURE == remove_blackborders(&ucrop,&Nxcrop,&Nycrop,u,dx,dy,Nx,Ny,nx,ny,ntranslations,sigma,vflag)) {
    printf("An error occured the image cropping process (internal module 'remove_blackborders')\n");
    free(dx); free(dy); 
    return EXIT_FAILURE; 
  }

  /*******************************/
  /* save (cropped) output image */
  /*******************************/
  if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_out,ucrop,Nxcrop,Nycrop,0,type,"w")) {
    printf("An error occured during the TIFF image saving process '%s'\n",fname_out);
    free(dx); free(dy); free(ucrop); 
    return EXIT_FAILURE;
  }

  /***************/
  /* free memory */
  /***************/
  free(dx);
  free(dy);
  free(ucrop); 

  return EXIT_SUCCESS;
}
