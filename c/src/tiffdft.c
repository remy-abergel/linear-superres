#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'

/* usage displayer */
static void display_usage()
{
  printf("\nCompute the 2D Fourier Transform (DFT) of a (monopage or multipage) TIFF image\n\n");
  printf("Usage: tiffdft [-ftype type] [-r real] [-i imag] [-m modulus] [-p phasis] [-l logmodulus] [-noshift] in\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -r real     : real part of DFT(in)\n");
  printf(" -i imag     : imaginary part of DFT(in)\n");
  printf(" -m modulus  : modulus of DFT(in)\n");
  printf(" -p phasis   : phasis (in radian) of DFT(in)\n");
  printf(" -l logmod   : log-modulus = log(1+|DFT(in)|)\n");
  printf(" -noshift    : do not center the frequency domain (put the (0,0) frequency\n");
  printf("               at the top-left corner of the image instead of its center)\n");
  printf(" in          : input TIFF (monopage or multipage) image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_real=NULL,*fname_imag=NULL,*fname_modulus=NULL,*fname_phasis=NULL,*fname_logmod=NULL,noshift=0;
  double **u=NULL,**re=NULL,**im=NULL,**rho=NULL,**theta=NULL,**logmod=NULL,val_re,val_im,val_rho;
  fftw_complex *dft_u=NULL;
  fftw_plan plan_u=NULL;
  int s,x,y,a,b,adr,adr2,err,k,type,nx,ny,nimages,half_nx,half_ny;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-noshift") == 0) noshift = 1;
    else if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-r") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -r requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_real = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-i") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -i requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_imag = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-m") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -m requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_modulus = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-p") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -p requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_phasis = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-l") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -l requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_logmod = argv[k+1]; k++; }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

  /*********************/
  /* check consistency */
  /*********************/
  if (NULL == fname_in) {
    display_usage();
    printf("Error: input 'in' is missing\n\n");
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

  if(!fname_real && !fname_imag && !fname_modulus && !fname_phasis && !fname_logmod) {
    display_usage();
    printf("Error: you must select at least one ouput option (-r, -i, -m, -p or -l)\n\n");
    return EXIT_FAILURE;    
  }

  /*************************/
  /* load input TIFF image */
  /*************************/
  if(EXIT_FAILURE == load_tiff_into_double(&u,&nx,&ny,&nimages,fname_in,(char)0)){
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }

  /*********************/
  /* memory allocation */
  /*********************/
  half_nx = nx/2; // integer division
  half_ny = ny/2; // integer division
  ASSERT_ALLOC(dft_u = (fftw_complex*) fftw_malloc ((1+half_nx)*ny*sizeof(fftw_complex)));

  if(fname_real) {
    ASSERT_ALLOC(re = (double**) malloc(nimages*sizeof(double*)));
    for(k=0;k<nimages;k++) {
      ASSERT_ALLOC(re[k] = (double*) malloc(nx*ny*sizeof(double)));
    }
  }
  
  if(fname_imag) {
    ASSERT_ALLOC(im = (double**) malloc(nimages*sizeof(double*)));
    for(k=0;k<nimages;k++) {
      ASSERT_ALLOC(im[k] = (double*) malloc(nx*ny*sizeof(double)));
    }
  }
  
  if(fname_modulus) {
    ASSERT_ALLOC(rho = (double**) malloc(nimages*sizeof(double*)));
    for(k=0;k<nimages;k++) {
      ASSERT_ALLOC(rho[k] = (double*) malloc(nx*ny*sizeof(double)));
    }
  }
  
  if(fname_phasis) {
    ASSERT_ALLOC(theta = (double**) malloc(nimages*sizeof(double*)));
    for(k=0;k<nimages;k++) {
      ASSERT_ALLOC(theta[k] = (double*) malloc(nx*ny*sizeof(double)));
    }
  }
  
  if(fname_logmod) {
    ASSERT_ALLOC(logmod = (double**) malloc(nimages*sizeof(double*)));
    for(k=0;k<nimages;k++) {
      ASSERT_ALLOC(logmod[k] = (double*) malloc(nx*ny*sizeof(double)));
    }
  }
    
  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  for(k=0;k<nimages;k++){

    // compute the DFT of u[k] (compute only half of the DFT coefficients)
    ASSERT_ALLOC(plan_u = fftw_plan_dft_r2c_2d(ny,nx,u[k],dft_u,FFTW_ESTIMATE));
    fftw_execute(plan_u);
    fftw_destroy_plan(plan_u);

    // compute output image(s) 
    for(y=0,b=-half_ny;b<=(ny-1)/2;b++,y++)
      for(x=0,a=-half_nx;a<=(nx-1)/2;a++,x++) {
	s = (a >= 0) ? 1 : -1;
	adr = ((s*b+ny)%ny)*(1+half_nx)+(s*a+nx)%nx; // location of frequency (|a|,b) in dft_u
	adr2 = noshift ? (ny+b)%ny*nx+(nx+a)%nx : y*nx+x; // location of frequency (a,b) in output image(s)
	val_re = dft_u[adr][0];
	val_im = s*dft_u[adr][1];
	val_rho = hypot(val_re,val_im); 
	if(re) re[k][adr2] = val_re;
	if(im) im[k][adr2] = val_im;
	if(rho) rho[k][adr2] = val_rho;
	if(theta) theta[k][adr2] = atan2(val_im,val_re);
	if(logmod) logmod[k][adr2] = log(1.+val_rho); 
      }

  }

  /*****************************/
  /* Save output TIFF image(s) */
  /*****************************/
  err = 0; 
  for(k=0;k<nimages;k++) {
    if(!err && re && EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_real,re[k],nx,ny,0,type,((k==0)?"w":"a"))) { printf("Error: failed to save output image '%s'\n",fname_real); err = 1; }
    if(!err && im && EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_imag,im[k],nx,ny,0,type,((k==0)?"w":"a"))) { printf("Error: failed to save output image '%s'\n",fname_imag); err = 1; }
    if(!err && rho && EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_modulus,rho[k],nx,ny,0,type,((k==0)?"w":"a"))) { printf("Error: failed to save output image '%s'\n",fname_modulus); err = 1; }
    if(!err && theta && EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_phasis,theta[k],nx,ny,0,type,((k==0)?"w":"a"))) { printf("Error: failed to save output image '%s'\n",fname_phasis); err = 1; }
    if(!err && logmod && EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_logmod,logmod[k],nx,ny,0,type,((k==0)?"w":"a"))) { printf("Error: failed to save output image '%s'\n",fname_logmod); err = 1; }
  }

  /***************/
  /* free memory */
  /***************/
  fftw_free(dft_u);
  for(k=0;k<nimages;k++) {
    free(u[k]);
    if(re) free(re[k]);
    if(im) free(im[k]);
    if(rho) free(rho[k]);
    if(theta) free(theta[k]);
    if(logmod) free(logmod[k]);
  }
  free(u);
  if(re) free(re);
  if(im) free(im);
  if(rho) free(rho);
  if(theta) free(theta);
  if(logmod) free(logmod);

  return err ? EXIT_FAILURE : EXIT_SUCCESS;
}
