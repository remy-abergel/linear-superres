#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_fftw_complex(fftw_complex***,int*,int*,int*,char*,char); // see source file 'tiffreadcomplex.c'
extern int save_monopage_striped_tiff_from_fftw_complex(char*,fftw_complex*,int,int,char,int,char*,char*,char*); // see source file 'tiffwritecomplex.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'
extern int getasciitranslations(double**,double**,int*,char*,char); // see source file 'ascii.c'
extern int complex_shannon_translation(fftw_complex*,fftw_complex*,double,double,int,int); // see source file 'tools_kernel.c'
extern int fuse_stack(double*,double*,fftw_complex**,int,int,int); // see source file 'tools_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nRegister a stack of TIFF images from a sequence of registration displacements.\n\n");
  printf("Usage: register-stack [-ftype type] [-r reg] [-a av] [-m med] [-v] in T\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -r reg      : output registered sequence (TIFF-stack)\n");
  printf(" -a av       : output the pixelwise average of the registered sequence (shift-and-add)\n");
  printf(" -m med      : output the pixelwise median of the registered sequence (shift-and-median)\n");
  printf(" -v          : enable verbose mode\n");
  printf(" in          : input sequence of low-resolution images (TIFF-stack)\n");
  printf(" T           : input sequence of displacements (ASCII format)\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_T=NULL,*fname_reg=NULL,*fname_av=NULL,*fname_med=NULL,vflag=0;
  double *dx=NULL,*dy=NULL;
  fftw_complex **u0=NULL,**reg=NULL;
  double *av=NULL,*med=NULL;
  int err,k,type,nx,ny,nimages,ntranslations;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-v") == 0) vflag = 1;
    else if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-r") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -r requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_reg = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-a") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -a requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_av = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-m") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -m requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_med = argv[k+1]; k++; }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_T) fname_T = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == fname_in) || (NULL == fname_T)) {
    display_usage();
    if (NULL == fname_in) printf("Error: input 'in' is missing\n\n");
    else printf("Error: input 'T' is missing\n\n");
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

  if(!fname_reg && !fname_av && !fname_med) { display_usage(); printf("Error: you must specify at least one output (options -r, -a or -m).\n\n"); return EXIT_FAILURE; }

  /*******************************************************************/
  /* load input TIFF-Stack image (sequence of low-resolution images) */
  /*******************************************************************/
  if(EXIT_FAILURE == load_tiff_into_fftw_complex(&u0,&nx,&ny,&nimages,fname_in,0)){
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
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
  if(fname_av) ASSERT_ALLOC(av = (double*) malloc (nx*ny*sizeof(double)));
  if(fname_med) ASSERT_ALLOC(med = (double*) malloc (nx*ny*sizeof(double)));
  ASSERT_ALLOC(reg = (fftw_complex**) malloc (nimages*sizeof(fftw_complex*)));
  for(k=0;k<nimages;k++) {
    ASSERT_ALLOC(reg[k] = (fftw_complex*) fftw_malloc (nx*ny*sizeof(fftw_complex)));
  }

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  if(vflag) {
    printf("Registration of the input sequence of images:\n");
    printf("=============================================\n");
    fflush(stdout);
  }
  for(k=0;k<nimages;k++) {
    complex_shannon_translation(reg[k],u0[k],dx[k],dy[k],nx,ny);
    if(vflag) printf("\33[2K\rimage %d/%d: done",k+1,nimages);
    fflush(stdout);
  }
  if(vflag) {
    printf("\n");
  }
  if(fname_av || fname_med) {
    fuse_stack(med,av,reg,nx,ny,nimages);
    if(vflag && fname_av) printf("compute the average of the registered stack (shift-and-add): done\n");
    if(vflag && fname_med) printf("compute the median of the registered stack (shift-and-median): done\n");
  }
  if(vflag) {
    printf("\n");
  }


  /********************************************/
  /* Save output (high-resolution) TIFF image */
  /********************************************/
  if(vflag) {
    printf("Save outputs TIFF images and free memory:\n");
    printf("=========================================\n");
  }

  // output reg: registered stack
  if(fname_reg) {
    for(k=0;k<nimages;k++) {
      if(EXIT_FAILURE == save_monopage_striped_tiff_from_fftw_complex(fname_reg,reg[k],nx,ny,0,type,((k==0)?"w":"a"),NULL,NULL)) {
	printf("Error: failed to save output image '%s'\n",fname_reg);
	for(k=0;k<nimages;k++) { fftw_free(u0[k]); fftw_free(reg[k]); }
	free(u0); free(reg); free(dx); free(dy);
	if(fname_med) free(med);
	if(fname_av) free(av);
	return EXIT_FAILURE;
      }
      if(vflag) {
	printf("\33[2K\rsave '%s' page %d/%d (precision %s): done",fname_reg,k+1,nimages,((NULL == datatype)?"float":datatype));
	fflush(stdout);
      }
    }
    if(vflag) printf("\n");
  }

  // output av (shift-and-add): average of the registered stack along the temporal axis
  if(fname_av) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_av,av,nx,ny,0,type,"w")) {
      printf("Error: failed to save output image '%s'\n",fname_av);
      for(k=0;k<nimages;k++) { fftw_free(u0[k]); fftw_free(reg[k]); }
      free(u0); free(reg); free(dx); free(dy); free(av);
      if(fname_med) free(med);
      return EXIT_FAILURE;
    }
    if(vflag) {
      printf("save '%s' (precision %s): done\n",fname_av,((NULL == datatype)?"float":datatype));
      fflush(stdout);
    }
  }

  // output med (shift-and-median): median of the registered stack along the temporal axis
  if(fname_med) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_med,med,nx,ny,0,type,"w")) {
      printf("Error: failed to save output image '%s'\n",fname_med);
      for(k=0;k<nimages;k++) { fftw_free(u0[k]); fftw_free(reg[k]); }
      free(u0); free(reg); free(dx); free(dy); free(med);
      if(fname_av) free(av);
      return EXIT_FAILURE;
    }
    if(vflag) {
      printf("save '%s' (precision %s): done\n",fname_med,((NULL == datatype)?"float":datatype));
      fflush(stdout);
    }
  }

  if(vflag) printf("\n");

  /***************/
  /* free memory */
  /***************/
  for(k=0;k<nimages;k++) { fftw_free(u0[k]); fftw_free(reg[k]); }
  free(u0);
  free(reg);
  free(dx);
  free(dy);
  if(fname_av) free(av);
  if(fname_med) free(med);

  return EXIT_SUCCESS;
}
