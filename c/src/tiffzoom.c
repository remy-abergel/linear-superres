#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_fftw_complex(fftw_complex***,int*,int*,int*,char*,char); // see source file 'tiffreadcomplex.c'
extern int save_monopage_striped_tiff_from_fftw_complex(char*,fftw_complex*,int,int,char,int,char*,char*,char*); // see source file 'tiffwritecomplex.c'
extern int complex_shannon_zooming(fftw_complex*,fftw_complex*,int,int,int,int); // see source file 'tools_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nZoom a (monopage or multipage) TIFF image using the Shannon interpolation.\n\n");
  printf("Usage: tiffzoom [-ftype type] [-zx zx] [-zy zy] [-z z] in out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -zx zx      : (default zx=2.) zooming factor along the horizontal axis (double >= 1.)\n");
  printf(" -zz zy      : (default zy=2.) zooming factor along the vertical axis (double >= 1.)\n");
  printf(" -z z        : set zx = zy = z\n");
  printf(" in          : input (monopage or multipage) TIFF image\n");
  printf(" out         : output (monopage or multipage) TIFF image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_out=NULL,*zx_value=NULL,*zy_value=NULL,*z_value=NULL;
  fftw_complex **in=NULL,**out=NULL;
  double zx,zy,z,zx_tmp,zy_tmp;
  int err,k,type,nx,ny,nimages,Nx,Ny;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-zx") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -zx requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	zx_value = argv[k+1]; err = sscanf(zx_value,"%lf",&zx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -zx\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-zy") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -zy requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	zy_value = argv[k+1]; err = sscanf(zy_value,"%lf",&zy); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -zy\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-z") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -z requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	z_value = argv[k+1]; err = sscanf(z_value,"%lf",&z); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -z\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /**********************/
  /* set default values */
  /**********************/
  if(NULL == zx_value) zx = 2.;
  if(NULL == zy_value) zy = 2.;

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == fname_in) || (NULL == fname_out)) {
    display_usage();
    if (NULL == fname_in) printf("Error: input 'in' is missing\n\n");
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

  if(zx_value && z_value) { display_usage(); printf("Error: options -z and -zx cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(zy_value && z_value) { display_usage(); printf("Error: options -z and -zy cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(z_value) zx = zy = z;

  if(zx < 1.){ display_usage(); printf("Error: input 'zx' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(zy < 1.){ display_usage(); printf("Error: input 'zy' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }

  /*******************************************************************/
  /* load input TIFF-Stack image (sequence of low-resolution images) */
  /*******************************************************************/
  if(EXIT_FAILURE == load_tiff_into_fftw_complex(&in,&nx,&ny,&nimages,fname_in,0)){
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }

  /*******************************************/
  /* retrieve & control I/O image dimensions */
  /*******************************************/
  Nx = (int)round(nx*zx);
  Ny = (int)round(ny*zy);
  zx_tmp = zx; zy_tmp = zy;
  zx = (double)Nx/(double)nx;
  zy = (double)Ny/(double)ny;
  if(zx_tmp != zx) printf("WARNING: zx has been changed into %g (instead of zx=%g) to allow integer width of the output image\n",zx,zx_tmp);
  if(zy_tmp != zy) printf("WARNING: zy has been changed into %g (instead of zy=%g) to allow integer height of the output image\n",zy,zy_tmp);

  /*********************/
  /* memory allocation */
  /*********************/
  ASSERT_ALLOC(out = (fftw_complex**) malloc (nimages*sizeof(fftw_complex*)));
  for(k=0;k<nimages;k++) { ASSERT_ALLOC(out[k] = (fftw_complex*) fftw_malloc (Nx*Ny*sizeof(fftw_complex))); }

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  for(k=0;k<nimages;k++) complex_shannon_zooming(out[k],in[k],nx,ny,Nx,Ny);

  /********************************************/
  /* Save output (high-resolution) TIFF image */
  /********************************************/
  for(k=0;k<nimages;k++) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_fftw_complex(fname_out,out[k],Nx,Ny,0,type,((k==0)?"w":"a"),NULL,NULL)) {
      printf("Error: failed to save output image '%s'\n",fname_out);
      for(k=0;k<nimages;k++) { fftw_free(in[k]); fftw_free(out[k]); }
      free(in); free(out);
      return EXIT_FAILURE;
    }
  }

  /***************/
  /* free memory */
  /***************/
  for(k=0;k<nimages;k++) { fftw_free(in[k]); fftw_free(out[k]); }
  free(in); free(out);

  return EXIT_SUCCESS;
}
