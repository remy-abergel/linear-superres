#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_monopage_tiff_into_fftw_complex(fftw_complex**,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_fftw_complex(char*,fftw_complex*,int,int,char,int,char*,char*,char*); // see source file 'tiffwritecomplex.c'
extern int getasciitranslations(double**,double**,int*,char*,char); // see source file 'ascii.c'
extern int gendataset(fftw_complex*, fftw_complex***, fftw_complex**, double*, double*, int*, int*, int*, int*, int, double*, double*, int, int, int, int, char, char); // see source file 'gendataset_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nCompute a realistic low-resolution stack from a high resolution image (avoiding periodic-like boundaries).\n\n");
  printf("Usage: gendataset [-ftype type] [-zx zx] [-zy zy] [-v] [-r ref] in T out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -cx cx      : (default 10) maximal number of columns to be discarded from\n");
  printf("               the input image during the first cropping step (int >= 0))\n");
  printf(" -cy cy      : (default 10) maximal number of rows to be discarded from\n");
  printf("               the input image during the first cropping step (int >= 0))\n");
  printf(" -p p        : (default 5) the final high-resolution ROI size will be equal\n");
  printf("               to (q/p) * the dimensions of the cropped HR image (int >= 1)\n");
  printf(" -q q        : (default 4) see above (int in [1,p-1])\n");
  printf(" -zx zx      : (default 2.) X-axis subsampling factor (double >= 1.)\n");
  printf(" -zy zy      : (default 2.) Y-axis subsampling factor (double >= 1.)\n");
  printf(" -v          : enable verbose mode\n");
  printf(" -r ref      : output high resolution reference image (cropping of in)\n");
  printf(" in          : input high-resolution TIFF image (large domain)\n");
  printf(" T           : input sequence of displacements (ASCII format)\n");
  printf(" out         : output realistic sequence of low-resolution images (TIFF-stack)\n");
  printf("               with restricted domain to get rid of periodization artifacts\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_T=NULL,*fname_out=NULL,*fname_ref=NULL,*zx_value=NULL,*zy_value=NULL,*p_value=NULL,*q_value=NULL,*cx_value=NULL,*cy_value=NULL,vflag=0;
  double *dx=NULL,*dy=NULL,zx=2.,zy=2.;
  fftw_complex **u0=NULL,*utmp=NULL,*ref=NULL;
  int err,k,type,nx,ny,nimages,ntranslations,Nx,Ny,cx=10,cy=10,p=5,q=4;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-v") == 0) vflag = 1;
    else if (strcmp(argv[k],"-ftype") == 0) {
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
    else if (strcmp(argv[k],"-p") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -p requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	p_value = argv[k+1]; err = sscanf(p_value,"%d",&p); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -p\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-q") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -q requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	q_value = argv[k+1]; err = sscanf(q_value,"%d",&q); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -q\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-cx") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -cx requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	cx_value = argv[k+1]; err = sscanf(cx_value,"%d",&cx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -cx\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-cy") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -cy requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	cy_value = argv[k+1]; err = sscanf(cy_value,"%d",&cy); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -cy\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-r") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -r requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_ref = argv[k+1]; k++; }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_T) fname_T = argv[k];
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

  if(zx < 1.) { display_usage(); printf("Error: input 'zx' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(zy < 1.) { display_usage(); printf("Error: input 'zy' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(cx < 0) { display_usage(); printf("Error: input 'cx' must be greater than or equal to zero.\n\n"); return EXIT_FAILURE; }
  if(cy < 0) { display_usage(); printf("Error: input 'cy' must be greater than or equal to zero.\n\n"); return EXIT_FAILURE; }
  if(p <= 0) { display_usage(); printf("Error: input 'p' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(q < 1 || q >= p) { display_usage(); printf("Error: input 'q' must be greater than in the range [1,p-1].\n\n"); return EXIT_FAILURE; }
  
  /*****************************************************************/
  /* load input TIFF image & perform dimensions consistency checks */
  /*****************************************************************/
  if(EXIT_FAILURE == load_monopage_tiff_into_fftw_complex(&utmp,&Nx,&Ny,fname_in,0)) {
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }
  
  if(Nx <= 1) {
    printf("Error: the width of the input image 'in' must be at least 2\n\n");
    fftw_free(utmp); 
    return EXIT_FAILURE; 
  } 
  if(Ny <= 1) {
    printf("Error: the height of the input image 'in' must be at least 2\n\n");
    fftw_free(utmp); 
    return EXIT_FAILURE; 
  } 
  
  /**************************************/
  /* load the sequence of displacements */
  /**************************************/
  getasciitranslations(&dx,&dy,&ntranslations,fname_T,vflag);
  nimages = ntranslations;
  
  /*********************************************************************/
  /* CORE OF THE ROUTINE : compute realistic low-resolution stack (u0) */
  /* from the provided input HR image                                  */
  /*********************************************************************/
  gendataset(utmp,&u0,&ref,dx,dy,&Nx,&Ny,&nx,&ny,nimages,&zx,&zy,cx,cy,p,q,fname_ref!=NULL,vflag);
  
  /**************************************/
  /* Save output (multipage) TIFF image */
  /**************************************/
  if(vflag) {
    printf("Save output images\n");
    printf("==================\n");
  }
  if(ref) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_fftw_complex(fname_ref,ref,Nx,Ny,0,type,"w",NULL,NULL)) {
      printf("Error: failed to save output image '%s' (precision %s)\n",fname_ref,((NULL == datatype)?"float":datatype));
      for(k=0;k<nimages;k++) fftw_free(u0[k]);
      fftw_free(ref); fftw_free(utmp);
      free(u0); free(dx); free(dy);
      return EXIT_FAILURE;
    }
    if(vflag) {
      printf("save reference high-resolution image '%s': done\n",fname_ref);
    }
  }
  for(k=0;k<nimages;k++) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_fftw_complex(fname_out,u0[k],nx,ny,0,type,((k==0)?"w":"a"),NULL,NULL)) {
      printf("Error: failed to save output image '%s'\n",fname_out);
      for(k=0;k<nimages;k++) fftw_free(u0[k]);
      if(ref) fftw_free(ref);
      fftw_free(utmp);
      free(u0); free(dx); free(dy);
      return EXIT_FAILURE;
    }
    if(vflag) {
      printf("\33[2K\rsave low-resolution stack '%s' page %d/%d (precision %s): done",fname_out,k+1,nimages,((NULL == datatype)?"float":datatype));
      fflush(stdout);
    }
  }
  if(vflag) printf("\n\n");
  
  /***************/
  /* free memory */
  /***************/
  for(k=0;k<nimages;k++) { fftw_free(u0[k]); }
  if(ref) fftw_free(ref);
  free(u0);
  fftw_free(utmp);
  free(dx);
  free(dy);
  
  return EXIT_SUCCESS;
}
