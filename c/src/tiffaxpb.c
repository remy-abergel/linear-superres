#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'

/* usage displayer */
static void display_usage()
{
  printf("\nGain/Offset correction to a (monopage or multipage) TIFF image: out=a*in+b.\n\n");
  printf("Usage: tiffaxpb [-ftype type] [-a a] [-b b] in out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -v          : display gain and offset values in the standard output\n");
  printf(" -i          : apply the inverse transformation (set out = (in-b)/a)\n");
  printf(" -n          : normalize pixel values from actual to [0,255]\n");
  printf(" [-a a]      : (default 1.) gain value\n");
  printf(" [-b b]      : (default 0.) offset value\n");
  printf(" in          : input TIFF (monopage or multipage) image\n");
  printf(" out         : output TIFF (same size as input) image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_out=NULL,*a_value=NULL,*b_value=NULL,vflag=0,nflag=0,iflag=0;
  double **u,a=1.,b=0.,m,M;
  int err,k,type,nx,ny,nimages,adr;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-i") == 0) iflag = 1;
    else if (strcmp(argv[k],"-v") == 0) vflag = 1; 
    else if (strcmp(argv[k],"-n") == 0) nflag = 1; 
    else if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-a") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -a requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	a_value = argv[k+1]; err = sscanf(a_value,"%lf",&a); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -a\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-b") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -b requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	b_value = argv[k+1]; err = sscanf(b_value,"%lf",&b); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -b\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

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

  if(iflag && nflag) {
    printf("Error: options -i and -n cannot be used simultaneously\n\n");
    return EXIT_FAILURE; 
  }

  if(iflag && a == 0.) {
    printf("Error: you must set a nonzero value for option -a when option -i is enabled\n\n");
    return EXIT_FAILURE; 
  }

  if(nflag && a_value) {
    printf("WARNING : option -a will be ignored because option -n is enabled\n");
  }

  if(nflag && b_value) {
    printf("WARNING : option -b will be ignored because option -n is enabled\n");
  }
    
  /*************************/
  /* load input TIFF image */
  /*************************/
  if(EXIT_FAILURE == load_tiff_into_double(&u,&nx,&ny,&nimages,fname_in,(char)0)){
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/

  // deal with option -n
  if(nflag) {
    // compute min (m) & max (M) of u
    for(m=M=u[0][0],k=0;k<nimages;k++)
      for(adr=0;adr<nx*ny;adr++) {
	m = fmin(m,u[k][adr]); 
	M = fmax(M,u[k][adr]);
      }
    // set gain & offset coefficients to map [m,M] to [0,255.]
    a = 255./(M-m); b = -m*a;
  }

  // deal with -v option
  if(vflag) printf("a = %.17e\nb = %.17e\n",a,b);
  
  // apply contrast change
  for(k=0;k<nimages;k++)
    for(adr=0;adr<nx*ny;adr++)
      u[k][adr] = iflag ? (u[k][adr]-b)/a : a*u[k][adr]+b;

  /**************************/
  /* Save output TIFF image */
  /**************************/
  for(k=0;k<nimages;k++) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_out,u[k],nx,ny,0,type,((k==0)?"w":"a"))) {
      printf("Error: failed to save output image '%s'\n",fname_out);
      for(adr=0;adr<nimages;adr++) free(u[adr]);
      free(u);
      return EXIT_FAILURE;
    }
  }

  /***************/
  /* free memory */
  /***************/
  for(adr=0;adr<nimages;adr++) free(u[adr]);
  free(u);

  return EXIT_SUCCESS;
}
