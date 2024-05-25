#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'
extern int sharpening(double*,double*,double,int,int); // see source file 'sharpening_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nImage sharpening by frequency amplification\n\n");
  printf("Usage: sharpening [-ftype type] [-l lambda] in out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -l lambda   : (default 5.) amplification parameter\n");
  printf(" in          : input (monopage or multipage) TIFF image\n");
  printf(" out         : output TIFF (same size as input) image\n\n");  
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_out=NULL,*l_value=NULL;
  double **u,lambda;
  int err,k,type,nx,ny,nimages,adr;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-l") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -l requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	l_value = argv[k+1]; err = sscanf(l_value,"%lf",&lambda); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -l\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /**********************/
  /* set default values */
  /**********************/
  if(NULL == l_value) lambda = 5.;

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
  for(k=0;k<nimages;k++) {
    sharpening(u[k],u[k],lambda,nx,ny); 
  }

  /***************************/
  /* Save output TIFF images */
  /***************************/

  // save output 'out'
  for(k=0;k<nimages;k++) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_out,u[k],nx,ny,0,type,((k==0)?"w":"a"))) {
      printf("Error: failed to save output image '%s'\n",fname_out);
      for(adr=0;adr<nimages;adr++) { free(u[adr]); }
      free(u);
      return EXIT_FAILURE;
    }
  }

  /***************/
  /* free memory */
  /***************/
  for(k=0;k<nimages;k++) { free(u[k]); }
  free(u);

  return EXIT_SUCCESS;
}
