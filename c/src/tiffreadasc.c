#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

// external modules
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'

/* usage displayer */
static void display_usage()
{
  printf("\nRead a monopage TIFF image in ascii format (row-major order).\n\n");
  printf("Usage: tiffreadasc [-ftype type] out width height\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" out         : output TIFF image (monopage)\n");
  printf(" witdh       : image width\n");
  printf(" height      : image height\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_out=NULL,*nx_value=NULL,*ny_value=NULL;
  double *u;
  int err,k,type,nx,ny,adr;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (NULL == fname_out) fname_out = argv[k];
    else if (NULL == nx_value) nx_value = argv[k];
    else if (NULL == ny_value) ny_value = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == fname_out) || (NULL == nx_value) || (NULL == ny_value)) {
    display_usage();
    if (NULL == fname_out) printf("Error: output 'out' is missing\n\n");
    else if (NULL == nx_value) printf("Error: input 'width' is missing\n\n");
    else printf("Error: input 'height' is missing\n\n");
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

  if(sscanf(nx_value,"%d",&nx) != 1) { display_usage(); printf("Error: could not retrieve properly input 'width'\n\n"); return EXIT_FAILURE; }
  if(sscanf(ny_value,"%d",&ny) != 1) { display_usage(); printf("Error: could not retrieve properly input 'height'\n\n"); return EXIT_FAILURE; }

  /*********************/
  /* memory allocation */
  /*********************/
  ASSERT_ALLOC(u = (double*) malloc (nx*ny * sizeof(double)));

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  for(adr=0;adr<nx*ny;adr++) {
    if(scanf("%lf",u+adr) != 1) {
      printf("\nError: some image graylevels are missing or could not be read properly.\n\n");
      return EXIT_FAILURE;
    }
  }

  /**************************/
  /* Save output TIFF image */
  /**************************/
  if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_out,u,nx,ny,0,type,"w")) {
    printf("Error: failed to save output image '%s'\n",fname_out);
    free(u);
    return EXIT_FAILURE;
  }

  /***************/
  /* free memory */
  /***************/
  free(u);

  return EXIT_SUCCESS;
}
