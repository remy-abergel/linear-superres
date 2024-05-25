#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'

/* usage displayer */
static void display_usage()
{
  printf("\nCreate a monopage or multipage TIFF image with constant gray level value.\n\n");
  printf("Usage: tiffaxpb [-ftype type] [-g g] [-L L] in width height out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -g g        : (default 0.) gray level value to set\n");
  printf(" -L L        : (default 1) number of TIFF pages to create\n");
  printf(" width       : width of the output image\n");
  printf(" height      : height of the output image\n");
  printf(" out         : output TIFF image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_out=NULL,*g_value=NULL,*L_value=NULL,*nx_value=NULL,*ny_value=NULL;
  double **u,g;
  int err,k,adr,type,nx,ny,L;
  
  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-g") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -g requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	g_value = argv[k+1]; err = sscanf(g_value,"%lf",&g); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -g\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-L") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -L requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	L_value = argv[k+1]; err = sscanf(L_value,"%d",&L); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -L\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == nx_value) {
      nx_value = argv[k]; err = sscanf(nx_value,"%d",&nx); 
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input argument 'width'\n\n"); return EXIT_FAILURE; }
    }
    else if (NULL == ny_value) {
      ny_value = argv[k]; err = sscanf(ny_value,"%d",&ny);
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input argument 'height'\n\n"); return EXIT_FAILURE; }
    }
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }
  
  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == nx_value) || (NULL == ny_value) || (NULL == fname_out)) {
    display_usage();
    if (NULL == nx_value) printf("Error: input 'width' is missing\n\n");
    else if (NULL == ny_value) printf("Error: input 'height' is missing\n\n");
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
  
  if(NULL == g_value) g = 0.; 
  if(NULL == L_value) L = 1; 

  /*********************/
  /* memory allocation */
  /*********************/
  ASSERT_ALLOC(u = (double**) malloc(L*sizeof(double*)));
  for(k=0;k<L;k++) { ASSERT_ALLOC(u[k] = (double*) malloc(nx*ny*sizeof(double))); }
  
  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  for(k=0;k<L;k++)
    for(adr=0;adr<nx*ny;adr++) u[k][adr] = g; 
  
  /**************************/
  /* Save output TIFF image */
  /**************************/
  for(k=0;k<L;k++) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_out,u[k],nx,ny,0,type,k==0?"w":"a")) {
      printf("Error: failed to save output image '%s'\n",fname_out);
      free(u);
      return EXIT_FAILURE;
    }
  }
  
  /***************/
  /* free memory */
  /***************/
  for(k=0;k<L;k++) free(u[k]); 
  free(u);
  
  return EXIT_SUCCESS;
}
