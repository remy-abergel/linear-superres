#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'

/* usage displayer */
static void display_usage()
{
  printf("\nPrint the graylevels of a TIFF image in ascii format (row major order).\n\n");
  printf("Usage: tiffprintasc [-s s] in\n\n");
  printf("   -s s : only print graylevels of slice with page index s\n");
  printf("          (the TIFF pages are numbered from 0)\n");
  printf("   in   : intput TIFF image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *fname_in=NULL,*s_value=NULL;
  double **u;
  int err,k,nx,ny,nimages,adr,s;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-s") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -s requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	s_value = argv[k+1]; err = sscanf(s_value,"%d",&s); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -s\n\n"); return EXIT_FAILURE; }
      }
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

  /*************************/
  /* load input TIFF image */
  /*************************/
  if(EXIT_FAILURE == load_tiff_into_double(&u,&nx,&ny,&nimages,fname_in,(char)0)){
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }

  /*******************************/
  /* printf the graylevel values */
  /*******************************/
  for(k=(s_value)?s:0;k<((s_value)?s+1:nimages);k++) {
    for(adr=0;adr<nx*ny;adr++) printf("%.17g\n",u[k][adr]);
  }

  /***************/
  /* free memory */
  /***************/
  for(adr=0;adr<nimages;adr++) free(u[adr]);
  free(u);

  return EXIT_SUCCESS;
}
