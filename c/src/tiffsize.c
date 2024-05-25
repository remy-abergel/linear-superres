#include <stdlib.h>
#include <stdio.h>
#include <tiffio.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'

/* usage displayer */
static void display_usage()
{
  printf("\nDisplay the dimensions (width, height, number of pages) of a TIFF image.\n\n");
  printf("Usage: tiffsize in \n\n");
  printf("   in      : input TIFF image\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *filename_in=NULL;
  TIFF *tiffp;
  unsigned int npages = 0;
  int width,height;

  /* main parser */
  if(argc < 2){ display_usage(); printf("Error: input 'in' is missing\n\n"); return EXIT_FAILURE; }
  else if(argc > 2){ display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  filename_in = argv[1];

  /*** kernel ***/

  // open the input TIFF file
  TIFFSetWarningHandler(NULL);
  tiffp = TIFFOpen(filename_in,"r");
  if (!tiffp) {
    printf("Failed to open TIFF File '%s'\n",filename_in);
    return EXIT_FAILURE;
  }

  // retrieve number of pages & TIFF dimensions //
  do ++npages; while (TIFFReadDirectory(tiffp));
  TIFFGetField(tiffp,TIFFTAG_IMAGEWIDTH,&width);
  TIFFGetField(tiffp,TIFFTAG_IMAGELENGTH,&height);

  // close TIFF image & printf results
  TIFFClose(tiffp);
  printf("width = %d\n",width);
  printf("height = %d\n",height);
  printf("number of pages = %d\n",npages);

  return EXIT_SUCCESS;
}
