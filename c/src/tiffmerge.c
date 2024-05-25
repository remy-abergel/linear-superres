#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'

static void display_usage()
{
  printf("\nMerge a set of TIFF images into a multipage TIFF image (TIFF-stack).\n\n");
  printf("Usage: tiffmerge [-ftype type] [-v] out in1 in2 ...\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -v          : enable verbose mode\n");
  printf(" out         : output TIFF stack\n");
  printf(" in1,in2,... : input (monopage or multipage) TIFF images\n\n");
}


int main(int argc, char **argv)
{
  char *filename_out = NULL,*datatype = NULL,vflag = 0;
  double **u;
  int i,id,id2,nx,ny,nz,nz2,type=0,page_id;

  /******************************************/
  /*            arguments parser            */
  /******************************************/

  /* main parser */
  for(i=1;i<argc;i++) {
    if (strcmp(argv[i],"-v") == 0) vflag = 1;
    else if (strcmp(argv[i],"-ftype") == 0) {
      if(i==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[i+1]; i++; }
    }
    else if (NULL == filename_out) { filename_out = argv[i]; break; }
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }
  nz = argc - i - 1;

  /* check consistency */
  if ((NULL == filename_out) || (nz < 1)) {
    display_usage();
    if (NULL == filename_out) printf("Error: output 'out' is missing\n\n");
    else printf("Error: at least one input image is required\n\n");
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

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  for(page_id=id=0;id<nz;id++) {
    if(EXIT_FAILURE == load_tiff_into_double(&u,&nx,&ny,&nz2,argv[i+1+id],0)) return EXIT_FAILURE;
    for(id2=0;id2<nz2;id2++,page_id++) {
      if(vflag) printf("Appending '%s' page (%d/%d):",argv[i+1+id],id2+1,nz2);
      if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(filename_out,u[id2],nx,ny,0,type,(page_id==0)?"w":"a")) return EXIT_FAILURE;
      if(vflag) printf("\33[2K\rAppending '%s' page (%d/%d): done\n",argv[i+1+id],id2+1,nz2);
      free(u[id2]);
    }
    free(u);
  }

  return EXIT_SUCCESS;
}
