#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_tiff_into_double(double***,int*,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'

static void display_usage()
{
  printf("\nExtract a subpart of a TIFF (monopage or multipage) image.\n\n");
  printf("Usage: tiffextract [-ftype type] [-i0 i0] [-i1 i1] [-r] in out X1 Y1 X2 Y2\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -i0 i0      : (default 0) starting index for multipage extraction\n");
  printf(" -i1 i1      : (default 0) ending index for multipage extraction\n");
  printf(" -r          : if set, X2 and Y2 must be the SIZE of the extracted region\n");
  printf(" in          : input image (TIFF)\n");
  printf(" out         : output image (TIFF)\n");
  printf(" X1          : upleft corner of the region to extract from input (x)\n");
  printf(" Y1          : upleft corner of the region to extract from input (y)\n");
  printf(" X2          : downright corner of the region to extract from input (x)\n");
  printf(" Y2          : downright corner of the region to extract from input (y)\n\n");
}


int main(int argc, char **argv)
{
  char *filename_in=NULL,*filename_out=NULL,*X1value=NULL,*X2value=NULL,*Y1value=NULL,*Y2value=NULL,*datatype=NULL,*i0value=NULL,*i1value=NULL,rflag=0;
  double **in=NULL,*out=NULL;
  int i,nx,ny,nz,X1,X2,Y1,Y2,i0,i1,x,y,err,nx_out,ny_out,type=0,vflag=0;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(i=1;i<argc;i++) {
    if (strcmp(argv[i],"-r") == 0) rflag = 1;
    else if (strcmp(argv[i],"-ftype") == 0) {
      if(i==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[i+1]; i++; }
    }
    else if (strcmp(argv[i],"-i0") == 0) {
      if(i==argc-1) { display_usage(); printf("Error: option -i0 requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	i0value = argv[i+1]; err = sscanf(i0value,"%d",&i0); i++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -i0\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[i],"-i1") == 0) {
      if(i==argc-1) { display_usage(); printf("Error: option -i0 requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	i1value = argv[i+1]; err = sscanf(i1value,"%d",&i1); i++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -i0\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == filename_in) filename_in = argv[i];
    else if (NULL == filename_out) filename_out = argv[i];
    else if (NULL == X1value) {
      X1value = argv[i]; err = sscanf(X1value,"%d",&X1);
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input 'X1'\n\n"); return EXIT_FAILURE; }
    }
    else if (NULL == Y1value) {
      Y1value = argv[i]; err = sscanf(Y1value,"%d",&Y1);
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input 'Y1'\n\n"); return EXIT_FAILURE; }
    }
    else if (NULL == X2value) {
      X2value = argv[i]; err = sscanf(X2value,"%d",&X2);
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input 'X2'\n\n"); return EXIT_FAILURE; }
    }
    else if (NULL == Y2value) {
      Y2value = argv[i]; err = sscanf(Y2value,"%d",&Y2);
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input 'Y2'\n\n"); return EXIT_FAILURE; }
    }
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /**********************/
  /* set default values */
  /**********************/
  if (NULL == i0value) i0 = 0.;
  if (NULL == i1value) i1 = 0.;
  if(i0>i1) { printf("Error: i0 must be smaller or equal to i1\n\n"); return EXIT_FAILURE; }
  if(i0<0) { printf("Error: i0 must be nonnegative\n\n"); return EXIT_FAILURE; }

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == filename_in) || (NULL == filename_out) || (NULL == X1value) || (NULL == X2value) || (NULL == Y1value) || (NULL == Y2value)) {
    display_usage();
    if (NULL == filename_in) printf("Error: input 'in' is missing\n\n");
    else if (NULL == filename_out) printf("Error: output 'out' is missing\n\n");
    else if (NULL == X1value) printf("Error: input 'X1' is missing\n\n");
    else if (NULL == X2value) printf("Error: input 'X2' is missing\n\n");
    else if (NULL == Y1value) printf("Error: input 'Y1' is missing\n\n");
    else printf("Error: input 'Y2' is missing\n\n");
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

  if ((!rflag)&&(X2<=X1)) { printf("\nError: X2 <= X1 is not allowed (except using option -r).\n\n"); return EXIT_FAILURE; }
  if ((!rflag)&&(Y2<=Y1)) { printf("\nError: X2 <= X1 is not allowed (except using option -r).\n\n"); return EXIT_FAILURE; }
  if (X1<0) { printf("\nError: input 'X1' must be nonnegative.\n\n"); return EXIT_FAILURE; }
  if (X2<0) { printf("\nError: input 'X2' must be nonnegative.\n\n"); return EXIT_FAILURE; }
  if (Y1<0) { printf("\nError: input 'Y1' must be nonnegative.\n\n"); return EXIT_FAILURE; }
  if (Y2<0) { printf("\nError: input 'Y2' must be nonnegative.\n\n"); return EXIT_FAILURE; }

  /********************************************************************/
  /* Open a TIFF-Stack image & get its graylevels in double precision */
  /********************************************************************/
  if(EXIT_FAILURE == load_tiff_into_double(&in,&nx,&ny,&nz,filename_in,vflag)) {
    printf("An error occured during the TIFF image reading process.\n");
    printf("Free memory and exit the program.\n\n");
    return EXIT_FAILURE;
  }
  if(i1>=nz) {
    printf("Error: inconsistent setting for i1, you entered i1=%d but you must select i1 < number of pages of the input image '%s' (%d pages detected)\n\n",i1,filename_in,nz);
    for(i=0;i<nz;i++) free(in[i]);
    free(in);
    return EXIT_FAILURE;
  }

  /******************************************/
  /* memory allocation for the output image */
  /******************************************/
  nx_out = (rflag) ? X2 : X2-X1+1;
  ny_out = (rflag) ? Y2 : Y2-Y1+1;
  ASSERT_ALLOC(out = (double*) calloc (nx_out*ny_out,sizeof(double)));

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  for(i=i0;i<=i1;i++){
    for(x=X1;x<=fmin(((rflag) ? X1+nx_out-1 : X2),nx-1);x++)
      for(y=Y1;y<=fmin(((rflag) ? Y1+ny_out-1 : Y2),ny-1);y++) {
	out[(y-Y1)*nx_out + (x-X1)] = in[i][y*nx+x];
      }
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(filename_out,out,nx_out,ny_out,vflag,type,((i==i0)?"w":"a"))) {
      printf("An error occured during the TIFF image saving process.\n");
      printf("Free memory and exit the program.\n\n");
      for(i0=0;i0<nz;i0++) free(in[i0]);
      free(in); free(out);
      return EXIT_FAILURE;
    }

  }

  /***************/
  /* free memory */
  /***************/
  for(i=0;i<nz;i++) free(in[i]);
  free(in); free(out);

  return EXIT_SUCCESS;
}
