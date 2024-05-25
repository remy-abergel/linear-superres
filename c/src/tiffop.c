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
  printf("\nPerform an elementary operation between two TIFF (monopage or multipage) images.\n\n");
  printf("Usage: tiffop [-t type] [-A A] [-a a] [-p] [-m] [-t] [-d] [-D] [-N] [-i] [-s] [-l] [-g] [-L] [-G] [-e] B out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -A A        : take TIFF image A as left term\n");
  printf(" -a a        : take constant a as left term\n");
  printf(" -p          : the A + B (plus) operator\n");
  printf(" -m          : the A - B (minus) operator\n");
  printf(" -t          : the A x B (times) operator\n");
  printf(" -d          : the A / B (divide) operator\n");
  printf(" -D          : the |A-B| (distance) operator\n");
  printf(" -N          : the hypot(A,B)=(A^2+B^2)^(1/2) (norm) operator\n");
  printf(" -i          : the inf(A,B) operator\n");
  printf(" -s          : the sup(A,B) operator\n");
  printf(" -l          : the A < B operator (result: 1=true, 0=false)\n");
  printf(" -g          : the A > B operator (result: 1=true, 0=false)\n");
  printf(" -L          : the A <= B operator (result: 1=true, 0=false)\n");
  printf(" -G          : the A >= B operator (result: 1=true, 0=false)\n");
  printf(" -e          : the A == B operator (result: 1=true, 0=false)\n");
  printf(" B           : input TIFF image B (right term)\n");
  printf(" out         : resulting TIFF image\n\n");
}


int main(int argc, char **argv)
{
  char *filename_A=NULL,*filename_B=NULL,*filename_out=NULL,*avalue=NULL,*datatype=NULL,aflag=0,Aflag=0,pflag=0,mflag=0,tflag=0,dflag=0,Dflag=0,Nflag=0,iflag=0,sflag=0,lflag=0,gflag=0,Lflag=0,Gflag=0,eflag=0;
  double **in_A=NULL,**in_B=NULL,*out=NULL,a,res,left;
  int i,id,adr,nx,ny,nz,nx_A,ny_A,nz_A,err,type=0,vflag=0;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(i=1;i<argc;i++) {
    if (strcmp(argv[i],"-p") == 0) pflag = 1;
    else if (strcmp(argv[i],"-m") == 0) mflag = 1;
    else if (strcmp(argv[i],"-t") == 0) tflag = 1;
    else if (strcmp(argv[i],"-d") == 0) dflag = 1;
    else if (strcmp(argv[i],"-D") == 0) Dflag = 1;
    else if (strcmp(argv[i],"-N") == 0) Nflag = 1;
    else if (strcmp(argv[i],"-i") == 0) iflag = 1;
    else if (strcmp(argv[i],"-s") == 0) sflag = 1;
    else if (strcmp(argv[i],"-l") == 0) lflag = 1;
    else if (strcmp(argv[i],"-g") == 0) gflag = 1;
    else if (strcmp(argv[i],"-L") == 0) Lflag = 1;
    else if (strcmp(argv[i],"-G") == 0) Gflag = 1;
    else if (strcmp(argv[i],"-e") == 0) eflag = 1;
    else if (strcmp(argv[i],"-ftype") == 0) {
      if(i==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[i+1]; i++; }
    }
    else if (strcmp(argv[i],"-A") == 0) {
      Aflag = 1;
      if(i==argc-1) { display_usage(); printf("Error: option -A requires an argument\n\n"); return EXIT_FAILURE; }
      else { filename_A = argv[i+1]; i++; }
    }
    else if (strcmp(argv[i],"-a") == 0) {
      aflag = 1;
      if(i==argc-1) { display_usage(); printf("Error: option -a requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	avalue = argv[i+1]; err = sscanf(avalue,"%lf",&a); i++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -a\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == filename_B) filename_B = argv[i];
    else if (NULL == filename_out) filename_out = argv[i];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == filename_B) || (NULL == filename_out)) {
    display_usage();
    if (NULL == filename_B) printf("Error: input 'B' is missing\n\n");
    else printf("Error: output 'out' is missing\n\n");
    return EXIT_FAILURE;
  }
  if(!Aflag && !aflag) {
    display_usage();
    printf("Error: please select at least one left term (-a or -A)\n\n");
  }
  if(Aflag && aflag) {
    display_usage();
    printf("Error: please select exactly one left term (-a or -A)\n\n");
  }
  if(1 != (pflag+mflag+tflag+dflag+Dflag+Nflag+iflag+sflag+lflag+gflag+Lflag+Gflag+eflag)){
    display_usage();
    printf("Error: please select exactly one of the operator options\n\n");
  }

  if((NULL == datatype)||(strcmp(datatype,"float")==0)) type = 0;
  else if (strcmp(datatype,"uint8")==0) type = 8;
  else if (strcmp(datatype,"int8")==0) type = -8;
  else if (strcmp(datatype,"uint16")==0) type = 16;
  else if (strcmp(datatype,"int16")==0) type = -16;
  else if (strcmp(datatype,"uint32")==0) type = 32;
  else if (strcmp(datatype,"int32")==0) type = -32;
  else { display_usage(); printf("Error: unrecognized value for option -ftype.\n\n"); return EXIT_FAILURE; }


  /**************************/
  /* Open input TIFF images */
  /**************************/
  if(EXIT_FAILURE == load_tiff_into_double(&in_B,&nx,&ny,&nz,filename_B,vflag)) {
    printf("An error occured during the TIFF image reading process '%s'.\n",filename_B);
    printf("Free memory and exit the program.\n\n");
    return EXIT_FAILURE;
  }
  if(Aflag){
    if(EXIT_FAILURE == load_tiff_into_double(&in_A,&nx_A,&ny_A,&nz_A,filename_A,vflag)) {
      printf("An error occured during the TIFF image reading process '%s'.\n",filename_A);
      printf("Free memory and exit the program.\n\n");
      return EXIT_FAILURE;
    }
    if((nx_A != nx) || (ny_A != ny) || (nz_A != nz)) {
      printf("Error: the input image A and B must have the same dimensions\n");
      for(id=0;id<nz;id++) free(in_B[id]);
      for(id=0;id<nz_A;id++) free(in_A[id]);
      free(in_A); free(in_B);
    }
  }

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  ASSERT_ALLOC(out = (double*) malloc (nx*ny*sizeof(double*)));
  for(id=0;id<nz;id++){
    for(adr=0;adr<nx*ny;adr++) {
      left = ((Aflag)?in_A[id][adr]:a);
      if      (pflag)    res = left + in_B[id][adr];
      else if (mflag)    res = left - in_B[id][adr];
      else if (tflag)    res = left * in_B[id][adr];
      else if (dflag)    res = left / in_B[id][adr];
      else if (Dflag)    res = fabs(left - in_B[id][adr]);
      else if (Nflag)    res = hypot(left,in_B[id][adr]);
      else if (iflag)    res = (left < in_B[id][adr] ? left : in_B[id][adr]);
      else if (sflag)    res = (left > in_B[id][adr] ? left : in_B[id][adr]);
      else if (lflag)    res = (left < in_B[id][adr] ? 1.0 : 0.0);
      else if (gflag)    res = (left > in_B[id][adr] ? 1.0 : 0.0);
      else if (Lflag)    res = (left <= in_B[id][adr] ? 1.0 : 0.0);
      else if (Gflag)    res = (left >= in_B[id][adr] ? 1.0 : 0.0);
      else if (eflag)    res = (left == in_B[id][adr] ? 1.0 : 0.0);
      out[adr] = res;
    }
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(filename_out,out,nx,ny,0,type,((id==0)?"w":"a"))) {
      printf("An error occured during the TIFF image saving process '%s'\n",filename_out);
      return EXIT_FAILURE;
    }
  }

  /***************/
  /* free memory */
  /***************/
  for(id=0;id<nz;id++) free(in_B[id]);
  if(Aflag) for(id=0;id<nz_A;id++) free(in_A[id]);
  if(Aflag) free(in_A);
  free(in_B);
  free(out);

  return EXIT_SUCCESS;
}
