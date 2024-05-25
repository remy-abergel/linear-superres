#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int load_monopage_tiff_into_double(double**,int*,int*,char*,char); // see source file 'tiffread.c'
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'
extern int metrics(double*,double*,double*,double*,double*,int,int,double*); // see source file 'tools_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nComputes the mean square difference between two TIFF images.\n\n");
  printf("Usage: tiffmse [-p peakval] u v\n\n");
  printf(" -p peakval    : set peak value (e.g., peakval=255 for 8-bits signals) for\n");
  printf("                 the PSNR (otherwise, peakval = max(u)-min(u) is used)\n");
  printf(" u             : first input (monopage) TIFF image\n");
  printf(" v             : second input (monopage) TIFF image\n");
  printf(" screen output : signal to noise ratio with respect to 'u' (SNR)\n");
  printf(" screen output : peak signal to noise ratio with respect to 'u' (PSNR)\n");
  printf(" screen output : mean square error between 'u' and 'v' (MSE)\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *fname_in1=NULL,*fname_in2=NULL,*p_value=NULL;
  double *in1,*in2,snr,psnr,mse,mrd,peakval;
  int err,k,type,nx,ny,nx2,ny2;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-p") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -p requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	p_value = argv[k+1]; err = sscanf(p_value,"%lf",&peakval); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -d\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_in1) fname_in1 = argv[k];
    else if (NULL == fname_in2) fname_in2 = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == fname_in1) || (NULL == fname_in2)) {
    display_usage();
    if (NULL == fname_in1) printf("Error: input 'in1' is missing\n\n");
    else printf("Error: input 'in2' is missing\n\n");
    return EXIT_FAILURE;
  }

  /*************************/
  /* load input TIFF image */
  /*************************/
  if(EXIT_FAILURE == load_monopage_tiff_into_double(&in1,&nx,&ny,fname_in1,(char)0)){
    printf("Error: failed to open input image '%s'\n",fname_in1);
    return EXIT_FAILURE;
  }
  if(EXIT_FAILURE == load_monopage_tiff_into_double(&in2,&nx2,&ny2,fname_in2,(char)0)){
    printf("Error: failed to open input image '%s'\n",fname_in2);
    free(in1);
    return EXIT_FAILURE;
  }
  if(nx!=nx2 || ny!=ny2) {
    display_usage();
    printf("Error: input images 'in1' and 'in2' should have the same dimensions\n\n");
    free(in1); free(in2);
    return EXIT_FAILURE;
  }

  /**********************/
  /* Core of the module */
  /**********************/
  metrics(&psnr,&snr,&mse,in1,in2,nx,ny,(p_value)?&peakval:NULL);
  printf("MSE = %.4e\n",mse);
  printf("SNR (dB) = %02.4g\n",snr);
  printf("PSNR (dB) = %02.4g\n",psnr);

  /***************/
  /* free memory */
  /***************/
  free(in1); free(in2);

  return EXIT_SUCCESS;
}
