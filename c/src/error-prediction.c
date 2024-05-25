#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int save_monopage_striped_tiff_from_double(char*,double*,int,int,char,int,char*); // see source file 'tiffwrite.c'
extern int getasciitranslations(double**,double**,int*,char*,char); // see source file 'ascii.c'
extern int error_prediction(double*,double*,double*,double*,double*,double*,int,int,int,int,int,double,double,char,char); // see source file 'error-prediction_kernel.c'

/* usage displayer */
static void display_usage()
{
  printf("\nPrediction of the super-resolution reconstruction quality.\n\n");
  printf("Usage: error-prediction [-ftype type] [-M M] [-zx zx] [-N N] [-zy zy] [-A A] [-s sigma] [-e eps] [-noshift] [-v] m n T\n\n");
  printf(" -ftype type   : (default float) datatype of the output TIFF images, possible\n");
  printf("                 choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf(" -M M          : (default 2 * width of in) set width of the ouptut image...\n");
  printf(" -zx zx        : ... or set X-axis super-resolution factor (double >= 1.)\n");
  printf(" -N N          : (default 2 * height of in) set height of the output image...\n");
  printf(" -zy zy        : ... or set Y-axis super-resolution factor (double >= 1.)\n");
  printf(" -A A          : output error amplification map (Fourier domain)\n");
  printf(" -s sigma      : noise level (standard deviation) of the sequence to process,\n");
  printf("                 this input is mandatory to display the predicted MSE and PSNR\n");
  printf(" -p peakval    : (default 255.) set the peak-value for the PSNR\n");
  printf(" -e eps        : (default 1e-9) threshold for the singular values in the\n");
  printf("                 pseudo-inversion routine\n");
  printf(" -noshift      : do not shift output A of half the image size (put the (0,0)\n");
  printf("                 frequency at the top-left corner of A instead of its center)\n");
  printf(" -v            : enable verbose mode\n"); 
  printf(" m             : width of the low-resolution image domain\n"); 
  printf(" n             : height of the low-resolution image domain\n"); 
  printf(" T             : input sequence of displacements (ASCII format, two-columns)\n");
  printf(" screen output : (with option -s only) predicted peak signal-to-noise ratio (PSNR)\n");
  printf(" screen output : (with option -s only) predicted mean square error (MSE)\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_T=NULL,*fname_A=NULL,*M_value=NULL,*N_value=NULL,*zx_value=NULL,*zy_value=NULL,*m_value=NULL,*n_value=NULL,*sigma_value=NULL,*peak_value=NULL,*eps_value=NULL,vflag=0,noshiftflag=0;
  double *dx=NULL,*dy=NULL,*A=NULL,*Asmall=NULL,zx,zy,zx_tmp,zy_tmp,sigma,mse,psnr,peakval,eps;
  int err,k,type,nx,ny,ntranslations,Nx,Ny,gx,gy,u,v,r,x,y,adr,Nxsmall,Nysmall,nxsmall,nysmall;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-v") == 0) vflag = 1;
    else if (strcmp(argv[k],"-noshift") == 0) noshiftflag = 1;
    else if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-A") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -A requires an argument\n\n"); return EXIT_FAILURE; }
      else { fname_A = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-s") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -s requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	sigma_value = argv[k+1]; err = sscanf(sigma_value,"%lf",&sigma); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -s\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-e") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -e requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	eps_value = argv[k+1]; err = sscanf(eps_value,"%lf",&eps); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -e\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-p") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -p requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	peak_value = argv[k+1]; err = sscanf(peak_value,"%lf",&peakval); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -p\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-M") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -M requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	M_value = argv[k+1]; err = sscanf(M_value,"%d",&Nx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -M\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-zx") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -zx requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	zx_value = argv[k+1]; err = sscanf(zx_value,"%lf",&zx); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -zx\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-N") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -N requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	N_value = argv[k+1]; err = sscanf(N_value,"%d",&Ny); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -N\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-zy") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -zy requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	zy_value = argv[k+1]; err = sscanf(zy_value,"%lf",&zy); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -zy\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == m_value) {
      m_value = argv[k]; err = sscanf(m_value,"%d",&nx); 
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input argument 'm'\n\n"); return EXIT_FAILURE; }
    }
    else if (NULL == n_value) {
      n_value = argv[k]; err = sscanf(n_value,"%d",&ny); 
      if(err != 1) { display_usage(); printf("Error: could not retrieve properly input argument 'n'\n\n"); return EXIT_FAILURE; }
    }
    else if (NULL == fname_T) fname_T = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /********************************************/
  /* check consistency and set default values */
  /********************************************/
  if(NULL != M_value && NULL != zx_value) { display_usage(); printf("Error: options -M and -zx cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != N_value && NULL != zy_value) { display_usage(); printf("Error: options -N and -zy cannot be used simultaneously\n\n"); return EXIT_FAILURE; }
  if(NULL != zx_value && zx < 1.) { display_usage(); printf("Error: input 'zx' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(NULL != zy_value && zy < 1.) { display_usage(); printf("Error: input 'zy' must be greater than or equal to one.\n\n"); return EXIT_FAILURE; }
  if(NULL == fname_A && NULL == sigma_value) { display_usage(); printf("Error: at least one option must be used among -A and -s.\n\n"); return EXIT_FAILURE; }

  if ((NULL == m_value) || (NULL == n_value) || (NULL == fname_T)) {
    display_usage();
    if (NULL == m_value) printf("Error: input 'm' is missing\n\n");
    else if (NULL == n_value) printf("Error: input 'n' is missing\n\n");
    else printf("Error: input 'T' is missing\n\n");
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

  if(NULL == fname_A && datatype) { printf("WARNING: option -ftype is useless without option -A\n"); };
  if(NULL == fname_A && noshiftflag) { printf("WARNING: option -noshift is useless without option -A\n"); };  
  if(NULL == sigma_value && peak_value) { printf("WARNING: option -p is useless without option -s\n"); };
  if(NULL == peak_value) peakval = 255.; 
  if(NULL == eps_value) eps = 1e-9; 

  /******************************************************************/
  /* retrieve dimensions of the high-resolution image domain (M, N) */
  /******************************************************************/
  if(M_value && Nx < nx) {
    display_usage();
    printf("Error: input 'M' must be greater than or equal to 'm' (M >= %d).\n\n",nx);
    return EXIT_FAILURE;
  }
  if(N_value && Ny < ny) {
    display_usage();
    printf("Error: input 'N' must be greater than or equal to 'n' (N >= %d).\n\n",ny);
    return EXIT_FAILURE;
  }

  if(NULL == zx_value && NULL == M_value) { zx = 2.; Nx = 2*nx; }
  else if(NULL == zx_value) zx = (double)Nx/(double)nx;
  else {
    Nx = (int)round(nx*zx);
    zx_tmp = zx;
    zx = (double)Nx/(double)nx;
    if(zx_tmp != zx) printf("WARNING: zx has been changed into %.17g (instead of zx=%g) to allow integer width for the high-resolution image domain\n",zx,zx_tmp);
  }

  if(NULL == zy_value && NULL == N_value) { zy = 2.; Ny = 2*ny; }
  else if(NULL == zy_value) zy = (double)Ny/(double)ny;
  else {
    Ny = (int)round(ny*zy);
    zy_tmp = zy;
    zy = (double)Ny/(double)ny;
    if(zy_tmp != zy) printf("WARNING: zy has been changed into %.17g (instead of zy=%g) to allow integer height for the high-resolution image domain\n",zy,zy_tmp);
  }

  if(vflag) {
    printf("\nRetrieve image domain dimensions\n");
    printf("================================\n");
    printf("width of the low-resolution domain: m = %d\n",nx);
    printf("height of the low-resolution domain: n = %d\n",ny);
    printf("width of the high-resolution domain: M = %d\n",Nx);
    printf("height of the high-resolution domain: N = %d\n",Ny);
    printf("horizontal super-resolution factor: zx = M/m = %g\n",zx);
    printf("vertical super-resolution factor: zy = N/n = %g\n",zy);
    printf("\n");
  }

  /**************************************/
  /* load the sequence of displacements */
  /**************************************/
  getasciitranslations(&dx,&dy,&ntranslations,fname_T,vflag);

  /*********************/
  /* memory allocation */
  /*********************/
  if(NULL != fname_A) { ASSERT_ALLOC(A = (double*) calloc (Nx*Ny,sizeof(double))); }
  
  /**********************/
  /* CORE OF THE MODULE */
  /**********************/

  /* compute gx = greatest common divisor of Nx and nx */
  u = Nx; v = nx; 
  while (v != 0) {
    r = u % v;
    u = v;
    v = r;
  }
  gx = u;
  
  /* compute gy = greatest common divisor of Ny and ny */
  u = Ny; v = ny; 
  while (v != 0) {
    r = u % v;
    u = v;
    v = r;
  }
  gy = u;

  /* compute domains HR & LR domain withsmallest dimensions (Nxsmall,Nysmall) & (nxsmall,nysmall) such as Nxsmall/nxsmall = Nx/nx & Nysmall/nysmall = Ny/nu */
  nxsmall = nx/gx;
  nysmall = ny/gy;
  Nxsmall = Nx/gx;
  Nysmall = Ny/gy; 

  /* memory allocation */
  ASSERT_ALLOC(Asmall = (double*) calloc (Nxsmall*Nysmall,sizeof(double)));
  
  /* perform error prediction (shift (0,0) frequency at the center of Asmall) */
  error_prediction(Asmall,&psnr,&mse,dx,dy,(sigma_value)?&sigma:NULL,nxsmall,nysmall,Nxsmall,Nysmall,ntranslations,peakval,eps,vflag,0);

  /* (optional) printf predicted PSNR and MSE */
  if(sigma_value) {
    if(vflag) printf("Predicted reconstruction quality\n================================\n");
    printf("predicted MSE = %.4e\n",mse);   
    printf("predicted PSNR (dB) = %02.4g\n",psnr);   
  }

  /* (optional) zoom Asmall by factor (gx,gy) using the nearest neighbor interpolation (+ unshift frequency domain if needed) */
  if(A) {
    if(noshiftflag) { // zoom + unshift 
      for(y=0;y<Ny;y++)
	for(x=0;x<Nx;x++)
	  A[((Ny+(y-Ny/2))%Ny)*Nx + (Nx+(x-Nx/2))%Nx] = Asmall[(y/gy)*Nxsmall+(x/gx)]; 
    }
    else { // zoom only
      for(adr=y=0;y<Ny;y++)
	for(x=0;x<Nx;x++,adr++)
	  A[adr] = Asmall[(y/gy)*Nxsmall+(x/gx)]; 
    }
  }

  /*********************************************************/
  /* if needed, save output (high-resolution) TIFF image A */
  /*********************************************************/
  if(fname_A && EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_A,A,Nx,Ny,0,type,"w")) {
    free(A); free(Asmall); free(dx); free(dy);
    return EXIT_FAILURE;
  }
  
  /***************/
  /* free memory */
  /***************/
  free(Asmall);  
  free(dx);
  free(dy);
  if(fname_A) free(A); 

  return EXIT_SUCCESS;
}
