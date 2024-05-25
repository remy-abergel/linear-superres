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
  printf("\nThreshold/normalize the pixel's gray-levels of a TIFF image.\n\n");
  printf("Usage: tiffthre [-ftype type] [-n] [-N] [-m m] [-M M] [-p p] [-q q] [-d d] [-a] [-l] in out\n\n");
  printf(" -ftype type : (default float) datatype of the output TIFF images, possible\n");
  printf("               choices are {uint8,int8,uint16,int16,uint32,int32,float}\n");
  printf("  -n          : to normalize pixel values from actual to [min,max]\n");
  printf("  -N          : to normalize pixel values from actual or [min,max] to [0,255]\n");
  printf("  -m min      : minimal pixel value\n");
  printf("  -M max      : maximal pixel value\n");
  printf("  -p p        : (value in [0,100]) to discard a ratio of p percent pixels with minimal values\n");
  printf("  -q q        : (value in [0,100]) to discard a ratio of q percent pixels with maximal values\n");
  printf("  -d d        : (value in [0,100]) to discard a ratio of d percent pixels with extremal values\n");
  printf("  -a          : to prevent thresholding (affine normalization only)\n");
  printf("  -l          : force linear normalization (preserve 0)\n");
  printf(" in          : input TIFF (monopage or multipage) image\n");
  printf(" out         : output thresholded/renormalized TIFF image\n\n");
}

/* comparison function for qsort */
int compare(const void *a,const void *b)
{
  double v;
  v = *((double *)a)-*((double *)b);
  return(v>0.?1:(v<0.?-1:0));
}

/* command line interface */
int main(int argc, char **argv)
{
  char *datatype=NULL,*fname_in=NULL,*fname_out=NULL,*m_value=NULL,*M_value=NULL,*p_value=NULL,*q_value=NULL,*d_value=NULL,nflag=0,Nflag=0,aflag=0,lflag=0;
  double **in,**out;
  int err,id,type,nx,ny,nimages,i,k,bestk,l,n,adr;
  double min,max,p,q,d,m0,m1,a,b,w,bestw,*v;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-n") == 0) nflag = 1;
    else if (strcmp(argv[k],"-N") == 0) Nflag = 1;
    else if (strcmp(argv[k],"-a") == 0) aflag = 1;
    else if (strcmp(argv[k],"-l") == 0) lflag = 1;
    else if (strcmp(argv[k],"-ftype") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -ftype requires an argument\n\n"); return EXIT_FAILURE; }
      else { datatype = argv[k+1]; k++; }
    }
    else if (strcmp(argv[k],"-m") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -m requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	m_value = argv[k+1]; err = sscanf(m_value,"%lf",&min); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -m\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-M") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -M requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	M_value = argv[k+1]; err = sscanf(M_value,"%lf",&max); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -M\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-p") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -p requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	p_value = argv[k+1]; err = sscanf(p_value,"%lf",&p); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -p\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-q") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -q requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	q_value = argv[k+1]; err = sscanf(q_value,"%lf",&q); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -q\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-d") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -d requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	d_value = argv[k+1]; err = sscanf(d_value,"%lf",&d); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -d\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
    }

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

  if ((!nflag && !Nflag) && (p_value || q_value || d_value || aflag)) {
    display_usage(); printf("Error: -p -q -d -a options useless for simple thresholding\n\n"); return EXIT_FAILURE;
  }
  if ((!nflag && !Nflag) && (!m_value && !M_value)) {
    display_usage(); printf("Error: Please specify at least min (option -m) or max (option -M) threshold value\n\n"); return EXIT_FAILURE;
  }
  if (nflag && Nflag) {
    display_usage(); printf("-n and -N options cannot be used together\n\n"); return EXIT_FAILURE;
  }
  if (aflag && !lflag && !p_value && !q_value && !d_value) {
    display_usage(); printf("-a option is useless without -p -q -d or -l option\n\n"); return EXIT_FAILURE;
  }
  if(d_value && (p_value || q_value)) {
    display_usage(); printf("You cannot combine -d with -p or -q option\n\n"); return EXIT_FAILURE;
  }
  if(p_value && (p < 0 || p > 100)) {
    display_usage(); printf("you must set a value in [0,100] for option -p\n\n"); return EXIT_FAILURE;
  }
  if(q_value && (q < 0 || q > 100)) {
    display_usage(); printf("you must set a value in [0,100] for option -q\n\n"); return EXIT_FAILURE;
  }
  if(d_value && (d < 0 || d > 100)) {
    display_usage(); printf("you must set a value in [0,100] for option -d\n\n"); return EXIT_FAILURE;
  }
  if(nflag && (!m_value || !M_value)) {
    display_usage(); printf("Please specify min (option -m) and max (option -M) values when using option -n\n\n"); return EXIT_FAILURE;
  }

  /*************************/
  /* load input TIFF image */
  /*************************/
  if(EXIT_FAILURE == load_tiff_into_double(&in,&nx,&ny,&nimages,fname_in,(char)0)){
    printf("Error: failed to open input image '%s'\n",fname_in);
    return EXIT_FAILURE;
  }
  n = nimages*nx*ny;

  /*********************/
  /* memory allocation */
  /*********************/
  ASSERT_ALLOC(out = (double**) malloc (nimages*sizeof(double*)));
  for(k=0;k<nimages;k++) {
    ASSERT_ALLOC(out[k] = (double*) malloc (nx*ny*sizeof(double)));
  }

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/
  if (!nflag && !Nflag) { // simple thresholding

    for(id=0;id<nimages;id++) {
      for (adr=0;adr<nx*ny;adr++) {
	w = in[id][adr];
	if (m_value && w<min) w=min;
	if (M_value && w>max) w=max;
	out[id][adr] = w;
      }
    }

  } else { // affine normalization and maybe thresholding

    if (!p_value && !q_value && !d_value) {

      // compute min (m0) and max (m1)
      m0 = m1 = in[0][0];
      for(id=0;id<nimages;id++) {
	for (adr=0;adr<nx*ny;adr++) {
	  w = in[id][adr];
	  if (w<m0) m0=w;
	  if (w>m1) m1=w;
	}
      }

    }
    else { // discard some values

      // sort values
      ASSERT_ALLOC(v = (double*)malloc(n*sizeof(double)));
      for(k=id=0;id<nimages;id++) {
	for (adr=0;adr<nx*ny;adr++,k++) {
	  v[k] = in[id][adr];
	}
      }
      qsort((void *)v,n,sizeof(double),&compare);
      m0 = v[0]; m1 = v[n-1];

      // compute (pseudo) min (m0) and max (m1)
      if (p_value) {
	k = (int)(0.01*p*(double)n);
	m0 = v[k<0?0:(k>=n?n-1:k)];
      }
      if (q_value) {
	k = n-1-(int)(0.01*q*(double)n);
	m1 = v[k<0?0:(k>=n?n-1:k)];
      }
      if (d_value) {
	l = (int)(0.01*d*(double)n);
        if (l<0) l=0;
	if (l>=n) l=n-1;
	for (k=0;k<=l;k++) {
	  w = v[k+n-l-1]-v[k];
	  if (k==0 || w<bestw) {
	    bestw = w;
	    bestk = k;
	  }
	}
	m0 = v[bestk];
	m1 = v[bestk+n-l-1];
      }

      free(v);

    }

    // compute affine transform u -> a*u+b
    if (nflag) {
      if (lflag) {
	a = max/(m1==0.?1.:m1);
	b = 0.;
	m0 = min; m1 = max; // for thresholding
      }
      else {
	a = (max-min)/(m1-m0==0.?1.:m1-m0);
	b = min-a*m0;
	m0 = min; m1 = max; // for thresholding
      }
    }
    else {
      if (m_value) m0=min;
      if (M_value) m1=max;
      if (lflag) m0=0.;
      a = 255./(m1-m0==0.?1.:m1-m0);
      b = -a*m0;
      m0 = 0.; m1 = 255.; // for thresholding
    }

    // apply normalization (and maybe thresholding)
    for(id=0;id<nimages;id++) {
      for (adr=0;adr<nx*ny;adr++) {
	w = in[id][adr]*a+b;
	if (!aflag) {
	  if (w<m0) w=m0;
	  if (w>m1) w=m1;
	}
	out[id][adr] = w;
      }
    }

  }

  /**************************/
  /* Save output TIFF image */
  /**************************/
  for(k=0;k<nimages;k++) {
    if(EXIT_FAILURE == save_monopage_striped_tiff_from_double(fname_out,out[k],nx,ny,0,type,((k==0)?"w":"a"))) {
      printf("Error: failed to save output image '%s'\n",fname_out);
      for(adr=0;adr<nimages;adr++) free(out[adr]);
      free(out);
      return EXIT_FAILURE;
    }
  }

  /***************/
  /* free memory */
  /***************/
  for(adr=0;adr<nimages;adr++) free(out[adr]);
  free(out);

  return EXIT_SUCCESS;
}
