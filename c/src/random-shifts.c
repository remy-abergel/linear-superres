#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* usage displayer */
static void display_usage()
{
  printf("\nGenerate a random sequence of 2D displacements (uniform distribution).\n\n");
  printf("Usage: random-shifts [-m m] [-M M] size out\n\n");
  printf(" -m m      : (default -0.5) minimal displacement\n");
  printf(" -M M      : (default  0.5) maximal displacement\n");
  printf(" -s seed   : use a specified seed for the random number generator\n");
  printf(" size      : number of 2D displacements to generate\n");
  printf(" out       : output sequence of displacements (ASCII format)\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *fname_out=NULL,*m_value=NULL,*M_value=NULL,*seed_value=NULL,*size_value=NULL;
  double dx,dy,m,M;
  long int seed;
  int size,k,err;
  FILE *fout=NULL;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-m") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -m requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	m_value = argv[k+1]; err = sscanf(m_value,"%lf",&m); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -m\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-M") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -M requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	M_value = argv[k+1]; err = sscanf(M_value,"%lf",&M); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -M\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-s") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -s requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	seed_value = argv[k+1]; err = sscanf(seed_value,"%ld",&seed); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -s\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == size_value) { size_value = argv[k]; err = sscanf(size_value,"%d",&size); }
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /**********************/
  /* set default values */
  /**********************/
  if(NULL == m_value) m = -.5;
  if(NULL == M_value) M = .5;

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == size_value) || (NULL == fname_out)) {
    display_usage();
    if (NULL == size_value) printf("Error: input 'size' is missing\n\n");
    else printf("Error: output 'out' is missing\n\n");
    return EXIT_FAILURE;
  }
  if(m>M) printf("WARNING: you selected m > M (m=%g, M=%g).\n",m,M);

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/

  // initialize the random number generator
  srand48((NULL!=seed_value) ? seed : (long int) time (NULL) + (long int) getpid());

  // fill output file
  if(NULL == (fout = fopen(fname_out,"w"))) {
    printf("Error: could not open file '%s'\n",fname_out);
    return EXIT_FAILURE;
  }
  for(k=0;k<size;k++) {
    dx = m + (M-m)*drand48();
    dy = m + (M-m)*drand48();
    fprintf(fout,"%.17e %.17e\n",dx,dy);
  }
  fclose(fout);
  fout = NULL;

  return EXIT_SUCCESS;
}
