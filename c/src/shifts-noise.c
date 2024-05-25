#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* external modules */
extern int getasciitranslations(double**,double**,int*,char*,char); // see source file 'ascii.c'

/* internal structure */
typedef struct valueandindex {
  double value;   // value
  int id;         // index
} *Valueandindex;

/* comparison function for qsort (sort the structures in ascending order of their 'value' attribute) */
static int compare_valueandindex(const void *a,const void *b)
{
  struct valueandindex const *pa = a;
  struct valueandindex const *pb = b;
  double diff = pa->value - pb->value;
  return ((diff < 0) ? -1 : ((diff > 0) ? 1 : 0));
}

/* usage displayer */
static void display_usage()
{
  printf("\nAdd two kind of Gaussian noise to a sequence of 2D displacements.\n\n");
  printf("Usage: shifts-noise [-G std1] [-g std2] [-n n] [-s seed] [-f] in out\n\n");
  printf(" -G std1   : (default 0.3) standard deviation of X1 and Y1\n");
  printf(" -g std2   : (default 0.01) standard deviation of X2 and Y2\n");
  printf(" -n n      : (default 0) number of 2D displacements perturbated by (X1,Y1),\n");
  printf("             the others being perturbated by (X2,Y2)\n");
  printf(" -s seed   : use a specified seed for the random number generator\n");
  printf(" -f        : force the (X1,Y1) perturbations to affect the n first 2D\n");
  printf("             displacements of the sequence (otherwise, those perturbations\n");
  printf("             are randomly affected to n entries of the displacement sequence)\n");
  printf(" in        : output sequence of displacements (ASCII format)\n");
  printf(" out       : output sequence of displacements (ASCII format)\n\n");
}

/* command line interface */
int main(int argc, char **argv)
{
  char *fname_in=NULL,*fname_out=NULL,*g_value=NULL,*G_value=NULL,*n_value=NULL,*seed_value=NULL,fflag=0;
  double *dx,*dy,std1,std2,std,a,b,tx,ty;
  long int seed;
  int n,size,k,err;
  struct valueandindex *id;
  FILE *fout=NULL;

  /*************************/
  /* parse input arguments */
  /*************************/
  for(k=1; k<argc; k++) {
    if (strcmp(argv[k],"-f") == 0) fflag = 1;
    else if (strcmp(argv[k],"-G") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -G requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	G_value = argv[k+1]; err = sscanf(G_value,"%lf",&std1); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -G\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-g") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -g requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	g_value = argv[k+1]; err = sscanf(g_value,"%lf",&std2); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -g\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-n") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -n requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	n_value = argv[k+1]; err = sscanf(n_value,"%d",&n); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -n\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (strcmp(argv[k],"-s") == 0) {
      if(k==argc-1) { display_usage(); printf("Error: option -s requires an argument\n\n"); return EXIT_FAILURE; }
      else {
	seed_value = argv[k+1]; err = sscanf(seed_value,"%ld",&seed); k++;
	if(err != 1) { display_usage(); printf("Error: could not retrieve properly the argument of option -s\n\n"); return EXIT_FAILURE; }
      }
    }
    else if (NULL == fname_in) fname_in = argv[k];
    else if (NULL == fname_out) fname_out = argv[k];
    else { display_usage(); printf("Error: could not retrieve properly input arguments.\n\n"); return EXIT_FAILURE; }
  }

  /**********************/
  /* set default values */
  /**********************/
  if(NULL == G_value) std1 = .3;
  if(NULL == g_value) std2 = .01;
  if(NULL == n_value) n = 0;

  /*********************/
  /* check consistency */
  /*********************/
  if ((NULL == fname_in) || (NULL == fname_out)) {
    display_usage();
    if (NULL == fname_in) printf("Error: input 'in' is missing\n\n");
    else printf("Error: output 'out' is missing\n\n");
    return EXIT_FAILURE;
  }
  if(std1 < 0) { printf("Error: input 'std1' must be nonnegative\n\n"); return EXIT_FAILURE; }
  if(std2 < 0) { printf("Error: input 'std2' must be nonnegative\n\n"); return EXIT_FAILURE; }

  /**************************************/
  /* load the sequence of displacements */
  /**************************************/
  getasciitranslations(&dx,&dy,&size,fname_in,(char)0);

  /*********************/
  /* memory allocation */
  /*********************/
  ASSERT_ALLOC(id = (struct valueandindex *) malloc (size*sizeof(struct valueandindex)));

  /**********************/
  /* CORE OF THE MODULE */
  /**********************/

  // initialize the random number generator
  srand48((NULL!=seed_value) ? seed : (long int) time (NULL) + (long int) getpid());

  // associate to each index k=0..size-1 a random value computed using drand48
  for(k=0;k<size;k++) { id[k].id = k; id[k].value = drand48(); };

  // if needed, shuffle indexes (sort the elements of id according to their 'value' attribute to generate a random permutation)
  if(!fflag) qsort((void*)id,size,sizeof(struct valueandindex),compare_valueandindex);

  // fill output file
  if(NULL == (fout = fopen(fname_out,"w"))) {
    printf("Error: could not open file '%s'\n",fname_out);
    free(id); free(dx); free(dy);
    return EXIT_FAILURE;
  }
  for(k=0;k<size;k++) {
    std = (k<n) ? std1 : std2;
    a = drand48();
    b = drand48();
    dx[id[k].id] = dx[id[k].id] + std*sqrt(-2.0*log(a))*cos(2.0*M_PI*b);
    a = drand48();
    b = drand48();
    dy[id[k].id] = dy[id[k].id] + std*sqrt(-2.0*log(a))*cos(2.0*M_PI*b);
  }
  for(k=0;k<size;k++) {
    fprintf(fout,"%.17e %.17e\n",dx[k],dy[k]);
  }
  fclose(fout);
  fout = NULL;

  /***************/
  /* free memory */
  /***************/
  free(id);
  free(dx);
  free(dy);

  return EXIT_SUCCESS;

}
