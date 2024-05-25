#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* internal modules */
int getasciitranslations(double**,double**,int*,char*,char);
int getasciiweights(double**,int*,char*,char);

/* getasciitranslations: retrieve a sequence translation vectors
                         stored in ASCII format.

   Format of the input ASCII file
   ==============================

   Each line of the file must

   + start with character '#': in this case, the line is a comment and
   will be ignored

   OR

   + contain exactly two floating points numbers, <dx> <dy>, separated
   with spaces (at least one space). <dx> denotes the horizontal
   component of the translation vector and <dy> denotes the vertical
   component of the translation vector.

   Important note: empty lines are not allowed (even at the end of the
   file).

   Example of admissible input ASCII file
   ======================================

   --------8<-------- FILE STARTS AT THE NEXT LINE --------8<--------
   # This files contains the two components, dx and dy, of the
   # translation vector associated to each image of the stack.
   #
   # <dx>     <dy>
   0.438744359656398   0.751267059305653
   0.381558457093008   0.255095115459269
   0.765516788149002   0.505957051665142
   0.795199901137063   0.699076722656686
   0.186872604554379   0.890903252535798
   0.489764395788231   0.959291425205444
   0.445586200710899   0.547215529963803
   0.646313010111265   0.138624442828679
   0.709364830858073   0.149294005559057
   0.754686681982361   0.257508254123736
   0.276025076998578   0.840717255983663
   0.679702676853675   0.254282178971531
   0.655098003973841   0.814284826068816
   0.162611735194631   0.243524968724989
   0.118997681558377   0.929263623187228
   0.498364051982143   0.349983765984809
   0.959743958516081   0.196595250431208
   0.340385726666133   0.251083857976031
   0.585267750979777   0.616044676146639
   0.223811939491137   0.473288848902729
   -------8<------- FILE ENDED AT THE PREVIOUS LINE -------8<--------

   Module parameters
   =================

   + filename: path to the file containing the sequence translation
               vectors stored in ASCII format

   + vflag: use vflag = 1 to enable verbose mode, and vflag != 1 to
            disable verbose mode.

   + size: on exit, *size is equal to the number of translations
           retrieved in the ASCII file.

   + dx: on entry, dx is a non allocated double**, on exit, (*dx) is a
         double array containing the (*size) horizontal displacements

   + dy: on entry, dy is a non allocated double**, on exit, (*dy) is a
         double array containing the (*size) vertical displacements

   Module Output 
   =============

   This module should return EXIT_SUCCESS if no error occured, or
   EXIT_FAILURE otherwise.

*/
int getasciitranslations(dx,dy,size,filename,vflag)
     double **dx,**dy;
     int *size;
     char *filename,vflag;
{
  int e,n,c,len=0,lenMax=0,nlines=0,ncomments=0;
  FILE *file;
  char *line;

  /* deal with verbose mode */
  if(vflag) {
    printf("Retrieve translation vectors\n");
    printf("============================\n");
  }

  /* open input file, count number of (nonempty) lines and compute the maximal length of a line */
  if(NULL == (file = fopen(filename,"r"))) {
    printf("Error: could not open file '%s'.\n",filename);
    return EXIT_FAILURE;
  }
  while((c = fgetc(file)) != EOF) {
    len++;
    if(c == '\n') {
      nlines++;
      if (len > lenMax) lenMax = len;
      len=0;
    }
  }
  rewind(file);
  lenMax++;
  if(NULL == (line = (char*) malloc (lenMax*sizeof(char)))) {
    printf("Error: Not enough memory.\n");
    fclose(file);
    return EXIT_FAILURE;
  }

  /* count the number of lines starting with character '#' (comments) */
  for (n=0;n<nlines;n++) {
    if(NULL == fgets(line,lenMax,file)) {
      printf("An error occured while reading the input file '%s'.\n",filename);
      fclose(file); free(line);
      return EXIT_FAILURE;
    }
    if(line[0]=='#') ncomments++;
  }
  rewind(file);

  /* allocate memory for the displacements sequences and retrieve data */
  *size = nlines - ncomments;
  if(NULL == ((*dx) = (double*) malloc ((*size)*sizeof(double)))) {
    printf("Error: Not enough memory.\n");
    fclose(file); free(line);
    return EXIT_FAILURE;
  }
  if(NULL == ((*dy) = (double*) malloc ((*size)*sizeof(double)))) {
    printf("Error: Not enough memory.\n");
    fclose(file); free(line);
    return EXIT_FAILURE;
  }
  for(n=c=0;c<nlines;c++) {
    if(NULL == fgets(line,lenMax,file)) {
      printf("An error occured while reading the input file '%s'.\n",filename);
      fclose(file); free(line); free(*dx); free(*dy);
      return EXIT_FAILURE;
    }
    if(line[0] != '#') {
      if (n < (*size)) e = sscanf(line,"%lf %lf",&((*dx)[n]),&((*dy)[n]));
      if ((e != 2)||(n>=(*size))) {
	printf("Error: could not correctly retrieve the input translations.\n");
	printf("Please check the format of your input translation, each line must\nstart with two double numbers ``<dx> <dy>'' (separated with space) OR\nstart with character '#' (comment). Notice that empty lines are not allowed!\n");
	fclose(file); free(line); free(*dx); free(*dy);
	return EXIT_FAILURE;
      }
      if(vflag) printf("image %d/%d: dx=%5.5g, dy=%5.5g\n",n+1,*size,(*dx)[n],(*dy)[n]);
      n++;
    }
  }
  fclose(file); free(line);

  if(vflag) printf("\n");

  return EXIT_SUCCESS;

}


/* getasciiweights: retrieve a sequence of weights stored in ASCII
                    format.

   Format of the input ASCII file
   ==============================

   Each line of the file must contain exactly one floating point number

   Example of admissible input ASCII file containing 20 weights
   ============================================================

   --------8<-------- FILE STARTS AT THE NEXT LINE --------8<--------
   9.57799194675232161e-04
   3.04322920006941562e-03
   3.15418067262123684e-03
   8.19003708841888620e-04
   2.79957490731773525e-03
   2.87483951519565666e-03
   3.34354897092811482e-03
   1.00490790193313580e-03
   3.17214484686116627e-03
   2.99460875231187659e-03
   2.53073719875250103e-04
   4.63674004372124969e-04
   3.22172167712223793e-03
   2.69724456523737655e-03
   3.24357349703752358e-03
   3.05444131550708245e-03
   3.20313563173085791e-03
   2.89236417433967738e-03
   2.92677171390602576e-03
   3.03672843741324954e-03
   -------8<------- FILE ENDED AT THE PREVIOUS LINE -------8<--------

   Module parameters
   =================

   + filename: path to the file containing the sequence translation
               vectors stored in ASCII format

   + vflag: use vflag = 1 to enable verbose mode, and vflag != 1 to
            disable verbose mode.

   + size: on exit, *size is equal to the number of weights retrieved
           in the ASCII file.

   + weights: on entry, weights is a non allocated double**, on exit,
              (*weights) is a double array containing the (*size) weight values

   Module Output 
   =============

   This module should return EXIT_SUCCESS if no error occured, or EXIT_FAILURE otherwise.

*/
int getasciiweights(weights,size,filename,vflag)
     double **weights;
     int *size;
     char *filename,vflag;
{
  int e,n,c,len=0,lenMax=0,nlines=0,ncomments=0;
  FILE *file;
  char *line;

  /* deal with verbose mode */
  if(vflag) {
    printf("Retrieve weights\n");
    printf("================\n");
  }

  /* open input file, count number of (nonempty) lines and compute the maximal length of a line */
  if(NULL == (file = fopen(filename,"r"))) {
    printf("Error: could not open file '%s'.\n",filename);
    return EXIT_FAILURE;
  }
  while((c = fgetc(file)) != EOF) {
    len++;
    if(c == '\n') {
      nlines++;
      if (len > lenMax) lenMax = len;
      len=0;
    }
  }
  rewind(file);
  lenMax++;
  if(NULL == (line = (char*) malloc (lenMax*sizeof(char)))) {
    printf("Error: Not enough memory.\n");
    fclose(file);
    return EXIT_FAILURE;
  }

  /* allocate memory for the weights sequence and retrieve data */
  *size = nlines;
  if(NULL == ((*weights) = (double*) malloc ((*size)*sizeof(double)))) {
    printf("Error: Not enough memory.\n");
    fclose(file); free(line);
    return EXIT_FAILURE;
  }
  for(n=c=0;c<nlines;c++) {
    if(NULL == fgets(line,lenMax,file)) {
      printf("An error occured while reading the input file '%s'.\n",filename);
      fclose(file); free(line); free(*weights);
      return EXIT_FAILURE;
    }
    if (n < (*size)) e = sscanf(line,"%lf",&((*weights)[n]));
      if ((e != 1)||(n>=(*size))) {
	printf("Error: could not correctly retrieve the input weights.\n");
	printf("Please check the format of your input weights, each line\n");
	printf("must contain one double number. Notice that empty lines\n");
	printf("are not allowed!\n");
	fclose(file); free(line); free(*weights);
	return EXIT_FAILURE;
      }
      if(vflag) printf("image %d/%d: weight=%-.17e\n",n+1,*size,(*weights)[n]);
      n++;
  }
  fclose(file); free(line);

  if(vflag) printf("\n");

  return EXIT_SUCCESS;
}
