#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) { printf("Not enough memory\n"); exit(EXIT_FAILURE); }

/* sharpening : sharpening by frequency amplification with a radial
 *              profile
 *
 * set DFT(out)(a,b) = 1+lambda*(1-exp(-rho(a,b)))
 * for all (a,b) in [-width/2,width/2) x [-height/2,height/2) n Z^2
 * denoting rho(a,b) = sqrt((2*a/width)^2+(2*b/height)^2)
 *
 * Module parameters 
 * =================
 *
 * + width : width of the input and output images
 *
 * + height : height of the input and output images
 *
 * + lambda : amplification parameter (see description above)
 *
 * + in : graylevel values of the input image stored in row-major
 *        order (double array with size width*height)
 *
 * + out : on entry -> a preallocated double array with size
 *         width*height, on exit -> graylevel values of the output
 *         (sharpened) image (row-major order), NB: calling this
 *         module with out = in is possible
 * 
 * Module output
 * =============
 * 
 * This module should return EXIT_SUCCESS
 *
 */
int sharpening(double *out, double *in, double lambda, int width, int height)
{
  fftw_complex *dft_in,*dft_out;
  fftw_plan plan_in,iplan_out;
  double nrm,val,rho;
  int adr,alf,bet;
  
  // memory allocation
  ASSERT_ALLOC(dft_in = (fftw_complex*) fftw_malloc ((1+width/2)*height * sizeof(fftw_complex)));
  ASSERT_ALLOC(dft_out = (fftw_complex*) fftw_malloc ((1+width/2)*height * sizeof(fftw_complex)));
  ASSERT_ALLOC(plan_in = fftw_plan_dft_r2c_2d(height,width,in,dft_in,FFTW_ESTIMATE)); // WARNING: choosing FFTW_MEASURE would destroy 'in'
  ASSERT_ALLOC(iplan_out = fftw_plan_dft_c2r_2d(height,width,dft_out,out,FFTW_ESTIMATE)); // WARNING: choosing FFTW_MEASURE would destroy 'out'
  
  // kernel
  fftw_execute(plan_in); 
  nrm = 1./((double)width*(double)height);
  for(alf=0; alf<=width/2; alf++)
    for(bet=-(height/2); bet<=(height-1)/2; bet++) {
      adr = ((height+bet)%height)*(1+width/2)+alf;
      rho = hypot(2.*(double)alf/(double)width,2.*(double)bet/(double)height); 
      val = 1.+lambda*(1.-exp(-rho));
      dft_out[adr][0] = nrm * dft_in[adr][0] * val;
      dft_out[adr][1] = nrm * dft_in[adr][1] * val;
    }
  fftw_execute(iplan_out);

  // free memory
  fftw_free(dft_in);
  fftw_free(dft_out);
  fftw_destroy_plan(plan_in);
  fftw_destroy_plan(iplan_out);

  return EXIT_SUCCESS;
}
