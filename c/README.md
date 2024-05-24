# Linear Super-Resolution through Translational Motion (C package)

## Installation of the C modules
### Dependencies

This package was designed with **as less dependencies as
possible**. Only the following libraries are required:

+ `libfftw3-dev`   : development packages for the FFTW library
+ `libtiff5-dev` : development packages for the Tag Image File Format
  library (TIFF)
+ `liblapacke-dev` : development packages for the Linear Algebra Package
  (C language API for LAPACK)

Under a Debian GNU/Linux distribution, one can easily get those
libraries using the following apt-get command:

```bash
sudo apt-get install libfftw3-dev libtiff5-dev liblapacke-dev
```

The remaining dependencies are the following standard libraries:
`stdio`, `stdlib`, `string`, `math`, `float`, `limits`,
`unistd`, `time`. Those libraries are already installed by default
on most systems.

### Compilation with gcc

All the C modules of this package can be installed at once by
executing the [`SETUP_C_MODULES`](SETUP_C_MODULES) executable file
(bash):

```bash
./SETUP_C_MODULES
```

When the compilation is successful, the following message is displayed
in the standard output:

```
****************************************************************
*                     compile main modules                     *
****************************************************************

  + compilation of module 'simulator': success
  + compilation of module 'stack-apodization': success
  + compilation of module 'remove-blackborders': success
  + compilation of module 'leastsquares-superres': success
  + compilation of module 'error-prediction': success
  + compilation of module 'irls': success
  + compilation of module 'luckyimaging': success
  + compilation of module 'sharpening': success

****************************************************************
*   compile secondary modules (tools for image manipulation)   *
****************************************************************

  + compilation of module 'register-stack': success
  + compilation of module 'random-shifts': success
  + compilation of module 'shifts-noise': success
  + compilation of module 'tiffmse': success
  + compilation of module 'tiffzoom': success
  + compilation of module 'tiffaddnoise': success
  + compilation of module 'tiffdft': success
  + compilation of module 'tiffcopy': success
  + compilation of module 'tiffconst': success
  + compilation of module 'tiffsqrt': success
  + compilation of module 'tiffaxpb': success
  + compilation of module 'tiffprintasc': success
  + compilation of module 'tiffreadasc': success
  + compilation of module 'tiffsize': success
  + compilation of module 'tiffmerge': success
  + compilation of module 'tiffextract': success
  + compilation of module 'tiffop': success
  + compilation of module 'tiffthre': success

``` 

Notice that you can also manually compile each module using gcc (the
compilation commands can be adapted to other C-compilers). To that
aim, place yourself into the [`src`](src) directory, and run the
following commands:

```bash
gcc -w -O3 simulator.c operators_kernel.c tiffreadcomplex.c tiffwritecomplex.c ascii.c -lm -ltiff -lfftw3 -o simulator
gcc -w -O3 stack-apodization.c apodization_kernel.c tiffread.c tiffwrite.c ascii.c -lm -ltiff -o stack-apodization
gcc -w -O3 remove-blackborders.c remove-blackborders_kernel.c tiffread.c tiffwrite.c ascii.c -lm -ltiff -o remove-blackborders
gcc -w -O3 leastsquares-superres.c leastsquares_kernel.c operators_kernel.c blockmatrix_kernel.c tiffreadcomplex.c tiffwritecomplex.c tiffwrite.c ascii.c -llapacke -lm -ltiff -lfftw3 -o leastsquares-superres
gcc -w -O3 error-prediction.c error-prediction_kernel.c ascii.c blockmatrix_kernel.c tiffwrite.c -llapacke -lm -ltiff -o error-prediction
gcc -w -O3 irls.c leastsquares_kernel.c operators_kernel.c blockmatrix_kernel.c tiffreadcomplex.c tiffwritecomplex.c ascii.c -llapacke -lm -ltiff -lfftw3 -o irls
gcc -w -O3 luckyimaging.c leastsquares_kernel.c operators_kernel.c blockmatrix_kernel.c tiffreadcomplex.c tiffwritecomplex.c ascii.c -llapacke -lm -ltiff -lfftw3 -o luckyimaging
gcc -w -O3 sharpening.c sharpening_kernel.c tiffread.c tiffwrite.c -lm -ltiff -lfftw3 -o sharpening
gcc -w -O3 register-stack.c tools_kernel.c tiffread.c tiffwrite.c tiffreadcomplex.c tiffwritecomplex.c ascii.c -lm -ltiff -lfftw3 -o register-stack
gcc -w -O3 random-shifts.c -o random-shifts
gcc -w -O3 shifts-noise.c ascii.c -lm -o shifts-noise
gcc -w -O3 tiffmse.c tools_kernel.c tiffread.c tiffwrite.c -lm -ltiff -lfftw3 -o tiffmse
gcc -w -O3 tiffzoom.c tools_kernel.c tiffreadcomplex.c tiffwritecomplex.c -lm -ltiff -lfftw3 -o tiffzoom
gcc -w -O3 tiffaddnoise.c tiffread.c tiffwrite.c -lm -ltiff -o tiffaddnoise
gcc -w -O3 tiffdft.c tiffread.c tiffwrite.c -lm -ltiff -lfftw3 -o tiffdft
gcc -w -O3 tiffcopy.c tiffread.c tiffwrite.c -lm -ltiff -lfftw3 -o tiffcopy
gcc -w -O3 tiffconst.c tiffread.c tiffwrite.c -ltiff -o tiffconst
gcc -w -O3 tiffsqrt.c tiffread.c tiffwrite.c -ltiff -lm -o tiffsqrt
gcc -w -O3 tiffaxpb.c tiffread.c tiffwrite.c -ltiff -lm -o tiffaxpb
gcc -w -O3 tiffprintasc.c tiffread.c -ltiff -o tiffprintasc
gcc -w -O3 tiffreadasc.c tiffwrite.c -ltiff -o tiffreadasc
gcc -w -O3 tiffsize.c tiffread.c -ltiff -o tiffsize
gcc -w -O3 tiffmerge.c tiffread.c tiffwrite.c -ltiff -o tiffmerge
gcc -w -O3 tiffextract.c tiffread.c tiffwrite.c -lm -ltiff -o tiffextract
gcc -w -O3 tiffop.c tiffread.c tiffwrite.c -lm -ltiff -o tiffop
gcc -w -O3 tiffthre.c tiffread.c tiffwrite.c -lm -ltiff -o tiffthre
```
When the installation is done, you should be able to find the executable files
listed below in the [`src`](src) directory of the package:

+ `simulator`
+ `stack-apodization`
+ `remove-blackborders`
+ `leastsquares-superres`
+ `error-prediction`
+ `irls`
+ `luckyimaging`
+ `sharpening`
+ `register-stack`
+ `random-shifts`
+ `shifts-noise`
+ `tiffmse`
+ `tiffzoom`
+ `tiffaddnoise`
+ `tiffdft`
+ `tiffcopy`
+ `tiffconst`
+ `tiffsqrt`
+ `tiffaxpb`
+ `tiffprintasc`
+ `tiffreadasc`
+ `tiffsize`
+ `tiffmerge`
+ `tiffextract`
+ `tiffop`
+ `tiffthre`

### Modules usage and (short) documentation
   
You can check the usage and see a (short) documentation of each module
by running the corresponding executable file without argument. For
instance, if you place yourself in the [`src`](src), the following
instruction

```bash
./leastsquares-superres
```
   
should result in the following message in your terminal:

```

Super-resolution using the least-squares estimator.

Usage: leastsquares-superres [-ftype type] [-M M] [-zx zx] [-N N] [-zy zy] [-e eps] [-v] in T out

 -ftype type : (default float) datatype of the output TIFF images, possible
               choices are {uint8,int8,uint16,int16,uint32,int32,float}
 -M M        : (default 2 * width of in) set width of the ouptut image...
 -zx zx      : ... or set X-axis super-resolution factor (double >= 1.)
 -N N        : (default 2 * height of in) set height of the output image...
 -zy zy      : ... or set Y-axis super-resolution factor (double >= 1.)
 -e eps      : (default 1e-9) threshold for the singular values in the
               pseudo-inversion routine
 -v          : enable verbose mode
 in          : input sequence of low-resolution images (TIFF-stack)
 T           : input sequence of displacements (ASCII format, two-columns)
 out         : output high-resolution TIFF image

Error: input 'in' is missing

```


### Optional additional libraries for data visualization

In the
[Examples](c#user-content-examples-reproduce-several-experiments-of-our-publication)
Section of this documentation, we shall use
[ImageJ](https://imagej.net/ij/) and
[Gnuplot](http://www.gnuplot.info/) for image displaying and graphing
from the command line interface. Under a Debian GNU/Linux
distribution, one can install ImageJ and Gnuplot using the following
apt-get command:

```bash
sudo apt-get install imagej gnuplot
```
   
Of course, the installation of ImageJ and Gnuplot is not mandatory
(they are only useful for data visualization from the command line
interface).

Now, you may want to jump to [practical
examples](c#examples-reproduce-several-experiments-of-the-companion-research-article)
(and reproduce some experiments presented the companion article) or to
have a closer look to the [source code organization and
content](c#user-content-package-description-and-organization).

## Package description and organization
### Organization

This package is organized as follows:

+ the [`data`](../data) directory contains some testing datasets,
+ the [`src`](src) directory contains all the source files (.c).

The source files can be separated into three categories named as
"kernel" (`KNL`), "Input/Output" (`I/O`), or "Command Line Interface"
(`CLI`):
	
+ `KNL` : those files contain all the core routines of the package,
they are not interactive and do not contain any main function.
	
+ `I/O` : those files provide routines for handling the
input/outputs (in particular, images in TIFF format and files in ASCII
format). Those files do not contain any main function.
	
+ `CLI` : those files correspond to the interactive modules, they
are in charge of the argument parsing operations and rely on the
routines from the `KNL` and `I/O` source files for further
operations.
	
Notice that several modules of this package, dedicated to basic image
manipulations in TIFF-Format, were mostly inspired from the MegaWave
image processing library.

### Description of the kernel source files of the package

We give below a synthetic description of the kernel source files of
the package.

| SOURCE FILE                                                        | TYPE  | CONTENT / PURPOSE                                                                                        | RELATED PARTS OF THE COMPANION ARTICLE                        |
|--------------------------------------------------------------------|-------|----------------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| [`operators_kernel.c`](src/operators_kernel.c)                     | `KNL` | routines related to the A operator and its adjoint                                                       | Section 2.2, Section 3, Section 6                             |
| [`apodization_kernel.c`](src/apodization_kernel.c)                 | `KNL` | routines related to the apodization of a sequence of low-resolution images                               | Section 2.3                                                   |
| [`blockmatrix_kernel.c`](src/blockmatrix_kernel.c)                 | `KNL` | routines dedicated to the computation of the block matrix defined in (37) and to matrix pseudo inversion | Section 3.3, pseudocode Algorithm 1                           |
| [`error-prediction_kernel.c`](src/error-prediction_kernel.c)       | `KNL` | routines dedicated to the prediction of the reconstruction error                                         | Section 4                                                     |
| [`leastsquares_kernel.c`](src/leastsquares_kernel.c)               | `KNL` | routines related to super-resolution using the least-squares, IRLS and lucky-imaging                     | Pseudocode algorithms 1 & 2, Section 6 (IRLS & lucky-imaging) |
| [`sharpening_kernel.c`](src/sharpening_kernel.c)                   | `KNL` | routines related to image sharpening                                                                     | Section 7                                                     |
| [`tools_kernel.c`](src/tools_kernel.c)                             | `KNL` | routines dedicated to several image manipulations (circular convolution, sub-pixel shift, zooming, etc.) | --                                                            |
| [`remove-blackborders_kernel.c`](src/remove-blackborders_kernel.c) | `KNL` | post-treatment tool for removing borders caused by apodization in the processed high-resolution images   | Section 2.3, Section 4 (Figure 6)                             |

More details about the relations between the content of the kernel
source files and the companion article are given below.


| NAME OF THE ROUTINE    | SOURCE FILE                                                  | RELATION WITH THE COMPANION ARTICLE                                                                                                               |
|------------------------|--------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| `apodization`          | [`apodization_kernel.c`](src/apodization_kernel.c)           | computes the apodization filter *gamma* involved in Equation (14) and (15) as well the apodized low-resolution sequence defined in (16)           |
| `direct_operator_dft`  | [`operators_kernel.c`](src/operators_kernel.c)               | computes the **Discrete Fourier Transform (DFT) of Aj(u)** (operator Aj applied to the high-resolution image u) using Equation (27)               |
| `adjoint_operator_dft` | [`operators_kernel.c`](src/operators_kernel.c)               | computes the *DFT of adj_Aj(v)* (adjoint of Aj applied to the low-resolution image v) using Equation (28)                                         |
| `compute_blockmatrix`  | [`leastsquares_kernel.c`](src/leastsquares_kernel.c)         | implements the pseudocode **Algorithm 1**                                                                                                         |
| `leastsquares`         | [`leastsquares_kernel.c`](src/leastsquares_kernel.c)         | implements the pseudocode **Algorithm 2**                                                                                                         |
| `irls`                 | [`leastsquares_kernel.c`](src/leastsquares_kernel.c)         | implements the **IRLS scheme** described in Equation (56)                                                                                         |
| `luckyimaging`         | [`leastsquares_kernel.c`](src/leastsquares_kernel.c)         | implements the **lucky-imaging procedure** described in Section 6                                                                                 |
| `error_prediction`     | [`error-prediction_kernel.c`](src/error-prediction_kernel.c) | compute the error amplification map A defined in Equation (41) and the predicted MSE and PSNR (46) associated to the least-squares reconstruction |
| `sharpening`           | [`sharpening_kernel.c`](src/sharpening_kernel.c)             | apply the sharpening enhancement post-processing to a super-resolved image using (58)                                                             |
	
### Description of the Command line interface source files of the package

| SOURCE FILE                                              | TYPE  | CONTENT / PURPOSE                                                             | RELATED PARTS OF THE COMPANION ARTICLE                                                            |
|----------------------------------------------------------|-------|-------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| [`simulator.c`](src/simulator.c)                         | `CLI` | Compute a stack of (shifted and subsampled) low-resolution images             | corresponds to the operator A defined in Equation (12)                                            |
| [`stack-apodization.c`](src/stack-apodization.c)         | `CLI` | Compute low/high resolution multiplicative apodization filters                | Section 2.3                                                                                       |
| [`leastsquares-superres.c`](src/leastsquares-superres.c) | `CLI` | Super-resolution using the least-squares estimator                            | Section 3                                                                                         |
| [`error-prediction.c`](src/error-prediction.c)           | `CLI` | Prediction of the super-resolution reconstruction quality                     | Section 4                                                                                         |
| [`irls.c`](src/irls.c)                                   | `CLI` | Iteratively Reweighted Least-Squares                                          | Section 6                                                                                         |
| [`luckyimaging.c`](src/luckyimaging.c)                   | `CLI` | Lucky-imaging procedure                                                       | Section 6                                                                                         |
| [`sharpening.c`](src/sharpening.c)                       | `CLI` | Image sharpening using a Wiener filter                                        | Section 7                                                                                         |
| [`remove-blackborders.c`](src/remove-blackborders.c)     | `CLI` | Crop borders caused by apodization in high-resolution images                  | Figure 6                                                                                          |
| [`register-stack.c`](src/register-stack.c)               | `CLI` | Register a stack of TIFF images from a sequence of registration displacements | Figure 11(b) or Figure 12(b) (*shift-and-add* & *shift-and-median*)                               |
| [`random-shifts.c`](src/random-shifts.c)                 | `CLI` | Generate a random sequence of 2D displacements (uniform distribution)         | sections 5.2 to 5.4 (simulation of random shift sequences)                                        |
| [`shifts-noise.c`](src/shifts-noise.c)                   | `CLI` | Add noise to a sequence of 2D displacements                                   | Section 5.4 and Section 6 (for corrupting a sequence of displacements using random perturbations) |
| [`tiffzoom.c`](src/tiffzoom.c)                           | `CLI` | Magnification of a TIFF image using the Shannon interpolation                 | Figure 10 (third row), Figure 14 (second column)                                                  |
| [`tiffaddnoise.c`](src/tiffaddnoise.c)                   | `CLI` | Add a white Gaussian noise to a TIFF image                                    | Section 5 (for adding noise to synthetic sequences)                                               |
| [`tiffdft.c`](src/tiffdft.c)                             | `CLI` | Compute the 2D Fourier Transform (DFT) of a TIFF image                        | Figure 7, Figure 10, Figure 12, Figure 14                                                         |
| [`tiffmse.c`](src/tiffmse.c)                             | `CLI` | Compute the distance between two TIFF images (SNR, MSE and PSNR metrics)      | --                                                                                                |
| [`tiffmerge.c`](src/tiffmerge.c)                         | `CLI` | Merge a set of TIFF images into a multipage TIFF image (TIFF-stack)           | --                                                                                                |
| [`tiffsqrt.c`](src/tiffsqrt.c)                           | `CLI` | Pixelwise square root of a TIFF image                                         | --                                                                                                |
| [`tiffsize.c`](src/tiffsize.c) (*)                       | `CLI` | Display the dimensions (width, height, number of pages) of a TIFF image       | --                                                                                                |
| [`tiffaxpb.c`](src/tiffaxpb.c) (*)                       | `CLI` | Gain/Offset correction to a TIFF image                                        | --                                                                                                |
| [`tiffconst.c`](src/tiffconst.c) (*)                     | `CLI` | Create TIFF image with constant gray level                                    | --                                                                                                |
| [`tiffcopy.c`](src/tiffcopy.c) (*)                       | `CLI` | Copy a TIFF image (useful for TIFF datatype conversion)                       | --                                                                                                |
| [`tiffextract.c`](src/tiffextract.c) (*)                 | `CLI` | Extract a subpart of a TIFF image                                             | --                                                                                                |
| [`tiffop.c`](src/tiffop.c) (*)                           | `CLI` | Perform an elementary operation between two TIFF images                       | --                                                                                                |
| [`tiffprintasc.c`](src/tiffprintasc.c) (*)               | `CLI` | Print the graylevels of a TIFF image in ascii format                          | --                                                                                                |
| [`tiffreadasc.c`](src/tiffreadasc.c) (*)                 | `CLI` | Read a monopage TIFF image in ascii format                                    | --                                                                                                |
| [`tiffthre.c`](src/tiffthre.c) (*)                       | `CLI` | Threshold/normalize the pixel's gray-levels of a TIFF image                   | --                                                                                                |

(*) modules adapted or inspired from the MegaWave image processing library [2].


### Description of the I/O source files of the package

| SOURCE FILE                                    | TYPE  | CONTENT / PURPOSE                                                                                                   |
|------------------------------------------------|-------|---------------------------------------------------------------------------------------------------------------------|
| [`tiffread.c`](src/tiffread.c)                 | `I/O` | read TIFF images in signed or unsigned 8, 16, 32 or float precision store the graylevels in double precision        |
| [`tiffwrite.c`](src/tiffwrite.c)               | `I/O` | write TIFF images in signed or unsigned 8, 16, 32 or float precision from its graylevels stored in double precision |
| [`tiffreadcomplex.c`](src/tiffreadcomplex.c)   | `I/O` | complex variant of tiffread.c (store the graylevels as the real part of a `fftw_complex`)                           |
| [`tiffwritecomplex.c`](src/tiffwritecomplex.c) | `I/O` | complex variant of tiffwrite.c (write as a TIFF image the real part of the `fftw_complex` graylevels)               |
| [`ascii.c`](src/ascii.c)                       | `I/O` | routines dedicated to manipulation of files in ASCII format                                                         |

## Examples (reproduce several experiments of the companion research article)
### Simulate realistic low-resolution sequences

We illustrate a procedure for synthetizing realistic sequences of
low-resolution images from a given high-resolution image and a
sequence of displacements.
	
Place yourself in the [`src`](src) directory of this package and
run the following bash commands:

```bash
# display a high-resolution image using ImageJ
imagej ../data/bridge.tif 

# generate a sequence containing 10 random shifts in [-5,5] x [-5,5]
# using the 'random-shifts' module
./random-shifts -m -5 -M 5 10 /tmp/shifts.txt

# compute a realistic sequence by cropping a non-realistic sequence
# generated by the 'simulator' module (the crop get rid of the
# unrealistic periodic-like boundaries)
./simulator -zx 3 -zy 3 ../data/bridge.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif
./tiffextract -r -i0 0 -i1 9 /tmp/u0_nonrealistic.tif /tmp/u0_realistic.tif 81 36 150 108

# extract the same area from the high-resolution image (remark the
# correspondance with the low-resolution cropping area: 243=zx*81,
# 108=zy*36, 450=zx*150 and 324=zy*108, with zx=3 and zy=3)
./tiffextract -r ../data/bridge.tif /tmp/ref.tif 243 108 450 324     

# display the realistic low-resolution sequence and the corresponding
# high-resolution image using ImageJ
imagej /tmp/u0_realistic.tif /tmp/ref.tif	
```

### Using apodization to avoid boundary artifacts (reproduce Figure 1 of the companion article)
	
The following experiment illustrates the importance of properly
dealing with periodization artifacts induced by the periodicity of the
Shannon interpolate when computing a super-resolved image using the
leastsquare-superres module.

Place yourself in the [`src`](src) directory of this package and
run the following commands:

```bash
# display the realistic low-resolution sequence (TIFF-Stack) using
# ImageJ (see the previous section for details about how such sequence
# can be generated)
imagej ../data/bridge_sequence_fig1.tif

# perform super-resolution over a precomputed realistic sequence of
# low-resolution images: we observe important artefacts in the
# super-resolved image
./leastsquares-superres -ftype uint8 -zx 3 -zy 3 ../data/bridge_sequence_fig1.tif ../data/shifts_fig1.txt /tmp/out1.tif
imagej /tmp/out1.tif 

# apodize the sequence using the 'stack-apodization' module
./stack-apodization -zx 3 -zy 3 -a /tmp/u0_apod.tif ../data/bridge_sequence_fig1.tif ../data/shifts_fig1.txt
imagej /tmp/u0_apod.tif 

# perform super-resolution over the apodized sequence (the super-resolved image
# does not exhibits artifact anymore)
./leastsquares-superres -ftype uint8 -zx 3 -zy 3 /tmp/u0_apod.tif ../data/shifts_fig1.txt /tmp/out2.tif
imagej /tmp/out2.tif
```

### Super-resolution using the least-squares (reproduce Figure 7 of the companion article)

The next experiment intents to perform super-resolution with factor 2
along both dimensions (zx=zy=2) from a sequence containing L = 20
noisy low-resolution images, with no perturbation on the
displacements.

Place yourself in the [`src`](src) directory of this package and
run the following bash commands:

```bash
# compute a realistic sequence containing L=20 low-resolution images
# corrupted by an additive Gaussian noise (stdev = 2)
./random-shifts -m -5 -M 5 20 /tmp/shifts.txt
./simulator -zx 2 -zy 2 ../data/bridge.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif
./tiffextract -r -i0 0 -i1 19 /tmp/u0_nonrealistic.tif /tmp/u0_realistic.tif 121 54 225 162
./tiffaddnoise -g 2 /tmp/u0_realistic.tif /tmp/u0_noisy.tif

# apodize this sequence to avoid boundary artifacts in further
# processings (see previous section for more details)
./stack-apodization -zx 2 -zy 2 -G /tmp/apod_hr.tif -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts.txt
imagej /tmp/u0_apod.tif 

# display the spectra (modulus of the DFT) in logarithmic scale of the
# low-resolution images (remark the aliasing in the frequency domain)
./tiffdft -l /tmp/dft_u0_apod.tif /tmp/u0_apod.tif 
imagej /tmp/dft_u0_apod.tif 

# perform super-resolution, display the result and its spectrum
./leastsquares-superres -zx 2 -zy 2 /tmp/u0_apod.tif /tmp/shifts.txt /tmp/out.tif
./tiffdft -l /tmp/dft_out.tif /tmp/out.tif 
imagej /tmp/out.tif /tmp/dft_out.tif 

# compare to the apodized reference image (compute PSNR & MSE metrics)
# note about the reference extraction: remark the correspondance with
# the low-resolution cropping area: 242=zx*121, 108=zy*54, 450=zx*225
# and 324=zy*162, with zx=zy=2)
./tiffextract -r ../data/bridge.tif /tmp/ref.tif 242 108 450 324
./tiffop -A /tmp/apod_hr.tif -t /tmp/ref.tif /tmp/ref_apod.tif
./tiffmse -p 255 /tmp/ref_apod.tif /tmp/out.tif
```

### Prediction of the reconstruction quality (reproduce Figure 4, 5, 6 of the companion article)

The following experiment illustrates how the reconstruction quality
provided by the least-square estimator can be blindly (i.e. without an
reference image) and efficiently predicted. It also illustrates the
non-stationarity of the reconstruction error in the Fourier domain.
	
**Remark**: in this experiment, to simplify the image synthesis
process, we generate a (non-realistic) sequence using the `simulator`
module and then we apodize it using the `stack-apodization`
module. This provides a realistic apodized sequence, in the sense that
the apodized sequence obtained in this way is visually
indistinguishable from the apodized sequence that we would have
obtained if we had generated a realistic sequence and apodized it
afterwards (because, by construction, the apodization filter vanishes
where the periodic-like artifact occur in the non-realistic simulated
sequence).
	
Place yourself in the [`src`](src) directory of this package and run
the following bash commands:

```bash
# configure experimental parameters (using bash variables)
sigma_noise=2; # noise standard deviation 
sigma_apod=1; # apodization parameter

# compute an apodized low-resolution sequence (L=10, zx=2.5, zy=2.4),
# use noise standard deviation sigma = 2
./tiffextract -r ../data/bridge.tif /tmp/ref_crop.tif 242 108 450 324
./random-shifts -m -5 -M 5 10 /tmp/shifts.txt
./simulator -zx 2.5 -zy 2.4 /tmp/ref_crop.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif 
./tiffaddnoise -g $sigma_noise /tmp/u0_nonrealistic.tif /tmp/u0_noisy.tif
./stack-apodization -s $sigma_apod -zx 2.5 -zy 2.4 -G /tmp/apod_hr.tif -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts.txt
./tiffop -A /tmp/ref_crop.tif -t /tmp/apod_hr.tif /tmp/ref_apod.tif
imagej /tmp/u0_apod.tif 

# use 'tiffsize' module + basic BASH commands to retrieve the
# dimensions of the low-resolution domain
m=`./tiffsize /tmp/u0_apod.tif | grep width | cut -d =  -f 2`
n=`./tiffsize /tmp/u0_apod.tif | grep height | cut -d =  -f 2`

# perform prediction of the least-square super-resolution
# reconstruction quality (MSE & PSNR) from the dimensions of the
# low-resolution image domai (m,n), the super-resolution factors
# (zx,zy), the displacement sequence and the noise level (sigma)
./error-prediction -p 255 -A /tmp/A.tif -s $sigma_noise -zx 2.5 -zy 2.4 $m $n /tmp/shifts.txt

# perform super-resolution reconstruction and compute observed MSE &
# PSNR
./leastsquares-superres -zx 2.5 -zy 2.4  /tmp/u0_apod.tif /tmp/shifts.txt /tmp/uls.tif
./tiffmse -p 255 /tmp/ref_apod.tif /tmp/uls.tif 

# we remark that the observed MSE & PSNR are a bit better than the
# predictions, this is due to apodization which is not taken into
# account in our theoretical study of the error. Indeed, the
# reconstruction error is very small in the areas affected by
# apodization (black borders), as we can check here (absdiff.tif
# = |uls.tif - ref_apod.tif| is almost zero in areas affected by
# apodization)
./tiffop -A /tmp/uls.tif -D /tmp/ref_apod.tif /tmp/absdiff.tif
imagej /tmp/absdiff.tif

# removing the black borders (that do not correspond to useful signal)
# from the apodized high-resolution images, the observed MSE and PSNR
# are even closer to their predicted values
./remove-blackborders -s $sigma_apod -zx 2.5 -zy 2.4 /tmp/shifts.txt /tmp/uls.tif /tmp/uls_crop.tif
./remove-blackborders -s $sigma_apod -zx 2.5 -zy 2.4 /tmp/shifts.txt /tmp/ref_apod.tif /tmp/ref_apod_crop.tif
./tiffmse -p 255 /tmp/uls_crop.tif /tmp/ref_apod_crop.tif

# compute the normalized reconstruction error in the Fourier domain
./tiffop -A /tmp/uls.tif -m /tmp/ref_apod.tif /tmp/diff.tif # diff.tif = uls.tif - ref_apod.tif 
./tiffdft -m /tmp/observed_error_fourier_domain.tif /tmp/diff.tif # compute |DFT(diff.tif)|
M=`./tiffsize /tmp/uls.tif | grep width | cut -d =  -f 2` # retrieve M = width of the high-resolution domain
N=`./tiffsize /tmp/uls.tif | grep height | cut -d =  -f 2`;  # retrieve N = height of the high-resolution domain
COF=`echo "1/sqrt($M*$N)" | bc -l` # compute normalization coefficient COF = 1/sqrt(M*N)
./tiffop -a $COF -t /tmp/observed_error_fourier_domain.tif /tmp/observed_normalized_error_fourier_domain.tif

# compare the observed normalized error in the Fourier domain to its
# expectation (we use tiffmerge module to stack the images in order to
# facilitate their comparison using imagej)
./tiffop -a $sigma_noise -t /tmp/A.tif /tmp/predicted_normalized_error_fourier_domain.tif 
./tiffmerge /tmp/mov.tif /tmp/observed_normalized_error_fourier_domain.tif /tmp/predicted_normalized_error_fourier_domain.tif
imagej /tmp/mov.tif # first image = observed normalized error, second image = expected normalized error (use keyboard left/right to switch)	
```

In the following script, we focus on the reconstruction error in both
Fourier & spatial domains, we reproduce (with larger image domains)
experiments similar to that described in Figure 4 and Figure 5.

```bash
# set size of the low and high resolution domains (use bash
# variables)
m=200; # width of the low-resolution domain
n=100; # height of the low-resolution domain
M=460; # width of the high-resolution domain 
N=220; # height of the high-resolution domain
COF=`echo "1/sqrt($M*$N)" | bc -l` # normalization coefficient COF = 1/sqrt(M*N)

# compute a pure Gaussian noise sequence (L=9 images, noise
# standard deviation = 1), a sequence of L=9 random
# displacements, and the associated reconstruction error (=
# least-squares operator applied to the pure noise sequence)
./tiffconst -L 9 $m $n /tmp/eps_L9.tif; 
./tiffaddnoise -g 1. /tmp/eps_L9.tif /tmp/eps_L9.tif 
./random-shifts -m -5 -M 5 9 /tmp/shifts_L9.txt
./leastsquares-superres -M $M -N $N /tmp/eps_L9.tif /tmp/shifts_L9.txt /tmp/eps_prime_L9.tif; 

# proceed similarly to compute a reconstruction error for L = 14
./tiffconst -L 14 $m $n /tmp/eps_L14.tif; 
./tiffaddnoise -g 1. /tmp/eps_L14.tif /tmp/eps_L14.tif 
./random-shifts -m -5 -M 5 14 /tmp/shifts_L14.txt
./leastsquares-superres -M $M -N $N /tmp/eps_L14.tif /tmp/shifts_L14.txt /tmp/eps_prime_L14.tif; 

# proceed similarly to compute a reconstruction error for L = 20
./tiffconst -L 20 $m $n /tmp/eps_L20.tif; 
./tiffaddnoise -g 1. /tmp/eps_L20.tif /tmp/eps_L20.tif 
./random-shifts -m -5 -M 5 20 /tmp/shifts_L20.txt
./leastsquares-superres -M $M -N $N /tmp/eps_L20.tif /tmp/shifts_L20.txt /tmp/eps_prime_L20.tif; 

# display the (normalized) reconstruction error in the Fourier
# domain and the associated map of error amplification
# coefficients (L = 9): reproduce a similar result as that 
# displayed in the first column of Figure 4 (use gray colormap 
# instead of false colors)
./tiffdft -m /tmp/eps_prime_L9_fourier.tif /tmp/eps_prime_L9.tif
./tiffop -a $COF -t /tmp/eps_prime_L9_fourier.tif /tmp/eps_prime_L9_fourier.tif 
./error-prediction -A /tmp/A_L9.tif -M $M -N $N $m $n /tmp/shifts_L9.txt
imagej /tmp/eps_prime_L9_fourier.tif /tmp/A_L9.tif

# display the (normalized) reconstruction error in the Fourier
# domain and the associated map of error amplification
# coefficients (L = 14): reproduce a similar result as that 
# displayed in the second column of Figure 4 (use gray colormap 
# instead of false colors)
./tiffdft -m /tmp/eps_prime_L14_fourier.tif /tmp/eps_prime_L14.tif
./tiffop -a $COF -t /tmp/eps_prime_L14_fourier.tif /tmp/eps_prime_L14_fourier.tif 
./error-prediction -A /tmp/A_L14.tif -M $M -N $N $m $n /tmp/shifts_L14.txt 
imagej /tmp/eps_prime_L14_fourier.tif /tmp/A_L14.tif

# display the (normalized) reconstruction error in the Fourier
# domain and the associated map of error amplification
# coefficients (L = 20): reproduce a similar result as that 
# displayed in the last column of Figure 4 (use gray colormap 
# instead of false colors)
./tiffdft -m /tmp/eps_prime_L20_fourier.tif /tmp/eps_prime_L20.tif
./tiffop -a $COF -t /tmp/eps_prime_L20_fourier.tif /tmp/eps_prime_L20_fourier.tif 
./error-prediction -A /tmp/A_L20.tif -M $M -N $N $m $n /tmp/shifts_L20.txt 
imagej /tmp/eps_prime_L20_fourier.tif /tmp/A_L20.tif

# display the reconstruction error in the spatial domain 
# (reproduce similar results as those displayed in the first 
# row of Figure 5)
imagej /tmp/eps_prime_L9.tif /tmp/eps_prime_L14.tif /tmp/eps_prime_L20.tif

# compute and display the empirical standard deviation of the
# reconstruction error in the spatial domain (this simulation takes
# several minutes) to reproduce similar results as those displayed 
# in the last row of Figure 5 (use gray colormap instead of false 
# colors)
./tiffconst -g 0 $M $N /tmp/emp_var_L9.tif
./tiffconst -g 0 $M $N /tmp/emp_var_L14.tif
./tiffconst -g 0 $M $N /tmp/emp_var_L20.tif
Nsimu=100; # number of random simulations to perform
for((id=1; id<=$Nsimu; id+=1))
do
	
    # compute a realization of the reconstruction error for L = 9 and
    # update empirical variance estimation 
    ./tiffconst -L 9 $m $n /tmp/eps_L9.tif;
    ./tiffaddnoise -g 1. /tmp/eps_L9.tif /tmp/eps_L9.tif;
    ./leastsquares-superres -M $M -N $N /tmp/eps_L9.tif /tmp/shifts_L9.txt /tmp/eps_prime_L9.tif;
    ./tiffop -A /tmp/eps_prime_L9.tif -t /tmp/eps_prime_L9.tif /tmp/a.tif;
    ./tiffop -A /tmp/emp_var_L9.tif -p /tmp/a.tif /tmp/emp_var_L9.tif;
	
    # compute a realization of the reconstruction error for L = 14 and
    # update empirical variance estimation
    ./tiffconst -L 14 $m $n /tmp/eps_L14.tif;
    ./tiffaddnoise -g 1. /tmp/eps_L14.tif /tmp/eps_L14.tif;
    ./leastsquares-superres -M $M -N $N /tmp/eps_L14.tif /tmp/shifts_L14.txt /tmp/eps_prime_L14.tif;
    ./tiffop -A /tmp/eps_prime_L14.tif -t /tmp/eps_prime_L14.tif /tmp/a.tif;
    ./tiffop -A /tmp/emp_var_L14.tif -p /tmp/a.tif /tmp/emp_var_L14.tif;
	
    # compute a realization of the reconstruction error for L = 20 and
    # update empirical variance estimation 
    ./tiffconst -L 20 $m $n /tmp/eps_L20.tif;
    ./tiffaddnoise -g 1. /tmp/eps_L20.tif /tmp/eps_L20.tif;
    ./leastsquares-superres -M $M -N $N /tmp/eps_L20.tif /tmp/shifts_L20.txt /tmp/eps_prime_L20.tif;
    ./tiffop -A /tmp/eps_prime_L20.tif -t /tmp/eps_prime_L20.tif /tmp/a.tif
    ./tiffop -A /tmp/emp_var_L20.tif -p /tmp/a.tif /tmp/emp_var_L20.tif
	
    echo simulation \#$id/$Nsimu: done
	
done
./tiffop -a `echo "1/$Nsimu" | bc -l` -t /tmp/emp_var_L9.tif /tmp/emp_var_L9.tif
./tiffop -a `echo "1/$Nsimu" | bc -l` -t /tmp/emp_var_L14.tif /tmp/emp_var_L14.tif
./tiffop -a `echo "1/$Nsimu" | bc -l` -t /tmp/emp_var_L20.tif /tmp/emp_var_L20.tif
./tiffsqrt /tmp/emp_var_L9.tif /tmp/emp_std_L9.tif
./tiffsqrt /tmp/emp_var_L14.tif /tmp/emp_std_L14.tif
./tiffsqrt /tmp/emp_var_L20.tif /tmp/emp_std_L20.tif
imagej /tmp/emp_std_L9.tif /tmp/emp_std_L14.tif /tmp/emp_std_L20.tif
```

Now, let us compare the accuracy of the PSNR prediction over
several random experiments (we reproduce an experiment similar to
that described in Figure 6, this experiment takes several minutes).

```bash
# configure experimental parameters (using bash variables)
sigma_noise=2; # noise standard deviation 
sigma_apod=1; # apodization parameter
peakvalue=255; # peak-value used for (observed/predicted) PSNR computation

# prepare reference image
./tiffextract -r ../data/bridge_large.tif /tmp/ref_crop.tif 242 108 450 324
M=`./tiffsize /tmp/ref_crop.tif | grep width | cut -d =  -f 2` # retrieve M = width of the high-resolution domain
N=`./tiffsize /tmp/ref_crop.tif | grep height | cut -d =  -f 2`;  # retrieve N = height of the high-resolution domain

#########################################################################################
# reproduce Figure 6 (a) : compute synthetic low-resolution stack of images, perform    #
# least-squares reconstruction (zx = zy = 2) and evaluate PSNR                          #
#########################################################################################

# set simulation parameters
L_list="4 5 6 8 10 12 14" # values of L to be tested
Nsimu=50; # number of simulations per value of L (in Figure 5 (a) we used Nsimu = 100)

# set size of the low-resolution image domain
m=`echo "$M/2" | bc -l | cut -d . -f 1` # compute m = floor(M/2) = width of the low-resolution image domain
n=`echo "$N/2" | bc -l | cut -d . -f 1` # compute n = floor(N/2) = height of the low-resolution image domain

# prepare gnuplot script for data plotting
echo "set size ratio .5" > /tmp/scatterplot.gnuplot 
echo "set xlabel 'predicted PSNR (dB)'" >> /tmp/scatterplot.gnuplot 
echo "set ylabel 'observed PSNR (dB)'" >> /tmp/scatterplot.gnuplot 
echo "set key bottom box 3" >> /tmp/scatterplot.gnuplot 
echo "set grid" >> /tmp/scatterplot.gnuplot 
echo -n "plot " >> /tmp/scatterplot.gnuplot 

# main loop 
for L in $L_list 
do 
	
    echo "# PREDICTED PSNR (dB) , OBSERVED PSNR (dB)" > /tmp/PSNR_PRED_VS_OBSERVED_L$L.txt
	
    for ((id=1; id<=$Nsimu; id+=1))
    do
		
        # generate apodized low-resolution sequence
        ./random-shifts -m -5 -M 5 $L /tmp/shifts.txt
        ./simulator -m $m -n $n /tmp/ref_crop.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif 
        ./tiffaddnoise -g $sigma_noise /tmp/u0_nonrealistic.tif /tmp/u0_noisy.tif
        ./stack-apodization -s $sigma_apod -M $M -N $N -G /tmp/apod_hr.tif -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts.txt
		
        # compute PSNR prediction
		psnr_pred=`./error-prediction -p $peakvalue -s $sigma_noise -M $M -N $N $m $n /tmp/shifts.txt | grep PSNR | cut -d = -f 2`
		echo -n $psnr_pred >> /tmp/PSNR_PRED_VS_OBSERVED_L$L.txt
        
        # compute reference image associated to the apodized sequence (and
        # remove black borders caused by apodization)
        ./tiffop -A /tmp/ref_crop.tif -t /tmp/apod_hr.tif /tmp/ref_apod.tif
        ./remove-blackborders -s $sigma_apod -m $m -n $n /tmp/shifts.txt /tmp/ref_apod.tif /tmp/ref_noborder.tif; 
        
        # perform least-squares reconstruction & remove black borders caused
        # by apodization
        ./leastsquares-superres -M $M -N $N /tmp/u0_apod.tif /tmp/shifts.txt /tmp/uls.tif
		./remove-blackborders -s $sigma_apod -m $m -n $n /tmp/shifts.txt /tmp/uls.tif /tmp/uls_noborder.tif; 
        
        # compute observed PSNR between the least-squares reconstruction and
		# the apodized reference image (remove black-borders)
		psnr_observed=`./tiffmse -p $peakvalue /tmp/ref_noborder.tif /tmp/uls_noborder.tif | grep PSNR | cut -d = -f 2`
		echo ",$psnr_observed" >> /tmp/PSNR_PRED_VS_OBSERVED_L$L.txt
		
		# display progression 
		echo L = $L simulation \#$id/$Nsimu : predicted PSNR = $psnr_pred dB, observed PSNR = $psnr_observed dB
		
	done
	
	# prepare gnuplot script for display
	echo -n "'/tmp/PSNR_PRED_VS_OBSERVED_L$L.txt' using 1:2 title 'L = $L', " >> /tmp/scatterplot.gnuplot
	
done

# display results using gnuplot 
gnuplot -p /tmp/scatterplot.gnuplot


#########################################################################################
# reproduce Figure 6 (b): compute synthetic low-resolution stack of images, perform     #
# least-squares reconstruction (zx = 2.1, zy = 2.3) and evaluate PSNR)                  #
#########################################################################################

# set simulation parameters
L_list="9 10 11 13 15 17 19" # values of L to be tested 
Nsimu=50; # number of simulations per value of L (in Figure 5 (a) we used Nsimu = 100)

# set size of the low-resolution image domain
m=`echo "$M/2.1" | bc -l | cut -d . -f 1` # compute m = width of the low-resolution image domain
n=`echo "$N/2.29" | bc -l | cut -d . -f 1` # compute n = height of the low-resolution image domain

# prepare gnuplot script for data plotting
echo "set size ratio .5" > /tmp/scatterplot.gnuplot 
echo "set xlabel 'predicted PSNR (dB)'" >> /tmp/scatterplot.gnuplot 
echo "set ylabel 'observed PSNR (dB)'" >> /tmp/scatterplot.gnuplot 
echo "set key bottom box 3" >> /tmp/scatterplot.gnuplot 
echo "set grid" >> /tmp/scatterplot.gnuplot 
echo -n "plot " >> /tmp/scatterplot.gnuplot 

# main loop 
for L in $L_list 
do 
    
    echo "# PREDICTED PSNR (dB) , OBSERVED PSNR (dB)" > /tmp/PSNR_PRED_VS_OBSERVED_L$L.txt
    
    for ((id=1; id<=$Nsimu; id+=1))
    do
        
        # generate apodized low-resolution sequence
        ./random-shifts -m -5 -M 5 $L /tmp/shifts.txt
        ./simulator -m $m -n $n /tmp/ref_crop.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif 
        ./tiffaddnoise -g $sigma_noise /tmp/u0_nonrealistic.tif /tmp/u0_noisy.tif
        ./stack-apodization -s $sigma_apod -M $M -N $N -G /tmp/apod_hr.tif -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts.txt
        
        # compute PSNR prediction
        psnr_pred=`./error-prediction -p $peakvalue -s $sigma_noise -M $M -N $N $m $n /tmp/shifts.txt | grep PSNR | cut -d = -f 2`
        echo -n $psnr_pred >> /tmp/PSNR_PRED_VS_OBSERVED_L$L.txt
        
        # compute reference image associated to the apodized sequence (and
        # remove black borders caused by apodization)
        ./tiffop -A /tmp/ref_crop.tif -t /tmp/apod_hr.tif /tmp/ref_apod.tif
        ./remove-blackborders -s $sigma_apod -m $m -n $n /tmp/shifts.txt /tmp/ref_apod.tif /tmp/ref_noborder.tif; 
        
        # perform least-squares reconstruction & remove black borders caused
        # by apodization
        ./leastsquares-superres -M $M -N $N /tmp/u0_apod.tif /tmp/shifts.txt /tmp/uls.tif
        ./remove-blackborders -s $sigma_apod -m $m -n $n /tmp/shifts.txt /tmp/uls.tif /tmp/uls_noborder.tif; 
        
        # compute observed PSNR between the least-squares reconstruction and
        # the apodized reference image (remove black-borders)
        psnr_observed=`./tiffmse -p $peakvalue /tmp/ref_noborder.tif /tmp/uls_noborder.tif | grep PSNR | cut -d = -f 2`
        echo ",$psnr_observed" >> /tmp/PSNR_PRED_VS_OBSERVED_L$L.txt
        
        # display progression 
        echo L = $L simulation \#$id/$Nsimu : predicted PSNR = $psnr_pred dB, observed PSNR = $psnr_observed dB
        
    done
    
    # prepare gnuplot script for display
    echo -n "'/tmp/PSNR_PRED_VS_OBSERVED_L$L.txt' using 1:2 title 'L = $L', " >> /tmp/scatterplot.gnuplot
    
done

# display results using gnuplot
gnuplot -p /tmp/scatterplot.gnuplot
```

### Influence of the displacement configuration over the quality of the reconstruction (reproduce Figure 9 of the companion article)
   
In the next experiment, we reproduce some results similar to that
displayed in the first row of Figure 9, showing how the displacement
configuration may affect the quality reconstruction (especially when L
is close to zx*zy).

Place yourself in the [`src`](src) directory of this package and run
the following bash commands:

```bash
# compute three sequences containing 4 low-resolution images using the
# same displacement sequences as that used to compute the images
# displayed in the first row of Figure 9
./simulator -zx 2 -zy 2 ../data/bridge.tif ../data/shifts_fig9a.txt /tmp/u0_nonrealistic1.tif
./simulator -zx 2 -zy 2 ../data/bridge.tif ../data/shifts_fig9b.txt /tmp/u0_nonrealistic2.tif
./simulator -zx 2 -zy 2 ../data/bridge.tif ../data/shifts_fig9c.txt /tmp/u0_nonrealistic3.tif
./tiffextract -r -i0 0 -i1 3 /tmp/u0_nonrealistic1.tif /tmp/u0_1.tif 121 54 225 162
./tiffextract -r -i0 0 -i1 3 /tmp/u0_nonrealistic2.tif /tmp/u0_2.tif 121 54 225 162
./tiffextract -r -i0 0 -i1 3 /tmp/u0_nonrealistic3.tif /tmp/u0_3.tif 121 54 225 162
./tiffaddnoise -g 2 /tmp/u0_1.tif /tmp/u0_1.tif
./tiffaddnoise -g 2 /tmp/u0_2.tif /tmp/u0_2.tif
./tiffaddnoise -g 2 /tmp/u0_3.tif /tmp/u0_3.tif
./stack-apodization -zx 2 -zy 2 -G /tmp/apod_hr1.tif -a /tmp/u0_1_apod.tif /tmp/u0_1.tif ../data/shifts_fig9a.txt
./stack-apodization -zx 2 -zy 2 -G /tmp/apod_hr2.tif -a /tmp/u0_2_apod.tif  /tmp/u0_2.tif ../data/shifts_fig9b.txt
./stack-apodization -zx 2 -zy 2 -G /tmp/apod_hr3.tif -a /tmp/u0_3_apod.tif  /tmp/u0_3.tif ../data/shifts_fig9c.txt
m=`./tiffsize /tmp/u0_1_apod.tif | grep width | cut -d =  -f 2` # retrieve m = width of the low-resolution domain
n=`./tiffsize /tmp/u0_1_apod.tif | grep height | cut -d =  -f 2` # retrieve n = height of the low-resolution domain

# compute the corresponding apodized high-resolution reference images
# (one reference image per apodization filter)
./tiffextract -r ../data/bridge.tif /tmp/ref.tif 242 108 450 324
./tiffop -A /tmp/apod_hr1.tif -t /tmp/ref.tif /tmp/ref_apod1.tif
./tiffop -A /tmp/apod_hr2.tif -t /tmp/ref.tif /tmp/ref_apod2.tif
./tiffop -A /tmp/apod_hr3.tif -t /tmp/ref.tif /tmp/ref_apod3.tif

# perform super-resolution reconstruction from each sequence of
# low-resolution images
./leastsquares-superres -ftype uint8 -zx 2 -zy 2 /tmp/u0_1_apod.tif ../data/shifts_fig9a.txt /tmp/out1.tif
./leastsquares-superres -ftype uint8 -zx 2 -zy 2 /tmp/u0_2_apod.tif ../data/shifts_fig9b.txt /tmp/out2.tif
./leastsquares-superres -ftype uint8 -zx 2 -zy 2 /tmp/u0_3_apod.tif ../data/shifts_fig9c.txt /tmp/out3.tif

# perform super-resolution reconstruction from each sequence of
# low-resolution images and compute the predicted PSNR for each
# reconstruction
imagej /tmp/out1.tif /tmp/out2.tif /tmp/out3.tif;
./error-prediction -p 255 -s 2 -zx 2 -zy 2 $m $n ../data/shifts_fig9a.txt
./error-prediction -p 255 -s 2 -zx 2 -zy 2 $m $n ../data/shifts_fig9b.txt
./error-prediction -p 255 -s 2 -zx 2 -zy 2 $m $n ../data/shifts_fig9c.txt
```
   
### Least-squares reconstruction using erroneous displacements (reproduce Figure 10 of the companion article)
   
In the following experiment, we perform least-squares reconstruction
from inexact sequences of displacements, which corresponds to a
similar experiment as that displayed in Figure 10.
   
Place yourself in the [`src`](src) directory of this package and
run the following bash commands:

```bash
# compute a synthetic sequence containing L=20 low-resolution images
# corrupted by an additive Gaussian noise (stdev = 2), the sequence
# is not apodized yet
./random-shifts -m -5 -M 5 20 /tmp/shifts.txt
./simulator -zx 2 -zy 2 ../data/bridge.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif
./tiffextract -r -i0 0 -i1 19 /tmp/u0_nonrealistic.tif /tmp/u0.tif 121 54 225 162
./tiffaddnoise -g 2 /tmp/u0.tif /tmp/u0_noisy.tif

# simulate erroneous displacement sequences (this models errors in the
# estimations of the displacements) by adding Gaussian noise with
# different stdev to the displacement sequence. Then, perform the
# reconstruction (stack apodization followed by super-resolution)

# displacement noise stdev = 0.01 (display output image and its spectrum)
./shifts-noise -g 0.01 /tmp/shifts.txt /tmp/shifts_sig0.01.txt
./stack-apodization -zx 2 -zy 2 -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts_sig0.01.txt
./leastsquares-superres -zx 2 -zy 2 /tmp/u0_apod.tif /tmp/shifts_sig0.01.txt /tmp/out_sig0.01.tif
./tiffdft -l /tmp/dft_out_sig0.01.tif /tmp/out_sig0.01.tif 
imagej /tmp/out_sig0.01.tif /tmp/dft_out_sig0.01.tif;

# displacement noise stdev = 0.05 (display output image and its spectrum)
./shifts-noise -g 0.05 /tmp/shifts.txt /tmp/shifts_sig0.05.txt
./stack-apodization -zx 2 -zy 2 -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts_sig0.05.txt
./leastsquares-superres -zx 2 -zy 2 /tmp/u0_apod.tif /tmp/shifts_sig0.05.txt /tmp/out_sig0.05.tif
./tiffdft -l /tmp/dft_out_sig0.05.tif /tmp/out_sig0.05.tif 
imagej /tmp/out_sig0.05.tif /tmp/dft_out_sig0.05.tif; 

# displacement noise stdev = 0.1 (display output image and its spectrum)
./shifts-noise -g 0.1 /tmp/shifts.txt /tmp/shifts_sig0.1.txt
./stack-apodization -zx 2 -zy 2 -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts_sig0.1.txt
./leastsquares-superres -zx 2 -zy 2 /tmp/u0_apod.tif /tmp/shifts_sig0.1.txt /tmp/out_sig0.1.tif
./tiffdft -l /tmp/dft_out_sig0.1.tif /tmp/out_sig0.1.tif 
imagej /tmp/out_sig0.1.tif /tmp/dft_out_sig0.1.tif;

# displacement noise stdev = 0.2 (display output image and its spectrum)
./shifts-noise -g 0.2 /tmp/shifts.txt /tmp/shifts_sig0.2.txt
./stack-apodization -zx 2 -zy 2 -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts_sig0.2.txt
./leastsquares-superres -zx 2 -zy 2 /tmp/u0_apod.tif /tmp/shifts_sig0.2.txt /tmp/out_sig0.2.tif
./tiffdft -l /tmp/dft_out_sig0.2.tif /tmp/out_sig0.2.tif 
imagej /tmp/out_sig0.2.tif /tmp/dft_out_sig0.2.tif;

# magnify (and crop) output images using the Shannon interpolation to
# better visualize the artifacts in the spatial domain
./tiffzoom -z 4 /tmp/out_sig0.01.tif /tmp/out_sig0.01_zoom.tif
./tiffzoom -z 4 /tmp/out_sig0.05.tif /tmp/out_sig0.05_zoom.tif
./tiffzoom -z 4 /tmp/out_sig0.1.tif /tmp/out_sig0.1_zoom.tif
./tiffzoom -z 4 /tmp/out_sig0.2.tif /tmp/out_sig0.2_zoom.tif
./tiffextract -r /tmp/out_sig0.01_zoom.tif /tmp/out_sig0.01_zoom.tif 1296  336  340  244
./tiffextract -r /tmp/out_sig0.05_zoom.tif /tmp/out_sig0.05_zoom.tif 1296  336  340  244
./tiffextract -r /tmp/out_sig0.1_zoom.tif /tmp/out_sig0.1_zoom.tif 1296  336  340  244
./tiffextract -r /tmp/out_sig0.2_zoom.tif /tmp/out_sig0.2_zoom.tif 1296  336  340  244
imagej /tmp/out_sig0.01_zoom.tif /tmp/out_sig0.05_zoom.tif /tmp/out_sig0.1_zoom.tif /tmp/out_sig0.2_zoom.tif
```

### Least-squares reconstruction over real data (reproduce Figure 12 of the companion article)

This experiments intents to perform super-resolution with factor 1.8
along both dimensions (zx=zy=1.8) from a real-life sequence containing
L = 500 thermal infrared images that we acquired by ourselves (the
associated sequence of displacements was estimated from the
low-resolution sequence using Keren's algorithm [3]). This experiments
partially reproduces Figure 12 of the companion article.

Place yourself in the [`src`](src) directory of this package and
run the following bash commands:

```bash
# normalize and apodize the FLIR sequence to avoid edge effects in the
# reconstruction
./tiffaxpb -n ../data/flirT640.tif /tmp/u0.tif
./stack-apodization -zx 1.8 -zy 1.8 -a /tmp/u0_apod.tif /tmp/u0.tif ../data/shifts_flirT640.txt

# extract the first image of the apodized low-resolution sequence
./tiffextract -i0 0 -i1 0 -r /tmp/u0_apod.tif /tmp/first.tif 0 0 200 200 

# compute the shift-and-median, i.e., the median of the registered
# low-resolution sequence
./register-stack -m /tmp/shiftandmedian.tif /tmp/u0_apod.tif ../data/shifts_flirT640.txt

# perform least-squares super-resolution
./leastsquares-superres -zx 1.8 -zy 1.8 /tmp/u0_apod.tif ../data/shifts_flirT640.txt /tmp/leastsquares.tif

# display and compare those three images (reproduce the first row of Figure 12)
imagej /tmp/first.tif /tmp/shiftandmedian.tif /tmp/leastsquares.tif 
	
# display and compare their DFT spectra (reproduce the last row of Figure 12)
./tiffdft -l /tmp/first_dft.tif /tmp/first.tif 
./tiffdft -l /tmp/shiftandmedian_dft.tif /tmp/shiftandmedian.tif 
./tiffdft -l /tmp/leastsquares_dft.tif /tmp/leastsquares.tif 
imagej /tmp/first_dft.tif /tmp/shiftandmedian_dft.tif /tmp/leastsquares_dft.tif
```

### IRLS & Lucky imaging (reproduce Figure 14 of the companion article)

We reproduce here a similar experiment to that proposed in Figure 14
of the companion article (again, the result will formally slightly
differ due to the random generation of the experiment
parameters). This experiment illustrates how the lucky-imaging
procedure can be used to remove outliers from the initial sequence in
order to improve the reconstruction quality.

Place yourself in the [`src`](src) directory of this package and
run the following bash commands:
   
```bash
# compute a synthetic sequence containing L=20 low-resolution images
# corrupted by an additive Gaussian noise (stdev = 2), the sequence
# is not apodized yet
./random-shifts -m -5 -M 5 20 /tmp/shifts.txt
./simulator -zx 2 -zy 2 ../data/bridge.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif
./tiffextract -r -i0 0 -i1 19 /tmp/u0_nonrealistic.tif /tmp/u0.tif 121 54 225 162
./tiffaddnoise -g 2 /tmp/u0.tif /tmp/u0_noisy.tif

# compute perturbated displacement sequence by adding a small amounts
# of Gaussian noise (stdev=0.01) to 14 of those displacements and a
# large amount of Gaussian noise (stdev=3) to the remaining 6
# displacements
./shifts-noise -G 0.3 -g 0.01 -n 6 /tmp/shifts.txt /tmp/shifts_perturbated.txt
./stack-apodization -zx 2 -zy 2 -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts_perturbated.txt
./leastsquares-superres -zx 2 -zy 2 /tmp/u0_apod.tif /tmp/shifts_perturbated.txt /tmp/out_ls.tif

# use the IRLS Algorithm to compute the minimizer of the l1-l2 residual
./irls -w /tmp/weights.txt -v /tmp/u0_apod.tif /tmp/shifts_perturbated.txt /tmp/out_l1l2.tif

# apply the lucky imaging procedure: remove from the sequence the 6
# images corresponding to the 6 smallest weights computed with the
# IRLS Algorithm and reconstruct the super-resolved image using the 14
# remaining low-resolution images
./luckyimaging -n 14 /tmp/u0_apod.tif /tmp/shifts_perturbated.txt /tmp/weights.txt /tmp/out_lucky.tif

# display the computed weights (ideally, the 6 smallest weights should
# correspond to the 6 low-resolution images associated with large
# perturbation on the displacement)
cat /tmp/weights.txt

# display and compare the images at the pixel scale
imagej /tmp/out_ls.tif /tmp/out_l1l2.tif /tmp/out_lucky.tif

# zoom and compare the images at the subpixel scale
./tiffzoom -z 4 /tmp/out_ls.tif /tmp/out_ls_zoom.tif
./tiffzoom -z 4 /tmp/out_l1l2.tif /tmp/out_l1l2_zoom.tif
./tiffzoom -z 4 /tmp/out_lucky.tif /tmp/out_lucky_zoom.tif
./tiffextract -r /tmp/out_ls_zoom.tif /tmp/out_ls_zoom.tif 1296 336 338 242
./tiffextract -r /tmp/out_l1l2_zoom.tif /tmp/out_l1l2_zoom.tif 1296 336 338 242
./tiffextract -r /tmp/out_lucky_zoom.tif /tmp/out_lucky_zoom.tif 1296 336 338 242
imagej /tmp/out_ls_zoom.tif /tmp/out_l1l2_zoom.tif /tmp/out_lucky_zoom.tif

# display and compare spectra 
./tiffdft -l /tmp/dft_ls.tif /tmp/out_ls.tif
./tiffdft -l /tmp/dft_l1l2.tif /tmp/out_l1l2.tif
./tiffdft -l /tmp/dft_lucky.tif /tmp/out_lucky.tif
./tiffaxpb -a -1.2 -b -24 /tmp/dft_ls.tif /tmp/dft_ls.tif
./tiffaxpb -a -1.2 -b -24 /tmp/dft_l1l2.tif /tmp/dft_l1l2.tif
./tiffaxpb -a -1.2 -b -24 /tmp/dft_lucky.tif /tmp/dft_lucky.tif
imagej /tmp/dft_ls.tif /tmp/dft_l1l2.tif /tmp/dft_lucky.tif
```

**Remarks**:

+ in the previous synthetic experiment, you can use the `-f` option of
  module shift-noise to force the 'n' large perturbations to affect
  the 'n' first 2D displacements of the sequence. Then, it will be
  easier to check wether those large perturbations are well identified
  by the IRLS procedure (i.e. the 'n' first weight components should
  correspond to the first 'n' weight component) or not.

+ without option `-n`, the `luckyimaging` module generates a movie
  instead of the single image u_lucky^{n}. The frame #j (for j =
  0,1,...) of this movie will contain the image u_lucky^{nj} (see
  module documentation for information about how each nj is being
  computed). This option is usefull when we have no idea of the number
  of image to keep. The proper number of images to keep can be then
  determined from a visual inspection of this movie (in the spatial
  and/or in the Fourier domain)
   
See another example below:

```bash
# compute a synthetic sequence containing L=20 low-resolution images
# corrupted by an additive Gaussian noise (stdev = 2), the sequence
# is not apodized yet
./random-shifts -m -5 -M 5 20 /tmp/shifts.txt
./simulator -zx 2 -zy 2 ../data/bridge.tif /tmp/shifts.txt /tmp/u0_nonrealistic.tif
./tiffextract -r -i0 0 -i1 19 /tmp/u0_nonrealistic.tif /tmp/u0.tif 121 54 225 162
./tiffaddnoise -g 2 /tmp/u0.tif /tmp/u0_noisy.tif
	
# compute perturbated displacement sequence by adding a large amount
# of Gaussian noise (stdev=3) to the 6 *FIRST* displacements (use
# option -f of the shift-noise module) and a small amounts of Gaussian
# noise (stdev=0.01) to the 14 remaining displacements
./shifts-noise -f -G 0.3 -g 0.01 -n 6 /tmp/shifts.txt /tmp/shifts_perturbated.txt
./stack-apodization -zx 2 -zy 2 -a /tmp/u0_apod.tif /tmp/u0_noisy.tif /tmp/shifts_perturbated.txt
./leastsquares-superres -zx 2 -zy 2 /tmp/u0_apod.tif /tmp/shifts_perturbated.txt /tmp/out_ls.tif
	
# use the IRLS Algorithm to compute the minimizer of the l1-l2 residual
./irls -w /tmp/weights.txt -v /tmp/u0_apod.tif /tmp/shifts_perturbated.txt /tmp/out_l1l2.tif 
	
# display the computed weights (ideally, the 6 smallest weights should
# correspond to the 6 low-resolution images associated with large
# perturbation on the displacement). Since we used the -f option of
# the 'shift-noise' module to force the large displacements to affect
# the 6 first low-resolution images, hopefully we should observe that
# the first 6 weights computed using the IRLS procedure correspond to
# the 6 smallest weights.
cat /tmp/weights.txt
	
# apply the lucky imaging procedure to compute a movie with frames
# corresponding to reconstructions u_lucky^{nj} with various values of
# nj (see module documentation for more details about the value of nj)
./luckyimaging -v /tmp/u0_apod.tif /tmp/shifts_perturbated.txt /tmp/weights.txt /tmp/out_lucky.tif 
	
# by a visual inspection of this movie in both the spatial and the
# frequency domain, we can try to determine a posteriori what was the
# most appropriate number of image to remove (ideally, we should
# observe disparition of artifacts in the image u_lucky^{6} which
# corresponds the the frame 7 of the ouptut movie when numbering the
# frames from 1 to 11).
./tiffdft -l /tmp/out_lucky_dft.tif /tmp/out_lucky.tif
imagej /tmp/out_lucky.tif /tmp/out_lucky_dft.tif
```

### Super-resolution and deconvolution (reproduce Figure 17 of the companion article)
	
The next experiment illustrates the benefit of applying frequency
enhancement procedure (i.e., the sharpening filter described in
Section 7 of the companion article) after the least-squares
super-resolution reconstruction process over a real data sequence (the
FLIR T640 thermal infrared image sequence). This experiments
reproduces Figure 17 of the companion article.
	
Place yourself in the [`src`](src) directory of this package and
run the following bash commands

```bash
# normalize & apodize the FLIR sequence to avoid edge effects in the
# reconstruction
./tiffaxpb -n ../data/flirT640.tif /tmp/u0.tif
./stack-apodization -zx 1.8 -zy 1.8 -a /tmp/u0_apod.tif /tmp/u0.tif ../data/shifts_flirT640.txt
	
# compute the shift-and-median, i.e., the median of the registered
# low-resolution sequence, and zoom it by a factor 1.8 along both
# dimensions (Shannon interpolation)
./register-stack -m /tmp/shiftandmedian.tif /tmp/u0_apod.tif ../data/shifts_flirT640.txt
./tiffzoom -z 1.8 /tmp/shiftandmedian.tif /tmp/shiftandmedian_zoomed.tif
	
# perform least-squares super-resolution
./leastsquares-superres -zx 1.8 -zy 1.8 /tmp/u0_apod.tif ../data/shifts_flirT640.txt /tmp/leastsquares.tif
	
# perform frequency amplification (sharpening)
./sharpening -l 5 /tmp/shiftandmedian_zoomed.tif /tmp/shiftandmedian_zoomed_sharpened.tif
./sharpening -l 5 /tmp/leastsquares.tif /tmp/leastsquares_sharpened.tif
	
# display and compare high-resolution images 
imagej /tmp/shiftandmedian_zoomed.tif /tmp/shiftandmedian_zoomed_sharpened.tif /tmp/leastsquares.tif /tmp/leastsquares_sharpened.tif
```
