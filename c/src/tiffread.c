#include <stdlib.h>
#include <stdio.h>
#include <tiffio.h>

/* preprocessor functions */
#define MIN(a,b) (((a)>(b)) ? (b) : (a))
#define ASSERT_ALLOC(cmd) if(NULL == (cmd)) \
                          { \
			    printf("Not enough memory\n"); \
			    return EXIT_FAILURE;	   \
			  }
#define RAISE_ERROR(message) TIFFClose(tiffp); printf("%s",message);


/* internal modules */
int load_tiff_into_double(double***,int*,int*,int*,char*,char);
int load_monopage_tiff_into_double(double**,int*,int*,char*,char);

static int _load_monopage_tiff_striped_uint8_into_double(double*,TIFF*);
static int _load_monopage_tiff_striped_int8_into_double(double*,TIFF*);
static int _load_monopage_tiff_striped_uint16_into_double(double*,TIFF*);
static int _load_monopage_tiff_striped_int16_into_double(double*,TIFF*);
static int _load_monopage_tiff_striped_uint32_into_double(double*,TIFF*);
static int _load_monopage_tiff_striped_int32_into_double(double*,TIFF*);
static int _load_monopage_tiff_striped_float_into_double(double*,TIFF*);

static int _load_monopage_tiff_tiled_uint8_into_double(double*,TIFF*,const uint32,const uint32);
static int _load_monopage_tiff_tiled_int8_into_double(double*,TIFF*,const uint32,const uint32);
static int _load_monopage_tiff_tiled_uint16_into_double(double*,TIFF*,const uint32,const uint32);
static int _load_monopage_tiff_tiled_int16_into_double(double*,TIFF*,const uint32,const uint32);
static int _load_monopage_tiff_tiled_uint32_into_double(double*,TIFF*,const uint32,const uint32);
static int _load_monopage_tiff_tiled_int32_into_double(double*,TIFF*,const uint32,const uint32);
static int _load_monopage_tiff_tiled_float_into_double(double*,TIFF*,const uint32,const uint32);

/* globale variables */
static uint16 bitspersample,sampleformat = SAMPLEFORMAT_UINT;
static uint32 nx,ny;


/* load the "n" images of a multipages TIFF (n is the number of pages,
   n=1 is allowed) and cast all values in double

   the images corresponding to the page index n=0...npages-1 is stored
   in (*output)[n], so that the value for any pixel position (x,y)
   such as

   x in {0,...,width-1}, y in {0,...,height-1},

   (*output)[n][y*width+x] = the graylevel value of the image with
                             page index n at position (x,y).

*/
int load_tiff_into_double(double ***output, int *width, int *height, int *npages, char *filename, char verbose)
{

  TIFF *tiffp;
  unsigned int nnpages = 0;
  uint16 compression_type,samplesperpixel,photo,config;
  uint32 tw,th;
  int page_id; // = page index = 0,1,...,npages-1

  /* consistency checks */

  // open the input TIFF file
  if(verbose) {
    printf("load input TIFF image '%s': ",filename);
    fflush(stdout);
  }
  TIFFSetWarningHandler(NULL);
  tiffp = TIFFOpen(filename,"r");
  if (!tiffp) {
    printf("Failed to open TIFF File '%s'\n",filename);
    return EXIT_FAILURE;
  }
  if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': in progress",filename); fflush(stdout); }

  // retrieve number of pages and check it is >= 1 //
  do ++nnpages; while (TIFFReadDirectory(tiffp));
  if (nnpages < 1) {
    if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': failure.\n",filename); fflush(stdout); }
    RAISE_ERROR("Number of page should be >= 1.\n");
    return EXIT_FAILURE;
  }

  // retrieve image dimension (last page) and output number of pages
  TIFFGetField(tiffp,TIFFTAG_IMAGEWIDTH,width);
  TIFFGetField(tiffp,TIFFTAG_IMAGELENGTH,height);
  *npages = nnpages;

  /* memory allocation four output */
  ASSERT_ALLOC((*output) = (double**) malloc ((*npages)*sizeof(double*)));

  /*** MAIN LOOP : retrieve each page individually ***/
  for (page_id=0;page_id<(*npages);page_id++) {

    // check consistency //
    TIFFSetDirectory(tiffp,page_id);
    TIFFGetField(tiffp,TIFFTAG_IMAGEWIDTH,&nx);
    TIFFGetField(tiffp,TIFFTAG_IMAGELENGTH,&ny);
    if((*width != nx)||(*height != ny)) {
      if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': failure.\n",filename); fflush(stdout); }
      RAISE_ERROR("Inconsistent size between the pages of the Multipage TIFF image (all pages must have the same width and height).\n");
      return EXIT_FAILURE;
    }
    TIFFGetFieldDefaulted(tiffp,TIFFTAG_SAMPLESPERPIXEL,&samplesperpixel);
    if (samplesperpixel != 1) {
      if (verbose) { printf("\33[2K\rLoad input TIFF image '%s' (page %d): failure.\n",filename,page_id); fflush(stdout); }
      RAISE_ERROR("Multichannel TIFF files are not unsupported (only multipage monochannel TIFF files are supported).\n");
      return EXIT_FAILURE;
    }
    TIFFGetFieldDefaulted(tiffp,TIFFTAG_SAMPLEFORMAT,&sampleformat);
    TIFFGetFieldDefaulted(tiffp,TIFFTAG_BITSPERSAMPLE,&bitspersample);
    TIFFGetFieldDefaulted(tiffp,TIFFTAG_COMPRESSION,&compression_type);

    if((compression_type != COMPRESSION_NONE)&&(compression_type != COMPRESSION_LZW)) {
      if (verbose) { printf("\33[2K\rLoad input TIFF image '%s' (page %d): failure.\n",filename,page_id); fflush(stdout); }
      printf("Unsupported TIFF compression type (%d). This reader only supports Uncompressed (compression type = %d) and LZW-compressed (compression_type =%d) images.",compression_type,COMPRESSION_NONE,COMPRESSION_LZW);
      RAISE_ERROR("Unsupported TIFF compression (only non-compressed or LZW compressed TIFF images are supported).");
      return EXIT_FAILURE;
    }

    if (!((sampleformat==SAMPLEFORMAT_UINT)||(sampleformat==SAMPLEFORMAT_INT)||(sampleformat==SAMPLEFORMAT_IEEEFP)||(sampleformat==SAMPLEFORMAT_VOID))) {
      if (verbose) { printf("\33[2K\rLoad input TIFF image '%s' (page %d): failure.\n",filename,page_id); fflush(stdout); }
      printf("Unsupported TIFF sample format (%d), this reader only supports (u)int and float images (sample format = %d, %d or %d).\nNote also that undefined data format are accepted (sample format=%d) but they are loaded as unsigned integer.",sampleformat,SAMPLEFORMAT_UINT,SAMPLEFORMAT_INT,SAMPLEFORMAT_IEEEFP,SAMPLEFORMAT_VOID);
      TIFFClose(tiffp);
      return EXIT_FAILURE;
    }

    if(sampleformat==SAMPLEFORMAT_VOID) sampleformat = SAMPLEFORMAT_UINT;


    /* memory allocation for (*output)[page_id] */
    ASSERT_ALLOC((*output)[page_id] = (double*) malloc (nx*ny*sizeof(double)));

    /* retrieve graylevels of the current page and  store them in double precision */
    TIFFGetField(tiffp,TIFFTAG_PLANARCONFIG,&config);
    TIFFGetField(tiffp,TIFFTAG_PHOTOMETRIC,&photo);

    if (TIFFIsTiled(tiffp)) { // tiled TIFF case
      TIFFGetField(tiffp,TIFFTAG_TILEWIDTH,&tw);
      TIFFGetField(tiffp,TIFFTAG_TILELENGTH,&th);
      switch (bitspersample) {
      case 8 :
	if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_tiled_uint8_into_double((*output)[page_id],tiffp,tw,th);
	else _load_monopage_tiff_tiled_int8_into_double((*output)[page_id],tiffp,tw,th);
	break;
      case 16 :
	if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_tiled_uint16_into_double((*output)[page_id],tiffp,tw,th);
	else _load_monopage_tiff_tiled_int16_into_double((*output)[page_id],tiffp,tw,th);
	break;
      case 32 :
	if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_tiled_uint32_into_double((*output)[page_id],tiffp,tw,th);
	else if (sampleformat==SAMPLEFORMAT_INT) _load_monopage_tiff_tiled_int32_into_double((*output)[page_id],tiffp,tw,th);
	else _load_monopage_tiff_tiled_float_into_double((*output)[page_id],tiffp,tw,th);
	break;
      default :
	if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': failure.\n",filename); fflush(stdout); }
	return EXIT_FAILURE;
      }
    }
    else { // striped (=non-tiled) TIFF case
      switch (bitspersample) {
      case 8 :
	if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_striped_uint8_into_double((*output)[page_id],tiffp);
	else _load_monopage_tiff_striped_int8_into_double((*output)[page_id],tiffp);
	break;
      case 16 :
	if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_striped_uint16_into_double((*output)[page_id],tiffp);
	else _load_monopage_tiff_striped_int16_into_double((*output)[page_id],tiffp);
	break;
      case 32 :
	if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_striped_uint32_into_double((*output)[page_id],tiffp);
	else if (sampleformat==SAMPLEFORMAT_INT) _load_monopage_tiff_striped_int32_into_double((*output)[page_id],tiffp);
	else _load_monopage_tiff_striped_float_into_double((*output)[page_id],tiffp);
	break;
      default :
	if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': failure.\n",filename); fflush(stdout); }
	return EXIT_FAILURE;
      }
    }

    /* deal with verbose mode */
    if(verbose) {
      printf("\33[2K\rload input TIFF image '%s' (page %d): done\n",filename,page_id);
      fflush(stdout);
      printf(" > width x height = %d x %d\n",nx,ny);
      printf(" > precision = %d bits per pixel %s\n",bitspersample,sampleformat==SAMPLEFORMAT_UINT ? "(unsigned)" : "(signed)");
      if (sampleformat==SAMPLEFORMAT_IEEEFP) printf(" > this image is a tiff-float image\n");
    }

  }

  TIFFClose(tiffp);

  return EXIT_SUCCESS;
}

/* load the monopage TIFF pixel values into a double pointer */
int load_monopage_tiff_into_double(double **output, int *width, int *height, char *filename, char verbose)
{

  TIFF *tiffp;
  unsigned int npages = 0;
  uint16 compression_type,samplesperpixel,photo,config;
  uint32 tw,th;

  /* consistency checks */

  // open the input TIFF file
  if(verbose) {
    printf("load input TIFF image '%s': ",filename);
    fflush(stdout);
  }
  TIFFSetWarningHandler(NULL);
  tiffp = TIFFOpen(filename,"r");
  if (!tiffp) {
    printf("Failed to open TIFF File '%s'\n",filename);
    return EXIT_FAILURE;
  }
  if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': in progress",filename); fflush(stdout); }

  // retrieve number of pages and check consistency //
  do ++npages; while (TIFFReadDirectory(tiffp));
  if(npages > 1) {
    if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': failure.\n",filename); fflush(stdout); }
    RAISE_ERROR("Input image is not monopage, you must load it using the 'load_tiff_into_double' module.\n");
    return EXIT_FAILURE;
  }

  // retrieve image dimension
  TIFFSetDirectory(tiffp,0);
  TIFFGetField(tiffp,TIFFTAG_IMAGEWIDTH,width);
  TIFFGetField(tiffp,TIFFTAG_IMAGELENGTH,height);
  nx = *width; ny = *height;

  /* allocate memory and fill output */
  ASSERT_ALLOC((*output) = (double*) malloc ((*width)*(*height) * sizeof(double)));

  TIFFGetFieldDefaulted(tiffp,TIFFTAG_SAMPLESPERPIXEL,&samplesperpixel);
  if (samplesperpixel != 1) {
    if (verbose) { printf("\33[2K\rLoad input TIFF image '%s' : failure.\n",filename); fflush(stdout); }
    RAISE_ERROR("Multichannel TIFF files are not unsupported (only monopage and monochannel TIFF files are supported).\n");
    return EXIT_FAILURE;
  }
  TIFFGetFieldDefaulted(tiffp,TIFFTAG_SAMPLEFORMAT,&sampleformat);
  TIFFGetFieldDefaulted(tiffp,TIFFTAG_BITSPERSAMPLE,&bitspersample);
  TIFFGetFieldDefaulted(tiffp,TIFFTAG_COMPRESSION,&compression_type);

  if((compression_type != COMPRESSION_NONE)&&(compression_type != COMPRESSION_LZW)) {
    if (verbose) { printf("\33[2K\rLoad input TIFF image '%s' : failure.\n",filename); fflush(stdout); }
    printf("Unsupported TIFF compression type (%d). This reader only supports Uncompressed (compression type = %d) and LZW-compressed (compression_type =%d) images.",compression_type,COMPRESSION_NONE,COMPRESSION_LZW);
    RAISE_ERROR("Unsupported TIFF compression (only non-compressed or LZW compressed TIFF images are supported).");
    return EXIT_FAILURE;
  }

  if (!((sampleformat==SAMPLEFORMAT_UINT)||(sampleformat==SAMPLEFORMAT_INT)||(sampleformat==SAMPLEFORMAT_IEEEFP)||(sampleformat==SAMPLEFORMAT_VOID))) {
    if (verbose) { printf("\33[2K\rLoad input TIFF image '%s' : failure.\n",filename); fflush(stdout); }
    printf("Unsupported TIFF sample format (%d), this reader only supports (u)int and float images (sample format = %d, %d or %d).\nNote also that undefined data format are accepted (sample format=%d) but they are loaded as unsigned integer.",sampleformat,SAMPLEFORMAT_UINT,SAMPLEFORMAT_INT,SAMPLEFORMAT_IEEEFP,SAMPLEFORMAT_VOID);
    TIFFClose(tiffp);
    return EXIT_FAILURE;
  }

  if(sampleformat==SAMPLEFORMAT_VOID) sampleformat = SAMPLEFORMAT_UINT;

  /* retrieve graylevels of the current page and  store them in double precision */
  TIFFGetField(tiffp,TIFFTAG_PLANARCONFIG,&config);
  TIFFGetField(tiffp,TIFFTAG_PHOTOMETRIC,&photo);

  if (TIFFIsTiled(tiffp)) { // tiled TIFF case
    TIFFGetField(tiffp,TIFFTAG_TILEWIDTH,&tw);
    TIFFGetField(tiffp,TIFFTAG_TILELENGTH,&th);
    switch (bitspersample) {
    case 8 :
      if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_tiled_uint8_into_double(*output,tiffp,tw,th);
      else _load_monopage_tiff_tiled_int8_into_double(*output,tiffp,tw,th);
      break;
    case 16 :
      if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_tiled_uint16_into_double(*output,tiffp,tw,th);
      else _load_monopage_tiff_tiled_int16_into_double(*output,tiffp,tw,th);
      break;
    case 32 :
      if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_tiled_uint32_into_double(*output,tiffp,tw,th);
      else if (sampleformat==SAMPLEFORMAT_INT) _load_monopage_tiff_tiled_int32_into_double(*output,tiffp,tw,th);
      else _load_monopage_tiff_tiled_float_into_double(*output,tiffp,tw,th);
      break;
    default :
      if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': failure.\n",filename); fflush(stdout); }
      return EXIT_FAILURE;
    }
  }
  else { // striped (=non-tiled) TIFF case
    switch (bitspersample) {
    case 8 :
      if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_striped_uint8_into_double(*output,tiffp);
      else _load_monopage_tiff_striped_int8_into_double(*output,tiffp);
      break;
    case 16 :
      if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_striped_uint16_into_double(*output,tiffp);
      else _load_monopage_tiff_striped_int16_into_double(*output,tiffp);
      break;
    case 32 :
      if (sampleformat==SAMPLEFORMAT_UINT) _load_monopage_tiff_striped_uint32_into_double(*output,tiffp);
      else if (sampleformat==SAMPLEFORMAT_INT) _load_monopage_tiff_striped_int32_into_double(*output,tiffp);
      else _load_monopage_tiff_striped_float_into_double(*output,tiffp);
      break;
    default :
      if (verbose) { printf("\33[2K\rLoad input TIFF image '%s': failure.\n",filename); fflush(stdout); }
      return EXIT_FAILURE;
    }
  }

  /* deal with verbose mode */
  if(verbose) {
    printf("\33[2K\rload input TIFF image '%s' : done\n",filename);
    fflush(stdout);
    printf(" > width x height = %d x %d\n",*width,*height);
    printf(" > precision = %d bits per pixel %s\n",bitspersample,sampleformat==SAMPLEFORMAT_UINT ? "(unsigned)" : "(signed)");
    if (sampleformat==SAMPLEFORMAT_IEEEFP) printf(" > this image is a tiff-float image\n");
  }

  TIFFClose(tiffp);

  return EXIT_SUCCESS;
}



/***************************************************************************************
 * single page TIFF loaders for:
 *
 *  - striped (=non-tiled) images
 *  - with sample format: uint8, int8, uint16, int16, uint32, int32 or float
 *
 ***************************************************************************************/

/* load a striped uint8 image into a pointer of double */
static int _load_monopage_tiff_striped_uint8_into_double(double *graylevels,TIFF *tiffp)
{
  unsigned int rr,cc;
  tstrip_t strip;
  uint32 row, rowsperstrip = (uint32)-1, nrow;
  uint8 *buf;

  ASSERT_ALLOC( buf = (uint8*)_TIFFmalloc(TIFFStripSize(tiffp)));

  TIFFGetField(tiffp,TIFFTAG_ROWSPERSTRIP,&rowsperstrip);

  for (row = 0; row<ny; row+= rowsperstrip) {
    nrow = (row+rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp, row, 0);
    if ((TIFFReadEncodedStrip(tiffp,strip,buf,-1))<0) { _TIFFfree(buf); RAISE_ERROR("Invalid strip found into the TIFF image\n"); }
    for (rr = 0;rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc)
	graylevels[(row+rr)*nx+cc] = (double) buf[rr*nx+cc];
  }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a striped int8 image into a pointer of double */
static int _load_monopage_tiff_striped_int8_into_double(double *graylevels,TIFF *tiffp)
{
  unsigned int rr,cc;
  tstrip_t strip;
  uint32 row, rowsperstrip = (uint32)-1, nrow;
  int8 *buf;

  ASSERT_ALLOC( buf = (int8*)_TIFFmalloc(TIFFStripSize(tiffp)));

  TIFFGetField(tiffp,TIFFTAG_ROWSPERSTRIP,&rowsperstrip);

  for (row = 0; row<ny; row+= rowsperstrip) {
    nrow = (row+rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp, row, 0);
    if ((TIFFReadEncodedStrip(tiffp,strip,buf,-1))<0) { _TIFFfree(buf); RAISE_ERROR("Invalid strip found into the TIFF image\n"); }
    for (rr = 0;rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc)
	graylevels[(row+rr)*nx+cc] = (double) buf[rr*nx+cc];
  }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a striped uint16 image into a pointer of double */
static int _load_monopage_tiff_striped_uint16_into_double(double *graylevels,TIFF *tiffp)
{
  unsigned int rr,cc;
  tstrip_t strip;
  uint32 row, rowsperstrip = (uint32)-1, nrow;
  uint16 *buf;

  ASSERT_ALLOC( buf = (uint16*)_TIFFmalloc(TIFFStripSize(tiffp)));

  TIFFGetField(tiffp,TIFFTAG_ROWSPERSTRIP,&rowsperstrip);

  for (row = 0; row<ny; row+= rowsperstrip) {
    nrow = (row+rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp, row, 0);
    if ((TIFFReadEncodedStrip(tiffp,strip,buf,-1))<0) { _TIFFfree(buf); RAISE_ERROR("Invalid strip found into the TIFF image\n"); }
    for (rr = 0;rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc)
	graylevels[(row+rr)*nx+cc] = (double) buf[rr*nx+cc];
  }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a striped int16 image into a pointer of double */
static int _load_monopage_tiff_striped_int16_into_double(double *graylevels,TIFF *tiffp)
{
  unsigned int rr,cc;
  tstrip_t strip;
  uint32 row, rowsperstrip = (uint32)-1, nrow;
  int16 *buf;

  ASSERT_ALLOC( buf = (int16*)_TIFFmalloc(TIFFStripSize(tiffp)));

  TIFFGetField(tiffp,TIFFTAG_ROWSPERSTRIP,&rowsperstrip);

  for (row = 0; row<ny; row+= rowsperstrip) {
    nrow = (row+rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp, row, 0);
    if ((TIFFReadEncodedStrip(tiffp,strip,buf,-1))<0) { _TIFFfree(buf); RAISE_ERROR("Invalid strip found into the TIFF image\n"); }
    for (rr = 0;rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc)
	graylevels[(row+rr)*nx+cc] = (double) buf[rr*nx+cc];
  }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a striped uint32 image into a pointer of double */
static int _load_monopage_tiff_striped_uint32_into_double(double *graylevels,TIFF *tiffp)
{
  unsigned int rr,cc;
  tstrip_t strip;
  uint32 row, rowsperstrip = (uint32)-1, nrow;
  uint32 *buf;

  ASSERT_ALLOC( buf = (uint32*)_TIFFmalloc(TIFFStripSize(tiffp)));

  TIFFGetField(tiffp,TIFFTAG_ROWSPERSTRIP,&rowsperstrip);

  for (row = 0; row<ny; row+= rowsperstrip) {
    nrow = (row+rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp, row, 0);
    if ((TIFFReadEncodedStrip(tiffp,strip,buf,-1))<0) { _TIFFfree(buf); RAISE_ERROR("Invalid strip found into the TIFF image\n"); }
    for (rr = 0;rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc)
	graylevels[(row+rr)*nx+cc] = (double) buf[rr*nx+cc];
  }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a striped int32 image into a pointer of double */
static int _load_monopage_tiff_striped_int32_into_double(double *graylevels,TIFF *tiffp)
{
  unsigned int rr,cc;
  tstrip_t strip;
  uint32 row, rowsperstrip = (uint32)-1, nrow;
  int32 *buf;

  ASSERT_ALLOC( buf = (int32*)_TIFFmalloc(TIFFStripSize(tiffp)));

  TIFFGetField(tiffp,TIFFTAG_ROWSPERSTRIP,&rowsperstrip);

  for (row = 0; row<ny; row+= rowsperstrip) {
    nrow = (row+rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp, row, 0);
    if ((TIFFReadEncodedStrip(tiffp,strip,buf,-1))<0) { _TIFFfree(buf); RAISE_ERROR("Invalid strip found into the TIFF image\n"); }
    for (rr = 0;rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc)
	graylevels[(row+rr)*nx+cc] = (double) buf[rr*nx+cc];
  }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a striped float image into a pointer of double */
static int _load_monopage_tiff_striped_float_into_double(double *graylevels,TIFF *tiffp)
{
  unsigned int rr,cc;
  tstrip_t strip;
  uint32 row, rowsperstrip = (uint32)-1, nrow;
  float *buf;

  ASSERT_ALLOC( buf = (float*)_TIFFmalloc(TIFFStripSize(tiffp)));

  TIFFGetField(tiffp,TIFFTAG_ROWSPERSTRIP,&rowsperstrip);

  for (row = 0; row<ny; row+= rowsperstrip) {
    nrow = (row+rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp, row, 0);
    if ((TIFFReadEncodedStrip(tiffp,strip,buf,-1))<0) { _TIFFfree(buf); RAISE_ERROR("Invalid strip found into the TIFF image\n"); }
    for (rr = 0;rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc)
	graylevels[(row+rr)*nx+cc] = (double) buf[rr*nx+cc];
  }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/***************************************************************************************
 * single page TIFF loaders for:
 *
 *  - tiled (=non-stripped) images
 *  - with sample format: uint8, int8, uint16, int16, uint32, int32 or float
 *
 ***************************************************************************************/

/* load a tiled uint8 image into a pointer of double */
static int _load_monopage_tiff_tiled_uint8_into_double(double *graylevels, TIFF *tiffp, const uint32 tw, const uint32 th)
{
  unsigned int row,col,rr,cc;
  uint8 *buf;

  ASSERT_ALLOC(buf = (uint8*)_TIFFmalloc(TIFFTileSize(tiffp)));

  for (row = 0; row<ny; row+=th)
    for (col = 0; col<nx; col+=tw) {
      if (TIFFReadTile(tiffp,buf,col,row,0,0)<0) { _TIFFfree(buf); RAISE_ERROR("Invalid tile found into the TIFF image\n"); }
      for (rr = row; rr<MIN(row+th,ny); ++rr)
	for (cc = col; cc<MIN(col+tw,nx); ++cc)
	  graylevels[rr*nx+cc] = (double) buf[(rr-row)*tw + (cc-col)];
    }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a tiled int8 image into a pointer of double */
static int _load_monopage_tiff_tiled_int8_into_double(double *graylevels, TIFF *tiffp, const uint32 tw, const uint32 th)
{
  unsigned int row,col,rr,cc;
  int8 *buf;

  ASSERT_ALLOC(buf = (int8*)_TIFFmalloc(TIFFTileSize(tiffp)));

  for (row = 0; row<ny; row+=th)
    for (col = 0; col<nx; col+=tw) {
      if (TIFFReadTile(tiffp,buf,col,row,0,0)<0) { _TIFFfree(buf); RAISE_ERROR("Invalid tile found into the TIFF image\n"); }
      for (rr = row; rr<MIN(row+th,ny); ++rr)
	for (cc = col; cc<MIN(col+tw,nx); ++cc)
	  graylevels[rr*nx+cc] = (double) buf[(rr-row)*tw + (cc-col)];
    }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a tiled uint16 image into a pointer of double */
static int _load_monopage_tiff_tiled_uint16_into_double(double *graylevels, TIFF *tiffp, const uint32 tw, const uint32 th)
{
  unsigned int row,col,rr,cc;
  uint16 *buf;

  ASSERT_ALLOC(buf = (uint16*)_TIFFmalloc(TIFFTileSize(tiffp)));

  for (row = 0; row<ny; row+=th)
    for (col = 0; col<nx; col+=tw) {
      if (TIFFReadTile(tiffp,buf,col,row,0,0)<0) { _TIFFfree(buf); RAISE_ERROR("Invalid tile found into the TIFF image\n"); }
      for (rr = row; rr<MIN(row+th,ny); ++rr)
	for (cc = col; cc<MIN(col+tw,nx); ++cc)
	  graylevels[rr*nx+cc] = (double) buf[(rr-row)*tw + (cc-col)];
    }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a tiled int16 image into a pointer of double */
static int _load_monopage_tiff_tiled_int16_into_double(double *graylevels, TIFF *tiffp, const uint32 tw, const uint32 th)
{
  unsigned int row,col,rr,cc;
  int16 *buf;

  ASSERT_ALLOC(buf = (int16*)_TIFFmalloc(TIFFTileSize(tiffp)));

  for (row = 0; row<ny; row+=th)
    for (col = 0; col<nx; col+=tw) {
      if (TIFFReadTile(tiffp,buf,col,row,0,0)<0) { _TIFFfree(buf); RAISE_ERROR("Invalid tile found into the TIFF image\n"); }
      for (rr = row; rr<MIN(row+th,ny); ++rr)
	for (cc = col; cc<MIN(col+tw,nx); ++cc)
	  graylevels[rr*nx+cc] = (double) buf[(rr-row)*tw + (cc-col)];
    }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a tiled uint32 image into a pointer of double */
static int _load_monopage_tiff_tiled_uint32_into_double(double *graylevels, TIFF *tiffp, const uint32 tw, const uint32 th)
{
  unsigned int row,col,rr,cc;
  uint32 *buf;

  ASSERT_ALLOC(buf = (uint32*)_TIFFmalloc(TIFFTileSize(tiffp)));

  for (row = 0; row<ny; row+=th)
    for (col = 0; col<nx; col+=tw) {
      if (TIFFReadTile(tiffp,buf,col,row,0,0)<0) { _TIFFfree(buf); RAISE_ERROR("Invalid tile found into the TIFF image\n"); }
      for (rr = row; rr<MIN(row+th,ny); ++rr)
	for (cc = col; cc<MIN(col+tw,nx); ++cc)
	  graylevels[rr*nx+cc] = (double) buf[(rr-row)*tw + (cc-col)];
    }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}

/* load a tiled int32 image into a pointer of double */
static int _load_monopage_tiff_tiled_int32_into_double(double *graylevels, TIFF *tiffp, const uint32 tw, const uint32 th)
{
  unsigned int row,col,rr,cc;
  int32 *buf;

  ASSERT_ALLOC(buf = (int32*)_TIFFmalloc(TIFFTileSize(tiffp)));

  for (row = 0; row<ny; row+=th)
    for (col = 0; col<nx; col+=tw) {
      if (TIFFReadTile(tiffp,buf,col,row,0,0)<0) { _TIFFfree(buf); RAISE_ERROR("Invalid tile found into the TIFF image\n"); }
      for (rr = row; rr<MIN(row+th,ny); ++rr)
	for (cc = col; cc<MIN(col+tw,nx); ++cc)
	  graylevels[rr*nx+cc] = (double) buf[(rr-row)*tw + (cc-col)];
    }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}


/* load a tiled float image into a pointer of double */
static int _load_monopage_tiff_tiled_float_into_double(double *graylevels, TIFF *tiffp, const uint32 tw, const uint32 th)
{
  unsigned int row,col,rr,cc;
  float *buf;

  ASSERT_ALLOC(buf = (float*)_TIFFmalloc(TIFFTileSize(tiffp)));

  for (row = 0; row<ny; row+=th)
    for (col = 0; col<nx; col+=tw) {
      if (TIFFReadTile(tiffp,buf,col,row,0,0)<0) { _TIFFfree(buf); RAISE_ERROR("Invalid tile found into the TIFF image\n"); }
      for (rr = row; rr<MIN(row+th,ny); ++rr)
	for (cc = col; cc<MIN(col+tw,nx); ++cc)
	  graylevels[rr*nx+cc] = (double) buf[(rr-row)*tw + (cc-col)];
    }
  _TIFFfree(buf);

  return EXIT_SUCCESS;
}
