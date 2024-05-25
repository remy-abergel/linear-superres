#include <stdlib.h>
#include <stdio.h>
#include <tiffio.h>
#include <float.h>
#include <limits.h>
#include <fftw3.h>

#define UINT8_MIN 0
#define UINT8_MAX UCHAR_MAX
#define INT8_MIN SCHAR_MIN
#define INT8_MAX SCHAR_MAX
#define UINT16_MIN 0
#define UINT16_MAX USHRT_MAX
#define INT16_MIN SHRT_MIN
#define INT16_MAX SHRT_MAX
#define UINT32_MIN 0
#define UINT32_MAX UINT_MAX
#define INT32_MIN INT_MIN
#define INT32_MAX INT_MAX

/* external modules */
int save_monopage_striped_tiff_from_fftw_complex(char*,fftw_complex*,int,int,char,int,char*,char*,char*);

/* internal modules */
static int _save_page_striped_tiff_uint8_from_fftw_complex(TIFF*,fftw_complex*,int,int,uint32);
static int _save_page_striped_tiff_int8_from_fftw_complex(TIFF*,fftw_complex*,int,int,uint32);
static int _save_page_striped_tiff_uint16_from_fftw_complex(TIFF*,fftw_complex*,int,int,uint32);
static int _save_page_striped_tiff_int16_from_fftw_complex(TIFF*,fftw_complex*,int,int,uint32);
static int _save_page_striped_tiff_uint32_from_fftw_complex(TIFF*,fftw_complex*,int,int,uint32);
static int _save_page_striped_tiff_int32_from_fftw_complex(TIFF*,fftw_complex*,int,int,uint32);
static int _save_page_striped_tiff_float_from_fftw_complex(TIFF*,fftw_complex*,int,int,uint32);

/* save the real part of the fftw_complex input signal (pixmap) as a monopage striped image */
int save_monopage_striped_tiff_from_fftw_complex(filename,pixmap,nx,ny,verbose,type,mode,pname,desc)
  char *filename,*mode,verbose,*pname,*desc;
     fftw_complex *pixmap;
     int nx,ny,type;
{
  int err;
  uint32 rowsperstrip;

  TIFF *tiffp = TIFFOpen(filename,mode);

  //TIFFSetDirectory(tiffp,0);
  if(pname) TIFFSetField(tiffp,TIFFTAG_PAGENAME,pname);
  if(desc) TIFFSetField(tiffp,TIFFTAG_IMAGEDESCRIPTION,desc);
  TIFFSetField(tiffp,TIFFTAG_IMAGEWIDTH,nx);
  TIFFSetField(tiffp,TIFFTAG_IMAGELENGTH,ny);
  TIFFSetField(tiffp,TIFFTAG_ORIENTATION,ORIENTATION_TOPLEFT);
  TIFFSetField(tiffp,TIFFTAG_SAMPLESPERPIXEL,(uint16)1);
  TIFFSetField(tiffp,TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
  TIFFSetField(tiffp,TIFFTAG_PHOTOMETRIC,PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tiffp,TIFFTAG_COMPRESSION,(uint16)1);
  rowsperstrip = TIFFDefaultStripSize(tiffp,(uint32)-1);
  TIFFSetField(tiffp,TIFFTAG_ROWSPERSTRIP,rowsperstrip);
  TIFFSetField(tiffp,TIFFTAG_FILLORDER,FILLORDER_MSB2LSB);

  switch (type) { // output TIFF precision (personal convention)
  case 0 : // TIFF float precision
    TIFFSetField(tiffp,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tiffp,TIFFTAG_BITSPERSAMPLE,(uint16)32);
    if(verbose) {
      printf("save '%s' image (precision float): in progress",filename);
      fflush(stdout);
    }
    if (tiffp) err = _save_page_striped_tiff_float_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip);
    else err = EXIT_FAILURE;
    if(verbose) {
      printf("\33[2K\rsave '%s' image (precision float): %s\n",filename,(EXIT_SUCCESS == err)?"done ":"failure");
      fflush(stdout);
    }
    break;
  case 8 : // TIFF uint8 precision
    TIFFSetField(tiffp,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_UINT);
    TIFFSetField(tiffp,TIFFTAG_BITSPERSAMPLE,(uint16)8);
    if(verbose) {
      printf("save '%s' image (precision uint8): in progress",filename);
      fflush(stdout);
    }
    if (tiffp) err = _save_page_striped_tiff_uint8_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip);
    else err = EXIT_FAILURE;
    if(verbose) {
      printf("\33[2K\rsave '%s' image (precision uint8): %s\n",filename,(EXIT_SUCCESS == err)?"done ":"failure");
      fflush(stdout);
    }
    break;
  case -8 : // TIFF int8 precision
    TIFFSetField(tiffp,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_INT);
    TIFFSetField(tiffp,TIFFTAG_BITSPERSAMPLE,(uint16)8);
    if(verbose) {
      printf("save '%s' image (precision int8): in progress",filename);
      fflush(stdout);
    }
    if (tiffp) err = _save_page_striped_tiff_int8_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip);
    else err = EXIT_FAILURE;
    if(verbose) {
      printf("\33[2K\rsave '%s' image (precision int8): %s\n",filename,(EXIT_SUCCESS == err)?"done ":"failure");
      fflush(stdout);
    }
    break;
  case 16 : // TIFF uint16 precision
    TIFFSetField(tiffp,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_UINT);
    TIFFSetField(tiffp,TIFFTAG_BITSPERSAMPLE,(uint16)16);
    if(verbose) {
      printf("save '%s' image (precision uint16): in progress",filename);
      fflush(stdout);
    }
    if (tiffp) err = _save_page_striped_tiff_uint16_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip);
    else err = EXIT_FAILURE;
    if(verbose) {
      printf("\33[2K\rsave '%s' image (precision uint16): %s\n",filename,(EXIT_SUCCESS == err)?"done ":"failure");
      fflush(stdout);
    }
    break;
  case -16 : // TIFF int16 precision
    TIFFSetField(tiffp,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_INT);
    TIFFSetField(tiffp,TIFFTAG_BITSPERSAMPLE,(uint16)16);
    if(verbose) {
      printf("save '%s' image (precision int16): in progress",filename);
      fflush(stdout);
    }
    if (tiffp) err = _save_page_striped_tiff_int16_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip);
    else err = EXIT_FAILURE;
    if(verbose) {
      printf("\33[2K\rsave '%s' image (precision int16): %s\n",filename,(EXIT_SUCCESS == err)?"done ":"failure");
      fflush(stdout);
    }
    break;
  case 32 : // TIFF uint32 precision
    TIFFSetField(tiffp,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_UINT);
    TIFFSetField(tiffp,TIFFTAG_BITSPERSAMPLE,(uint16)32);
    if(verbose) {
      printf("save '%s' image (precision uint32): in progress",filename);
      fflush(stdout);
    }
    if (tiffp) err = _save_page_striped_tiff_uint32_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip);
    else err = EXIT_FAILURE;
    if(verbose) {
      printf("\33[2K\rsave '%s' image (precision uint32): %s\n",filename,(EXIT_SUCCESS == err)?"done ":"failure");
      fflush(stdout);
    }
    break;
  case -32 : // TIFF int32 precision
    TIFFSetField(tiffp,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_INT);
    TIFFSetField(tiffp,TIFFTAG_BITSPERSAMPLE,(uint16)32);
    if(verbose) {
      printf("save '%s' image (precision int32): in progress",filename);
      fflush(stdout);
    }
    if (tiffp) err = _save_page_striped_tiff_int32_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip);
    else err = EXIT_FAILURE;
    if(verbose) {
      printf("\33[2K\rsave '%s' image (precision int32): %s\n",filename,(EXIT_SUCCESS == err)?"done ":"failure");
      fflush(stdout);
    }
    break;
  default : // unrecognized precision
    if(verbose) printf("used unrecognized type in internal module 'save_monopage_striped_tiff_from_double'.\n");
    err = EXIT_FAILURE;
  }

  TIFFClose(tiffp);
  return err;
}



/* save the real part of the input fftw_complex signal (pixmap) as a monopage striped uint8 TIFF image */
static int _save_page_striped_tiff_uint8_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip)
     TIFF *tiffp;
     fftw_complex *pixmap;
     int nx,ny;
     uint32 rowsperstrip;
{
  unsigned int row,rr,cc;
  uint32 nrow;
  tstrip_t strip;
  tsize_t i;
  uint8 *buf;
  double val;
  int nout=0;

  buf = (uint8*) _TIFFmalloc(TIFFStripSize(tiffp));
  if(!buf) { printf("Not enough memory\n"); return EXIT_FAILURE; }

  for (row = 0; row<ny; row+=rowsperstrip){
    nrow = (row + rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp,row,0);
    i = 0;
    for (rr = 0; rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc) {
	val = pixmap[(row+rr)*nx + cc][0];
	if(val > UINT8_MAX) { val = UINT8_MAX; nout++; }
	if(val < UINT8_MIN) { val = UINT8_MIN; nout++; }
	buf[i++] = (uint8) val;
      }
    if (TIFFWriteEncodedStrip(tiffp,strip,buf,i*sizeof(uint8))<0) return EXIT_FAILURE;
  }
  _TIFFfree(buf);

  TIFFWriteDirectory(tiffp);

  if(nout > 0) printf("WARNING: %d graylevels were out of [%d,%d]\n",nout,UINT8_MIN,UINT8_MAX);

  return EXIT_SUCCESS;
}


/* save the real part of the input fftw_complex signal (pixmap) as a monopage striped int8 TIFF image */
static int _save_page_striped_tiff_int8_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip)
     TIFF *tiffp;
     fftw_complex *pixmap;
     int nx,ny;
     uint32 rowsperstrip;
{
  unsigned int row,rr,cc;
  uint32 nrow;
  tstrip_t strip;
  tsize_t i;
  int8 *buf;
  double val;
  int nout=0;

  buf = (int8*) _TIFFmalloc(TIFFStripSize(tiffp));
  if(!buf) { printf("Not enough memory\n"); return EXIT_FAILURE; }

  for (row = 0; row<ny; row+=rowsperstrip){
    nrow = (row + rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp,row,0);
    i = 0;
    for (rr = 0; rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc) {
	val = pixmap[(row+rr)*nx + cc][0];
	if(val > INT8_MAX) { val = INT8_MAX; nout++; }
	if(val < INT8_MIN) { val = INT8_MIN; nout++; }
	buf[i++] = (int8) val;
      }
    if (TIFFWriteEncodedStrip(tiffp,strip,buf,i*sizeof(int8))<0) return EXIT_FAILURE;
  }
  _TIFFfree(buf);

  TIFFWriteDirectory(tiffp);

  if(nout > 0) printf("WARNING: %d graylevels were out of [%d,%d]\n",nout,INT8_MIN,INT8_MAX);

  return EXIT_SUCCESS;
}


/* save the real part of the input fftw_complex signal (pixmap) as a monopage striped uint16 TIFF image */
static int _save_page_striped_tiff_uint16_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip)
     TIFF *tiffp;
     fftw_complex *pixmap;
     int nx,ny;
     uint32 rowsperstrip;
{
  unsigned int row,rr,cc;
  uint32 nrow;
  tstrip_t strip;
  tsize_t i;
  uint16 *buf;
  double val;
  int nout=0;

  buf = (uint16*) _TIFFmalloc(TIFFStripSize(tiffp));
  if(!buf) { printf("Not enough memory\n"); return EXIT_FAILURE; }

  for (row = 0; row<ny; row+=rowsperstrip){
    nrow = (row + rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp,row,0);
    i = 0;
    for (rr = 0; rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc) {
	val = pixmap[(row+rr)*nx + cc][0];
	if(val > UINT16_MAX) { val = UINT16_MAX; nout++; }
	if(val < UINT16_MIN) { val = UINT16_MIN; nout++; }
	buf[i++] = (uint16) val;
      }
    if (TIFFWriteEncodedStrip(tiffp,strip,buf,i*sizeof(uint16))<0) return EXIT_FAILURE;
  }
  _TIFFfree(buf);

  TIFFWriteDirectory(tiffp);

  if(nout > 0) printf("WARNING: %d graylevels were out of [%d,%d]\n",nout,UINT16_MIN,UINT16_MAX);

  return EXIT_SUCCESS;
}


/* save the real part of the input fftw_complex signal (pixmap) as a monopage striped int16 TIFF image */
static int _save_page_striped_tiff_int16_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip)
     TIFF *tiffp;
     fftw_complex *pixmap;
     int nx,ny;
     uint32 rowsperstrip;
{
  unsigned int row,rr,cc;
  uint32 nrow;
  tstrip_t strip;
  tsize_t i;
  int16 *buf;
  double val;
  int nout=0;

  buf = (int16*) _TIFFmalloc(TIFFStripSize(tiffp));
  if(!buf) { printf("Not enough memory\n"); return EXIT_FAILURE; }

  for (row = 0; row<ny; row+=rowsperstrip){
    nrow = (row + rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp,row,0);
    i = 0;
    for (rr = 0; rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc) {
	val = pixmap[(row+rr)*nx + cc][0];
	if(val > INT16_MAX) { val = INT16_MAX; nout++; }
	if(val < INT16_MIN) { val = INT16_MIN; nout++; }
	buf[i++] = (int16) val;
      }
    if (TIFFWriteEncodedStrip(tiffp,strip,buf,i*sizeof(int16))<0) return EXIT_FAILURE;
  }
  _TIFFfree(buf);

  TIFFWriteDirectory(tiffp);

  if(nout > 0) printf("WARNING: %d graylevels were out of [%d,%d]\n",nout,INT16_MIN,INT16_MAX);

  return EXIT_SUCCESS;
}


/* save the real part of the input fftw_complex signal (pixmap) as a monopage striped uint32 TIFF image */
static int _save_page_striped_tiff_uint32_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip)
     TIFF *tiffp;
     fftw_complex *pixmap;
     int nx,ny;
     uint32 rowsperstrip;
{
  unsigned int row,rr,cc;
  uint32 nrow;
  tstrip_t strip;
  tsize_t i;
  uint32 *buf;
  double val;
  int nout=0;

  buf = (uint32*) _TIFFmalloc(TIFFStripSize(tiffp));
  if(!buf) { printf("Not enough memory\n"); return EXIT_FAILURE; }

  for (row = 0; row<ny; row+=rowsperstrip){
    nrow = (row + rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp,row,0);
    i = 0;
    for (rr = 0; rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc) {
	val = pixmap[(row+rr)*nx + cc][0];
	if(val > UINT32_MAX) { val = UINT32_MAX; nout++; }
	if(val < UINT32_MIN) { val = UINT32_MIN; nout++; }
	buf[i++] = (uint32) val;
      }
    if (TIFFWriteEncodedStrip(tiffp,strip,buf,i*sizeof(uint32))<0) return EXIT_FAILURE;
  }
  _TIFFfree(buf);

  TIFFWriteDirectory(tiffp);

  if(nout > 0) printf("WARNING: %d graylevels were out of [%d,%u]\n",nout,UINT32_MIN,UINT32_MAX);

  return EXIT_SUCCESS;
}


/* save the real part of the input fftw_complex signal (pixmap) as a monopage striped int32 TIFF image */
static int _save_page_striped_tiff_int32_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip)
     TIFF *tiffp;
     fftw_complex *pixmap;
     int nx,ny;
     uint32 rowsperstrip;
{
  unsigned int row,rr,cc;
  uint32 nrow;
  tstrip_t strip;
  tsize_t i;
  int32 *buf;
  double val;
  int nout=0;

  buf = (int32*) _TIFFmalloc(TIFFStripSize(tiffp));
  if(!buf) { printf("Not enough memory\n"); return EXIT_FAILURE; }

  for (row = 0; row<ny; row+=rowsperstrip){
    nrow = (row + rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp,row,0);
    i = 0;
    for (rr = 0; rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc) {
	val = pixmap[(row+rr)*nx + cc][0];
	if(val > INT32_MAX) { val = INT32_MAX; nout++; }
	if(val < INT32_MIN) { val = INT32_MIN; nout++; }
	buf[i++] = (int32) val;
      }
    if (TIFFWriteEncodedStrip(tiffp,strip,buf,i*sizeof(int32))<0) return EXIT_FAILURE;
  }
  _TIFFfree(buf);

  TIFFWriteDirectory(tiffp);

  if(nout > 0) printf("WARNING: %d graylevels were out of [%d,%d]\n",nout,INT32_MIN,INT32_MAX);

  return EXIT_SUCCESS;
}


/* save the real part of the input fftw_complex signal (pixmap) as a monopage striped float TIFF image */
static int _save_page_striped_tiff_float_from_fftw_complex(tiffp,pixmap,nx,ny,rowsperstrip)
     TIFF *tiffp;
     fftw_complex *pixmap;
     int nx,ny;
     uint32 rowsperstrip;
{
  unsigned int row,rr,cc;
  uint32 nrow;
  tstrip_t strip;
  tsize_t i;
  float *buf;
  double val;
  int nout=0;

  buf = (float*) _TIFFmalloc(TIFFStripSize(tiffp));
  if(!buf) { printf("Not enough memory\n"); return EXIT_FAILURE; }

  for (row = 0; row<ny; row+=rowsperstrip){
    nrow = (row + rowsperstrip>ny?ny-row:rowsperstrip);
    strip = TIFFComputeStrip(tiffp,row,0);
    i = 0;
    for (rr = 0; rr<nrow; ++rr)
      for (cc = 0; cc<nx; ++cc) {
	val = pixmap[(row+rr)*nx + cc][0];
	if(val > FLT_MAX) { val = FLT_MAX; nout++; }
	if(val < -FLT_MAX) { val = -FLT_MIN; nout++; }
	buf[i++] = (float) val;
      }
    if (TIFFWriteEncodedStrip(tiffp,strip,buf,i*sizeof(float))<0) return EXIT_FAILURE;
  }
  _TIFFfree(buf);

  TIFFWriteDirectory(tiffp);

  if(nout > 0) printf("WARNING: %d graylevels were out of [%g,%g]\n",nout,-FLT_MAX,FLT_MAX);

  return EXIT_SUCCESS;
}
