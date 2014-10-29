#include <cstdio>
#include <cstdlib>
#include "pcdecoder.h"
#include "rcqsmodel.h"
#include "front.h"
#include "fpzip.h"
#include "codec.h"
#include "read.h"

#if FPZIP_FP == FPZIP_FP_FAST || FPZIP_FP == FPZIP_FP_SAFE
// decompress 3D array at specified precision using floating-point arithmetic
template <typename T, unsigned bits>
static void
decompress3d(
  RCdecoder* rd,   // entropy decoder
  T*         data, // flattened 3D array to decompress to
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz    // number of z samples
)
{
  // initialize decompressor
  typedef PCmap<T, bits> MAP;
  RCmodel* rm = new RCqsmodel(false, PCdecoder<T, MAP>::symbols);
  PCdecoder<T, MAP>* fd = new PCdecoder<T, MAP>(rd, &rm);
  FRONT<T> f(nx, ny);

  // decode difference between predicted (p) and actual (a) value
  unsigned x, y, z;
  for (z = 0, f.advance(0, 0, 1); z < nz; z++)
    for (y = 0, f.advance(0, 1, 0); y < ny; y++)
      for (x = 0, f.advance(1, 0, 0); x < nx; x++) {
        #if FPZIP_FP == FPZIP_FP_SAFE
        volatile T p = f(1, 1, 1);
        p += f(1, 0, 0);
        p -= f(0, 1, 1);
        p += f(0, 1, 0);
        p -= f(1, 0, 1);
        p += f(0, 0, 1);
        p -= f(1, 1, 0);
        #else
        T p = f(1, 0, 0) - f(0, 1, 1) +
              f(0, 1, 0) - f(1, 0, 1) +
              f(0, 0, 1) - f(1, 1, 0) +
              f(1, 1, 1);
        #endif
        T a = fd->decode(p);
        *data++ = a;
        f.push(a);
      }

  delete fd;
  delete rm;
}
#elif FPZIP_FP == FPZIP_FP_EMUL
#include "fpe.h"
// decompress 3D array at specified precision using floating-point emulation
template <typename T, unsigned bits>
static void
decompress3d(
  RCdecoder* rd,   // entropy encoder
  T*         data, // flattened 3D array to decompress to
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz    // number of z samples
)
{
  // initialize decompressor
  typedef PCmap<T, bits> MAP;
  typedef FPE<T> FLOAT;
  RCmodel* rm = new RCqsmodel(false, PCdecoder<T, MAP>::symbols);
  PCdecoder<T, MAP>* fd = new PCdecoder<T, MAP>(rd, &rm);
  FRONT<FLOAT> f(nx, ny);

  // decode difference between predicted (p) and actual (a) value
  unsigned x, y, z;
  for (z = 0, f.advance(0, 0, 1); z < nz; z++)
    for (y = 0, f.advance(0, 1, 0); y < ny; y++)
      for (x = 0, f.advance(1, 0, 0); x < nx; x++) {
        FLOAT p = f(1, 0, 0) - f(0, 1, 1) +
                  f(0, 1, 0) - f(1, 0, 1) +
                  f(0, 0, 1) - f(1, 1, 0) +
                  f(1, 1, 1);
        T a = fd->decode(T(p));
        *data++ = a;
        f.push(a);
      }
                                                                                
  delete fd;
  delete rm;
}
#else // FPZIP_FP_INT
// decompress 3D array at specified precision using integer arithmetic
template <typename T, unsigned bits>
static void
decompress3d(
  RCdecoder* rd,   // entropy decoder
  T*         data, // flattened 3D array to decompress to
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz    // number of z samples
)
{
  // initialize decompressor
  typedef PCmap<T, bits> TMAP;
  typedef typename TMAP::RANGE U;
  typedef PCmap<U, bits, U> UMAP;
  RCmodel* rm = new RCqsmodel(false, PCdecoder<U, UMAP>::symbols);
  PCdecoder<U, UMAP>* fd = new PCdecoder<U, UMAP>(rd, &rm);
  TMAP map;
  FRONT<U> f(nx, ny, map.forward(0));

  // decode difference between predicted (p) and actual (a) value
  unsigned x, y, z;
  for (z = 0, f.advance(0, 0, 1); z < nz; z++)
    for (y = 0, f.advance(0, 1, 0); y < ny; y++)
      for (x = 0, f.advance(1, 0, 0); x < nx; x++) {
        U p = f(1, 0, 0) - f(0, 1, 1) +
              f(0, 1, 0) - f(1, 0, 1) +
              f(0, 0, 1) - f(1, 1, 0) +
              f(1, 1, 1);
        U a = fd->decode(p);
        *data++ = map.inverse(a);
        f.push(a);
      }

  delete fd;
  delete rm;
}
#endif

// decompress p-bit float, 2p-bit double
#define decompress_case(p)\
  case subsize(T, p):\
    decompress3d<T, subsize(T, p)>(rd, data, nx, ny, nz);\
    break

// decompress 4D array
template <typename T>
static void
decompress4d(
  RCdecoder* rd,   // entropy decoder
  T*         data, // flattened 4D array to decompress to
  int*       prec, // per field precision
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz,   // number of z samples
  unsigned   nf    // number of fields
)
{
  // decompress one field at a time
  for (unsigned i = 0; i < nf; i++) {
    int bits = rd->template decode <unsigned>(32);
    if (prec)
      prec[i] = bits;
    switch (bits) {
      decompress_case( 2);
      decompress_case( 3);
      decompress_case( 4);
      decompress_case( 5);
      decompress_case( 6);
      decompress_case( 7);
      decompress_case( 8);
      decompress_case( 9);
      decompress_case(10);
      decompress_case(11);
      decompress_case(12);
      decompress_case(13);
      decompress_case(14);
      decompress_case(15);
      decompress_case(16);
      decompress_case(17);
      decompress_case(18);
      decompress_case(19);
      decompress_case(20);
      decompress_case(21);
      decompress_case(22);
      decompress_case(23);
      decompress_case(24);
      decompress_case(25);
      decompress_case(26);
      decompress_case(27);
      decompress_case(28);
      decompress_case(29);
      decompress_case(30);
      decompress_case(31);
      decompress_case(32);
      default:
        fprintf(stderr, "fpzip: invalid precision %d in file\n", bits);
        abort();
        break;
    }
    data += nx * ny * nz;
  }
}

void
read_header(
  RCdecoder* rd,
  unsigned   nx,
  unsigned   ny,
  unsigned   nz,
  unsigned   nf,
  int        dp
)
{
  // magic
  if (rd->decode<unsigned>(8) != 'f' ||
      rd->decode<unsigned>(8) != 'p' ||
      rd->decode<unsigned>(8) != 'z' ||
      rd->decode<unsigned>(8) != '\0') {
    fprintf(stderr, "fpzip: not an fpz stream\n");
    abort();
  }

  // format version
  if (rd->decode<unsigned>(16) != FPZ_MAJ_VERSION ||
      rd->decode<unsigned>(16) != FPZ_MIN_VERSION) {
    fprintf(stderr, "fpzip: format version not supported\n");
    abort();
  }

  // array dimensions
  if (rd->decode<unsigned>(32) != nf ||
      rd->decode<unsigned>(32) != nz ||
      rd->decode<unsigned>(32) != ny ||
      rd->decode<unsigned>(32) != nx) {
    fprintf(stderr, "fpzip: array dimensions do not match\n");
    abort();
  }

  // single or double precision
  if (rd->decode() != !!dp) {
    fprintf(stderr, "fpzip: floating-point type does not match\n");
    abort();
  }
}

static void
fpzip_stream_read(
  RCdecoder* rd,   // entropy decoder
  void*      data, // array to read
  int*       prec, // per field bits of precision
  int        dp,   // double precision array if nonzero
  unsigned   nx,   // number of x samples
  unsigned   ny,   // number of y samples
  unsigned   nz,   // number of z samples
  unsigned   nf    // number of fields
)
{
  rd->init();
  read_header(rd, nx, ny, nz, nf, dp);
  if (dp)
    decompress4d(rd, (double*)data, prec, nx, ny, nz, nf);
  else{
    decompress4d(rd, (float*)data, prec, nx, ny, nz, nf);
  }
}

// read and decompress a single or double precision 4D array from file
unsigned
fpzip_file_read(
  FILE*       file, // binary input stream
  void*       data, // array to read
  int*        prec, // per field bits of precision
  int         dp,   // double precision array if nonzero
  unsigned    nx,   // number of x samples
  unsigned    ny,   // number of y samples
  unsigned    nz,   // number of z samples
  unsigned    nf    // number of fields
)
{
  RCfiledecoder* rd = new RCfiledecoder(file);
  fpzip_stream_read(rd, data, prec, dp, nx, ny, nz, nf);
  unsigned bytes = rd->error ? 0 : rd->bytes();
  delete rd;
  return bytes;
}

// read and decompress a single or double precision 4D array from file
unsigned
fpzip_memory_read(
  const void* buffer, // pointer to compressed data
  void*       data,   // array to read
  int*        prec,   // per field bits of precision
  int         dp,     // double precision array if nonzero
  unsigned    nx,     // number of x samples
  unsigned    ny,     // number of y samples
  unsigned    nz,     // number of z samples
  unsigned    nf      // number of fields
)
{
  RCmemdecoder* rd = new RCmemdecoder(buffer);
  fpzip_stream_read(rd, data, prec, dp, nx, ny, nz, nf);
  unsigned bytes = rd->error ? 0 : rd->bytes();
  delete rd;
  return bytes;
}
unsigned
fpzip_memory_read_py(
  float* buffer, // pointer to compressed data
  float*       data,   // array to read
  int*        prec,   // per field bits of precision
  int         dp,     // double precision array if nonzero
  unsigned    nx,     // number of x samples
  unsigned    ny,     // number of y samples
  unsigned    nz,     // number of z samples
  unsigned    nf      // number of fields
)
{
  if (!fpzip_memory_read((const void *)buffer, (void *)data,prec,dp, nx, ny,nz, nf)) {
    fprintf(stderr, "fpzip: cannot read file from memory\n");
    abort();
  }
  return 1; 
}
void 
fpzip_memory_read_f(
  void* buffer,
  void*   data,
  int*    prec,
  const int*     dp,
  const int*     nx,
  const int*     ny,
  const int*     nz,
  const int*     nf
)
{
  if (!fpzip_memory_read(buffer, data,prec,*dp, *nx, *ny,*nz, *nf)) {
    fprintf(stderr, "fpzip: cannot read file from memory\n");
    abort();
  }
}
void
fpzip_memory_read_f_(
  void* buffer,
  void*   data,
  int*    prec,
  const int*     dp,
  const int*     nx,
  const int*     ny,
  const int*     nz,
  const int*     nf
)
{
  fpzip_memory_read_f(buffer,data,prec,dp,nx,ny,nz,nf);
}
// wrappers for fortran calls
void
fpzip_file_read_f(
  const char* path, // path to input file
  void*       data, // array to read
  int*        prec, // per field bits of precision
  const int*  dp,   // double precision array if nonzero
  const int*  nx,   // number of x samples
  const int*  ny,   // number of y samples
  const int*  nz,   // number of z samples
  const int*  nf    // number of fields
)
{
  FILE* file = fopen(path, "rb");
  if (!file) {
    fprintf(stderr, "fpzip: cannot open file %s\n", path);
    abort();
  }
  if (!fpzip_file_read(file, data, prec, *dp, *nx, *ny, *nz, *nf)) {
    fprintf(stderr, "fpzip: cannot read file %s\n", path);
    abort();
  }
  fclose(file);
}

void
fpzip_file_read_f_(
  const char* path, // path to input file
  void*       data, // array to read
  int*        prec, // per field bits of precision
  const int*  dp,   // double precision array if nonzero
  const int*  nx,   // number of x samples
  const int*  ny,   // number of y samples
  const int*  nz,   // number of z samples
  const int*  nf    // number of fields
)
{
  fpzip_file_read_f(path, data, prec, dp, nx, ny, nz, nf);
}
