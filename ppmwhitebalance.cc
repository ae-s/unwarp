// ppmwhitebalance.cc
// Copyright (C) 2013 Michael Rose

// ============================================================================
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ============================================================================

// This program can be used to correct the white balance of an image
// with a calibrating photo made from a gray card, based on a calculated
// grid of correction values.
// The primary goal is to normalize the colors of pictures from book pages
// or other material using a digital camera instead of a scanner.
// The setup for several images using the same calibration shot should
// remain as fixed as possible (lights, shadows, camera settings).
// All pictures read and written by this program must be in raw PPM format.
// The program ist testet with Linux but should also run on other platforms
// with only minor modifications. No special libraries are used.
// To compile use something like:
//
// $ g++ -o ppmwhitebalance ppmwhitebalance.cc
//
// A brief help text is printed, when wrong command syntax is used or simply
// by typing:
//
// $ ppmwhitebalance -h
//
// The user uses a gray card fitting the whole area of the image, e.g.
// the gray card should cover the picture area of interest.
// The user puts this gray card at the same position as the pages to be
// color corrected and takes a digital photo.
// By using conversion tools like 'jpegtopnm' she creates a file,
// e.g. 'calibration.ppm'. From this file, the command
//
// $ ppmwhitebalance calibration.ppm > calibration.bin
//
// extracts the information needed to normalize the color of the book sheets.
// To adjust the detection process, the user can use the following options:
//
// -gb        Correct only brightness values not individual RGB-channels.
// -g1        Maximal interpolation order is linear.
// -gc  80    Changes the desired RGB gray value of the gray card to 80.
//            Take mean value of calibration picture if zero.
// -gm 4.0    Maximal Multiplicator for the RGB channels.
// -gd 4.0    Maximal Divisor       for the RGB channels.
// -nx 10     Number of grid points in x-direction (default 1 if nx=ny=0).
// -ny 10     Number of grid points in y-direction (default 1 if nx=ny=0).
//
// If one of the last two options is zero, the value is choosen
// approximately by the aspect ratio of the image.
//
// If 'page01.ppm' is one of the images of the camera to be color
// adjusted, this can be achieved by:
//
// $ ppmwhitebalance -d calibration.bin page01.ppm > page01_enhanced.ppm
//
// If a book is digitized page after page, and the environmental conditions
// change, a new calibration picture shot should be taken from time to time.
//
// If the border colors of the calibration image 'calibration.ppm' changes
// rapidly, the interpolation in the neighborhood can be disturbed.
// Either change the calibration image manually by a graphic editor,
// cut the problematic errors in all pictures before color calibration or
// avoid such a setup completely.

// ============================================================================

// Include files; no other libraries needed than 'libgpp':

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <new>

#define VERSION "1.0"

// ----------------------------------------------------------------------------

// Example to compile successfully under Visual C++ 6.0 with
// cl -GX ppmwhitebalance.cpp

#ifdef _WIN32
  typedef unsigned __int8   u_int8_t;
  typedef unsigned __int16  u_int16_t;
  typedef unsigned __int32  u_int32_t;
  #include <io.h>
  #include <fcntl.h>
  #define binmode(fh) _setmode(_fileno(fh),_O_BINARY)
#else
  #define binmode(fh)
#endif

// ----------------------------------------------------------------------------

// Global types:

typedef u_int8_t  pixel_t;  // Type for RGB components of a pixel.
typedef u_int16_t color_t;  // Type for (r,g,b)-components of color grid.

typedef unsigned long ulong;

// ----------------------------------------------------------------------------

// Global variables:

const  int  nGlbBuf = 4096;
static char glbBuf[nGlbBuf];  // Used for different purposes.

// ----------------------------------------------------------------------------

// Howto obtain gray values:

#define PTR_TO_GRAY(value,ptr) \
  value += 0.299 * double(*ptr++); \
  value += 0.587 * double(*ptr++); \
  value += 0.114 * double(*ptr++)

// ============================================================================

// Global parameters/options:

struct Parameter {

  // Program name:
  const char *prgName;

  // Suppression of normal messages:
  int quiet;

  // File names:
  const char *calibPicName;     // Input image with gray card.
  const char *calibTextName;    // Calibration factors in text format.
  const char *calibColorName;   // Final color grid.
  const char *sourcePicName;    // Skewed input image.
  const char *destPicName;      // Unwarped output image.

  // Desired gray RGB value of gray card:
  int grayCard;

  // Should adjustment only act on brightness:
  int colorBrightness;

  // Maximal multiplicator for RGB values:
  double colorFactor;

  // Maximal divisor for RGB values:
  double colorDivisor;

  // Maximal linear interpolation:
  int gridLinear;

  // Number of grid points in x-direction:
  int gridNX;

  // Number of grid points in y-direction:
  int gridNY;

  Parameter();
  void Usage(int exitCode);
  void operator()(int argc, char *argv[]);
};

static Parameter param;  // Program parameters are global for ease of use.

// ============================================================================
// General common functions:
// ----------------------------------------------------------------------------

// Error message with program termination:

static void Error(const char *format, ...) {
  va_list ap;
  va_start(ap, format);
  vsprintf(glbBuf, format, ap);  // Hopefully this fits (long filenames)!
  va_end(ap);
  throw glbBuf;
}

// Normal messages, which can be suppressed:

static void Print(const char *format, ...) {
  if (!param.quiet) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
} }

// ============================================================================

// Onedimensional constant/linear/cubic interpolation:
//
// 'f[4]' contains up to four interpolation coefficients, depending on
// the fraction between the aequidistant grid points. 'm' is the maximum
// possible grid index (between 0 and 'm'), 'o' the row/column offset
// to adjust the grid pointer for adjacent rows/columns. 'r = 1' indicates
// that the weighting process should start one step before the current
// row/column and 'n' denotes the valid coefficients 'f[0 ... n]'.
// 'o2' is another offset suitable to create system matrices.
// The cubic interpolation is valid for 'm >= 2'.
// The interpolation changes automatically to linear order, if 'm < 2'
// and to constant order if 'm == 0'.

struct Cubic1 {
  double f[4];
  long   order, m, o, o2, n, r;
  void SetSize(long m) {
    this->m = m;
    order   = m < 2 ? m : (param.gridLinear ? 1 : 3); }
  Cubic1(long m, long o) { SetSize(m);  this->o = o; }
  Cubic1(long m, long o, long o2) { SetSize(m); this->o = o; this->o2 = o2; }
  void operator()(long& i, double t1);
};

// ----------------------------------------------------------------------------

// Initializes interpolation at location 'i + t1', '0 <= t1 <= 1'.
// The interpolation is C1-continuous if 'order == 3':

void Cubic1::operator()(long& i, double t1) {
  double t0;
  switch (order) {
    case 0:   f[n = r = 0] = 1.0;                   break;
    case 1:   f[r = 0] = 1.0 - t1;  f[n = 1] = t1;  break;
    default:
      t0 = 1.0 - t1;
      if (i) {
        if (i == m - 1) {  // Highest grid interval:
          f[    0] = -0.5 * t0 * t1;
          f[    1] =       t0 * (2.0 - t0);
          f[n = 2] = 0.5 * t1 * (2.0 - t0);
        } else {  // General grid interval:
          f[    0] = -0.5 * t0 * t0 * t1;
          f[    1] = -t0 * (((1.5 * t1) - 1.0) * t1 - 1.0);
          f[    2] = -t1 * (((1.5 * t0) - 1.0) * t0 - 1.0);
          f[n = 3] = -0.5 * t1 * t1 * t0;
        }
        r = 1;
      } else {  // Lowest grid interval:
        f[r = 0] = 0.5 * t0 * (2.0 - t1);
        f[    1] =       t1 * (2.0 - t1);
        f[n = 2] = -0.5 * t0 * t1; } }
}

// ----------------------------------------------------------------------------

// Twodimensional constant/linear/cubic interpolation:

template<class ityp>
struct Cubic2 {
  Cubic1 cx,cy;  // Onedimensional interpolation along x-axis and y-axis.
  Cubic2(long mx, long ox, long my, long oy) : cx(mx,ox), cy(my,oy) { }
  Cubic2(long mx, long ox, long ox2, long my, long oy, long oy2) :
    cx(mx,ox,ox2), cy(my,oy,oy2) { }
  double operator()(ityp *raw);
  void AssembleRhs(double scale, double *rhs);
  void AssembleMat(double *mat);
};

// ----------------------------------------------------------------------------

// Returns interpolated value for previously specified grid point
// in the cell, whose lower left pointer is given by 'raw':
//
// The interpolation on a rectangular grid has
// the tensor product property.

template<class ityp>
double Cubic2<ityp>::operator()(ityp *raw) {
  double  value = 0.0;
  long    ix,iy;
  ityp   *r;
  if (cx.r)  raw -= cx.o;
  if (cy.r)  raw -= cy.o;
  for (ix=0; ix <= cx.n; ++ix, raw += cx.o) {
    for (iy=0, r=raw; iy <= cy.n; ++iy, r += cy.o) {
      value += double(*r) * (cx.f[ix] * cy.f[iy]); } }
  return value;
}

// ----------------------------------------------------------------------------

// Assembles right hand side vector "G' * rhs":

template <class ityp>
void Cubic2<ityp>::AssembleRhs(double scale, double *rhs) {
  long   ix, iy;
  double *p, v;
  if (cx.r) rhs -= cx.o;
  if (cy.r) rhs -= cy.o;
  for (ix=0; ix <= cx.n; ++ix, rhs += cx.o) {
    v = scale * cx.f[ix];
    for (iy=0, p=rhs; iy <= cy.n; ++iy, p += cy.o) {
      *p += v * cy.f[iy]; } }
}

// ----------------------------------------------------------------------------

// Assembles full or band system matrix "G' * G":

template <class ityp>
void Cubic2<ityp>::AssembleMat(double *mat) {
  long   ix, iy, ix2, iy2;
  double *p1, *p2, *p3, v1, v2, v3;
  if (cx.r) mat -= cx.o + cx.o2;
  if (cy.r) mat -= cy.o + cy.o2;
  for (ix=0; ix <= cx.n; ++ix, mat += cx.o) {
    v1 = cx.f[ix];
    for (iy=0, p1=mat; iy <= cy.n; ++iy, p1 += cy.o) {
      v2   = v1 * cy.f[iy];
      for (ix2=0, p2=p1; ix2 <= cx.n; ++ix2, p2 += cx.o2) {
        v3 = v2 * cx.f[ix2];
        for (iy2=0, p3=p2; iy2 <= cy.n; ++iy2, p3 += cy.o2) {
          *p3 += v3 * cy.f[iy2]; } } } }
}

// ============================================================================

// Type for image data:

struct Picture {
  long     width,   // Image width
           height,  // Image height
           size;    // 3 * width * height, each pixel has RGB-tripel.
  pixel_t *pixel;   // 3 * width * height RGB bytes.

  Picture()  { pixel = 0;  Reset(); }
  ~Picture() {             Reset(); }

  void    Reset(long width=0, long height=0);
  void    Read(const char *filename, FILE *fh=0);
  void    Write(const char *filename, FILE *fh=0);
  pixel_t GetMeanValue();
};

// ----------------------------------------------------------------------------

// Resizes the image array by given size parameters:

void Picture::Reset(long width, long height) {
  if (pixel)  delete pixel;
  this->width  = width;
  this->height = height;
  size         = 3 * width * height;
  pixel        = 0;
  if (size)  pixel = new pixel_t[size];
}

// ----------------------------------------------------------------------------

// Read raw PPM image file:

void Picture::Read(const char *filename, FILE *fh) {
  int mode = 0;
  if (fh != stdin)  fh = fopen(filename, "rb");
  else              binmode(fh);  // Windows.
  if (!fh)  Error("Couldn't read '%s' as raw PPM image file", filename);
  try {
    while (mode < 3) {
      if (!fgets(glbBuf,nGlbBuf,fh))  throw 1;
      if (mode) {
        if (glbBuf[0] != '#') {
          if (mode == 1) {
            if (sscanf(glbBuf, "%ld %ld\n", &width, &height) != 2 ||
                width <= 0 || height <= 0 ||
                width > 25000 || height > 25000)  throw 1;
            mode = 2;
          } else {
            if (strcmp(glbBuf, "255\n"))  throw 1;
            mode = 3;
      } } }
      else if (strcmp(glbBuf,"P6\n"))  throw 1;  else  mode = 1; }
    Reset(width, height);
    if (fread(pixel, sizeof(pixel_t), size, fh) != ulong(size))  throw 2;
    if (fgetc(fh) != EOF)  throw 3;
    if (fh != stdin) fclose(fh); }
  catch (int nmr) {
    const char *msg;
    if (fh != stdin)  fclose(fh);
    switch (nmr) {
      case 1:   msg = "has wrong preamble";  break;
      case 2:   msg = "is too small";        break;
      default:  msg = "is too big"; }
    Error("PPM image file '%s' %s", filename, msg); }
}

// ----------------------------------------------------------------------------

// Write raw PPM image file:

void Picture::Write(const char *filename, FILE *fh) {
  if (fh != stdout)  fh = fopen(filename, "wb");
  else               binmode(fh);  // Windows.
  if (!fh)  Error("Couldn't write '%s' as raw PPM file", filename);
  fprintf(fh,"P6\n# CREATOR: %s\n%ld %ld\n255\n",param.prgName,width,height);
  fwrite(pixel, sizeof(pixel_t), size, fh);
  if (fh != stdout)  fclose(fh);
}

// ----------------------------------------------------------------------------

pixel_t Picture::GetMeanValue() {
  pixel_t *ptr   = pixel;
  long     nn    = width * height, n = nn;
  double   value = double(n >> 1);
  while (n--) { PTR_TO_GRAY(value,ptr); }
  value /= double(nn);
  if      (value < 16.0)   value = 16.0;
  else if (value > 240.0)  value = 240.0;
  return pixel_t(value);
}

// ============================================================================

// System modul to solve linear equation with band matrix.
// Matrix is of size 'n x n' with lower and upper bandwidth 'b'.
// The offset between valid entries on adjacent columns on the
// same row is 'o=2*b', the adjacent diagonal elements have offset 'o + 1':

struct System {
  long    n, b, o;
  double *sol, *rhs, *mat;

  System()  { sol = 0;  Reset(); }
  ~System() {           Reset(); }

  double *GetBaseMat(long i) { return mat + (o + 1) * i; }
  double *GetBaseRhs(long i) { return rhs + i;           }

  void Reset(long n=0, long b=0);
  void ClearRhs();
  void ClearMat();
  void Cholesky();
  void Solve();
};

// ----------------------------------------------------------------------------

void System::Reset(long n, long b) {
  long size;
  if (sol)  delete sol;
  this->n = n;
  this->b = b;
  o       = b << 1;
  size    = (o + 3) * n;
  sol     = rhs = mat = 0;
  if (size) {
    sol = new double[size];
    rhs = sol + n;
    mat = rhs + n;
    for (long i=0; i < size; ++i)  sol[i] = 0.0; }
}

// ----------------------------------------------------------------------------

void System::ClearRhs() {
  long    s = n;
  double *p = rhs;
  while (s--)  *p++ = 0.0; }

// ----------------------------------------------------------------------------

void System::ClearMat() {
  long    s = (o + 1) * n;
  double *p = mat;
  while (s--)  *p++ = 0.0; }

// ----------------------------------------------------------------------------

// Cholesky "C' * C = M" decomposition of band matrix (band is preserved):

void System::Cholesky() {
  double *pi, *pj, *pji, *pkj, *pki, s;
  long    i, j, k, m, mo;
  for (i=0, mo=-o, pi=mat; i < n; ++i, pi += o + 1) {
    if (i <= b) { m = i;  mo += o; } else m = b;
    for (j = 0, pji=pi-m, pj=pji-mo; j <= m; ++j, ++pji, pj += o+1) {
      for (k=0, pkj=pj, pki=pji, s=*pji; k < j; ++k)  s -= *--pkj * *--pki;
      if (j == m) {
        if (s < 1e-14)  Error("System matrix not positive definite");
        *pji = sqrt(s); }
      else  *pji = s / *pj; } }
}

// ----------------------------------------------------------------------------

// Solve "C' * C * sol = rhs":

void System::Solve() {
  double *ps, *pi, *pji, *pjs, s;
  long    i, j, m;
  for (i=0; i < n; ++i) sol[i] = rhs[i];
  // Solve for transposed cholesky factor:
  for (i=0, pi=mat, ps=sol; i < n; ++i, pi += o + 1) {
    m = i <= b ? i : b;
    for (j=0, pji=pi, pjs=ps, s=*ps; j < m; ++j)  s -= *--pji * *--pjs;
    *ps++ = s / *pi; }
  // Solve for cholesky factor:
  for (i=0, pi=mat+(n-1)*(o+1), ps=sol+(n-1); i < n; ++i, pi -= o + 1) {
    m = i <= b ? i : b;
    for (j=0, pji=pi+o, pjs=ps, s=*ps; j < m; ++j, pji += o)
      s -= *pji * *++pjs;
    *ps-- = s / *pi; }
}

// ============================================================================

// Conversion modul:

struct Converter {
  Picture *pic;
  double  *map;
  System   sys;
  long     pn1, po1, pn2, po2;
  long     gn1, gn2, go2, gb;
  int      swap;

  Converter()  { map = 0;  Reset(); }
  ~Converter() {           Reset(); }

  void Reset(Picture* pic=0, long nx=0, long ny=0);
  void SetMap(int channel = -1);
  void BuildtAndSolveSystem(int doMat=0);
  void ReadSolution(color_t *rgbMap);
  void GetRgbMap(color_t *rgbMap);

};

// ----------------------------------------------------------------------------

void Converter::Reset(Picture *pic, long nx, long ny) {
  long size;
  if (map)  delete map;
  this->pic = pic;
  map       = 0;
  pn1       = po1 = pn2 = po2 = gn1 = gn2 = go2 = gb = 0;
  swap      = 0;
  if (pic) {
    pn1 = pic->width;
    po1 = 1;
    pn2 = pic->height;
    po2 = pn1;
    gn1 = nx;
    gn2 = ny;
    go2 = nx;
    if (gn1 > gn2) {
      po1 = pn1;  pn1 = pn2;  pn2 = po1; po1 = po2;  po2 = 1;
      go2 = gn1;  gn1 = gn2;  gn2 = go2; go2 = gn1;
      swap = 1; }
    gb   = (param.gridLinear ? 1 : 3) * (go2 + 1);
    size = pn1 * pn2;
    map  = new double[size];
    for (long i=0; i < size; ++i)  map[i] = 0.0;
    sys.Reset(gn1 * gn2, gb); }
  else  sys.Reset();
}

// ----------------------------------------------------------------------------

// Creates the map for the RGB-channel (0,1,2) or for gray values (-1):

void Converter::SetMap(int channel) {
  pixel_t *pixel  = pic->pixel;
  double  *ptr    = map;
  long     n      = pn1 * pn2;
  double   orig;
  while (n--) {
    if (channel < 0) { orig = 0; PTR_TO_GRAY(orig,pixel);          }
    else             { orig = double(pixel[channel]);  pixel += 3; }
    *ptr++ = orig / 255.0; }
}

// ----------------------------------------------------------------------------

// Creates the band matrix and right hand side for the selected channel:

void Converter::BuildtAndSolveSystem(int doMat) {
  long    gm1 = gn1 - 1, gmm1 = gm1 ? gm1 : 1;
  long    gm2 = gn2 - 1, gmm2 = gm2 ? gm2 : 1;
  double  scale1 = pn1 > 1 ? double(gm1) / double(pn1 - 1) : 1.0;
  double  scale2 = pn2 > 1 ? double(gm2) / double(pn2 - 1) : 1.0;
  double *src1, *src2, ti1, ti2, t1, t2;
  long    i1, i2, k1, k2, k;
  Cubic2<double>  cub(gm1, 1, 2*gb, gm2, go2, go2*(2*gb));
  if (doMat)  sys.ClearMat();
  else        sys.ClearRhs();
  for (i1=0, src1=map; i1 < pn1; ++i1, src1 += po1) {
    t1  = scale1 * double(i1);
    k1  = long(ti1 = floor(t1));
    t1 -= ti1;
    if (k1 >= gmm1) { k1 = gmm1 - 1;  t1 = 1.0; }
    cub.cx(k1, t1);
    for (i2=0, src2=src1; i2 < pn2; ++i2, src2 += po2) {
      t2  = scale2 * double(i2);
      k2  = long(ti2 = floor(t2));
      t2 -= ti2;
      if (k2 >= gmm2) { k2 = gmm2 - 1;  t2 = 1.0; }
      cub.cy(k2, t2);
      k = k2 * go2 + k1;
      if (doMat)  cub.AssembleMat(sys.GetBaseMat(k));
      else        cub.AssembleRhs(*src2, sys.GetBaseRhs(k)); } }
  if (doMat) {
    Print("Cholesky factorization ...\n");
    sys.Cholesky(); }
  else sys.Solve();
}

// ----------------------------------------------------------------------------

// Solve the band matrix system:

void Converter::ReadSolution(color_t *rgbMap) {
  int     brightness = param.colorBrightness;
  long    gnx, gox, gny, goy, ix, iy;
  double *px, *py, z;
  if (swap) { gnx = gn2;  gox = go2;  gny = gn1;  goy = 1;   }
  else      { gnx = gn1;  gox = 1;    gny = gn2;  goy = go2; }
  for (iy=0, py=sys.sol; iy < gny; ++iy, py += goy) {
    for (ix=0, px=py; ix < gnx; ++ix, px += gox) {
      z = *px;
      z = z < 0.0 ? 0.0 : z > 1.0 ? 1.0 : z;
      z = 65535.0 * z + 0.5;
      *rgbMap = color_t(z);
      if (brightness)  rgbMap[2] = rgbMap[1] = rgbMap[0];
      rgbMap += 3; } }
}

// ----------------------------------------------------------------------------

// Combines all actions to obtain color grid:

void Converter::GetRgbMap(color_t *rgbMap) {
  int channel;
  Print("Build system matrix ...\n");
  BuildtAndSolveSystem(1);
  if (param.colorBrightness) {
    // Set factors for brightness only:
    Print("Calibrate brightness ...\n");
    SetMap();
    BuildtAndSolveSystem();
    ReadSolution(rgbMap); }
  else {
    for (channel = 0; channel < 3; ++channel) {
      // Set factors for each pixel of selected channel:
      Print("Calibrate channel %d ...\n", channel);
      SetMap(channel);
      BuildtAndSolveSystem();
      ReadSolution(rgbMap++); } }
}

// ============================================================================

// Type for grid of color points:

struct Grid {
  long       nx,      // Point indices in x-direction: 0 ... nx-1
             ny,      // Point indices in y-direction: 0 ... ny-1
             size;    // 3 * nx * ny
  int        linear;  // Maximum interpolation order is linear
  pixel_t    target;  // Target RGB-value for gray map.
  u_int32_t  smin;    // Minimal scaling factor.
  u_int32_t  smax;    // Maximal scaling factor.
  color_t   *rgbMap;  // normalized RGB factors.

  Grid()  { rgbMap = 0;  Reset(); }
  ~Grid() {              Reset(); }

  void Reset(long nx=0, long ny=0);
  void TextRead(const char *filename, FILE *fh=0);
  void TextWrite(const char *filename, FILE *fh=0);
  void Read(const char *filename, FILE *fh=0);
  void Write(const char *filename, FILE *fh=0);
  void Calibrate(Picture& pic);
  void Convert(Picture& src, Picture& dst);
};

// ----------------------------------------------------------------------------

void Grid::Reset(long nx, long ny) {
  size = 3 * nx * ny;
  if (nx > 1000 || ny > 1000 || size > 30000)
    Error("Color grid size %ldx%ld is too high", nx, ny);
  if (rgbMap)  delete rgbMap;
  this->nx = nx;
  this->ny = ny;
  rgbMap   = 0;
  linear   = 0;
  target   = 0;
  smin     = smax = 0.0;
  if (size) {
    rgbMap = new color_t[size];
    for (long i=0; i < size; ++i)  rgbMap[i] = 0; }
}

// ----------------------------------------------------------------------------

// Reads the color calibration grid from text file:

void Grid::TextRead(const char *filename, FILE *fh) {
  int      z, zz, r, g, b;
  long     k, n, ix, iy;
  color_t *rgb;
  if (fh != stdin)  fh = fopen(filename, "r");
  if (!fh)  Error("Couldn't read color calibration file '%s'", filename);
  try {
    if (fscanf(fh, "Colorcalibrationgrid %ld %ld\n", &ix, &iy) != 2 ||
        ix <= 0 || iy <= 0 || ix > 1000 || iy > 1000 ||
        ix * iy > 10000)  throw 0;
    Reset(ix, iy);
    if (fscanf(fh, "Onlylinear %d\n", &linear) != 1 ||
        linear < 0 || linear > 1)  throw 0;
    if (fscanf(fh, "Target %d\n", &z) != 1 || z < 16 || z > 240)  throw 0;
    target = pixel_t(z);
    if (fscanf(fh, "Factor %x %x ", &z, &zz) != 2 ||
        z < 4096 || z > 1048576 || zz < 4096 || zz > 1048576)  throw 0;
    smin = u_int32_t(z);
    smax = u_int32_t(zz);
    while (!feof(fh))  if (fgetc(fh) == '\n')  break;
    rgb = rgbMap;
    n   = nx * ny;
    for (k=0; k < n; ++k) {
      if (fscanf(fh, "%ld\t%ld\t%x\t%x\t%x ", &ix, &iy, &r, &g, &b) != 5 ||
          ix < 0 || ix >= nx || iy < 0 || iy >= ny ||
          r < 0 || r > 65535 ||
          g < 0 || g > 65535 ||
          b < 0 || b > 65535)  throw 1;
      *rgb++ = color_t(r);
      *rgb++ = color_t(g);
      *rgb++ = color_t(b);
      while (!feof(fh))  if (fgetc(fh) == '\n')  break; }
    if (fh != stdin)  fclose(fh); }
  catch (int nmr) {
    const char *msg = nmr ? "has wrong point syntax" : "has wrong preamble";
    if (fh != stdin)  fclose(fh);
    Error("Calibration color file '%s' %s", filename, msg); }
  // Set global parameter for interpolation restriction:
  param.gridLinear = linear;
}

// ----------------------------------------------------------------------------

// Writes the color calibration grid to text file:

void Grid::TextWrite(const char *filename, FILE *fh) {
  color_t *rgb = rgbMap;
  if (fh != stdout)  fh = fopen(filename, "w");
  if (!fh)  Error("Couldn't write color calibration file '%s'", filename);
  fprintf(fh, "Colorcalibrationgrid %ld %ld\n", nx, ny);
  fprintf(fh, "Onlylinear %d\n", linear);
  fprintf(fh, "Target %d\n", target);
  fprintf(fh, "Factor %06x %06x   (%2.2lf  %2.2lf)\n",
          smin, smax, double(smin) / 65536.0, double(smax) / 65536.0);
  for (long iy=0; iy < ny; ++iy) {
    for (long ix=0; ix < nx; ++ix, rgb += 3) {
      fprintf(fh,"%ld\t%ld\t%04x\t%04x\t%04x   (%2.2lf  %2.2lf  %2.2lf)\n",
              ix, iy, rgb[0], rgb[1], rgb[2],
              double(rgb[0]) / 65535.0,
              double(rgb[1]) / 65535.0,
              double(rgb[2]) / 65535.0); } }
  if (fh != stdout)  fclose(fh);
}

// ----------------------------------------------------------------------------

// Reads the color calibration grid from binary file:

void Grid::Read(const char *filename, FILE *fh) {
  u_int32_t magic;
  if (fh != stdin)  fh = fopen(filename, "rb");
  else              binmode(fh);  // Windows.
  if (!fh ||
      fread(&magic, sizeof(u_int32_t), 1, fh) != 1 ||
      magic != 0x7f4418bd ||
      fread(&nx, sizeof(long), 1, fh) != 1 ||
      nx < 1 || nx > 1000 ||
      fread(&ny, sizeof(long), 1, fh) != 1 ||
      ny < 1 || ny > 1000 ||
      fread(&linear, sizeof(int), 1, fh) != 1 ||
      linear < 0 || linear > 1 ||
      fread(&target, sizeof(pixel_t), 1, fh) != 1 ||
      target < 16 || target > 240 ||
      fread(&smin, sizeof(u_int32_t), 1, fh) != 1 ||
      smin < 4096 || smin > 1048576 ||
      fread(&smax, sizeof(u_int32_t), 1, fh) != 1 ||
      smax < 4096 || smax > 1048576 ||
      smin > smax)  magic = 0;
  if (magic) {
    size = 3 * nx * ny;
    if (size <= 30000) {
      rgbMap = new color_t[size];
      if (fread(rgbMap, sizeof(color_t), size, fh) != ulong(size))
        magic = 0; }
    else magic = 0; }
  if (fh && fh != stdin)  fclose(fh);
  if (!magic) Error("Couldn't read color grid file '%s'", filename);
  // Set global parameter for interpolation restriction:
  param.gridLinear = linear;
}

// ----------------------------------------------------------------------------

// Writes the color calibration grid to binary file:

void Grid::Write(const char *filename, FILE *fh) {
  u_int32_t magic = 0x7f4418bd;
  if (fh != stdout)  fh = fopen(filename, "wb");
  else               binmode(fh);  // Windows.
  if (!fh ||
      fwrite(&magic, sizeof(u_int32_t), 1, fh) != 1 ||
      fwrite(&nx, sizeof(long), 1, fh)         != 1 ||
      fwrite(&ny, sizeof(long), 1, fh)         != 1 ||
      fwrite(&linear, sizeof(int), 1, fh)      != 1 ||
      fwrite(&target, sizeof(pixel_t), 1, fh)  != 1 ||
      fwrite(&smin, sizeof(u_int32_t), 1, fh)  != 1 ||
      fwrite(&smax, sizeof(u_int32_t), 1, fh)  != 1 ||
      fwrite(rgbMap, sizeof(color_t), size, fh) != ulong(size))  magic = 0;
  if (fh && fh != stdout)  fclose(fh);
  if (!magic) Error("Couldn't write color grid file '%s'", filename);
}

// ----------------------------------------------------------------------------

// Create color calibration from gray card image:

void Grid::Calibrate(Picture& pic) {
  long nx = param.gridNX, ny = param.gridNY;
  Converter converter;
  if      (nx <= 0)  nx = (ny * pic.width + (pic.height >> 1)) / pic.height;
  else if (ny <= 0)  ny = (nx * pic.height + (pic.width >> 1)) / pic.width;
  if (nx <= 0)  nx = 1;
  if (ny <= 0)  ny = 1;
  if (nx > pic.width)   nx = pic.width;
  if (ny > pic.height)  ny = pic.height;
  Print("Calibration grid: %ld x %ld\n", nx, ny);
  if (!param.grayCard) {
    param.grayCard = int(pic.GetMeanValue());
    Print("Automatic graycard value: %d\n", param.grayCard); }
  Reset(nx, ny);
  converter.Reset(&pic, nx, ny);
  converter.GetRgbMap(rgbMap);
  // Save important parameters for later conversion:
  linear = param.gridLinear;
  target = pixel_t(param.grayCard);
  smin   = u_int32_t(65536.0 / param.colorDivisor + 0.5);
  smax   = u_int32_t(65536.0 * param.colorFactor  + 0.5);
}

// ----------------------------------------------------------------------------

// Applies the color enhancement to image:

void Grid::Convert(Picture& src, Picture& dst) {
  long     width = src.width,  height = src.height;
  long     mx = nx - 1,  mmx = mx ? mx : 1, ox = 3;
  long     my = ny - 1,  mmy = my ? my : 1, oy = 3 * nx;
  double   scalex = width  > 1 ? double(mx) / double(width  - 1) : 1.0;
  double   scaley = height > 1 ? double(my) / double(height - 1) : 1.0;
  double   gray   = double(target) /   255.0;
  double   div    = double(smin)   / 65536.0;
  double   fac    = double(smax)   / 65536.0;
  long     ix, iy, kx, ky;
  int      channel;
  double   tix, tiy, tx, ty, z;
  pixel_t *ps, *pd;
  color_t *rgb1, *rgb2;
  Cubic2<color_t> cub(mx, ox, my, oy);
  dst.Reset(width, height);
  for (iy=0, ps=src.pixel, pd=dst.pixel; iy < height; ++iy) {
    ty  = scaley * double(iy);
    ky  = long(tiy = floor(ty));
    ty -= tiy;
    if (ky >= mmy) { ky = mmy-1;  ty = 1.0; }
    cub.cy(ky, ty);
    rgb2 = rgbMap + (oy * ky);
    for (ix=0; ix < width; ++ix) {
      tx  = scalex * double(ix);
      kx  = long(tix = floor(tx));
      tx -= tix;
      if (kx >= mmx) { kx=mmx-1;  tx = 1.0; }
      cub.cx(kx, tx);
      rgb1 = rgb2 + (ox * kx);
      for (channel=0; channel < 3; ++channel) {
        z = cub(rgb1++) / 65535.0;
        z = fac * z <= gray ? fac : div * z >= gray ? div : gray / z;
        z = double(*ps++) * z;
        if      (z < 0.0)    z = 0.0;
        else if (z > 255.0)  z = 255.0;
        *pd++ = pixel_t(z); } } }
}

// ============================================================================

// Main routine:

void Main() {
  int     hasData = 0;
  Grid    grid;
  if (param.calibPicName) {
    Picture pic;
    if (*param.calibPicName)  pic.Read(param.calibPicName);
    else                      pic.Read("stdin", stdin);
    grid.Calibrate(pic);
    if (param.calibTextName)  grid.TextWrite(param.calibTextName);
    hasData = 1; }
  if (!hasData && param.calibTextName) {
    grid.TextRead(param.calibTextName);
    hasData = 1; }
  if (param.calibColorName || param.sourcePicName) {
    if (hasData) {
      if (param.calibColorName) {
        if (*param.calibColorName)  grid.Write(param.calibColorName);
        else                        grid.Write("stdout", stdout); } }
    else if (param.calibColorName && *param.calibColorName) {
      grid.Read(param.calibColorName);
      hasData = 1; }
    if (hasData && param.sourcePicName && param.destPicName) {
      Picture src, dst;
      if (*param.sourcePicName)  src.Read(param.sourcePicName);
      else                       src.Read("stdin", stdin);
      grid.Convert(src, dst);
      if (*param.destPicName)  dst.Write(param.destPicName);
      else                     dst.Write("stdout", stdout); } }
}

// ============================================================================

// Standard values for global parameters:

Parameter::Parameter() {
  prgName         = "ppmwhitebalance";
  quiet           = 0;
  calibPicName    = 0;
  calibTextName   = 0;
  calibColorName  = 0;
  sourcePicName   = 0;
  destPicName     = 0;
  grayCard        = 128;
  colorBrightness = 0;
  colorFactor     = 1.5;
  colorDivisor    = 1.5;
  gridLinear      = 0;
  gridNX          = 0;
  gridNY          = 0;
}

// ----------------------------------------------------------------------------

// Usage/help information:

void Parameter::Usage(int exitCode) {
  if (exitCode < 0) {
    fprintf(stdout,
            "%s " VERSION "\n\n"
            "Copyright (C) 2013 Michael Rose\n"
            "License GPLv3+: GNU GPL version 3 or later"
            " <http://gnu.org/licenses/gpl.html>\n"
            "This is free software: you are free to change"
            " and redistribute it.\n"
            "There is NO WARRANTY, to the extent permitted by law.\n\n",
            prgName);
    exit(0);
  }
  fprintf(stderr,
    "Usage: %s [options] [--] [inpname or stdin]\n\n"
    "Options:\n"
    "  --version     Print program version.\n"
    "  -h            Print program usage.\n"
    "  -q            Suppress normal program messages.\n"
    "  -c            Enforce calibration mode.\n"
    "  -cc (inpname) Set input PPM picture with gray card calibration image.\n"
    "  -cp           Set file name for textual calibration factors.\n"
    "  -d  (stdout)  Set file name for binary color grid.\n"
    "  -i  (inpname) Set input file name.\n"
    "  -o  (stdout)  Set ouput file name.\n"
    "  -gb           Only change brightness of color.\n"
    "  -g1           Maximum order of interpolation is linear.\n"
    "  -gc (128)     Set the desired RGB gray card value.\n"
    "                A zero value indicates picture mean gray value.\n"
    "  -gm (1.5)     Set maximal multiplicator for RGB values.\n"
    "  -gd (1.5)     Set maximal divisor       for RGB values.\n"
    "  -nx (0)       Number of grid points in x-direction (10 if nx=ny=0).\n"
    "  -ny (0)       Number of grid points in y-direction.\n"
    "                If nx or ny is zero, the image aspect ratio is used.\n\n"
    "Simple calibration:  ppwhitebalance calib.ppm > color.bin\n"
    "Simple usage:        ppwhitebalance -d color.bin src.ppm > dst.ppm\n",
    prgName);
  exit(exitCode);
}

// ----------------------------------------------------------------------------

// Sets global parameters/options from the command line:

void Parameter::operator()(int argc,char *argv[]) {
  const char *inpName = 0;
  const char **doName = 0;
  int doOpt = 1, doCalib = 0, *doNumber = 0, argi;
  double *doReal = 0;
  char   *end;
  for (argi = 1; argi < argc; ++argi) {
    char *arg = argv[argi];
    if (!*arg) Usage(1);
    if (doName) { *doName = arg;  doName = 0; }
    else if (doNumber) {
      long  z = strtol(arg, &end, 10);
      if ((doNumber == &grayCard && z && (z < 16 || z >  240)) ||
          (doNumber != &grayCard && (z < 1  || z > 1000)))  Usage(1);
      *doNumber = int(z);
      doNumber  = 0; }
    else if (doReal) {
      double z = strtod(arg, &end);
      if (z < 1.0 || z > 16.0)  Usage(1);
      *doReal = z;
      doReal  = 0; }
    else if (doOpt && *arg == '-') {
      if      (!strcmp(arg, "--version")) Usage(-1);
      else if (!strcmp(arg, "--"))        doOpt = 0;
      else if (!strcmp(arg, "-h"))        Usage(0);
      else if (!strcmp(arg, "-q"))        quiet           = 1;
      else if (!strcmp(arg, "-c"))        doCalib         = 1;
      else if (!strcmp(arg, "-cc"))       doName          = &calibPicName;
      else if (!strcmp(arg, "-cp"))       doName          = &calibTextName;
      else if (!strcmp(arg, "-d"))        doName          = &calibColorName;
      else if (!strcmp(arg, "-i"))        doName          = &sourcePicName;
      else if (!strcmp(arg, "-o"))        doName          = &destPicName;
      else if (!strcmp(arg, "-gb"))       colorBrightness = 1;
      else if (!strcmp(arg, "-g1"))       gridLinear      = 1;
      else if (!strcmp(arg, "-gc"))       doNumber        = &grayCard;
      else if (!strcmp(arg, "-gm"))       doReal          = &colorFactor;
      else if (!strcmp(arg, "-gd"))       doReal          = &colorDivisor;
      else if (!strcmp(arg, "-nx"))       doNumber        = &gridNX;
      else if (!strcmp(arg, "-ny"))       doNumber        = &gridNY;
      else  Usage(1);
    } else {
      if (argi+1 != argc || !*arg)  Usage(1);
      inpName = arg;
  } }
  if (doName || doNumber || doReal)  Usage(1);
  // Check options and set missing file names for calibration/convert mode:
  argi = 0;
  if (calibPicName)    argi++;
  if (calibTextName)   argi++;
  if (calibColorName)  argi++;
  if (argi != 1)  doCalib = 1;
  if (doCalib) {
    if (sourcePicName || destPicName)  Usage(1);
    if (calibPicName) { if (inpName)  Usage(1); }
    else if (inpName) { calibPicName = inpName; }
    if (!calibPicName && !calibTextName) calibPicName = "";
    if (!calibColorName && (!calibPicName || !calibTextName))
      calibColorName = "";
  } else {
    if (!sourcePicName)  sourcePicName = inpName ? inpName : "";
    else if (inpName)    Usage(1);
    if (!destPicName)    destPicName = "";
  }
  if (!gridNX && !gridNY)  gridNX = gridNY = 1;
}

// ============================================================================

// Main program to supply command line args, catch errors and
// call the main routine:

int main(int argc, char *argv[]) {
  try {
    try {
      param(argc, argv);
      Main();
    }
    // Catch exception from 'new' operator:
    catch (std::bad_alloc& ba) { throw "Out of memory"; }
  }
  catch (const char *msg) {
    fprintf(stderr, "!!! Error in %s:\n!!! %s!\n", param.prgName, msg);
    return 1;
  }
  return 0;
}
