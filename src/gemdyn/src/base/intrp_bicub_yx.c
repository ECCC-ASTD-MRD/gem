/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2019  Division de Recherche en Prevision Numerique
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#if ! defined(F_TEST)

#if defined(NEVER_TRUE)
 rm -f intrp_bicub_yx intrp_bicub_yx.o intrp_bicub_yx_f.F90
 ln -sf intrp_bicub_yx.c intrp_bicub_yx_f.F90

 to build (Intel compilers) :
 s.cc -O2 -DTIMING -c -march=core-avx2 intrp_bicub_yx.c
 s.f90 -no-wrap-margin -DF_TEST -o intrp_bicub_yx intrp_bicub_yx_f.F90 intrp_bicub_yx.o
 s.cc -O2 -DTIMING -c -march=core2 intrp_bicub_yx.c
 s.f90 -no-wrap-margin -DF_TEST -o intrp_bicub_yx0 intrp_bicub_yx_f.F90 intrp_bicub_yx.o

 to build (GNU compilers) :
 gcc -O3 -DTIMING -c -mavx2 -mfma intrp_bicub_yx.c
 gfortran -DF_TEST -o intrp_bicub_yx intrp_bicub_yx_f.F90 intrp_bicub_yx.o -ffree-line-length-none
 gcc -O3 -DTIMING -c -msse3 intrp_bicub_yx.c
 gfortran -DF_TEST -o intrp_bicub_yx0 intrp_bicub_yx_f.F90 intrp_bicub_yx.o -ffree-line-length-none

 ./intrp_bicub_yx
 ./intrp_bicub_yx0
#endif

static double cp167 =  0.166666666666666667E0;
static double cm167 = -0.166666666666666667E0;
static double cp5 =  .5;
static double cm5 = -.5;
static double one = 1.0;
static double two = 2.0;

#if defined(TIMING)
#include <stdio.h>
#include <stdint.h>
#pragma weak rdtsc_=rdtsc
uint64_t rdtsc_(void);
uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__)
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  static uint64_t time0;
  return time0++;
#endif
}
#endif

#if defined(__x86_64__)
#include <immintrin.h>
// use separate multiply and add instructions if fused multiply-add not available
#if defined(__AVX__)
#if ! defined(__FMA__)
#define _mm256_fmadd_ps(a,b,c) _mm256_add_ps(_mm256_mul_ps(a,b),c)
#define _mm256_fmadd_pd(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b),c)
#define _mm_fmadd_ps(a,b,c) _mm_add_ps(_mm_mul_ps(a,b),c)
#define _mm_fmadd_pd(a,b,c) _mm_add_pd(_mm_mul_pd(a,b),c)
#endif
#endif
#endif


// default values mean practically "no limit"
static int qminx = -999999999;  // limits for the 'q' grid, in "origin 0"
static int qmaxx =  999999999;
static int qminy = -999999999;
static int qmaxy =  999999999;

static int uminx = -999999999;  // limits for the 'u' grid
static int umaxx =  999999999;
static int uminy = -999999999;
static int umaxy =  999999999;

static int vminx = -999999999;  // limits for the 'v' grid
static int vmaxx =  999999999;
static int vminy = -999999999;
static int vmaxy =  999999999;

static int qoffx = 0;
static int qoffy = 0;
static int uoffx = 0;
static int uoffy = 0;
static int voffx = 0;
static int voffy = 0;

// set offsets for q grid, replicate q values for u and v grids
void set_intrp_bicub_off_yx(int qoffsetx, int qoffsety){
  qoffx = qoffsetx;
  qoffy = qoffsety;
  uoffx = qoffsetx;
  uoffy = qoffsety;
  voffx = qoffsetx;
  voffy = qoffsety;
}

// set offsets for q, u, and v grids
void set_intrp_bicub_off_quv(int qoffsetx, int qoffsety, int uoffsetx, int uoffsety, int voffsetx, int voffsety){
  qoffx = qoffsetx;
  qoffy = qoffsety;
  uoffx = uoffsetx;
  uoffy = uoffsety;
  voffx = voffsetx;
  voffy = voffsety;
}

// set bounds for 'q' grid only (scalars)
void set_intrp_bicub_yx(int qxmin, int qxmax, int qymin, int qymax){
  qminx = qxmin - 1 ;   // qxmin is in "origin 1", qminx is in "origin 0"
  qmaxx = qxmax - 1 ;   // bounds for scalars, see set_intrp_bicub_quv
  qminy = qymin - 1 ;
  qmaxy = qymax - 1 ;
}
// set bounds for 'q' (scalars), 'u', and 'v' grids
void set_intrp_bicub_quv(int qxmin, int qxmax, int qymin, int qymax,
                         int uxmin, int uxmax, int uymin, int uymax,
                         int vxmin, int vxmax, int vymin, int vymax){
  // qxmin is in "origin 1", qminx is in "origin 0", same for other limits
  // uxmin is in "origin 1", uminx is in "origin 0", same for other limits
  // vxmin is in "origin 1", vminx is in "origin 0", same for other limits
  // bounds for scalars
  qminx = qxmin - 1 ;  // min index along x
  qmaxx = qxmax - 1 ;  // max index along x
  qminy = qymin - 1 ;  // min index along y
  qmaxy = qymax - 1 ;  // max index along y
  // bounds for u
  uminx = uxmin - 1 ;  // min index along x
  umaxx = uxmax - 1 ;  // max index along x
  uminy = uymin - 1 ;  // min index along y
  umaxy = uymax - 1 ;  // max index along y
  // bounds for v
  vminx = vxmin - 1 ;  // min index along x
  vmaxx = vxmax - 1 ;  // max index along x
  vminy = vymin - 1 ;  // min index along y
  vmaxy = vymax - 1 ;  // max index along y
}

/*
           interpolate a column from a 3D source array f, put results in array r (ru, rv for winds)

   f       3D Fortran indexed source array, real*4, dimension(mini:maxi,minj:maxj,nk)
           ni = (maxi-mini-1)
           ninj = (maxi-mini-1)*(maxj-minj-1)
   ni      distance between f(i,j,k) and f(i,j+1,k) (row size)
   ninj    distance between f(i,j,k) and f(i,j,k+1) (plane size)
   r,ru,rv dimension(nk), result array(s)
   nk      number of 2D planes
   xx      i coordinate in i j fractional index space of desired column ( xx(l) is assumed >= mini+1 )
   xxu     i coordinate in i j fractional index space of desired column ( xxu(l) is assumed >= mini+1 )
   xxv     i coordinate in i j fractional index space of desired column ( xxv(l) is assumed >= mini+1 )
   yy      j coordinate in i j fractional index space of desired column ( yy(l) is assumed >= minj+1 )
   yyu     j coordinate in i j fractional index space of desired column ( yyu(l) is assumed >= minj+1 )
   yyv     j coordinate in i j fractional index space of desired column ( yyv(l) is assumed >= minj+1 )
   s       rotation matrix for the wind components (intrp_bicub_yx_vec, intrp_bicub_yx_uv)

   f is assumed to point to the address of f(1,1,1)  (Fortran indexing)
   xx, yy positions are in fractional indes space (origin 1)
   the value at (1.0, 1.0) is f(1,1,1)

   xx = 2.5, yy = 2.5 would be the center ot the square formed by
   f(2,2,k) f(3,2,k) f(2,3,k) f(3.3.k)  (where 1 <= k <= nk)

   x and y coordinates are assumed to be within bounds set via a call to set_intrp_bicub_yx or set_intrp_bicub_quv
   3 sets of bounds are used, one for scalars, one for u and one for v

   to call from FORTRAN, the following interface is used

   interface                                                                                      !InTf!
    subroutine set_intrp_bicub_off_yx(qx,qy) bind(C,name='set_intrp_bicub_off_yx')                !InTf!
      integer, intent(IN), value :: qx,qy                                                         !InTf!
    end subroutine set_intrp_bicub_off_yx                                                         !InTf!
    subroutine set_intrp_bicub_off_quv(qx,qy,ux,uy,vx,vy) bind(C,name='set_intrp_bicub_off_quv')  !InTf!
      integer, intent(IN), value :: qx,qy,ux,uy,vx,vy                                             !InTf!
    end subroutine set_intrp_bicub_off_quv                                                        !InTf!
    subroutine set_intrp_bicub_yx(x1,x2,y1,y2) bind(C,name='set_intrp_bicub_yx')                  !InTf!
      integer, intent(IN), value :: x1,x2,y1,y2                                                   !InTf!
    end subroutine set_intrp_bicub_yx                                                             !InTf!
    subroutine set_intrp_bicub_quv(x1,x2,y1,y2,u1,u2,u3,u4,v1,v2,v3,v4) bind(C,name='set_intrp_bicub_quv')   !InTf!
      integer, intent(IN), value :: x1,x2,y1,y2,u1,u2,u3,u4,v1,v2,v3,v4                           !InTf!
    end subroutine set_intrp_bicub_quv                                                            !InTf!
    subroutine intrp_bicub_yx_s(f, r, ni, ninj, nk, x, y) bind(C,name='intrp_bicub_yx_s')         !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_s                                                               !InTf!
    subroutine intrp_bicub_yx_s_mono(f, r, ni, ninj, nk, x, y) bind(C,name='intrp_bicub_yx_s_mono') !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_s_mono                                                          !InTf!
    subroutine intrp_bicub_yx_vec(u,v,ru,rv,ni,ninj,nk,xu,yu,xv,yv,s) bind(C,name='intrp_bicub_yx_vec')  !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: u, v                                             !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: ru, rv                                          !InTf!
      real(C_DOUBLE), intent(IN), value :: xu, yu, xv, yv                                         !InTf!
      real(C_DOUBLE), intent(IN), dimension(*) :: s                                               !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_vec                                                             !InTf!
    subroutine intrp_bicub_yx_uv(u,v,r1,r2,ni,ninj,nk,x,y,s) bind(C,name='intrp_bicub_yx_uv')     !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: u, v                                             !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r1, r2                                          !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      real(C_DOUBLE), intent(IN), dimension(*) :: s                                               !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk                                           !InTf!
    end subroutine intrp_bicub_yx_uv                                                              !InTf!
   end interface                                                                                  !InTf!

 */
// compute 2D X Y cubic interpolation coefficients (Origin 0 internally)
// assuming CONSTANT spacing along x, CONSTANT spacing along y
//
// xx    : position along x in index space of the interpolation target, input (origin 1)
// yy    : position along x in index space of the interpolation target, input (origin 1)
// wx    : output array, 4 elements, coefficients for cubic interpolation along x
// wy    : output array, 4 elements, coefficients for cubic interpolation along y
// qminx, qmaxx : limits along x used as box corner constraints
// qminy, qmaxy : limits along y used as box corner constraints
// offsetx  : offset applied to xx (global to local optional mapping) (usually 0)
// offsety  : offset applied to yy (global to local optional mapping) (usually 0)
//
// function result : displacement from point(1,1) in array to lower left corner of 4 x 4 subarray
//                   used for the bicubic interpolation
//
// inline version, used by multiple routines, the optimizer will perform the appropriate inlining
//
static inline int bicub_coeffs_inline(double xx, double yy, int ni, double *wx, double *wy,
                                  int qminx, int qmaxx, int qminy, int qmaxy, int offsetx, int offsety){
  double x, y;
  int ix, iy;
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;      // ix > xx if xx is negative
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;      // iy > yy if yy is negative
  ix = (ix < qminx) ? qminx : ix;        // apply boundary limits
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = xx - 2 - ix;                      // xx and yy are in "ORIGIN 1"
  y  = yy - 2 - iy;
  ix = ix + offsetx;
  iy = iy + offsety;

  wx[0] = cm167*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp167*x*(x+one)*(x-one);

  wy[0] = cm167*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp167*y*(y+one)*(y-one);

  return ix + iy * ni;               // point to lower left corner of 4 x 4 subarray
}

// version for scalar variables
void intrp_bicub_yx_s(float *f, float *r, int ni, int ninj, int nk, double xx, double yy){
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt ;
  __m128  frt ;
  __m128d ft0, ft1;
#else
  double fd0[4], fd1[4], fd2[4], fd3[4] ;
#endif
  double  wx[4], wy[4] ;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  double x, y;
  int ix, iy;
#if defined(NO_INLINE)
  x = xx - 1.0 ; y = yy - 1.0;                      // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  ix = (ix < qminx) ? qminx : ix;       // apply boundary limits (set by set_intrp_bicub_quv or set_intrp_bicub_yx)
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;                  // point to lower left corner of 4 x 4 subarray

  wx[0] = cm167*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp167*x*(x+one)*(x-one);

  wy[0] = cm167*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp167*y*(y+one)*(y-one);
#else
  f += bicub_coeffs_inline(xx, yy, ni, wx, wy, qminx, qmaxx, qminy, qmaxy, qoffx, qoffy);  // use inline function
#endif
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  fwx = _mm256_loadu_pd(wx) ;             // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;          // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  // fetch 4 rows, level 1
  fd0 = _mm256_cvtps_pd(_mm_loadu_ps(f    )) ;     // row 1 : f[i:i+3 , j   ,k]
  fd1 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni )) ;     // row 2 : f[i:i+3 , j+1 ,k]
  fd2 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni2)) ;     // row 3 : f[i:i+3 , j+2 ,k]
  fd3 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni3)) ;     // row 4 : f[i:i+3 , j+3 ,k]
  for(k=0 ; k<nk-1 ; k++){    // interpolation along J, level k, prefetch 4 rows for level k+1
    f+= ninj;                 // next plane

    fdt = _mm256_mul_pd(fd0,fwy0) ;                // row * coefficient
    fd0 = _mm256_cvtps_pd(_mm_loadu_ps(f    )) ;   // row 1 : f[i:i+3 , j   ,k]

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;          // sum of row * coefficient
    fd1 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni )) ;   // row 2 : f[i:i+3 , j+1 ,k]

    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;          // sum of row * coefficient
    fd2 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni2)) ;   // row 3 : f[i:i+3 , j+2 ,k]

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;          // sum of row * coefficient
    fd3 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni3)) ;   // row 4 : f[i:i+3 , j+3 ,k]

    // interpolation along I: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    _mm_store_ss(r,frt) ;                   // store float result(k)
    r ++;
  }
  // interpolation along J , level nk
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row * coefficient
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // interpolation along I: multiply by coefficients along x , then sum elements
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  _mm_store_ss(r,frt) ;                   // store float result(nk)
#else
  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    f+= ninj;
    r ++;
  }
#endif
}

// interpolate u and v with rotation, u and v are on different grids
// hence the 2 sets of positions (xxu,yyu) and (xxv,yyv)
// xxu, yyu, xxv, yyv in fractional index space, same as intrp_bicub_yx_s
// ru, rv : rotated U and V components of the wind at the desired point
// for rotation matrix s, see the GEM model
void intrp_bicub_yx_vec(float *u, float *v, float *ru, float *rv, int ni, int ninj, int nk,
                        double xxu, double yyu, double xxv, double yyv, double s[4]){
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  __m256d fwyu0, fwyu1, fwyu2, fwyu3, fwyv0, fwyv1, fwyv2, fwyv3, fdtu, fdtv, fwxu, fwxv;
  __m256d fd0, fd1, fd2, fd3;
  __m128d ft0, ft1;
  __m128 frt;
#else
  double ud0[4];
  double vd0[4];
#endif
  double  wxu[4], wyu[4], wxv[4], wyv[4] ;
  double  xu, yu, xv, yv, tu, tv;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ixu, iyu, ixv, iyv;
#if defined(NO_INLINE)
  xu = xxu - 1.0 ; yu = yyu - 1.0; // xxu and yyu are in "ORIGIN 1"
  ixu = xxu ; if(ixu > xxu) ixu = ixu -1 ; ixu = ixu - 2;
  iyu = yyu ; if(iyu > yyu) iyu = iyu -1 ; iyu = iyu - 2;
  ixu = (ixu < uminx) ? uminx : ixu;       // apply boundary limits (set by set_intrp_bicub_quv or set_intrp_bicub_yx)
  ixu = (ixu > umaxx) ? umaxx : ixu;
  iyu = (iyu < uminy) ? uminy : iyu;
  iyu = (iyu > umaxy) ? umaxy : iyu;
  xu  = xu - 1 - ixu;
  yu  = yu - 1 - iyu;
  u = u + ixu + iyu * ni;

  xv = xxv - 1.0 ; yv = yyv - 1.0; // xxv and yyv are in "ORIGIN 1"
  ixv = xxv ; if(ixv > xxv) ixv = ixv -1 ; ixv = ixv - 2;
  iyv = yyv ; if(iyv > yyv) iyv = iyv -1 ; iyv = iyv - 2;
  ixv = (ixv < vminx) ? vminx : ixv;       // apply boundary limits (set by set_intrp_bicub_quv or set_intrp_bicub_yx)
  ixv = (ixv > vmaxx) ? vmaxx : ixv;
  iyv = (iyv < vminy) ? vminy : iyv;
  iyv = (iyv > vmaxy) ? vmaxy : iyv;
  xv  = xv - 1 - ixv;
  yv  = yv - 1 - iyv;
  v = v + ixv + iyv * ni;

  wxu[0] = cm167*xu*(xu-one)*(xu-two);       // polynomial coefficients along i
  wxu[1] = cp5*(xu+one)*(xu-one)*(xu-two);
  wxu[2] = cm5*xu*(xu+one)*(xu-two);
  wxu[3] = cp167*xu*(xu+one)*(xu-one);

  wxv[0] = cm167*xv*(xv-one)*(xv-two);       // polynomial coefficients along i
  wxv[1] = cp5*(xv+one)*(xv-one)*(xv-two);
  wxv[2] = cm5*xv*(xv+one)*(xv-two);
  wxv[3] = cp167*xv*(xv+one)*(xv-one);

  wyu[0] = cm167*yu*(yu-one)*(yu-two);       // polynomial coefficients along j
  wyu[1] = cp5*(yu+one)*(yu-one)*(yu-two);
  wyu[2] = cm5*yu*(yu+one)*(yu-two);
  wyu[3] = cp167*yu*(yu+one)*(yu-one);

  wyv[0] = cm167*yv*(yv-one)*(yv-two);       // polynomial coefficients along j
  wyv[1] = cp5*(yv+one)*(yv-one)*(yv-two);
  wyv[2] = cm5*yv*(yv+one)*(yv-two);
  wyv[3] = cp167*yv*(yv+one)*(yv-one);
#else
  u += bicub_coeffs_inline(xxu, yyu, ni, wxu, wyu, uminx, umaxx, uminy, umaxy, uoffx, uoffy);  // use inline function
  v += bicub_coeffs_inline(xxv, yyv, ni, wxv, wyv, vminx, vmaxx, vminy, vmaxy, voffx, voffy);  // use inline function
#endif
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  fwyu0 = _mm256_set1_pd(wyu[0]) ;           // scalar * vector not available, promote scalar to vector
  fwyu1 = _mm256_set1_pd(wyu[1]) ;
  fwyu2 = _mm256_set1_pd(wyu[2]) ;
  fwyu3 = _mm256_set1_pd(wyu[3]) ;
  fwyv0 = _mm256_set1_pd(wyv[0]) ;           // scalar * vector not available, promote scalar to vector
  fwyv1 = _mm256_set1_pd(wyv[1]) ;
  fwyv2 = _mm256_set1_pd(wyv[2]) ;
  fwyv3 = _mm256_set1_pd(wyv[3]) ;
  fwxu  = _mm256_loadu_pd(wxu) ;             // vector of coefficients along i
  fwxv  = _mm256_loadu_pd(wxv) ;             // vector of coefficients along i
  for(k=0 ; k<nk ; k++){
    // interpolation along J for level k
    fd0  = _mm256_cvtps_pd( _mm_loadu_ps(u    ) ) ;   // row 1 : f[i:i+3 , j   ,k]
    fd1  = _mm256_cvtps_pd( _mm_loadu_ps(v    ) ) ;
    fdtu = _mm256_mul_pd(fd0,fwyu0) ;                 // row[j] * coefficient[0]
    fdtv = _mm256_mul_pd(fd1,fwyv0) ;

    fd0  = _mm256_cvtps_pd( _mm_loadu_ps(u+ni ) ) ;   // row 2 : f[i:i+3 , j+1 ,k]
    fd1  = _mm256_cvtps_pd( _mm_loadu_ps(v+ni ) ) ;
    fdtu = _mm256_fmadd_pd(fd0,fwyu1,fdtu) ;          // add row[j+1] * coefficient[1]
    fdtv = _mm256_fmadd_pd(fd1,fwyv1,fdtv) ;

    fd0  = _mm256_cvtps_pd( _mm_loadu_ps(u+ni2) ) ;   // row 3 : f[i:i+3 , j+2 ,k]
    fd1  = _mm256_cvtps_pd( _mm_loadu_ps(v+ni2) ) ;
    fdtu = _mm256_fmadd_pd(fd0,fwyu2,fdtu) ;          // add row[j+2] * coefficient[2]
    fdtv = _mm256_fmadd_pd(fd1,fwyv2,fdtv) ;

    fd0  = _mm256_cvtps_pd( _mm_loadu_ps(u+ni3) ) ;   // row 4 : f[i:i+3 , j+3 ,k]
    fd1  = _mm256_cvtps_pd( _mm_loadu_ps(v+ni3) ) ;
    fdtu = _mm256_fmadd_pd(fd0,fwyu3,fdtu) ;          // add row[j+3] * coefficient[3]
    fdtv = _mm256_fmadd_pd(fd1,fwyv3,fdtv) ;

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdtu = _mm256_mul_pd(fdtu,fwxu) ;
    fdtv = _mm256_mul_pd(fdtv,fwxv) ;

    ft1 = _mm256_extractf128_pd(fdtu,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdtu,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;              // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;          // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;              // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    _mm_store_sd(&tu,ft0) ;                  // store scalar float

    ft1 = _mm256_extractf128_pd(fdtv,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdtv,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;              // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;          // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;              // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    _mm_store_sd(&tv,ft0) ;                  // store scalar float

    ru[k] = tu*s[0] + tv*s[1] ;              // rotate wind components
    rv[k] = tu*s[2] + tv*s[3] ;
    u+= ninj;
    v+= ninj;
  }
#else
  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      ud0[i] = ( u[i]*wyu[0] + u[i+ni]*wyu[1] + u[i+ni2]*wyu[2] + u[i+ni3]*wyu[3] ) * wxu[i];
      vd0[i] = ( v[i]*wyv[0] + v[i+ni]*wyv[1] + v[i+ni2]*wyv[2] + v[i+ni3]*wyv[3] ) * wxv[i];
    }
    tu = ud0[0] + ud0[1] + ud0[2] + ud0[3];
    tv = vd0[0] + vd0[1] + vd0[2] + vd0[3];
    ru[k] = tu*s[0] + tv*s[1] ;              // rotate wind components
    rv[k] = tu*s[2] + tv*s[3] ;
    u+= ninj;
    v+= ninj;
  }
#endif
}

// interpolate u and v with rotation, u and v on same grid, at position (xx,yy)
// xx, yy in fractional index space, same as intrp_bicub_yx_s
// ru, rv : rotated U and V components of the wind at the desired point
// for rotation matrix s, see the GEM model
void intrp_bicub_yx_uv(float *u, float *v, float *ru, float *rv, int ni, int ninj, int nk, double xx, double yy, double s[4]){
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  __m256d ud0, ud1, ud2, ud3, vd0, vd1, vd2, vd3, fwx, fwy0, fwy1, fwy2, fwy3, udt, vdt ;
  __m128d ft0, ft1;
  __m128  frt;
#else
  double ud0[4];
  double vd0[4];
#endif
  double  wx[4], wy[4] ;
  double  x, y;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;
  double tu, tv;
#if defined(NO_INLINE)
  x = xx - 1.0 ; y = yy - 1.0;                      // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  ix = (ix < qminx) ? qminx : ix;       // apply boundary limits (set by set_intrp_bicub_quv or set_intrp_bicub_yx)
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  u = u + ix + iy * ni;                  // point to lower left corner of 4 x 4 subarray
  v = v + ix + iy * ni;

  wx[0] = cm167*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp167*x*(x+one)*(x-one);

  wy[0] = cm167*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp167*y*(y+one)*(y-one);
#else
  i = bicub_coeffs_inline(xx, yy, ni, wx, wy, qminx, qmaxx, qminy, qmaxy, qoffx, qoffy);  // use inline function
  u += i; v+= i;
#endif
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;

  ud0 = _mm256_cvtps_pd(_mm_loadu_ps(u    )) ;  // row 0 : u[i:i+3 , j   ,k]
  ud1 = _mm256_cvtps_pd(_mm_loadu_ps(u+ni )) ;  // row 1 : u[i:i+3 , j+1 ,k]
  ud2 = _mm256_cvtps_pd(_mm_loadu_ps(u+ni2)) ;  // row 2 : u[i:i+3 , j+2 ,k]
  ud3 = _mm256_cvtps_pd(_mm_loadu_ps(u+ni3)) ;  // row 3 : u[i:i+3 , j+3 ,k]

  vd0 = _mm256_cvtps_pd(_mm_loadu_ps(v    )) ;  // row 0 : v[i:i+3 , j   ,k]
  vd1 = _mm256_cvtps_pd(_mm_loadu_ps(v+ni )) ;  // row 1 : v[i:i+3 , j+1 ,k]
  vd2 = _mm256_cvtps_pd(_mm_loadu_ps(v+ni2)) ;  // row 2 : v[i:i+3 , j+2 ,k]
  vd3 = _mm256_cvtps_pd(_mm_loadu_ps(v+ni3)) ;  // row 3 : v[i:i+3 , j+3 ,k]

  for(k=0 ; k<nk-1 ; k++){
    u += ninj; v += ninj;
    // interpolation along y, rows j to J=3, prefetch next level
    udt = _mm256_mul_pd(ud0,fwy0) ;               // sum of row[j] * coefficient[j] ( 0 <= j < 4)
    vdt = _mm256_mul_pd(vd0,fwy0) ;
    ud0 = _mm256_cvtps_pd(_mm_loadu_ps(u    )) ;  // row 0 : u[i:i+3 , j   ,k+1]
    vd0 = _mm256_cvtps_pd(_mm_loadu_ps(v    )) ;  // row 0 : v[i:i+3 , j   ,k+1]

    udt = _mm256_fmadd_pd(ud1,fwy1,udt) ;
    vdt = _mm256_fmadd_pd(vd1,fwy1,vdt) ;
    ud1 = _mm256_cvtps_pd(_mm_loadu_ps(u+ni )) ;  // row 1 : u[i:i+3 , j+1 ,k+1]
    vd1 = _mm256_cvtps_pd(_mm_loadu_ps(v+ni )) ;  // row 1 : v[i:i+3 , j+1 ,k+1]

    udt = _mm256_fmadd_pd(ud2,fwy2,udt) ;
    vdt = _mm256_fmadd_pd(vd2,fwy2,vdt) ;
    ud2 = _mm256_cvtps_pd(_mm_loadu_ps(u+ni2)) ;  // row 2 : u[i:i+3 , j+2 ,k+1]
    vd2 = _mm256_cvtps_pd(_mm_loadu_ps(v+ni2)) ;  // row 2 : v[i:i+3 , j+2 ,k+1]

    udt = _mm256_fmadd_pd(ud3,fwy3,udt) ;
    vdt = _mm256_fmadd_pd(vd3,fwy3,vdt) ;
    ud3 = _mm256_cvtps_pd(_mm_loadu_ps(u+ni3)) ;  // row 3 : u[i:i+3 , j+3 ,k+1]
    vd3 = _mm256_cvtps_pd(_mm_loadu_ps(v+ni3)) ;  // row 3 : v[i:i+3 , j+3 ,k+1]

    // interpolation along x: multiply by coefficients along x , then sum elements (using vector folding)
    udt = _mm256_mul_pd(udt,fwx) ;
    vdt = _mm256_mul_pd(vdt,fwx) ;

    ft1 = _mm256_extractf128_pd(udt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(udt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    _mm_store_sd(&tu,ft0) ;                  // store double

    ft1 = _mm256_extractf128_pd(vdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(vdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    _mm_store_sd(&tv,ft0) ;                  // store double

    ru[k] = tu*s[0] + tv*s[1] ;              // rotate wind components
    rv[k] = tu*s[2] + tv*s[3] ;
  }
  udt = _mm256_mul_pd(ud0,fwy0) ;               // sum of row[j] * coefficient[j] ( 0 <= j < 4)
  vdt = _mm256_mul_pd(vd0,fwy0) ;
  udt = _mm256_fmadd_pd(ud1,fwy1,udt) ;
  vdt = _mm256_fmadd_pd(vd1,fwy1,vdt) ;
  udt = _mm256_fmadd_pd(ud2,fwy2,udt) ;
  vdt = _mm256_fmadd_pd(vd2,fwy2,vdt) ;
  udt = _mm256_fmadd_pd(ud3,fwy3,udt) ;
  vdt = _mm256_fmadd_pd(vd3,fwy3,vdt) ;
  udt = _mm256_mul_pd(udt,fwx) ;
  vdt = _mm256_mul_pd(vdt,fwx) ;

  ft1 = _mm256_extractf128_pd(udt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm256_extractf128_pd(udt,0) ;    // fdt[0]                      fdt[1]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
  _mm_store_sd(&tu,ft0) ;                  // store double

  ft1 = _mm256_extractf128_pd(vdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm256_extractf128_pd(vdt,0) ;    // fdt[0]                      fdt[1]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
  _mm_store_sd(&tv,ft0) ;                  // store double

  ru[k] = tu*s[0] + tv*s[1] ;              // rotate wind components
  rv[k] = tu*s[2] + tv*s[3] ;
#else
  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      ud0[i] = ( u[i]*wy[0] + u[i+ni]*wy[1] + u[i+ni2]*wy[2] + u[i+ni3]*wy[3] ) * wx[i];
      vd0[i] = ( v[i]*wy[0] + v[i+ni]*wy[1] + v[i+ni2]*wy[2] + v[i+ni3]*wy[3] ) * wx[i];
    }
    tu = ud0[0] + ud0[1] + ud0[2] + ud0[3];
    tv = vd0[0] + vd0[1] + vd0[2] + vd0[3];
    ru[k] = tu*s[0] + tv*s[1] ;              // rotate wind components
    rv[k] = tu*s[2] + tv*s[3] ;
    u  += ninj;
    v  += ninj;
  }
#endif
}

// "monotonic" version, the result must be between the max and min values of the 2x2 center
// of the 4x4 subarray used to interpolate
// SIMD version now consistent with  non SIMD version
void intrp_bicub_yx_s_mono(float *f, float *r, int ni, int ninj, int nk, double xx, double yy){
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt, fmi, fma ;
  __m128  fr0, fr1, fr2, fr3, frt ; ;       // frt is used as a scalar, fr0->fr3 are 128 bit aliases for fd0->fd3
  __m128d ft0, ft1, rmi, rma ;
  double dd0[4], dd1[4], dd2[4], dd3[4] ;
#else
  double fd0[4], fd1[4], fd2[4], fd3[4] ;
#endif
  double  wx[4], wy[4] ;
  double  x, y, tm1, tm2;
  float minval, maxval;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;
#if defined(NO_INLINE)
  x = xx - 1.0 ; y = yy - 1.0;                      // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  ix = (ix < qminx) ? qminx : ix;       // apply boundary limits (set by set_intrp_bicub_quv or set_intrp_bicub_yx)
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;                  // lower left corner of 4 x 4 subarray used for (un)centered interpolation

  wx[0] = cm167*x*(x-one)*(x-two);       // polynomial coefficients along i (x)
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp167*x*(x+one)*(x-one);

  wy[0] = cm167*y*(y-one)*(y-two);       // polynomial coefficients along j (y)
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp167*y*(y+one)*(y-one);
#else
  f += bicub_coeffs_inline(xx, yy, ni, wx, wy, qminx, qmaxx, qminy, qmaxy, qoffx, qoffy);  // use inline function
#endif
#if defined(__AVX__) && defined(__x86_64__) && !defined(NO_SIMD)
  fwx  = _mm256_loadu_pd(wx) ;            // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;          // scalar * vector not available,
  fwy1 = _mm256_set1_pd(wy[1]) ;          // promote scalars to vectors
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  frt  = _mm_xor_ps(frt,frt) ;            // set frt to zero
  // prefetch and promote to double 4 rows, level 1
  fd0 = _mm256_cvtps_pd(_mm_loadu_ps(f)) ;           // row 1 : f[i:i+3 , j   ,k]
  fd1 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni)) ;        // row 2 : f[i:i+3 , j+1 ,k]
  fd2 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni2)) ;       // row 3 : f[i:i+3 , j+2 ,k]
  fd3 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni3)) ;       // row 4 : f[i:i+3 , j+3 ,k]

  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    // interpolation along J, level k,
    // prefetch and promote to double 4 rows, level k+1
    fma = _mm256_max_pd(fd2,fd1);                    // max of rows 1 and 2
    fmi = _mm256_min_pd(fd2,fd1);                    // min of rows 1 and 2

    fdt = _mm256_mul_pd(fd0,fwy0) ;                  // sum of row[j] * coefficient[j]
    fd0 = _mm256_cvtps_pd( _mm_loadu_ps(f    )) ;    // row 0 : f[i:i+3 , j   ,k+1]

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fd1 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni )) ;     // row 2 : f[i:i+3 , j+1 ,k+1]

    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fd2 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni2)) ;     // row 3 : f[i:i+3 , j+2 ,k+1]

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fd3 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni3)) ;     // row 4 : f[i:i+3 , j+3 ,k]

    // get maximum of some fma vector elements and minimum of some fmi vector elements
    // elements 0 and 3 ignored, max/min of elements 1 and 2
    rma = _mm_permute_pd( _mm256_extractf128_pd(fma,0), 0x3); // duplicate element 1 of low part [1, 1]
    rmi = _mm_permute_pd( _mm256_extractf128_pd(fmi,0), 0x3); // duplicate element 1 of low part [1, 1]
    rma = _mm_max_sd(rma,_mm_permute_pd(_mm256_extractf128_pd(fma,1),0x0)) ;  // max with low element of high part [1, 1] [2, 3]
    rmi = _mm_min_sd(rmi,_mm_permute_pd(_mm256_extractf128_pd(fmi,1),0x0)) ;  // min with low element of high part [1, 1] [2, 3]
    // rma/rmi[0] = max/min( fma[1], fma [2]),   rma/rmi[1] will be ognored

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]

    ft0 = _mm_max_sd(ft0,rmi);              // apply monotonic limits
    ft0 = _mm_min_sd(ft0,rma);
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert result to float

    _mm_store_ss(r,frt) ;                   // store float
    r ++;
  }
  // interpolation along j , level nk
  fma = _mm256_max_pd(fd2,fd1);
  fmi = _mm256_min_pd(fd2,fd1);
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // get maximum of some fma vector elements and minimum of some fmi vector elements
  // elements 0 and 3 ignored, max/min of elements 1 and 2
  rma = _mm_permute_pd( _mm256_extractf128_pd(fma,0), 0x3); // duplicate element 1 of low part
  rmi = _mm_permute_pd( _mm256_extractf128_pd(fmi,0), 0x3); // duplicate element 1 of low part
  rma = _mm_max_sd(rma,_mm_permute_pd(_mm256_extractf128_pd(fma,1),0x0)) ;  // max with low element of high part
  rmi = _mm_min_sd(rmi,_mm_permute_pd(_mm256_extractf128_pd(fmi,1),0x0)) ;  // min with low element of high part
  // rma/rmi[0] = max/min( fma[1], fma [2]),   rma/rmi[1] will be ognored

  // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]

  ft0 = _mm_max_sd(ft0,rmi);              // apply monotonic limits
  ft0 = _mm_min_sd(ft0,rma);
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert result to float

  _mm_store_ss(r,frt) ;                   // store float
#else
  for(k=0 ; k<nk ; k++){
    minval = f[1];
    maxval = minval;
    for(i=1 ; i<3 ; i++){   // compute monotonic limits, center 2x2 points
      minval = (f[i+ni ] < minval) ? f[i+ni ] : minval;
      minval = (f[i+ni2] < minval) ? f[i+ni2] : minval;
      maxval = (f[i+ni ] > maxval) ? f[i+ni ] : maxval;
      maxval = (f[i+ni2] > maxval) ? f[i+ni2] : maxval;
    }
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    r[0] = (r[0] < minval ) ? minval : r[0] ;  // apply monotonic limits
    r[0] = (r[0] > maxval ) ? maxval : r[0] ;
    f += ninj;
    r ++;
  }
#endif
}
#endif

#if defined(C_TEST)

// this test is no longer functional, it will need to be seriously revised

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define NI 65
#define NJ 27
#define NK 81
#define NP 20

int main(int argc,char **argv){
  float f[NK][NJ][NI] ;
  float r[NK][NP];
  double x[NP], y[NP];
  int i, j, k;
  uint64_t t1, t2;

  for(i=0 ; i<NI ; i++){
    for(j=0 ; j<NJ ; j++){
      for(k=0 ; k<NK ; k++){
	f[k][j][i] = i + j + 2;
      }
    }
  }
  for(i=0 ; i<NP ; i++){
    x[i] = 2 + i ;
    y[i] = 2 + i ;
  }
  t1 = rdtsc();
  for(i=0 ; i<NP ; i++){
    intrp_bicub_yx(&f[0][0][0], &r[0][i], NI, NI*NJ, NK, NP, x[i], y[i]) ;
  }
  t2 = rdtsc();
  k = t2 - t1;
  printf(" r = %f %f\n",r[0][0],r[0][1]);
  printf("time = %d clocks for %d values, %d flops\n",k,NP*NK,NP*NK*35);
  t1 = rdtsc();
  for(i=0 ; i<NP ; i++){
    intrp_bicub_yx_2(&f[0][0][0], &r[0][i], NI, NI*NJ, NK, NP, x[i], y[i]) ;
  }
  t2 = rdtsc();
  k = t2 - t1;
  printf(" r = %f %f\n",r[0][0],r[0][1]);
  printf("time = %d clocks for %d values, %d flops\n",k,NP*NK,NP*NK*35);
}
#endif

#if defined(F_TEST)
program test_interp   !! embedded Fortran test program
  use ISO_C_BINDING
  implicit none
  integer, parameter :: NI=65
  integer, parameter :: NJ=27
  integer, parameter :: NK=81
  integer, parameter :: NP=4
  integer, parameter :: HX=2
  integer, parameter :: HY=2
  integer, parameter :: NR=10
  integer, parameter :: FLOPS1=NP*NK*35+40   ! scalars, 1 position
  integer, parameter :: FLOPS2=NP*NK*70+80   ! vectors, 2 positions
  integer, parameter :: FLOPS3=NP*NK*70+40   ! vectors, 1 position
  real(C_FLOAT), dimension(1-HX:NI+HX , 1-HY:NJ+HY , NK) :: f, f2
  real(C_FLOAT), dimension(NK,NP*2) :: r, r2, r1
  real(C_DOUBLE), dimension(NP) :: x, y, xmin, xmax, xmin2, xmax2
  real(C_DOUBLE) :: xi, yi
  real(C_DOUBLE), dimension(4) :: s
  integer :: i, j, k, ix, iy
  integer :: i0, j0
  integer*8, external :: rdtsc
  integer*8 :: t1, t2, tmg1(NR), tmg2(nr), tmg3(nr), tmg4(nr)
  integer :: nidim, ninjdim
  include 'intrp_bicub_yx.inc'

#define FXY(A,B,C) (1.3*(A)**3 + 1.4*(A)**2 + (A)*1.5 + 2.3*(B)**3 + 2.4*(B)**2 + (B)*2.5 + (C))

  call set_intrp_bicub_quv(1,NI-3,1,NJ-3,1,NI-3,1,NJ-3,1,NI-3,1,NJ-3)
  s = [1.0, 0.0, 0.0, 1.0]

  r = 9999.99
  do k = 1 , NK
    do j = 1-HY , NJ+HY
      do i = 1-HX , NI+HX
!        f(i,j,k) = i + j  + k
        f(i,j,k) = FXY(i*1.0 , j*1.0 , k*1.0)
      enddo
    enddo
!    print *,f(1,1,k),f(2,2,k)
  enddo
  f2 = f + 1
  do i = 1 , NP
    x(i) = i + 1.0 !- .999
    y(i) = i + 1.0 !- .999
  enddo
  x(NP) = ni - 1!.01
  y(NP) = nj - 1!.01
  nidim = NI + 2*HX
  ninjdim = nidim * (NJ + HY*2)
!  print *,'nidim=',nidim,' , ninjdim=',ninjdim
!  print *,'x=',x,' y=',y
!  print *,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),2),f(nint(x(2)),nint(y(2)),2)
!  print 101,loc(f(1,1,1)), loc(r(1,1))
101 format(2Z17)
!  print *,f(1,1,1), r(1,1)
  r  = 0.0
  r2 = 0.0
  do i = 1 , NP
    call intrp_bicub_yx_s(f(1,1,1), r(1,i), nidim, ninjdim, NK, x(i), y(i))  ! prime the timinng pump !
  enddo
  do j = 1, NR
    t1 = rdtsc()
    do i = 1 , NP
!       ix = x(i)
!       ix = max(1,ix-1)
!       ix = min(ix,NI-3)
!       xi = x(i) - (ix + 1)
!       iy = y(i)
!       iy = max(1,iy-1)
!       iy = min(iy,NJ-3)
!       yi = y(i) - (iy - 1)
!       call intrp_bicub_yx_s(f(ix,iy,1), r(i,1), nidim, ninjdim, NK, xi, yi)
      call intrp_bicub_yx_s(f(1,1,1), r(1,i), nidim, ninjdim, NK, x(i), y(i))
    enddo
    t2 = rdtsc()
    tmg1(j) = t2 - t1
!    print *,'time=',t2-t1,' cycles for',NP*NK*35,' values'

    t1 = rdtsc()
    do i = 1 , NP
      call intrp_bicub_yx_vec( f(1,1,1), f2(1,1,1), r(1,i), r2(1,i), nidim, ninjdim, NK, x(i), y(i), x(1+mod(i+1,NP)), y(1+mod(i+1,NP)), s )
    enddo
    t2 = rdtsc()
    tmg3(j) = t2 - t1

    t1 = rdtsc()
    do i = 1 , NP
      call intrp_bicub_yx_uv( f(1,1,1), f2(1,1,1), r(1,i), r2(1,i), nidim, ninjdim, NK, x(i), y(i), s )
    enddo
    t2 = rdtsc()
    tmg4(j) = t2 - t1

    t1 = rdtsc()
    do i = 1 , NP
      call intrp_bicub_yx_s_mono( f(1,1,1), r(1,i), nidim, ninjdim, NK, x(i), y(i) )
    enddo
    t2 = rdtsc()
    tmg2(j) = t2 - t1
!    print *,'time=',t2-t1,' cycles for',NP*NK*35,' values'
  enddo
  print *,"Gflops (@ 3.7 GHz TSC)"
  print 100,'direct   =',FLOPS1*3.7/tmg1
  print 100,'mono     =',FLOPS1*3.7/tmg2
  print 100,'vec 2pos =',FLOPS2*3.7/tmg3
  print 100,'vec 1pos =',FLOPS3*3.7/tmg4
  print 1000,'flops    =',FLOPS1,FLOPS1,FLOPS2,FLOPS3
  print 1000,'Bytes    =',NP*NK*16*4
  print 100,'GBytes/s =',NP*NK*16*4*3.7/tmg1
100 format(A,40F8.3)
1000 format(A,40I6)
!  print 102,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),NK),f(nint(x(2)),nint(y(2)),NK)
!  print 102,x(1)+y(1)+1,x(2)+y(2)+1,x(1)+y(1)+NK,x(2)+y(2)+NK
  print *,"X coordinates, X dimension :",NI
  print 102,x(:)
  print *,"Y coordinates, Y dimension :",NJ
  print 102,y(:)
  print *,"F matrix (1)"
  do j = 6, 1-hy, -1
    print 102,f(1-hx:6,j,1)
  enddo
  print *,"F matrix (nk)"
  do j = 6, 1-hy, -1
    print 102,f(1-hx:6,j,NK)
  enddo
  r = -1.0
  print *,"MONO: limit (min, max)"
  do i = 1 , NP
    i0 = x(i)
    i0 = max(1,min(ni-3,i0-1))
    j0 = y(i)
    j0 = max(1,min(nj-3,j0-1))
    call intrp_bicub_yx_s(f(1,1,1), r(1,i), nidim, ninjdim, NK, x(i), y(i))
!   print *,'DBG', i0, j0
!   print 102, x(i), y(i), r(1,i), r(NK,i)
    call intrp_bicub_yx_s_mono( f(1,1,1), r(1,i), nidim, ninjdim, NK, x(i), y(i) )
!   print 102,r(1,i),f(i0+1:i0+2,j0+1:j0+2,1)
!   print 102,r(NK,i),f(i0+1:i0+2,j0+1:j0+2,NK)
    xmin(i)  = minval(f(i0+1:i0+2,j0+1:j0+2,1))
    xmax(i)  = maxval(f(i0+1:i0+2,j0+1:j0+2,1))
    xmin2(i) = minval(f(i0+1:i0+2,j0+1:j0+2,NK))
    xmax2(i) = maxval(f(i0+1:i0+2,j0+1:j0+2,NK))
  enddo
  print 102,xmin(:) , xmin2(:)
  print 102,xmax(:) , xmax2(:)
  print *,"MONO: expected"
  print 102 , max( xmin(:),min( xmax(:),FXY(x(:),y(:),1 ))), &
              max(xmin2(:),min(xmax2(:),FXY(x(:),y(:),NK)))
!  print 102,x(:)+y(:)+1, x(:)+y(:)+NK
  print *," got"
  print 102,r(1,1:NP),r(NK,1:NP)
  print *," delta"
  print 102,(r(1,1:NP)-FXY(x(:),y(:),1)),(r(NK,1:NP)-FXY(x(:),y(:),NK))
  print 103,(r(1,1:NP)-FXY(x(:),y(:), 1))  / FXY(x(:),y(:), 1) , &
            (r(NK,1:NP)-FXY(x(:),y(:),NK))  / FXY(x(:),y(:),NK)

  do i = 1 , NP
    ix = x(i)
    ix = max(1,ix-1)
    ix = min(ix,NI-3)
    xi = x(i) - (ix + 1)
    iy = y(i)
    iy = max(1,iy-1)
    iy = min(iy,NJ-3)
    yi = y(i) - (iy + 1)
!    print *,"DEBUG: x(i), ix, xi, Y(i), iy , y(i)",x(i),ix,xi,y(i),iy,yi
    call intrp_bicub_yx_s( f(1,1,1), r(1,i), nidim, ninjdim, NK, x(i), y(i) )
    call intrp_bicub_yx_uv( f(1,1,1), f2(1,1,1), r1(1,i), r2(1,i), nidim, ninjdim, NK, x(i), y(i), s )
  enddo
!  print 102,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),NK),f(nint(x(2)),nint(y(2)),NK)
!  print 102,x(1)+y(1)+1,x(2)+y(2)+1,x(1)+y(1)+NK,x(2)+y(2)+NK
  print *,"DIRECT: expected"
  print 102,FXY(x(:),y(:),1), FXY(x(:),y(:),NK)
  print *," got"
  print 102,r(1,1:NP),r(NK,1:NP)
  print *," delta"
  print 102,(r(1,1:NP)-FXY(x(:),y(:),1)),(r(NK,1:NP)-FXY(x(:),y(:),NK))
  print 103,(r(1,1:NP)-FXY(x(:),y(:), 1))  / FXY(x(:),y(:), 1) , &
            (r(NK,1:NP)-FXY(x(:),y(:),NK))  / FXY(x(:),y(:),NK)
  print *," gotuv (second line : first line + 1.0)"
  print 102,r1(1,1:NP),r1(NK,1:NP)
  print 102,r2(1,1:NP),r2(NK,1:NP)
  r1 = 999.999
  r2 = 999.999
  do i = 1 , NP  ! x, y positions are rotated right by 1
    call intrp_bicub_yx_vec( f(1,1,1), f2(1,1,1), r1(1,i), r2(1,i), nidim, ninjdim, NK, x(1+mod(i,NP)), y(1+mod(i,NP)), x(i), y(i), s )
  enddo
  print *," vec, first line rotated left by one position with respect to direct"
  print 102,r1(1,1:NP),r1(NK,1:NP)
  print 102,r2(1,1:NP),r2(NK,1:NP)

102 format(16F15.6)
103 format(16E15.2)
end

#endif
