#if ! defined(F_TEST)
/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2017  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
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
#if defined(NEVER_TRUE)
 to test (Intel compilers) :
 rm -f intrp_bicub_yx intrp_bicub_yx.o intrp_bicub_yx_f.F90
 s.cc -O2 -DTIMING -c -march=core-avx2 intrp_bicub_yx.c
 ln -sf intrp_bicub_yx.c intrp_bicub_yx_f.F90
 s.f90 -no-wrap-margin -DF_TEST -o intrp_bicub_yx intrp_bicub_yx_f.F90 intrp_bicub_yx.o
 ./intrp_bicub_yx
#endif

static float cp133 =  0.166666666666666667E0;
static float cm133 = -0.166666666666666667E0;
static float cp5 =  .5;
static float cm5 = -.5;
static float one = 1.0;
static float two = 2.0;

#if defined(TIMING)
#include <stdio.h>
#include <stdint.h>
#pragma weak rdtsc_=rdtsc
uint64_t rdtsc_(void);
uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}
#endif
#include <immintrin.h>

// use separate multiply and add instructions if fused multiply-add not available
#if defined(__AVX__) && defined(__x86_64__)
#if ! defined(__FMA__)
#define _mm256_fmadd_ps(a,b,c) _mm256_add_ps(_mm256_mul_ps(a,b),c)
#define _mm256_fmadd_pd(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b),c)
#define _mm_fmadd_ps(a,b,c) _mm_add_ps(_mm_mul_ps(a,b),c)
#define _mm_fmadd_pd(a,b,c) _mm_add_pd(_mm_mul_pd(a,b),c)
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
/*
interface
  subroutine set_intrp_bicub_yx(x1,x2,y1,y2) bind(C,name='set_intrp_bicub_yx')
  integer, intent(IN), value :: x1,x2,y1,y2
  end subroutine set_intrp_bicub_yx
  subroutine set_intrp_bicub_quv(qx1,qx2,qy1,qy2,ux1,ux2,uy1,uy2,vx1,vx2,vy1,vy2) bind(C,name='set_intrp_bicub_quv')
  integer, intent(IN), value :: qx1,qx2,qy1,qy2,ux1,ux2,uy1,uy2,vx1,vx2,vy1,vy2
  end subroutine set_intrp_bicub_quv
end interface
 */

// set limits for 'q' grid only
void set_intrp_bicub_yx(int qxmin, int qxmax, int qymin, int qymax){
  qminx = qxmin - 1 ;   // qxmin is in "origin 1", qminx is in "origin 0"
  qmaxx = qxmax - 1 ;
  qminy = qymin - 1 ;
  qmaxy = qymax - 1 ;
}
// set limits for 'q', 'u', and 'v' grids
void set_intrp_bicub_quv(int qxmin, int qxmax, int qymin, int qymax,
                         int uxmin, int uxmax, int uymin, int uymax,
                         int vxmin, int vxmax, int vymin, int vymax){
  qminx = qxmin - 1 ;  // qxmin is in "origin 1", qminx is in "origin 0"
  qmaxx = qxmax - 1 ;
  qminy = qymin - 1 ;
  qmaxy = qymax - 1 ;
  uminx = uxmin - 1 ;
  umaxx = uxmax - 1 ;
  uminy = uymin - 1 ;
  umaxy = uymax - 1 ;
  vminx = vxmin - 1 ;
  vmaxx = vxmax - 1 ;
  vminy = vymin - 1 ;
  vmaxy = vymax - 1 ;
}
/*
   interpolate a column from a 3D source array f, put results in array r

   INPUT:

   f       3D Fortran indexed source array, real*4, dimension(mini:maxi,minj:maxj,nk)
           ni = (maxi-mini-1)
           ninj = (maxi-mini-1)*(maxj-minj-1)
   ni      distance between f(i,j,k) and f(i,j+1,k)
   ninj    distance between f(i,j,k) and f(i,j,k+1)
   nk      number of 2D planes
   np      distance between 2 kevels in result array r
   x       fractional index along x with respect to submatrix(2,2)
           (centered case, 0<=x < 1.0) (uncentered case : -1.0 <= x <= 2.0)
   y       fractional index along y with respect to submatrix(2,2)
           (centered case, 0<=y < 1.0) (uncentered case : -1.0 <= y <= 2.0)

   OUTPUT:

   r       2D array, fortran dimension(np,nk) ( will store r(1,:) )

   this routine expects the address of f(ix,iy), lower left corner of the 4x4 submatrix used for interpolation
*/

void intrp_bicub_yx_s(float *f, float *r, int ni, int ninj, int nk, double xx, double yy){
#if defined(__AVX__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt ;
  __m128  fr0, fr1, fr2, fr3, frt ;
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

  x = xx - 1.0 ; y = yy - 1.0;                      // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  ix = (ix < qminx) ? qminx : ix;
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;                  // point to lower left corner of 4 x 4 subarray

  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  // fetch 4 rows, level 1
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    // interpolation along J level k, prefetch 4 rows for level k+1
    fdt = _mm256_mul_pd(fd0,fwy0) ;         // sum of row[j] * coefficient[j]
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    _mm_store_ss(r,frt) ;                   // store float
    r ++;
  }
  // interpolation along j , level nk
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // interpolation along i: multiply by coefficients along x , then sum elements
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  _mm_store_ss(r,frt) ;                   // store float
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
void intrp_nk2d_yx_s_mono(float *f, float *r, int ni, int ninj, int nk, double x, double y){
#if defined(__AVX__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt ;
  __m128  fr0, fr1, fr2, fr3, frt, fmi, fma ;
  __m128d ft0, ft1;
#else
  double fd0[4], fd1[4], fd2[4], fd3[4] ;
  float minval, maxval;
#endif
  double  wx[4], wy[4] ;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  float *m = f;    // pointer used for the min/max constraint

  if(x >= 0.0) m++;
  if(x >  1.0) m++;
  if(y >= 0.0) m+=ni;
  if(y >  1.0) m+=ni;  // m should now point to the lower left corner of the 2 by 2 limit matrix

  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  // fetch 4 rows, level 1
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    fma = _mm_loadu_ps(m);                // shuffle 0,1,0,1 might be needed
    fma = _mm_permute_ps(fma,0x44);       // [0] [1] [0] [1]
    fmi = fma;
    frt = _mm_loadu_ps(m+ni);             // shuffle 0,1,0,1 might be needed
    frt = _mm_permute_ps(frt,0x44);       // [0] [1] [0] [1]
    fma = _mm_max_ps(frt,fma);
    fmi = _mm_min_ps(frt,fmi);
    m+= ninj;
    // interpolation along J level k, prefetch 4 rows for level k+1
    fdt = _mm256_mul_pd(fd0,fwy0) ;         // sum of row[j] * coefficient[j]
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    frt = _mm_permute_ps(fma,0x55);         // [1] [1] [1] [1]
    fma = _mm_max_ss(frt,fma);              // max of first 2 values in fma
    frt = _mm_permute_ps(fmi,0x55);         // [1] [1] [1] [1]
    fmi = _mm_min_ss(frt,fmi);              // min of first 2 values in fma

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    frt = _mm_max_ss(frt,fmi) ;             // apply max(value,minval)
    frt = _mm_min_ss(frt,fma) ;             // apply min(value,maxval)
    _mm_store_ss(r,frt) ;                   // store float
    r ++;
  }
  // interpolation along j , level nk
  fma = _mm_loadu_ps(m);
  fma = _mm_permute_ps(fma,0x44);       // [0] [1] [0] [1]
  fmi = fma;
  frt = _mm_loadu_ps(m+ni);
  frt = _mm_permute_ps(frt,0x44);       // [0] [1] [0] [1]
  fma = _mm_max_ps(frt,fma);
  fmi = _mm_min_ps(frt,fmi);
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;

  frt = _mm_permute_ps(fma,0x55);         // [1] [1] [1] [1]
  fma = _mm_max_ss(frt,fma);              // max of first 2 values in fma
  frt = _mm_permute_ps(fmi,0x55);         // [1] [1] [1] [1]
  fmi = _mm_min_ss(frt,fmi);              // min of first 2 values in fma

  // interpolation along i: multiply by coefficients along x , then sum elements
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  frt = _mm_max_ss(frt,fmi) ;             // apply max(value,minval)
  frt = _mm_min_ss(frt,fma) ;             // apply min(value,maxval)
  _mm_store_ss(r,frt) ;                   // store float
#else
  for(k=0 ; k<nk ; k++){
    minval = (m[0]    < m[1]  ) ? m[0]    : m[1]   ;
    minval = (m[ni  ] < minval) ? m[ni]   : minval ;
    minval = (m[ni+1] < minval) ? m[ni+1] : minval ;
    maxval = (m[0]    > m[1]  ) ? m[0]    : m[1] ;
    maxval = (m[ni  ] < maxval) ? m[ni]   : maxval ;
    maxval = (m[ni+1] < maxval) ? m[ni+1] : maxval ;
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    r[0] = (maxval < r[0]) ? maxval : r[0] ;
    r[0] = (minval > r[0]) ? minval : r[0] ;
    f+= ninj;
    r ++;
  }
#endif
}
/*
           interpolate a column from a 3D source array f, put results in array r

   f       3D Fortran indexed source array, real*4, dimension(mini:maxi,minj:maxj,nk)
           ni = (maxi-mini-1)
           ninj = (maxi-mini-1)*(maxj-minj-1)
   ni      distance between f(i,j,k) and f(i,j+1,k)
   ninj    distance between f(i,j,k) and f(i,j,k+1)
   r       2D array, fortran dimension(np,nk)
   nk      number of 2D planes
   xx      i coordinate in i j fractional index space of desired column ( xx(l) is assumed >= mini+1 )
   xxu     i coordinate in i j fractional index space of desired column ( xxu(l) is assumed >= mini+1 )
   xxv     i coordinate in i j fractional index space of desired column ( xxv(l) is assumed >= mini+1 )
   yy      j coordinate in i j fractional index space of desired column ( yy(l) is assumed >= minj+1 )
   yyu     j coordinate in i j fractional index space of desired column ( yyu(l) is assumed >= minj+1 )
   yyv     j coordinate in i j fractional index space of desired column ( yyv(l) is assumed >= minj+1 )

   f is assumed to point to the address of f(1,1,1)  (Fortran indexing)
   xx, yy positions are in fractional indes space (origin 1)
   the value at (1.0, 1.0) is f(1,1,1)

   xx = 2.5, yy = 2.5 would be the center ot the square formed by
   f(2,2,k) f(3,2,k) f(2,3,k) f(3.3.k)  (where 1 <= k <= nk)

   to call from FORTRAN, the following interface is used

    subroutine intrp_bicub_yx(f, r, ni, ninj, nk, np, x, y)      bind(C,name='intrp_bicub_yx')
    subroutine intrp_bicub_yx_s_mono(f, r, ni, ninj, nk, np, x, y) bind(C,name='intrp_bicub_yx_s_mono')
      real(C_FLOAT), dimension(*), intent(IN) :: f
      real(C_FLOAT), dimension(*), intent(OUT) :: r
      real(C_DOUBLE), intent(IN), value :: x, y
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np
 */

// interpolate u and v with rotation, u and v are on different grids
// hence the 2 sets of positions (xxu,yyu) and (xxv,yyv)
// xxu, yyu, xxv, yyv in fractional index space, same as intrp_bicub_yx_s
void intrp_bicub_yx_vec(float *u, float *v, float *ru, float *rv, int ni, int ninj, int nk,
                        double xxu, double yyu, double xxv, double yyv, double s[4]){
#if defined(__AVX__) && defined(__x86_64__)
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

  xu = xxu - 1.0 ; yu = yyu - 1.0; // xxu and yyu are in "ORIGIN 1"
  ixu = xxu ; if(ixu > xxu) ixu = ixu -1 ; ixu = ixu - 2;
  iyu = yyu ; if(iyu > yyu) iyu = iyu -1 ; iyu = iyu - 2;
  ixu = (ixu < uminx) ? uminx : ixu;
  ixu = (ixu > umaxx) ? umaxx : ixu;
  iyu = (iyu < uminy) ? uminy : iyu;
  iyu = (iyu > umaxy) ? umaxy : iyu;
  xu  = xu - 1 - ixu;
  yu  = yu - 1 - iyu;
  u = u + ixu + iyu * ni;

  xv = xxv - 1.0 ; yv = yyv - 1.0; // xxv and yyv are in "ORIGIN 1"
  ixv = xxv ; if(ixv > xxv) ixv = ixv -1 ; ixv = ixv - 2;
  iyv = yyv ; if(iyv > yyv) iyv = iyv -1 ; iyv = iyv - 2;
  ixv = (ixv < vminx) ? vminx : ixv;
  ixv = (ixv > vmaxx) ? vmaxx : ixv;
  iyv = (iyv < vminy) ? vminy : iyv;
  iyv = (iyv > vmaxy) ? vmaxy : iyv;
  xv  = xv - 1 - ixv;
  yv  = yv - 1 - iyv;
  v = v + ixv + iyv * ni;

  wxu[0] = cm133*xu*(xu-one)*(xu-two);       // polynomial coefficients along i
  wxu[1] = cp5*(xu+one)*(xu-one)*(xu-two);
  wxu[2] = cm5*xu*(xu+one)*(xu-two);
  wxu[3] = cp133*xu*(xu+one)*(xu-one);

  wxv[0] = cm133*xv*(xv-one)*(xv-two);       // polynomial coefficients along i
  wxv[1] = cp5*(xv+one)*(xv-one)*(xv-two);
  wxv[2] = cm5*xv*(xv+one)*(xv-two);
  wxv[3] = cp133*xv*(xv+one)*(xv-one);

  wyu[0] = cm133*yu*(yu-one)*(yu-two);       // polynomial coefficients along j
  wyu[1] = cp5*(yu+one)*(yu-one)*(yu-two);
  wyu[2] = cm5*yu*(yu+one)*(yu-two);
  wyu[3] = cp133*yu*(yu+one)*(yu-one);

  wyv[0] = cm133*yv*(yv-one)*(yv-two);       // polynomial coefficients along j
  wyv[1] = cp5*(yv+one)*(yv-one)*(yv-two);
  wyv[2] = cm5*yv*(yv+one)*(yv-two);
  wyv[3] = cp133*yv*(yv+one)*(yv-one);

#if defined(__AVX__) && defined(__x86_64__)
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
void intrp_bicub_yx_uv(float *u, float *v, float *ru, float *rv, int ni, int ninj, int nk, double xx, double yy, double s[4]){
#if defined(__AVX__) && defined(__x86_64__)
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

  x = xx - 1.0 ; y = yy - 1.0;                      // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  ix = (ix < qminx) ? qminx : ix;
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  u = u + ix + iy * ni;                  // point to lower left corner of 4 x 4 subarray
  v = v + ix + iy * ni;

  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX__) && defined(__x86_64__)
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
//     frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    _mm_store_sd(&tu,ft0) ;                  // store float

    ft1 = _mm256_extractf128_pd(vdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(vdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
//     frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    _mm_store_sd(&tv,ft0) ;                  // store float

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
//   frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  _mm_store_sd(&tu,ft0) ;                  // store float

  ft1 = _mm256_extractf128_pd(vdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm256_extractf128_pd(vdt,0) ;    // fdt[0]                      fdt[1]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
//   frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  _mm_store_sd(&tv,ft0) ;                  // store float

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

void intrp_bicub_yx(float *f, float *r, int ni, int ninj, int nk, int np, double xx, double yy){
#if defined(__AVX__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt ;
  __m128  fr0, fr1, fr2, fr3, frt ;
  __m128d ft0, ft1;
#else
  double fd0[4], fd1[4], fd2[4], fd3[4] ;
#endif
  double  wx[4], wy[4] ;
  double  x, y;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;
// printf("DEBUG: f = %p, r = %p \n",f,r);
// printf("DEBUG: f[0] = %f, r[0] = %f, ni = %d, ninj = %d, nk = %d, np = %d, xx = %f, yy = %f\n",f[0],r[0],ni,ninj,nk,np,xx,yy);
  x = xx - 1.0 ; y = yy - 1.0; // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;   // xx and yy are in "ORIGIN 1"
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  ix = (ix < qminx) ? qminx : ix;
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;
//   printf("DEBUG: ix=%d, iy=%d, ni=%d\n",ix, iy, ni);
// printf("DEBUG: f[0] = %f, ix = %d, iy = %d, x = %f, y = %f\n",f[0],ix,iy,x,y);
  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  // fetch 4 rows, level 1
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    // interpolation along J level k, prefetch 4 rows for level k+1
    fdt = _mm256_mul_pd(fd0,fwy0) ;         // sum of row[j] * coefficient[j]
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    _mm_store_ss(r,frt) ;                   // store float
    r += np;
  }
  // interpolation along j , level nk
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // interpolation along i: multiply by coefficients along x , then sum elements
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  _mm_store_ss(r,frt) ;                   // store float
#else
  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    f+= ninj;
    r += np;
  }
#endif
}

void intrp_bicub_yx_s_mono(float *f, float *r, int ni, int ninj, int nk, double xx, double yy){
#if defined(__AVX__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt, fmi, fma ;
  __m128  fr0, fr1, fr2, fr3, frt ; ;       // frt is used as a scalar, fr0->fr3 are 128 bit aliases for fd0->fd3
  __m128d ft0, ft1, rmi, rma ;
  double dd0[4], dd1[4], dd2[4], dd3[4] ;
#else
  double fd0[4], fd1[4], fd2[4], fd3[4] ;
#endif
  double  wx[4], wy[4] ;
  double  x, y;
  float minval, maxval;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;

  x = xx - 1.0 ; y = yy - 1.0;                      // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  ix = (ix < qminx) ? qminx : ix;       // apply bounary limits (set by set_intrp_bicub_quv or set_intrp_bicub_yx)
  ix = (ix > qmaxx) ? qmaxx : ix;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;                  // lower left corner of 4 x 4 subarray used for (un)centered interpolation

  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i (x)
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j (y)
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX__) && defined(__x86_64__)
  fwx  = _mm256_loadu_pd(wx) ;            // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;          // scalar * vector not available,
  fwy1 = _mm256_set1_pd(wy[1]) ;          // promote scalars to vectors
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  frt  = _mm_xor_ps(frt,frt) ;            // set frt to zero
  // prefetch and promote to double 4 rows, level 1
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    // interpolation along J, level k,
    // prefetch and promote to double 4 rows, level k+1
    fma = _mm256_max_pd(fd1,fd0);
    fmi = _mm256_min_pd(fd1,fd0);
    fdt = _mm256_mul_pd(fd0,fwy0) ;          // sum of row[j] * coefficient[j]
    fd0 = _mm256_cvtps_pd( _mm_loadu_ps(f    )) ;    // row 0 : f[i:i+3 , j   ,k+1]

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fd1 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni )) ;     // row 2 : f[i:i+3 , j+1 ,k+1]

    fma = _mm256_max_pd(fd2,fma);
    fmi = _mm256_min_pd(fd2,fmi);
    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fd2 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni2)) ;     // row 3 : f[i:i+3 , j+2 ,k+1]

    fma = _mm256_max_pd(fd3,fma);
    fmi = _mm256_min_pd(fd3,fmi);
    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fd3 = _mm256_cvtps_pd(_mm_loadu_ps(f+ni3)) ;     // row 4 : f[i:i+3 , j+3 ,k]

    // get maximum of fma vector elements and minimum of fmi vector elements
    rma = _mm_max_pd(_mm256_extractf128_pd(fma,0),_mm256_extractf128_pd(fma,1)) ;
    rmi = _mm_min_pd(_mm256_extractf128_pd(fmi,0),_mm256_extractf128_pd(fmi,1)) ;
    rma = _mm_max_sd(rma,_mm_permute_pd(rma,0x3)) ;
    rmi = _mm_min_sd(rmi,_mm_permute_pd(rmi,0x3)) ;

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
  fma = _mm256_max_pd(fd1,fd0);
  fmi = _mm256_min_pd(fd1,fd0);
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fma = _mm256_max_pd(fd2,fma);
  fmi = _mm256_min_pd(fd2,fmi);
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fma = _mm256_max_pd(fd3,fma);
  fmi = _mm256_min_pd(fd3,fmi);
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // get maximum of fma vector elements and minimum of fmi vector elements
  rma = _mm_max_pd(_mm256_extractf128_pd(fma,0),_mm256_extractf128_pd(fma,1)) ;
  rmi = _mm_min_pd(_mm256_extractf128_pd(fmi,0),_mm256_extractf128_pd(fmi,1)) ;
  rma = _mm_max_sd(rma,_mm_permute_pd(rma,0x3)) ;
  rmi = _mm_min_sd(rmi,_mm_permute_pd(rmi,0x3)) ;

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
    minval = f[0];
    maxval = minval;
    for(i=0 ; i<4 ; i++){   // compute monotonic limits
      minval = (f[i    ] < minval) ? f[i    ] : minval;
      minval = (f[i+ni ] < minval) ? f[i+ni ] : minval;
      minval = (f[i+ni2] < minval) ? f[i+ni2] : minval;
      minval = (f[i+ni3] < minval) ? f[i+ni3] : minval;
      maxval = (f[i    ] > maxval) ? f[i    ] : maxval;
      maxval = (f[i+ni ] > maxval) ? f[i+ni ] : maxval;
      maxval = (f[i+ni2] > maxval) ? f[i+ni2] : maxval;
      maxval = (f[i+ni3] > maxval) ? f[i+ni3] : maxval;
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

#if defined(C_TEST)
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

#else
#if defined(NEVER_TRUE)
  interface                                                                                       !InTf!
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
  end interface                                                                                   !InTf!
#endif
program test_interp
  use ISO_C_BINDING
  implicit none
  integer, parameter :: NI=65
  integer, parameter :: NJ=27
  integer, parameter :: NK=81
  integer, parameter :: NP=4
  integer, parameter :: HX=2
  integer, parameter :: HY=2
  integer, parameter :: NR=10
  real(C_FLOAT), dimension(1-HX:NI+HX , 1-HY:NJ+HY , NK) :: f, f2
  real(C_FLOAT), dimension(NK,NP*2) :: r, r2, r1
  real(C_DOUBLE), dimension(NP) :: x, y, xmin, xmax, xmin2, xmax2
  real(C_DOUBLE) :: xi, yi
  real(C_DOUBLE), dimension(4) :: s
  integer :: i, j, k, ix, iy
  integer :: i0, j0
  integer(kind=INT64), external :: rdtsc
  integer(kind=INT64) :: t1, t2, tmg1(NR), tmg2(nr), tmg3(nr), tmg4(nr)
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
      end do
    end do
!    print *,f(1,1,k),f(2,2,k)
  end do
  f2 = f + 1
  do i = 1 , NP
    x(i) = i - .999 + 1.0
    y(i) = i - .999 + 1.0
  end do
  x(NP) = ni - .01
  y(NP) = nj - .01
  nidim = NI + 2*HX
  ninjdim = nidim * (NJ + HY*2)
!  print *,'nidim=',nidim,' , ninjdim=',ninjdim
!  print *,'x=',x,' y=',y
!  print *,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),2),f(nint(x(2)),nint(y(2)),2)
!  print 101,loc(f(1,1,1)), loc(r(1,1))
101 format(2Z17)
!  print *,f(1,1,1), r(1,1)
  do j = 1, NR
    t1 = rdtsc()
    do i = 1 , NP
      ix = x(i)
      ix = max(1,ix-1)
      ix = min(ix,NI-3)
      xi = x(i) - (ix + 1)
      iy = y(i)
      iy = max(1,iy-1)
      iy = min(iy,NJ-3)
      yi = y(i) - (iy - 1)
      call intrp_bicub_yx_s(f(ix,iy,1), r(i,1), nidim, ninjdim, NK, xi, yi)
    end do
    t2 = rdtsc()
    tmg1(j) = t2 - t1
!    print *,'time=',t2-t1,' cycles for',NP*NK*35,' values'

    t1 = rdtsc()
    do i = 1 , NP
      call intrp_bicub_yx_vec( f(1,1,1), f2(1,1,1), r(1,i), r2(1,i), nidim, ninjdim, NK, x(i), y(i), x(1+mod(i+1,NP)), y(1+mod(i+1,NP)), s )
    end do
    t2 = rdtsc()
    tmg3(j) = t2 - t1

    t1 = rdtsc()
    do i = 1 , NP
      call intrp_bicub_yx_uv( f(1,1,1), f2(1,1,1), r(1,i), r2(1,i), nidim, ninjdim, NK, x(i), y(i), s )
    end do
    t2 = rdtsc()
    tmg4(j) = t2 - t1

    t1 = rdtsc()
    do i = 1 , NP
      call intrp_bicub_yx_s_mono( f(1,1,1), r(1,i), nidim, ninjdim, NK, x(i), y(i) )
    end do
    t2 = rdtsc()
    tmg2(j) = t2 - t1
!    print *,'time=',t2-t1,' cycles for',NP*NK*35,' values'
  end do
  print 100,'direct   =',tmg1
  print 100,'mono     =',tmg2
  print 100,'vec 2pos =',tmg3
  print 100,'vec 1pos =',tmg4
  print 100,'flops    =',NP*NK*35+40,NP*NK*70+80
  print 100,'Bytes    =',NP*NK*16*4
100 format(A,40I6)
!  print 102,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),NK),f(nint(x(2)),nint(y(2)),NK)
!  print 102,x(1)+y(1)+1,x(2)+y(2)+1,x(1)+y(1)+NK,x(2)+y(2)+NK
  print *,'X coordinates'
  print 102,x(:)
  print *,'Y coordinates'
  print 102,y(:)
  print *,'F matrix (1)'
  do j = 6, 1-hy, -1
    print 102,f(1-hx:6,j,1)
  end do
  print *,'F matrix (nk)'
  do j = 6, 1-hy, -1
    print 102,f(1-hx:6,j,NK)
  end do
  print *,'MONO: limit (min, max)'
  do i = 1 , NP
    i0 = x(i)
    i0 = max(1,min(ni-3,i0-1))
    j0 = y(i)
    j0 = max(1,min(nj-3,j0-1))
!    print *, i0, j0, x(i), y(i)
    xmin(i)  = minval(f(i0:i0+3,j0:j0+3,1))
    xmax(i)  = maxval(f(i0:i0+3,j0:j0+3,1))
    xmin2(i) = minval(f(i0:i0+3,j0:j0+3,NK))
    xmax2(i) = maxval(f(i0:i0+3,j0:j0+3,NK))
  end do
  print 102,xmin(:) , xmin2(:)
  print 102,xmax(:) , xmax2(:)
  print *,"MONO: expected"
  print 102 , max( xmin(:),min( xmax(:),FXY(x(:),y(:),1 ))), &
              max(xmin2(:),min(xmax2(:),FXY(x(:),y(:),NK)))
!  print 102,x(:)+y(:)+1, x(:)+y(:)+NK
  print *,' got'
  print 102,r(1,1:NP),r(NK,1:NP)
  print *,' delta'
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
  end do
!  print 102,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),NK),f(nint(x(2)),nint(y(2)),NK)
!  print 102,x(1)+y(1)+1,x(2)+y(2)+1,x(1)+y(1)+NK,x(2)+y(2)+NK
  print *,"DIRECT: expected"
  print 102,FXY(x(:),y(:),1), FXY(x(:),y(:),NK)
  print *,' got'
  print 102,r(1,1:NP),r(NK,1:NP)
  print *,' gotuv'
  print 102,r1(1,1:NP),r1(NK,1:NP)
  print 102,r2(1,1:NP),r2(NK,1:NP)
  print *,' delta'
  print 102,(r(1,1:NP)-FXY(x(:),y(:),1)),(r(NK,1:NP)-FXY(x(:),y(:),NK))
  print 103,(r(1,1:NP)-FXY(x(:),y(:), 1))  / FXY(x(:),y(:), 1) , &
            (r(NK,1:NP)-FXY(x(:),y(:),NK))  / FXY(x(:),y(:),NK)
  r1 = 999.999
  r2 = 999.999
  do i = 1 , NP
    call intrp_bicub_yx_vec( f(1,1,1), f2(1,1,1), r1(1,i), r2(1,i), nidim, ninjdim, NK, x(1+mod(i,NP)), y(1+mod(i,NP)), x(i), y(i), s )
  end do
  print *,' vec'
  print 102,r1(1,1:NP),r1(NK,1:NP)
  print 102,r2(1,1:NP),r2(NK,1:NP)

102 format(16F15.6)
103 format(16E15.2)
end
#endif
