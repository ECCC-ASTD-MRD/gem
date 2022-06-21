/* RKL - RPN kernel Library for C and FORTRAN programming
 * Copyright (C) 2019  Recherche en Prevision Numerique
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

//****P* librkl/RPN kernel library interpolators
// DESCRIPTION
//
// Set of 3D interpolators,
// C and Fortran callable (C name and Fortran name are the same)
//
// interpolate the value at point (px, py, pz) fron a 3 Dimensional real(float) array
// Fortran : real, dimension(dimi,dimj,dimk) :: array
// C       : float array[dimk][dimj][dimi]
// (Fortran storage order is expected)
//
// the interpolation routines expect to receive the address of point
// array(i[+offseti],1[+offsetj],1)
// and will interpolate to get the value at point (px,py,pz)
// px, py, pz are expressed in "index" space, i.e.
// array(1.5,  1,  1) is halfway between array(1,1,1) and array(2,1,1)
// array(1  ,1.5,  1) is halfway between array(1,1,1) and array(1,2,1)
// array(1  ,1  ,2.5) is halfway between array(1,1,2) and array(1,1,3)
// array(2.5,3.5,4.5) is halfway between array(2,3,4) and array(3,4,5)
//
// constant spacing between array points is assumed for the first 2 dimensions.
//
// along the third dimension, non constant spacing is assumed, the positions along
// said third dimension are supplied via an opaque descriptor created with a call
// to the Vsearch_setup family of routines (returning a C pointer to said descriptor)
//
// 4 point lagrange polynomials are used for cubic interpolations. interpolation is always
// cubic along the first 2 dimensions, while along the third dimension linear interpolation
// is used in the first 2 (pz <= 2)  and the last 2 (pz >= nk-1) intervals, and cubic
// interpolation is used in the other intervals (2.0 < pz < nk - 1)
//
// there is a also monotonic interpolation option.
// in that case, 4 results are returned
// - normal interpolation result
// - tri-linear interpolation result
// - min and max values for the inner 2x2x2 set of values used for tri-linear interpolation
//
//     Name                        Purpose
//
// Vsearch_setup             : simple setup
// Vsearch_setup_plus        : Vsearch_setup + x and y direction offsets
// Tricublin_zyx1_n          : simple interpolator, source is a single variable in a 3D array
// Tricublin_zyx1_n_m        : source is a group of m  3D arrays
// Tricublin_zyx3_n          : source is a 3 variable ( dimension(3, ....) ) 3D array
// Tricublin_mono_zyx_n      : monotonic version of Tricublin_zyx1_n
// Tricublin_mono_zyx_n_m    : monotonic version of Tricublin_zyx1_n_m
//
// Tricublin_zyx1_n_m, Tricublin_zyx3_n, Tricublin_mono_zyx_n_m save computing cycles by
// reusing the computed interpolation coefficients for multiple arrays. The interpolation
// coefficients are dependent upon px, py, pz.
//
// the horizontal (first 2 dimensions) interpolation is bicubic.
// (an extra linear interpolation is performed in the monotonic case)
//
// the vertical (3rd dimension) interpolation is normally cubic, except when the target
// lies between the first 2 or the last 2 levels, in which case it is linear.
//
// the monotonic interpolators provides 4 outputs
// - the result of the centered tri-cubic or bi-cubic+linear interpolation (4x4x4 or 4x4x2)
// - the result of a  centered tri-linear interpolation (2x2x2)
// - the minimum of the 2x2x2 sub-cube used for tri-linear interpolation
// - the maximum of the 2x2x2 sub-cube used for tri-linear interpolation
//
// the interpolation routines assume that the address of the source array is the point
// px = 1.0, py = 1.0, pz = 1.0 which is why f(1,1,1) is normally passed to these routines.
// if this is not the case, offseti and offsetj are used to compensate
//
// EXAMPLES
//
//  use ISO_C_BINDING
//  include "tricublin_f90.inc"
//  type(C_PTR) :: lv, mv                            ! C pointer opaque descriptor object
//  real(kind=REAL64), dimension(NK) :: levels                  ! list of positions for the third dimension
//  real, dimension(-2:NI+4,-1:NJ+3,NK), target   :: f, f1, f2, f3   ! single valued source arrays
//  real, dimension(3,-2:NI+4,-1:NJ+3,NK), target :: f3              ! triple valued source array
//  type(C_PTR), dimension(3) :: f123
//  real, dimension(NPOINTS) :: d, l, mi, ma         ! normal, linear, min, max
//  real, dimension(3,NPOINTS) :: d3, l3, mi3, ma3   ! normal, linear, min, max
//  real, dimension(3,NPOINTS) :: pxpypz             ! (px,py,pz) triplets
//  ......                                           ! levels(1:NK) = ....
//  lv = vsearch_setup(levels, NK, NI+7, NJ+5)       ! NI+7, NJ+5 to account for halos
//  mv = vsearch_setup_plus(levels, NK, NI+7, NJ+5, 3, 2) ! passing f instead of f(1,1,1)
//  ......                                           ! pxpypz(1:3,1:NPOINTS) = ...., [1.0,1.0,1.0] points to  to f(1,1,1)
//  call tricublin_zyx1_n(d ,f(1,1,1)   ,pxpypz,lv,NPOINTS)      ! single interpolation
//  call tricublin_zyx1_n(d ,f          ,pxpypz,mv,NPOINTS)      ! with offset compensation (halos)
//  call tricublin_zyx3_n(d3,f3(1,1,1,1),pxpypz,lv,NPOINTS)      ! triple interpolation
//  f123(1) = C_LOC(f1(1,1,1))
//  f123(2) = C_LOC(f2(1,1,1))
//  f123(3) = C_LOC(f3(1,1,1))
//  call tricublin_zyx1_n_m(d3,f123,pxpypz,lv,NPOINTS,3)         ! 3 single interpolations
//  call tricublin_mono_zyx_n(d,l,mi,ma,f(1,1,1),pxpypz,lv,NPOINTS)      ! monotonic interpolation, 1 array
//  call tricublin_mono_zyx_n_m(d3,l3,mi3,ma3,f123,pxpypz,lv,NPOINTS,3)  ! monotonic interpolation, 3 arrays
//****

#define MAX(a,b) ((a) > (b)) ? (a) : (b)
#define MIN(a,b) ((a) < (b)) ? (a) : (b)

#if defined(__AVX2__) && defined(__x86_64__)
#include <immintrin.h>
#endif

#include <stdint.h>
#include <stdlib.h>

#if defined(TIMING)
#include <sys/time.h>
uint64_t Nanocycles(void) {
// uint64_t rdtscp_(void) {   // version "in order" avec "serialization"
//   uint32_t lo, hi;
//   __asm__ volatile ("rdtscp"
//       : /* outputs */ "=a" (lo), "=d" (hi)
//       : /* no inputs */
//       : /* clobbers */ "%rcx");
// //   __asm__ volatile ("mfence");
//   return (uint64_t)lo | (((uint64_t)hi) << 32);
  struct timeval t;
  uint64_t i;
  int status;
  status = gettimeofday(&t,NULL);
  i = t.tv_sec;
  i = i * 1000000;
  i = i + t.tv_usec;
  i = i * 1000;      // nanoseconds
  return i;
}
#endif

#if defined(__AVX2__) && defined(__x86_64__)
static double cp133 =  1.0/3.0;
#endif
static double cp167 =  1.0/6.0;
static double cm167 = -1.0/6.0;
static double cm5 = -0.5;
static double cp5 = 0.5;
static double one = 1.0;
static double two = 2.0;

typedef struct{
  float px;    // position along x in index space
  float py;    // position along y in index space
  float pz;    // position along z in index space
//   float z;     // absolute position along z
} pxpypz;

typedef struct{
  double *z;       // table for levels (nk doubles)
  double *ocz;     // pointer to inverse of coefficient denominators [4*nk doubles]
  uint32_t ni;     // distance between rows
  uint32_t nj;     // number of rows
  uint32_t nk;     // number of levels (max 255)
  uint32_t nij;    // ni*nj, distance between levels
  uint32_t offi;   // offset to add to i position (e.g. global to local grid remapping)
  uint32_t offj;   // offset to add to j position (e.g. global to local grid remapping)
}ztab;

static void cubic_coeff_denom(double *t, double *f, int n){
  int i, j;
#if defined(SIMD)
  __m256d cp1, cm1;
  static double onep[4] = { 1.0,  1.0,  1.0,  1.0};
  static double onem[4] = {-1.0, -1.0, -1.0, -1.0};
  cp1 = _mm256_loadu_pd(onep);
  cm1 = _mm256_loadu_pd(onem);
#endif
  t += 4;    // will store no coefficients for first interval (no point below)
  j = 0;
#if defined(SIMD)
  {
    __m256d a, b, c, d, ab, ac, ad, bc, bd, cd, v1, v2, v3, v4, t1, t2, t3, t4;
    for( ; j<n-6 ; j+=4){
      a  = _mm256_loadu_pd(f+0);
      b  = _mm256_loadu_pd(f+1);
      c  = _mm256_loadu_pd(f+2);
      d  = _mm256_loadu_pd(f+3);
      ab = _mm256_sub_pd(a,b);
      ac = _mm256_sub_pd(a,c);
      ad = _mm256_sub_pd(a,d);
      bc = _mm256_sub_pd(b,c);
      bd = _mm256_sub_pd(b,d);
      cd = _mm256_sub_pd(c,d);
      v1 = _mm256_mul_pd( ab ,  _mm256_mul_pd( ac , ad) );
      v2 = _mm256_mul_pd( ab ,  _mm256_mul_pd( bc , bd) );
      v3 = _mm256_mul_pd( ac ,  _mm256_mul_pd( bc , cd) );
      v4 = _mm256_mul_pd( ad ,  _mm256_mul_pd( bd , cd) );
      v1 = _mm256_div_pd(cp1 , v1);     // a1 a2 a3 a4   coefficient a for points 1 2 3 4
      v2 = _mm256_div_pd(cm1 , v2);     // b1 b2 b3 b4   coefficient b for points 1 2 3 4
      v3 = _mm256_div_pd(cp1 , v3);     // c1 c2 c3 c4   coefficient c for points 1 2 3 4
      v4 = _mm256_div_pd(cm1 , v4);     // d1 d2 d3 d4   coefficient d for points 1 2 3 4
      t1 = _mm256_unpacklo_pd(v1,v2);   // a1 b1 a3 b3
      t2 = _mm256_unpackhi_pd(v1,v2);   // a2 b2 a4 b4
      t3 = _mm256_unpacklo_pd(v3,v4);   // c1 d1 c3 d3
      t4 = _mm256_unpackhi_pd(v3,v4);   // c2 d2 c4 d4
      _mm_storeu_pd(t+ 0, _mm256_extractf128_pd(t1,0));   // a1 b1  coefficients a b c d for point 1
      _mm_storeu_pd(t+ 2, _mm256_extractf128_pd(t3,0));   // c1 d1
      _mm_storeu_pd(t+ 4, _mm256_extractf128_pd(t2,0));   // a2 b2  coefficients a b c d for point 2
      _mm_storeu_pd(t+ 6, _mm256_extractf128_pd(t4,0));   // c2 d2
      _mm_storeu_pd(t+ 8, _mm256_extractf128_pd(t1,1));   // a3 b3  coefficients a b c d for point 3
      _mm_storeu_pd(t+10, _mm256_extractf128_pd(t3,1));   // c3 d3
      _mm_storeu_pd(t+12, _mm256_extractf128_pd(t2,1));   // a4 b4  coefficients a b c d for point 4
      _mm_storeu_pd(t+14, _mm256_extractf128_pd(t4,1));   // c4 d4
      f++;
      t += 16;
    }
  }
#endif
  {
  double a, b, c, d, ab, ac, ad, bc, bd, cd;
    for( ; j<n-3 ; j++){     // will store no coefficients for the last 2 intervals
      a = f[0];
      b = f[1];
      c = f[2];
      d = f[3];
      ab = a - b;
      ac = a - c;
      ad = a - d;
      bc = b - c;
      bd = b - d;
      cd = c - d;
      t[0]  =  1.0/(ab*ac*ad);
      t[1]  = -1.0/(ab*bc*bd);
      t[2]  =  1.0/(ac*bc*cd);
      t[3]  = -1.0/(ad*bd*cd);
      t += 4;
      f++;
    }
  }
}

// triple product used for Lagrange cubic polynomials coefficients in the non constant case
#define TRIPRD(x,a,b,c) ((x-a)*(x-b)*(x-c))

// inverse denominators r from positions a b c d
static inline void denominators(double *r, double a, double b, double c, double d){
  r[0] = 1.0 / TRIPRD(a,b,c,d);
  r[1] = 1.0 / TRIPRD(b,a,c,d);
  r[2] = 1.0 / TRIPRD(c,a,b,d);
  r[3] = 1.0 / TRIPRD(d,a,b,c);
}

#if defined(NEVER_TO_BE_TRUE)
  interface                                                                      !InTf!
//****f* librkl/vsearch_setup  (Fortran version)
// Synopsis
//    the function returns an opaque object of type C_PTR, to be passed to  the interpolators
//
//    levels      : real(kind=REAL64) array of levels, dimension(nk)
//    nk          : number of levels
//    ni, nj      : used to compute indexing into source array(s) for interpolators
//
//    ni is used to compute the distance between 2 array elements for which index j differs by 1
//    ni*nj is used to compute the distance between 2 array elements for which index k differs by 1
//
//    when vsearch_setup is used, positions will be assumed to be in local index space
//    (equivalent to offseti = 0 and offsetj = 0 in vsearch_setup_plus)
//
//    see interpolators for more information
//
// ARGUMENTS
    function vsearch_setup(levels, nk, ni, nj) result (ztab) bind(C,name='Vsearch_setup')     !InTf!
      import :: C_PTR, C_DOUBLE, C_INT                                           !InTf!
      real(C_DOUBLE), dimension(nk), intent(IN) :: levels                        !InTf!
      integer(C_INT), intent(IN), value :: nk, ni, nj                            !InTf!
      type(C_PTR) :: ztab                                                        !InTf!
//****
    end function vsearch_setup                                                   !InTf!
//****f* librkl/vsearch_setup_plus  (Fortran version)
// Synopsis
//    the function returns an opaque object of type C_PTR, to be passed to  the interpolators
//
//    levels      : real(kind=REAL64) array of levels, dimension(nk)
//    nk          : number of levels
//    ni, nj      : used to compute indexing into source array(s) for interpolators
//    offseti     : used to convert global position along i into local array position
//    offsetj     : used to convert global position along j into local array position
//
//    ni is used to compute the distance between 2 array elements for which index j differs by 1
//    ni*nj is used to compute the distance between 2 array elements for which index k differs by 1
//    effective px = px + offseti
//    effective py = py + offsetj
//
//    see interpolators for more information
//
// ARGUMENTS
    function vsearch_setup_plus(levels, nk, ni, nj, offseti, offsetj) result (ztab) bind(C,name='Vsearch_setup_plus')     !InTf!
      import :: C_PTR, C_DOUBLE, C_INT                                           !InTf!
      real(C_DOUBLE), dimension(nk), intent(IN) :: levels                        !InTf!
      integer(C_INT), intent(IN), value :: nk, ni, nj                            !InTf!
      integer(C_INT), intent(IN), value :: offseti, offsetj                      !InTf!
      type(C_PTR) :: ztab                                                        !InTf!
//****
    end function vsearch_setup_plus                                              !InTf!
  end interface                                                                  !InTf!
#endif
// allocate lookup table set and
// return pointer to filled table set
// targets are expected to be positive, and monotonically increasing
//****f* librkl/Vsearch_setup
// Synopsis
//
// ARGUMENTS
ztab *Vsearch_setup(double *targets, int nk, int ni, int nj)
//****
{
  int i;
  ztab *lv;
  double pad;
  int nkl = MAX(nk,4);

  lv = malloc(sizeof(ztab));
  if(lv == NULL) return NULL ;

  lv->z = malloc(nk * sizeof(double));
  if(NULL == lv->z){      // malloc failed
    free(lv);             // deallocate lv
    return NULL;
  }
  for(i=0 ; i<nk  ; i++) { lv->z[i] = targets[i] ; }   // z coordinate table

// for(i=0 ; i<nk  ; i++) { printf("%12.7f ",lv->z[i]) ; } ; printf("\n");
  // denominators for Lagrange cubic polynomials coefficients ( entries 1 to n-2 make sense )
  // entries 0 and n-1 are fudged
  lv->ocz = (double *) malloc(4 * nkl * sizeof(double));
  if(NULL == lv->ocz){    // malloc failed
    free(lv->z);          // deallocate z
    free(lv);             // deallocate lv
    return NULL;
  }
  if(nk >= 4){             // protect against nk < 4
    denominators( &(lv->ocz[0]) , targets[0], targets[1], targets[2], targets[3]);                     // level 0 coeffs are normally not used
    for(i=1 ; i<nk - 2   ; i++) {
      denominators( &(lv->ocz[4*i]) , lv->z[i-1], lv->z[i  ], lv->z[i+1], lv->z[i+2]);
    }
    denominators( &(lv->ocz[4*(nk-2)]) , targets[nk-4], targets[nk-3], targets[nk-2], targets[nk-1]);  // level nk-2 coeffs are normally not used
    denominators( &(lv->ocz[4*(nk-1)]) , targets[nk-4], targets[nk-3], targets[nk-2], targets[nk-1]);  // level nk-1 coeffs are normally not used
  }else{
    for(i=0 ; i<nkl*4 ; i++) lv->ocz[i] = 1.0;            // not used, set to 1.0 to avoid Floating Point errors later on
  }

  lv->ni = ni;           // nb of points along x
  lv->nj = nj;           // nb of points along y
  lv->nk = nk;           // nb of points along z
  lv->nij = ni*nj;       // ni * nj
  lv->offi = 0;
  lv->offj = 0;
  return lv;             // return pointer to filled table
}
//****f* librkl/Vsearch_setup_plus
// Synopsis
//
// ARGUMENTS
ztab *Vsearch_setup_plus(double *targets, int nk, int ni, int nj, int offseti, int offsetj)
//****
{
  ztab *lv;
  lv = Vsearch_setup(targets, nk, ni, nj);
  if(lv != NULL) {
    lv->offi = offseti;
    lv->offj = offsetj;
  }
  return lv;             // return pointer to filled table
}

static inline int Vcoef_pxyz4_inline(double *cxyz, int *offset, float px8, float py8, float pz8, ztab *lv){
  int ix, iy, iz, ijk, zlinear;
  double pxy[2], *base, *pos;
  double zza, zzb, zzc, zzd, zzab, zzcd, dz, px, py, pz;
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m256i v0, v1, vc, t, ttop, tbot, ttemp;
  __m128d vxy, vc1, vc5, vc3, vp6, xxmx, xxpx, xxm1, x5p1, x6m3, x5m1, x6m6, dtmp;
  __m128d vr0, vr1, vr2, vr3;
  __m128i itmp;
  int m0, m1;
#else
  int i, j;
#endif
    px = px8 ;            // fractional index positions along x, y, z (float to double)
    py = py8 ;
    pz = pz8 ;
// printf("px py pz = %f %f %f\n",px,py,pz);
    ix = px ;
    px = px - ix;                   // px is now deltax (fractional part of px)
    ix = ix + lv->offi;             // global to local grid remapping
    ijk = ix - 2;                   // x displacement (elements), ix assumed to always be >1 and < ni-1
// printf("ix = %d, ijk = %d\n",ix,ijk);
    cxyz[12] = 0.0;
    cxyz[13] = 1.0 - px;             // linear interpolation coefficients along x
    cxyz[14] = px;
    cxyz[15] = 0.0;

    iy = py ;
    py = py - iy;                   // py is now deltay (fractional part of py)
    iy = iy + lv->offj;             // global to local grid remapping
    ijk = ijk + (iy - 2) * lv->ni;  // add y displacement (rows), ix assumed to always be >1 and < nj-1
// printf("iy = %d, ijk = %d\n",iy,ijk);
    cxyz[16] = 1.0 - py;             // linear interpolation coefficients along y
    cxyz[17] = py;

    iz = pz ;
    if(iz<1) iz = 1;
    if(iz>lv->nk-1) iz = lv->nk-1;  // iz < 1 or iz > nk-1 will result in linear extrapolation
    dz = pz - iz;                   // dz is now "fractional" part of pz  (may be <0 or >1 if extrapolating)
    ijk = ijk + (iz -1) * lv->nij;  // add z displacement (2D planes)
    cxyz[18] = 1.0 - dz;             // linear interpolation coefficients along z
    cxyz[19] = dz;

    iz--;                           // iz needs to be in "origin 0" (C index from Fortran index)
    zlinear = (iz - 1) | (lv->nk - 3 - iz);
    zlinear >>= 31;                  // nonzero only if iz < 1 or iz > nk -3 (top and bottom intervals)
    if(! zlinear) ijk = ijk - lv->nij;  // not the linear case, go down one 2D plane to get lower left corner of 4x4x4 cube
    *offset = ijk;
// printf("iz = %d, ijk = %d\n",iz,ijk);
    // now we can compute the coefficients along z using iz and dz
    if(zlinear){
      cxyz[ 8] = 0.0;                    // coefficients for linear interpolation along z
      cxyz[ 9] = 1.0 - dz;
      cxyz[10] = dz;
      cxyz[11] = 0.0;
    }else{
      base  = &(lv->ocz[4*iz]);  // precomputed inverses of denominators
      pos   = &(lv->z[iz]);
      pz  = dz * pos[1] + (1.0 - dz) * pos[0];   // pz is now an absolute position
      zza = pz - pos[-1] ; zzb = pz - pos[0] ; zzc = pz - pos[1] ; zzd = pz - pos[2] ;
      zzab = zza * zzb ; zzcd = zzc * zzd;
      // printf("target = %8.5f, dz = %8.5f, levels = %8.5f %8.5f\n",*PZ,dz,pos[0],pos[1]);
      // printf("target = %8.5f, base = %8.5f %8.5f %8.5f %8.5f\n",pz,base[0],base[1],base[2],base[3]);
      // printf("dz     = %8.5f, zz   = %8.5f %8.5f %8.5f %8.5f\n",dz,zza,zzb,zzc,zzd);
      // printf("iz     = %d, lv->odz[iz] = %8.5f\n",iz,lv->odz[iz]);
      cxyz[ 8] = zzb * zzcd * base[0];   //   cxyz[16] = TRIPRD(pz,pos[1],pos[2],pos[3]) * base[0];
      cxyz[ 9] = zza * zzcd * base[1];   //   cxyz[17] = TRIPRD(pz,pos[0],pos[2],pos[3]) * base[1];
      cxyz[10] = zzd * zzab * base[2];   //   cxyz[18] = TRIPRD(pz,pos[0],pos[1],pos[3]) * base[2];
      cxyz[11] = zzc * zzab * base[3];   //   cxyz[19] = TRIPRD(pz,pos[0],pos[1],pos[2]) * base[3];
    }

// printf("pz,dz,iz,nk,linear = %8.5f %8.5f %d %d %d\n",pz,dz,iz,lv->nk,linear);
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
    pxy[0] = px;                           // vector of length 2 to do x and y in one shot
    pxy[1] = py;
    vxy  = _mm_loadu_pd(pxy);              // coefficients for x and y will be interleaved

    vc1  = _mm_set1_pd(one);               //  one    (1.0)
    vc5  = _mm_set1_pd(cp5);               //  cp5    (0.5)
    vc3  = _mm_set1_pd(cp133);             //  cp133  (1/3)
    vp6  = _mm_set1_pd(cp167);             //  cp167  (1/6)

    // alternate formula to compute coefficients taking advantage of independent FMAs
    //   pa = cm167*xy*(xy-one)*(xy-two)     = xy*(xy-one) * (-cp167)*(xy-two)  =  (xy*xy  - xy)   * (-(cp167*xy) + cp133)
    //   pb = cp5*(xy+one)*(xy-one)*(xy-two) = cp5*(xy-two) * (xy+one)*(xy-one) =  (cp5*xy - one)  * (xy*xy       - one)
    //   pc = cm5*xy*(xy+one)*(xy-two)       = xy*(xy+one) * (-cp5)*(xy-two)    =  (xy*xy  + xy)   * (-(cp5*xy)   + one)
    //   pd = cp167*xy*(xy+one)*(xy-one)     = xy*(xy+one) * cp167*(xy-one)     =  (xy*xy  + xy)   * (cp167*xy    - cp167)
    // STEP1 : 7 independent terms using 3 different FMAs  a*b+c, a*b-c, (-a*b)+c
    xxmx = _mm_fmsub_pd(vxy,vxy,vxy);      //  xy*xy       - xy
    x6m3 = _mm_fnmadd_pd(vxy,vp6,vc3);     //  -(cp167*xy) + cp133
    x5p1 = _mm_fmsub_pd(vc5,vxy,vc1);      //  cp5*xy      - one
    xxm1 = _mm_fmsub_pd(vxy,vxy,vc1);      //  xy*xy       - one
    xxpx = _mm_fmadd_pd(vxy,vxy,vxy);      //  xy*xy       + xy
    x5m1 = _mm_fnmadd_pd(vxy,vc5,vc1);     //  -(cp5*xy)   + one
    x6m6 = _mm_fmsub_pd(vxy,vp6,vp6);      //  cp167*xy    - cp167

    // STEP2 : multiply STEP 1 terms (4 independent operations)
    vr0  = _mm_mul_pd(xxmx,x6m3);          // coefficients for x and y are interleaved
    vr1  = _mm_mul_pd(x5p1,xxm1);
    vr2  = _mm_mul_pd(xxpx,x5m1);
    vr3  = _mm_mul_pd(xxpx,x6m6);

    // final unshuffle to separate even terms (px) and odd terms (py) before storing them (independent operations)
    _mm_storeu_pd(cxyz   ,_mm_unpacklo_pd(vr0,vr1));  // cxyz[ 0: 1] = cx[0], cx[1]
    _mm_storeu_pd(cxyz+ 2,_mm_unpacklo_pd(vr2,vr3));  // cxyz[ 2: 3] = cx[2], cx[3]
    _mm_storeu_pd(cxyz+ 8,_mm_unpackhi_pd(vr0,vr1));  // cxyz[ 8: 9] = cy[0], cy[1]
    _mm_storeu_pd(cxyz+10,_mm_unpackhi_pd(vr2,vr3));  // cxyz[10:11] = cy[2], cy[3]

#else
//     cxyz[ 0] = (px*px  - px)   * (-(cp167*px) + cp133);
//     cxyz[ 1] = (cp5*px - one)  * (px*px       - one);
//     cxyz[ 2] = (px*px  + px)   * (-(cp5*px)   + one);
//     cxyz[ 3] = (px*px  + px)   * (cp167*px    - cp167);
//
//     cxyz[ 8] = (py*py  - py)   * (-(cp167*py) + cp133);
//     cxyz[ 9] = (cp5*py - one)  * (py*py       - one);
//     cxyz[10] = (py*py  + py)   * (-(cp5*py)   + one);
//     cxyz[11] = (py*py  + py)   * (cp167*py    - cp167);

    cxyz[ 0] = cm167*px*(px-one)*(px-two);        // coefficients for cubic interpolation along x
    cxyz[ 1] = cp5*(px+one)*(px-one)*(px-two);
    cxyz[ 2] = cm5*px*(px+one)*(px-two);
    cxyz[ 3] = cp167*px*(px+one)*(px-one);

    cxyz[ 4] = cm167*py*(py-one)*(py-two);        // coefficients for cubic interpolation along y
    cxyz[ 5] = cp5*(py+one)*(py-one)*(py-two);
    cxyz[ 6] = cm5*py*(py+one)*(py-two);
    cxyz[ 7] = cp167*py*(py+one)*(py-one);
#endif
  return zlinear;  // linear / cubic flag
}

// Fortran dimensions: f1(NI,NJ,NK)
// NI      : length of a line
// NINJ    : length of a plane
// zlinear : zero if cubic interpolation along z, non zero if linear interpolation along z
// f1      : address of lower left corner of the 4 x 4 x 4 or 4 x 4 x 2 interpolation box
// interpolation is done along z, then along y, then along x
// values are interpolated using the cubic/linear polynomial coefficients
// pxyz(8) and pxyz(11) both = 0.0 when linear interpolation along z
// pxyz(9) = 1.0 - dz, pxyz(10) = dz are expected in the linear case
// in the z linear case, planes 0 and 1 are the same as are planes 2 and 3
// static inline void Tricublin_zyxf_beta(float *d, float *f1, double *px, double *py, double *pz, int NI, int NINJ){
static inline void Tricublin_zyxf1_inline(float *d, float *f1, double *pxyz, int NI, int NINJ, int zlinear){
  int ni = NI;
  int ninj = NINJ;
  int ninjl;    // ninj (cubic along z) or 0 (linear along z)
  int ni2, ni3;
  double *px = pxyz;
  double *py = pxyz+4;
  double *pz = pxyz+8;
#if defined(__AVX2__) && defined(__x86_64__)
  float *s1;
  int32_t *l = (int32_t *)d;
  __m256d cz0, cz1, cz2, cz3, cy, cx;
  __m256d za;
  __m256d ya;
  __m128d va;
  __m128  ta;
#else
  double va4[4], vb4[4], vc4[4], vd4[4], dst[4];
  int i, ninj2, ninj3;
#endif

  ni2 = ni + ni;
  ni3 = ni2 + ni;
  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be set to 0
  if(zlinear) ninjl = 0;

#if defined(__AVX2__) && defined(__x86_64__)
  // ==== interpolation along Z, vector length is 12 (3 vectors of length 4 per plane) ====
  cz0 = _mm256_broadcast_sd(pz  );   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm256_broadcast_sd(pz+1);
  cz2 = _mm256_broadcast_sd(pz+2);
  cz3 = _mm256_broadcast_sd(pz+3);

  s1 = f1;      // row 0, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_mul_pd(za,cy);

  s1 = f1 + ni;        // row 1, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+1);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_fmadd_pd(za,cy,ya);

  s1 = f1 + ni2;     // row 2, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+2);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_fmadd_pd(za,cy,ya);

  s1 = f1 + ni3;     // row 3, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+3);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_fmadd_pd(za,cy,ya);

  // interpolation along x
  cx = _mm256_loadu_pd(px) ;     // make sure to use "unaligned" load to avoid segfault
  ya = _mm256_mul_pd(ya, cx );   // multiply by px
  // use folding to sum 4 terms
  va = _mm_add_pd( _mm256_extractf128_pd(ya,0) , _mm256_extractf128_pd(ya,1) ); va = _mm_hadd_pd(va,va);
  ta = _mm_cvtpd_ps(va); d[0] = _mm_cvtss_f32(ta) ;// l[0] = _mm_extract_epi32((__m128i) ta, 0);  // convert results to float and store
#else
  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  // field 1
  for (i=0 ; i<4 ; i++){
    va4[i] = f1[i    ]*pz[0] + f1[i    +ninjl]*pz[1] +  f1[i    +ninj2]*pz[2] + f1[i    +ninj3]*pz[3];
    vb4[i] = f1[i+ni ]*pz[0] + f1[i+ni +ninjl]*pz[1] +  f1[i+ni +ninj2]*pz[2] + f1[i+ni +ninj3]*pz[3];
    vc4[i] = f1[i+ni2]*pz[0] + f1[i+ni2+ninjl]*pz[1] +  f1[i+ni2+ninj2]*pz[2] + f1[i+ni2+ninj3]*pz[3];
    vd4[i] = f1[i+ni3]*pz[0] + f1[i+ni3+ninjl]*pz[1] +  f1[i+ni3+ninj2]*pz[2] + f1[i+ni3+ninj3]*pz[3];
    dst[i] = va4[i]*py[0] + vb4[i]*py[1] + vc4[i]*py[2] + vd4[i]*py[3];
  }
  d[0] = dst[0]*px[0] + dst[1]*px[1] + dst[2]*px[2] + dst[3]*px[3];
#endif
}

// Fortran dimensions: d(n) , f1(*), pxyz(3,n), lv(*)
// d       : results
// f1      : source array (dimensions found in opaque object lv)
//           first index along each dimension assumed to be 1
// pxyz    : positions along x, y, z , in "index" space
// lv      : opaque object containing vertical description (see Vsearch_setup)
// n       : number of points
#if defined(NEVER_TO_BE_TRUE)
  interface                                                                      !InTf!
//****f* librkl/tricublin_zyx1_n  (Fortran version)
// Synopsis
//    f         real input array of dimension(xmin:xmax,ymin,ymax,nk)
//              f(1,1,1) MUST be passed to this function
//              xmax MUST be xmin+ni-1, ymax MUST be ymin+nj-1  (see setup routines)
//    pxyz      real array, dimension(3,n), contains the target interpolation positions
//              positions are in global index space
//              pxyz(1,i) is the position along x of point i (i = 1..n)
//              pxyz(2,i) is the position along y of point i
//              pxyz(3,i) is the position along z of point i
//    d         real array, dimension(n), d(i) is the interpolated value at point i
//    lv        opaque object obtained from the setup routines
//    n         number of points to interpolate
// ARGUMENTS
    subroutine tricublin_zyx1_n(d,f,pxyz,lv,n) bind(C,name='Tricublin_zyx1_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(n), intent(OUT)   :: d                                     !InTf!
      real, dimension(*), intent(IN)    :: f                                     !InTf!
      real, dimension(3,n), intent(IN)  :: pxyz                                  !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n                                     !InTf!
//****
    end subroutine tricublin_zyx1_n                                              !InTf!
//****f* librkl/tricublin_zyx1_n_m  (Fortran version)
// ARGUMENTS
    subroutine tricublin_zyx1_n_m(d,f,pxyz,lv,n,m) bind(C,name='Tricublin_zyx1_n_m')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d                                     !InTf!
      type(C_PTR), dimension(*), intent(IN)    :: f                              !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n, m                                  !InTf!
//****
    end subroutine tricublin_zyx1_n_m                                            !InTf!
//****f* librkl/tricublin_zyx1_m_n  (Fortran version)
// ARGUMENTS
    subroutine tricublin_zyx1_m_n(d,f,pxyz,lv,n,m) bind(C,name='Tricublin_zyx1_m_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d                                     !InTf!
      type(C_PTR), dimension(*), intent(IN)    :: f                              !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n, m                                  !InTf!
//****
    end subroutine tricublin_zyx1_m_n                                            !InTf!
void Tricublin_zyx1_p(float **dp, float **fs, pxpypz *pxyz,  ztab *lv, int n, int m)
//****f* librkl/tricublin_zyx1_p  (Fortran version)
// ARGUMENTS
    subroutine tricublin_zyx1_p(d,f,pxyz,lv,n,m) bind(C,name='Tricublin_zyx1_p')   !InTf!
      import :: C_PTR                                                            !InTf!
      type(C_PTR), dimension(*), intent(IN)    :: f, d                           !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n, m                                  !InTf!
//****
    end subroutine tricublin_zyx1_p                                            !InTf!
//****f* librkl/tricublin_mono_zyx_n  (Fortran version)
// ARGUMENTS
    subroutine tricublin_mono_zyx_n(d,l,mi,ma,f,pxyz,lv,n) bind(C,name='Tricublin_mono_zyx_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d, l, mi, ma                          !InTf!
      real, dimension(*), intent(IN)    :: f                                     !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n                                     !InTf!
//****
    end subroutine tricublin_mono_zyx_n                                          !InTf!
//****f* librkl/tricublin_mono_zyx_n_m  (Fortran version)
// ARGUMENTS
    subroutine tricublin_mono_zyx_n_m(d,l,mi,ma,f,pxyz,lv,n,m) bind(C,name='Tricublin_mono_zyx_n_m')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d, l, mi, ma                          !InTf!
      type(C_PTR), dimension(*), intent(IN)    :: f                              !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n, m                                  !InTf!
//****
    end subroutine tricublin_mono_zyx_n_m                                        !InTf!
//****f* librkl/tricublin_mono_zyx_m_n  (Fortran version)
// ARGUMENTS
    subroutine tricublin_mono_zyx_m_n(d,l,mi,ma,f,pxyz,lv,n,m) bind(C,name='Tricublin_mono_zyx_m_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d, l, mi, ma                          !InTf!
      type(C_PTR), dimension(*), intent(IN)    :: f                              !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n, m                                  !InTf!
//****
    end subroutine tricublin_mono_zyx_m_n                                        !InTf!
//****f* librkl/tricublin_zyx3_n  (Fortran version)
// ARGUMENTS
    subroutine tricublin_zyx3_n(d,f,pxyz,lv,n) bind(C,name='Tricublin_zyx3_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d                                     !InTf!
      real, dimension(*), intent(IN)    :: f                                     !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n                                     !InTf!
//****
    end subroutine tricublin_zyx3_n                                              !InTf!
  end interface                                                                  !InTf!
#endif

// process n points
// for each point 1 value from ixyz, 3 values from pxyz are used (1 element)
// it is ASSUMED that along X and Y interpolation will always be CUBIC
// ( 2 <= px < "ni"-1 ) ( 2 <= py < "nj"-1 )
// pz < 2 and pz >= nk - 1 will induce linear interpolation or extrapolation
//****f* librkl/Tricublin_zyx1_n
// Synopsis
//
// ARGUMENTS
void Tricublin_zyx1_n(float *d, float *f1, pxpypz *pxyz,  ztab *lv, int n)
//****
{
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    Tricublin_zyxf1_inline(d, f1 + ixyz, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
    d++;         // next result
    pxyz += 1;   // next set of positions
  }
}

//****f* librkl/Tricublin_zyx1_n_m
// Synopsis
//
// ARGUMENTS
void Tricublin_zyx1_n_m(float *d, float **fs, pxpypz *pxyz,  ztab *lv, int n, int m)
//****
{  // multiple field version with sources pointer lists, d[n][m]  ( Fortran dimension(m,n))
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int i;
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    for(i=0 ; i< m ; i++){
      Tricublin_zyxf1_inline(d, fs[i] + ixyz, cxyz, lv->ni, lv->nij, zlinear);       // interpolate
      d++;   // next result
    }
    pxyz += 1;   // next set of positions
  }
}

//****f* librkl/Tricublin_zyx1_m_n
// Synopsis
//
// ARGUMENTS
void Tricublin_zyx1_m_n(float *dm, float **fs, pxpypz *pxyz,  ztab *lv, int n, int m)
//****
{  // multiple field version with sources pointer lists, d[m][n]  ( Fortran dimension(n,m))
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int i, j;
  int n0 = n;
  float *d0;
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  for(j=0 ; j<n ; j++){                                                           // loop over n points
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    d0 = dm;                                                                      // base address for results
    for(i=0 ; i< m ; i++){
      Tricublin_zyxf1_inline(d0, fs[i] + ixyz, cxyz, lv->ni, lv->nij, zlinear);       // interpolate
//       printf("DEBUG: i = %d, n= %d, offset = %d, value = %f \n",i,n,d0-dm,*d0);
      d0 = d0 + n0;                                                                   // next result
    }
    dm++;        // next set of results
    pxyz += 1;   // next set of positions
//     return;
  }
}

//****f* librkl/Tricublin_zyx1_p
// Synopsis
//
// ARGUMENTS
void Tricublin_zyx1_p(float **dp, float **fs, pxpypz *pxyz,  ztab *lv, int n, int m)
//****
{  // multiple field version with results and sources pointer lists
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int i, j;
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  for(j=0 ; j<n ; j++){                                                                  // loop over n points
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);         // compute coefficients
    for(i=0 ; i< m ; i++){                                                               // loop over m variables
      Tricublin_zyxf1_inline(dp[i] + j, fs[i] + ixyz, cxyz, lv->ni, lv->nij, zlinear);   // interpolate
    }
    pxyz += 1;   // next set of positions
  }
}

static inline void Tricublin_zyx_mm_d_inline(float *d, float *lin, float *min, float *max, float *f1, double *cxyz, int NI, int NINJ, int zlinear){
  int ni = NI;
  int ninj = NINJ;
  int ninjl;    // ninj (cubic along z) or 0 (linear along z)
  int ni2, ni3;
  double *px, *py, *pz;
#if defined(__AVX2__) && defined(__x86_64__)
  float *s1;
//   int32_t *l = (int32_t *)d;
  __m256d cz0, cz1, cz2, cz3, cy, cx, cl, cyl;
  __m256d za, zt;
  __m256d mi, ma;  // min, max
  __m256d ya;
  __m256d vzl, vyl; // accumulators for linear interpolation
  __m128d va, vmi, vma;
  __m128  ta;
  double dst[4];
#else
  double va4[4], vb4[4], vc4[4], vd4[4], dst[4], dsl[4];
  int i, ninj2, ninj3;
  double ma, mi;
#endif

  px = cxyz;
  py = cxyz + 4;
  pz = cxyz + 8;
  ni2 = ni + ni;
  ni3 = ni2 + ni;
  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be set to 0
  if(zlinear) ninjl = 0;

#if defined(__AVX2__) && defined(__x86_64__)
  // ==== interpolation along Z, vector length is 4 (1 vector1 of length 4 per plane) ====
  cz0 = _mm256_broadcast_sd(pz  );   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm256_broadcast_sd(pz+1);
  cz2 = _mm256_broadcast_sd(pz+2);
  cz3 = _mm256_broadcast_sd(pz+3);

  s1 = f1;      // row 0, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py);

  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 0
  za  = _mm256_mul_pd(zt, cz0);
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 1
  mi = zt ; ma = zt;                            // initialize min, max to first values
  za  = _mm256_fmadd_pd(zt, cz1, za);
  s1 += ninj;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 2
  mi = _mm256_min_pd(mi, zt); ma = _mm256_max_pd(ma, zt);
  za  = _mm256_fmadd_pd(zt, cz2, za);
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 3
  za  = _mm256_fmadd_pd(zt, cz3, za);

  ya  = _mm256_mul_pd(za,cy);  // row 0

  s1 = f1 + ni;        // row 1, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+1);
  cyl = _mm256_broadcast_sd(cxyz+16);

  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 0
  za  = _mm256_mul_pd(zt, cz0);
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 1
  mi = _mm256_min_pd(mi, zt); ma = _mm256_max_pd(ma, zt);
  za  = _mm256_fmadd_pd(zt, cz1, za);
  cl  = _mm256_broadcast_sd(cxyz+18);           // linear coefficient along z
  vzl = _mm256_mul_pd(zt, cl);
  s1 += ninj;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 2
  mi = _mm256_min_pd(mi, zt); ma = _mm256_max_pd(ma, zt);
  za  = _mm256_fmadd_pd(zt, cz2, za);
  cl  = _mm256_broadcast_sd(cxyz+19);           // linear coefficient along z
  vzl = _mm256_fmadd_pd(zt, cl, vzl);  // linear interpolation along z
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 3
  za  = _mm256_fmadd_pd(zt, cz3, za);

  ya  = _mm256_fmadd_pd(za,cy,ya);  // + row 1
  vyl = _mm256_mul_pd(vzl, cyl);    // linear interpolation along y

  s1 = f1 + ni2;     // row 2, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+2);
  cyl = _mm256_broadcast_sd(cxyz+17);

  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 0
  za  = _mm256_mul_pd(zt, cz0);
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 1
  mi = _mm256_min_pd(mi, zt); ma = _mm256_max_pd(ma, zt);
  za  = _mm256_fmadd_pd(zt, cz1, za);
  cl  = _mm256_broadcast_sd(cxyz+18);           // linear coefficient along z
  vzl = _mm256_mul_pd(zt, cl);
  s1 += ninj;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 2
  mi = _mm256_min_pd(mi, zt); ma = _mm256_max_pd(ma, zt);
  za  = _mm256_fmadd_pd(zt, cz2, za);
  cl  = _mm256_broadcast_sd(cxyz+19);           // linear coefficient along z
  vzl = _mm256_fmadd_pd(zt, cl, vzl);  // linear interpolation along z
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 3
  za  = _mm256_fmadd_pd(zt, cz3, za);

  ya  = _mm256_fmadd_pd(za,cy,ya);  // + row 2
  vyl = _mm256_fmadd_pd(vzl, cyl, vyl);    // linear interpolation along y

  s1 = f1 + ni3;     // row 3, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+3);

  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 0
  za  = _mm256_mul_pd(zt, cz0);
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 1
  mi = _mm256_min_pd(mi, zt); ma = _mm256_max_pd(ma, zt);
  za  = _mm256_fmadd_pd(zt, cz1, za);
  s1 += ninj;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 2
  mi = _mm256_min_pd(mi, zt); ma = _mm256_max_pd(ma, zt);
  za  = _mm256_fmadd_pd(zt, cz2, za);
  s1 += ninjl;
  zt = _mm256_cvtps_pd(_mm_loadu_ps(s1));       // plane 3
  za  = _mm256_fmadd_pd(zt, cz3, za);

  ya  = _mm256_fmadd_pd(za,cy,ya);  // + row 3
  // fold min and max from 4 values to 1 (old version foldins [0] [1] [2] [3]
//   vmi = _mm_min_pd(_mm256_extractf128_pd(mi,0), _mm256_extractf128_pd(mi,1));  // min max ([0] [1] , [2] [3])
//   vma = _mm_max_pd(_mm256_extractf128_pd(ma,0), _mm256_extractf128_pd(ma,1));
//   vmi = _mm_min_pd(vmi, _mm_shuffle_pd(vmi, vmi, 1));   // min max ([0] [1] , [1] [0])
//   vma = _mm_max_pd(vma, _mm_shuffle_pd(vma, vma, 1));
  // only fold to get min max ([1] [2])
  vmi = _mm256_extractf128_pd(mi,0);
  vma = _mm256_extractf128_pd(ma,0);
  vmi = _mm_shuffle_pd(vmi, vmi, 1); // shuffle to get [1] [0]
  vma = _mm_shuffle_pd(vma, vma, 1); // shuffle to get [1] [0]
  vmi = _mm_min_pd(vmi, _mm256_extractf128_pd(mi,1));   // min ([1] [0] , [2] [3])
  vma = _mm_max_pd(vma, _mm256_extractf128_pd(ma,1));   // max ([1] [0] , [2] [3])
//   l = (int32_t *)min;
//   l[0] = _mm_extract_epi32((__m128i) _mm_cvtpd_ps(vmi), 0);  // convert min to float and store
  min[0] = _mm_cvtss_f32( _mm_cvtpd_ps(vmi) );          // min ([1] [2]), convert to float and store
//   l = (int32_t *)max;
//   l[0] = _mm_extract_epi32((__m128i) _mm_cvtpd_ps(vma), 0);  // convert max to float and store
  max[0] = _mm_cvtss_f32( _mm_cvtpd_ps(vma) );          // max ([1] [2]), convert to float and store

  // interpolation along x (linear)
  cx = _mm256_loadu_pd(cxyz+12);    // linear interpolation coefficients
  vyl = _mm256_mul_pd(vyl,cx);
  va = _mm_add_pd( _mm256_extractf128_pd(vyl,0) , _mm256_extractf128_pd(vyl,1) );   // [0] [1] + [2] [3]
  va = _mm_hadd_pd(va,va);                                                          // [0]+[2] + [1]+[3]
//   l = (int32_t *)lin;
//   ta = _mm_cvtpd_ps(va); l[0] = _mm_extract_epi32((__m128i) ta, 0);  // convert linear interp result to float and store
  lin[0] = _mm_cvtss_f32( _mm_cvtpd_ps(va) ) ;   // convert to float and store
  // interpolation along x (cubic)
  cx = _mm256_loadu_pd(px) ;     // make sure to use "unaligned" load to avoid segfault
  ya = _mm256_mul_pd(ya, cx );   // multiply by px
  // use folding to sum 4 terms
  va = _mm_add_pd( _mm256_extractf128_pd(ya,0) , _mm256_extractf128_pd(ya,1) ); va = _mm_hadd_pd(va,va);
//   l = (int32_t *)d;
//   ta = _mm_cvtpd_ps(va); l[0] = _mm_extract_epi32((__m128i) ta, 0);  // convert results to float and store
  d[0] = _mm_cvtss_f32( _mm_cvtpd_ps(va) ) ;   // convert to float and store
#else
  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  // field 1
  // TODO: add min max code
  for (i=0 ; i<4 ; i++){   // tricubic or bicubic/linear interpolation
    va4[i] = f1[i    ]*pz[0] + f1[i    +ninjl]*pz[1] +  f1[i    +ninj2]*pz[2] + f1[i    +ninj3]*pz[3];
    vb4[i] = f1[i+ni ]*pz[0] + f1[i+ni +ninjl]*pz[1] +  f1[i+ni +ninj2]*pz[2] + f1[i+ni +ninj3]*pz[3];
    vc4[i] = f1[i+ni2]*pz[0] + f1[i+ni2+ninjl]*pz[1] +  f1[i+ni2+ninj2]*pz[2] + f1[i+ni2+ninj3]*pz[3];
    vd4[i] = f1[i+ni3]*pz[0] + f1[i+ni3+ninjl]*pz[1] +  f1[i+ni3+ninj2]*pz[2] + f1[i+ni3+ninj3]*pz[3];
    dst[i] = va4[i]*py[0] + vb4[i]*py[1] + vc4[i]*py[2] + vd4[i]*py[3];
  }
  d[0] = dst[0]*px[0] + dst[1]*px[1] + dst[2]*px[2] + dst[3]*px[3];
  for (i=1 ; i<3 ; i++){   // linear interpolation
    vb4[i] = f1[i+ni +ninjl]*cxyz[18] + f1[i+ni +ninj2]*cxyz[19];
    vc4[i] = f1[i+ni2+ninjl]*cxyz[18] + f1[i+ni2+ninj2]*cxyz[19];
    dsl[i] = vb4[i]*cxyz[16] + vc4[i]*cxyz[17];
  }
  lin[0] = dsl[1]*cxyz[13] + dsl[2]*cxyz[14];
  ma = f1[1 + ni + ninjl] ; mi = ma;   // point [1,1,1] of 2x2x2 inner box
  for (i=1 ; i<3 ; i++){              // min max of 2x2x2 inner box
    ma = MAX(ma , f1[i + ni  + ninjl]); ma = MAX(ma , f1[i + ni  + ninj2]);
    mi = MIN(mi , f1[i + ni  + ninjl]); mi = MIN(mi , f1[i + ni  + ninj2]);

    ma = MAX(ma , f1[i + ni2 + ninjl]); ma = MAX(ma , f1[i + ni2 + ninj2]);
    mi = MIN(mi , f1[i + ni2 + ninjl]); mi = MIN(mi , f1[i + ni2 + ninj2]);
  }
  *max = ma ; *min = mi;
#endif
}

//****f* librkl/Tricublin_mono_zyx_n
// Synopsis
//
// ARGUMENTS
void Tricublin_mono_zyx_n(float *d, float *l, float *mi, float *ma, float *f, pxpypz *pxyz,  ztab *lv, int n)
//****
{
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    Tricublin_zyx_mm_d_inline(d, l, mi, ma, f + ixyz, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
    d++;         // next result
    l++;
    mi++;
    ma++;
    pxyz += 1;   // next set of positions
  }
}

//****f* librkl/Tricublin_mono_zyx_n_m
// Synopsis
//
// ARGUMENTS
void Tricublin_mono_zyx_n_m(float *d, float *l, float *mi, float *ma, float **fs, pxpypz *pxyz,  ztab *lv, int n, int m)
//****
{
  int i;
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    for(i=0 ; i<m ; i++){
      Tricublin_zyx_mm_d_inline(d, l, mi, ma, fs[i] + ixyz, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
      d++;         // next result
      l++;
      mi++;
      ma++;
    }
    pxyz += 1;   // next set of positions
  }
}

//****f* librkl/Tricublin_mono_zyx_m_n
// Synopsis
//
// ARGUMENTS
void Tricublin_mono_zyx_m_n(float *dm, float *lm, float *mim, float *mam, float **fs, pxpypz *pxyz,  ztab *lv, int n, int m)
//****
{
  int i;
  float *d0, *l0, *mi0, *ma0;
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int n0 = n;
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    d0 = dm;
    l0 = lm;
    mi0 = mim;
    ma0 = mam;
    for(i=0 ; i<m ; i++){
      Tricublin_zyx_mm_d_inline(d0, l0, mi0, ma0, fs[i] + ixyz, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
      d0  += n0;         // next result
      l0  += n0;
      mi0 += n0;
      ma0 += n0;
    }
    dm++;         // next set of results
    lm++;
    mim++;
    mam++;
    pxyz += 1;   // next set of positions
  }
}
// Fortran dimensions: d(3) , f(3,NI,NJ,NK)
// NI   : size of a 1D line
// NINJ : size of a 2D plane
// f : address of lower left corner of the 4 x 4 x 4 or 4 x 4 x 2 box
// interpolation is done along z, then along y, then along x
// 3 values are interpolated using the same cubic/linear polynomial coefficients
// pz(2) = 1.0 - dz, pz(3) = dz are expected in the linear along z case (0.0 <= dz <= 1.0)
// in the z linear case, planes 0 and 1 are the same as are planes 2 and 3
// note for future expansion:
// should both interpolations have to be done, planes 1 and 2 can be used for the z linear case,
// and planes 0, 1, 2, 3 for the z cubic case
// in that case another mechanism will have to be used to signal the z linear case
static inline void Tricublin_zyxf3_inline(float *d, float *f, double *pxyz, int NI, int NINJ, int zlinear){
// static inline void Tricublin_zyx3f_beta(float *d, float *f, double *px, double *py, double *pz, int NI, int NINJ){
  int ni = 3*NI;
  int ninj = 3*NINJ;
  double *px = pxyz;
  double *py = pxyz+4;
  double *pz = pxyz+8;
  int ninjl;    // ninj (cubic along z) or 0 (linear along z)
  float *s = f;
  double dst[13];
//   int64_t *L = (int64_t *)d;
  int ni2, ni3;

#if defined(__AVX2__) && defined(__x86_64__)
  __m256d cz0, cz1, cz2, cz3, cy, cx;
  __m256d x0;
  __m256d za, zb, zc;
  __m256d ya, yb, yc;
  __m128  vd;
  int32_t *l = (int32_t *)d;
#else
  double va4[12], vb4[12], vc4[12], vd4[12];
  int ninj2, ninj3, i;
#endif

  ni2 = ni + ni;
  ni3 = ni2 + ni;
  dst[12] = 0;
  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be set to 0
  if(zlinear) ninjl = 0;

#if defined(__AVX2__) && defined(__x86_64__)
  // ==== interpolation along Z, vector length is 12 (3 vectors of length 4 per plane) ====
  cz0 = _mm256_broadcast_sd(pz  );   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm256_broadcast_sd(pz+1);
  cz2 = _mm256_broadcast_sd(pz+2);
  cz3 = _mm256_broadcast_sd(pz+3);

  s = f;                   // row 0, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )), cz0);           // plane 0
  zb  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)), cz0);
  zc  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)), cz0);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz1,za);       // plane 1
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz1,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz1,zc);
  s += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz2,za);       // plane 2
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz2,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz2,zc);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz3,za);       // plane 3
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz3,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz3,zc);

  ya  = _mm256_mul_pd(za,cy);
  yb  = _mm256_mul_pd(zb,cy);
  yc  = _mm256_mul_pd(zc,cy);

  s = f + ni;              // row 1, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+1);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )), cz0);           // plane 0
  zb  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)), cz0);
  zc  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)), cz0);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz1,za);       // plane 1
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz1,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz1,zc);
  s += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz2,za);       // plane 2
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz2,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz2,zc);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz3,za);       // plane 3
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz3,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz3,zc);

  ya  = _mm256_fmadd_pd(za,cy,ya);
  yb  = _mm256_fmadd_pd(zb,cy,yb);
  yc  = _mm256_fmadd_pd(zc,cy,yc);

  s = f + ni2;            // row 2, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+2);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )), cz0);           // plane 0
  zb  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)), cz0);
  zc  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)), cz0);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz1,za);       // plane 1
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz1,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz1,zc);
  s += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz2,za);       // plane 2
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz2,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz2,zc);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz3,za);       // plane 3
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz3,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz3,zc);

  ya  = _mm256_fmadd_pd(za,cy,ya);
  yb  = _mm256_fmadd_pd(zb,cy,yb);
  yc  = _mm256_fmadd_pd(zc,cy,yc);

  s = f + ni3;            // row 3, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+3);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )), cz0);           // plane 0
  zb  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)), cz0);
  zc  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)), cz0);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz1,za);       // plane 1
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz1,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz1,zc);
  s += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz2,za);       // plane 2
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz2,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz2,zc);
  s += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s  )),cz3,za);       // plane 3
  zb  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+4)),cz3,zb);
  zc  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s+8)),cz3,zc);

  ya  = _mm256_fmadd_pd(za,cy,ya);
  yb  = _mm256_fmadd_pd(zb,cy,yb);
  yc  = _mm256_fmadd_pd(zc,cy,yc);

  // interpolation along x, split 3 vectors of 4 elements into 4 vectors of 3 elements
  // dst : x(0,0) x(1,0) x(2,0) x(0,1)  x(1,1) x(2,1) x(0,2) x(1,2)   x(2,2) x(0,3) x(1,3) x(2,3)    0
          _mm256_storeu_pd(dst,ya);     _mm256_storeu_pd(dst+4,yb);   _mm256_storeu_pd(dst+8,yc);
  // cheating loads : only the first 3 values are of any interest
  //    x(0,0) x(1,0) x(2,0) x(0,1)  i.e. ya
  cx = _mm256_broadcast_sd(px) ;
  x0 = _mm256_mul_pd(ya, cx );                           // column 0
  //    x(0,1) x(1,1) x(2,1) x(0,2)
  cx = _mm256_broadcast_sd(px+1) ;
  x0 = _mm256_fmadd_pd(_mm256_loadu_pd(dst+3), cx ,x0);  // column 1
  //    x(0,2) x(1,2) x(2,2) x(0,3)
  cx = _mm256_broadcast_sd(px+2) ;
  x0 = _mm256_fmadd_pd(_mm256_loadu_pd(dst+6), cx ,x0);  // column 2
  //    x(0,3) x(1,3) x(2,3) 0
  cx = _mm256_broadcast_sd(px+3) ;
  x0 = _mm256_fmadd_pd(_mm256_loadu_pd(dst+9), cx ,x0);  // column 3

  vd =  _mm256_cvtpd_ps(x0);                  // back to float precision
  l[0] = _mm_extract_epi32((__m128i) vd, 0);  // extract element 0 and store
  l[1] = _mm_extract_epi32((__m128i) vd, 1);  // extract element 1 and store
  l[2] = _mm_extract_epi32((__m128i) vd, 2);  // extract element 2 and store
#else
  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  for (i=0 ; i<12 ; i++){
    va4[i] = s[i    ]*pz[0] + s[i    +ninjl]*pz[1] +  s[i    +ninj2]*pz[2] + s[i    +ninj3]*pz[3];
    vb4[i] = s[i+ni ]*pz[0] + s[i+ni +ninjl]*pz[1] +  s[i+ni +ninj2]*pz[2] + s[i+ni +ninj3]*pz[3];
    vc4[i] = s[i+ni2]*pz[0] + s[i+ni2+ninjl]*pz[1] +  s[i+ni2+ninj2]*pz[2] + s[i+ni2+ninj3]*pz[3];
    vd4[i] = s[i+ni3]*pz[0] + s[i+ni3+ninjl]*pz[1] +  s[i+ni3+ninj2]*pz[2] + s[i+ni3+ninj3]*pz[3];
    dst[i] = va4[i]*py[0] + vb4[i]*py[1] + vc4[i]*py[2] + vd4[i]*py[3];
  }
  d[0] = dst[0]*px[0] + dst[3]*px[1] + dst[6]*px[2] + dst[ 9]*px[3];
  d[1] = dst[1]*px[0] + dst[4]*px[1] + dst[7]*px[2] + dst[10]*px[3];
  d[2] = dst[2]*px[0] + dst[5]*px[1] + dst[8]*px[2] + dst[11]*px[3];
#endif
}

// process n triplets
// for each triplet 1 value from ixyz, 3 values from pxyz are used (1 element)
// it is ASSUMED that along X and Y interpolation will always be CUBIC
// ( 2 <= px < "ni"-1 ) ( 2 <= py < "nj"-1 )
// pz < 2 and pz >= nk - 1 will induce linear interpolation or extrapolation
//****f* librkl/Tricublin_zyx3_n
// Synopsis
//
// ARGUMENTS
void Tricublin_zyx3_n(float *d, float *f3, pxpypz *pxyz,  ztab *lv, int n)
//****
{
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
/*printf*/("%12.7f %12.7f %12.7f\n",pxyz->px, pxyz->py, pxyz->pz);
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    Tricublin_zyxf3_inline(d, f3 + ixyz*3, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
    d+=3;         // next result
    pxyz += 1;   // next set of positions
  }
}

#if defined(SELF_TEST)
// x, y, z are in delta form (0 <= delta <= 1.0)
static void Tricubic_coeffs_d(double *px, double *py, double *pz, double x, double y, double z){

  pz[0] = cm167*z*(z-one)*(z-two);        // coefficients for interpolation along z
  pz[1] = cp5*(z+one)*(z-one)*(z-two);
  pz[2] = cm5*z*(z+one)*(z-two);
  pz[3] = cp167*z*(z+one)*(z-one);
  pz[4] = 0.0;
  pz[5] = 1.0 - z;
  pz[6] = z;
  pz[7] = 0.0;

  py[0] = cm167*y*(y-one)*(y-two);        // coefficients for interpolation along y
  py[1] = cp5*(y+one)*(y-one)*(y-two);
  py[2] = cm5*y*(y+one)*(y-two);
  py[3] = cp167*y*(y+one)*(y-one);
  py[4] = 0.0;
  py[5] = 1.0 - y;
  py[6] = y;
  py[7] = 0.0;

  px[0] = cm167*x*(x-one)*(x-two);        // coefficients for interpolation along x
  px[1] = cp5*(x+one)*(x-one)*(x-two);
  px[2] = cm5*x*(x+one)*(x-two);
  px[3] = cp167*x*(x+one)*(x-one);
  px[4] = 0.0;
  px[5] = 1.0 - x;
  px[6] = x;
  px[7] = 0.0;
}

#define NI 300
#define NJ 200
#define NK 85
#include <stdio.h>
#include <sys/time.h>

float fx(float x) { x = x*.125 ; return (x*x*x*1.011 + x*x*1.022 + x*1.01 + .123); }
float fy(float y) { y = y*.125 ; return (y*y*y*1.044 + y*y*1.055 + y*1.03 + .234); }
float fz(float z) { return (z*1.123 + .345); }
float fxyz(float x, float y, float z) { return fx(x) * fy(y) * fz(z) ; }

int main(int argc, char **argv){
  float array[3*NI*NJ*NK];
  float a1[NI*NJ*NK];
  float a2[NI*NJ*NK];
  float a3[NI*NJ*NK];
  float dest[3*NI*NJ*NK];
  float *destp;
  float mi[3];
  float ma[3];
  float lin[3];
  double pz[16], px[16], py[16], pxyz[12];
  float avg=0.0;
  float dx, dy, dz;
  float expected1, expected2;
  struct timeval t1, t2;
  long long tm1, tm2;
  int i,j,k,ijk;
  float *p = array;
  float *p1 = a1;
  float *p2 = a2;
  float *p3 = a3;
  int nijk = NI*NJ*NK - 5*NI*NJ;
  int II, JJ, KK;
  pxpypz posxyz[NI*NJ*NK];
  double expected[NI*NJ*NK];
  pxpypz *ptrxyz;
  double *answer = expected;
  double posz[NK];
  ztab *lv;
  float e0, e1, e2;
  double e8;
  uint64_t tt0, tt1;
  int npts, exact;

  ptrxyz = posxyz;
  for (k=0 ; k<NK ; k++ ){
    posz[k] = k + 1 ;
    for (j=0 ; j<NJ ; j++ ){
      for (i=0 ; i<NI ; i++ ){
//         p[0] = i + j + k + 3;
        p[0] = fxyz(i+1.0,j+1.0,k+1.0);
        p[1] = p[0] + 10;
        p[2] = p[0] + 100;
        *p1++ = p[0];
        *p2++ = p[1];
        *p3++ = p[2];
        p += 3;
      }
    }
  }

  II = NI/2;
  JJ = NJ/2;
  KK = NK/2;
  dz = .333; posxyz[0].pz = dz ; posxyz[0].pz = posxyz[0].pz + KK;
  dx = .555; posxyz[0].px = dx + II;
  dy = .777; posxyz[0].py = dy + JJ;

  lv = Vsearch_setup(posz, NK, NI, NJ);

//   Tricubic_coeffs_d(pxyz, pxyz+4, pxyz+8, dx, dy, dz);

  dest[0] = -1.0 ; dest[1] = -1.0 ; dest[2] = -1.0 ;
  Tricublin_zyx1_n(&dest[0], &a1[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[1], &a2[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[2], &a3[0], posxyz,  lv, 1);
//   expected2 = II + dx + JJ + dy + KK + dz ;
//   expected2 = fxyz(II + dx, JJ + dy, KK + dz);
  e0 = fxyz(II + dx, JJ + dy, KK + dz) ; e1 = e0+10 ; e2 = e0+100;
  printf("Cubic case  :\ngot %12.7g, %12.7g, %12.7g, result/expected = %10.8f %10.8f %10.8f\n",dest[0],dest[1],dest[2],dest[0]/e0,dest[1]/e1,dest[2]/e2);

  dest[0] = -1.0 ; dest[1] = -1.0 ; dest[2] = -1.0 ;
  KK = 1;
  dz = .25; posxyz[0].pz = dz ; posxyz[0].pz = posxyz[0].pz + KK;
  Tricublin_zyx1_n(&dest[0], &a1[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[1], &a2[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[2], &a3[0], posxyz,  lv, 1);
//   expected2 = fxyz(II + dx, JJ + dy, KK + dz);
  e0 = fxyz(II + dx, JJ + dy, KK + dz) ; e1 = e0+10 ; e2 = e0+100;
  printf("Linear bottom case  :\ngot %12.7g, %12.7g, %12.7g, result/expected = %10.8f %10.8f %10.8f\n",dest[0],dest[1],dest[2],dest[0]/e0,dest[1]/e1,dest[2]/e2);

  dest[0] = -1.0 ; dest[1] = -1.0 ; dest[2] = -1.0 ;
  KK = 0;
  dz = .75; posxyz[0].pz = dz ; posxyz[0].pz = posxyz[0].pz + KK;
  Tricublin_zyx1_n(&dest[0], &a1[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[1], &a2[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[2], &a3[0], posxyz,  lv, 1);
//   expected2 = fxyz(II + dx, JJ + dy, KK + dz);
  e0 = fxyz(II + dx, JJ + dy, KK + dz) ; e1 = e0+10 ; e2 = e0+100;
  printf("Linear bot extrap case  :\ngot %12.7g, %12.7g, %12.7g, result/expected = %10.8f %10.8f %10.8f\n",dest[0],dest[1],dest[2],dest[0]/e0,dest[1]/e1,dest[2]/e2);

  dest[0] = -1.0 ; dest[1] = -1.0 ; dest[2] = -1.0 ;
  KK = NK - 1;
  dz = .75; posxyz[0].pz = dz ; posxyz[0].pz = posxyz[0].pz + KK;
  Tricublin_zyx1_n(&dest[0], &a1[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[1], &a2[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[2], &a3[0], posxyz,  lv, 1);
//   expected2 = fxyz(II + dx, JJ + dy, KK + dz);
  e0 = fxyz(II + dx, JJ + dy, KK + dz) ; e1 = e0+10 ; e2 = e0+100;
  printf("Linear top case  :\ngot %12.7g, %12.7g, %12.7g, result/expected = %10.8f %10.8f %10.8f\n",dest[0],dest[1],dest[2],dest[0]/e0,dest[1]/e1,dest[2]/e2);

  dest[0] = -1.0 ; dest[1] = -1.0 ; dest[2] = -1.0 ;
  KK = NK;
  dz = .25; posxyz[0].pz = dz ; posxyz[0].pz = posxyz[0].pz + KK;
  Tricublin_zyx1_n(&dest[0], &a1[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[1], &a2[0], posxyz,  lv, 1);
  Tricublin_zyx1_n(&dest[2], &a3[0], posxyz,  lv, 1);
//   expected2 = fxyz(II + dx, JJ + dy, KK + dz);
  e0 = fxyz(II + dx, JJ + dy, KK + dz) ; e1 = e0+10 ; e2 = e0+100;
  printf("Linear top extrap case  :\ngot %12.7g, %12.7g, %12.7g, result/expected = %10.8f %10.8f %10.8f\n",dest[0],dest[1],dest[2],dest[0]/e0,dest[1]/e1,dest[2]/e2);

  answer = expected;
  for (k=0 ; k<NK ; k++ ){
    for (j=1 ; j<NJ-2 ; j++ ){
      for (i=1 ; i<NI-2 ; i++ ){
        ptrxyz->px = i+dx+1;
        ptrxyz->py = j+dy+1;
        ptrxyz->pz = k+dz+1;
        *answer = fxyz(ptrxyz->px,ptrxyz->py,ptrxyz->pz);
        ptrxyz++;
        answer++;
      }
    }
  }
  npts = NK * (NI-3) * (NJ-3); // npts /=2.295; npts = 2167503;
  printf("npts = %d, ninjnk = %d\n",npts,NK*NI*NJ);
  for(i=0 ; i<npts ; i++) dest[i] = -1.0;
  tt0 = rdtscp_();
  Tricublin_zyx1_n(dest,a1,posxyz,lv,npts);
  tt1 = rdtscp_();
  tt1 = tt1 - tt0;
  printf("cycles per interpolation = %f\n",tt1 / (1.0*npts));

//   npts = 10;
  e0 = 0.0;
  e8 = 0.0;
  exact = 0;
  for(i=0 ; i<npts ; i++){
    if(i%300000 == 0) printf("expected %12.7g, got %12.7g\n",expected[i],dest[i]);
    e1 = 1.0 - (dest[i]/expected[i]);
    if(e1 < 0) e1 = -e1;
    if(e1 == 0)exact++;
    e8 += e1;
    if(e1 > e0) e0 = e1;
  }
  e1 = npts ; e1 = (exact / e1) * 100.0;
  printf("worst case error : %12.7g, avg error : %12.7g, exact : %8.2f percent\n",e0,e8/npts,e1);
  printf("END\n");
}
#endif
