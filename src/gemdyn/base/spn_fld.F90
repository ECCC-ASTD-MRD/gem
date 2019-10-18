!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END --------------------------------

!*s/r spn_fld - doing forward 2-D FFT, applying a filter, doing backward FFT
!             - and applying nudging tendency

      subroutine spn_fld ( F_Minx, F_Maxx, F_Miny, F_Maxy, F_Njl,   &
                           F_Minz, F_Maxz, F_Nk, F_Nkl, F_Gni, F_Gnj, &
                           F_Minij, F_Maxij, F_nij, F_nij0,      &
                           F_npex1, F_npey1, Fld_S )

      use cstv
      use gem_options
      use spn_options
      use glb_ld
      use gmm_itf_mod
      use gmm_vt1
      use gmm_nest
      use spn_work_mod
      use step_options
      use tdpack
      use glb_pil
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer  F_Minx, F_Maxx, F_Miny, F_Maxy, F_Njl
      integer  F_Minz, F_Maxz, F_Nk, F_Nkl, F_Gni, F_Gnj
      integer  F_Minij, F_Maxij, F_nij, F_nij0
      integer  F_npex1, F_npey1
      character (len=1) Fld_S

!author
!     Minwei Qian (CCRD) & Bernard Dugas, Syed Husain  (MRB)  - summer 2015
!
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_Minx       I    - minimum index on X (ldnh_minx)
! F_Maxx       I    - maximum index on X (ldnh_maxx)
! F_Miny       I    - minimum index on Y (ldnh_miny)
! F_Maxy       I    - maximum index on Y (ldnh_maxy)
! F_Njl        I    - number of points on local PEy for J (ldnh_nj)
! F_Minz       I    - minimum index on local PEx for K (trp_12smin)
! F_Maxz       I    - maximum index on local PEx for K (trp_12smax)
! F_Nk         I    - G_nk points in Z direction globally
! F_Nkl        I    - number of points on local PEx for K (trp_12sn)
! F_Gni        I    - number of points in X direction globally (G_ni)
! F_Gnj        I    - number of points in Y direction globally (G_nj)
! F_Minij      I    - minimum index on local PEy for I (trp_22min)
! F_Maxij      I    - maximum index on local PEy for I (trp_22max)
! F_nij        I    - number of points on local PEy for I (trp_22n)
! F_nij0       I    - global offset of the first I element on PEy
! F_npex1      I    - number of processors in X
! F_npey1      I    - number of processors in Y
! Fld_S        I    - name of variable to treat (either of 't','u','v')


      external rpn_comm_transpose

      real(kind=REAL64) fdwfft(F_Miny:F_Maxy,F_Minz :F_Maxz ,F_Gni+2+F_npex1)
      real(kind=REAL64)   fdg2(F_Minz:F_Maxz,F_Minij:F_Maxij,F_Gnj+2+F_npey1)
      real(kind=REAL64)  pri

      integer gmmstat
      type(gmm_metadata):: metadata
      real, dimension(:,:,:), pointer :: fld3d=>null(), fld_nest3d=>null()
      integer i,  j, k
      integer ii,jj,kk

      integer no_steps, tmdt
      real spn_wt
!
!----------------------------------------------------------------------
!
      tmdt    = int(Cstv_dt_8)
      no_steps= Step_nesdt/tmdt
      spn_wt = 1.0

      if (Spn_weight_L) then
        spn_wt = sqrt((cos(pi_8*(float(Lctl_step)/float(no_steps))))**2)**Spn_wt_pwr
      end if

      if (Fld_S == 't') then

         gmmstat = gmm_getmeta (gmmk_tt1_s, metadata)
         gmmstat = gmm_get (gmmk_tt1_s, fld3d, metadata)
         gmmstat = gmm_get (gmmk_nest_t_s, fld_nest3d, metadata)

         Ldiff3D (F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk)= &
              fld_nest3d(F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk) - &
                   fld3d(F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk)

!     call glbstat ( fld3d,'TTA',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj,1,F_nk )
!     call glbstat ( fld_nest3d,'NTTA',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj,1,F_nk )

      end if

      if (Fld_S == 'u') then

         gmmstat = gmm_getmeta (gmmk_ut1_s, metadata)
         gmmstat = gmm_get (gmmk_ut1_s, fld3d, metadata)
         gmmstat = gmm_get (gmmk_nest_u_s, fld_nest3d, metadata)

         Ldiff3D (F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk)= &
              fld_nest3d(F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk) - &
                   fld3d(F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk)

!     call glbstat ( fld3d,'UT1',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni-1,1,G_nj,1,F_nk )
!     call glbstat ( fld_nest3d,'NUT1',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni-1,1,G_nj,1,F_nk )

!     call glbstat ( fld3d,'UTA',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj,1,F_nk )
!     call glbstat ( fld_nest3d,'NUTA',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj,1,F_nk )

      end if

      if (Fld_S == 'v') then

         gmmstat = gmm_getmeta (gmmk_vt1_s, metadata)
         gmmstat = gmm_get (gmmk_vt1_s, fld3d, metadata)
         gmmstat = gmm_get (gmmk_nest_v_s, fld_nest3d, metadata)

         Ldiff3D (F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk)= &
              fld_nest3d(F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk) - &
                   fld3d(F_Minx:F_Maxx,F_Miny:F_Maxy,1:F_Nk)
!     call glbstat ( fld3d,'VT1',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj-1,1,F_nk )
!     call glbstat ( fld_nest3d,'NVT1',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj-1,1,F_nk )
!     call glbstat ( fld3d,'VTA',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj,1,F_nk )
!     call glbstat ( fld_nest3d,'NVTA',"spnfld",F_minx,F_maxx,F_miny,F_maxy, &
!                     1,F_nk, 1,G_ni,1,G_nj,1,F_nk )

      end if


! do transpose from (i,j,k) to (j,k,i)
      call rpn_comm_transpose                         &
           ( Ldiff3D, F_Minx, F_Maxx, F_Gni, (F_Maxy-F_Miny+1), &
                      F_Minz, F_Maxz, F_Nk,  fdwfft, 1, 2 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! projection ( wfft = x transposed * g )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$omp parallel
!$omp do
      do i= 1,F_Gni               ! trimming in X
         do k= F_Minz, F_Nkl
            do j= F_Njl+1-pil_n,F_Maxy
!               fdwfft(F_Njl+1-pil_n:F_Maxy,F_Minz:F_Nkl,i) = 0.
               fdwfft(j,k,i) = 0.
            end do
         end do
         do k= F_Minz, F_Nkl
            do j= F_Miny, pil_s
!               fdwfft(F_Miny:pil_s,F_Minz:F_Nkl,i) = 0.
               fdwfft(j,k,i) = 0.
            end do
         end do
         do k= F_Nkl+1,F_Maxz
            do j= F_Miny,F_Maxy
!               fdwfft(F_Miny:F_Maxy,F_Nkl+1:F_Maxz,i) = 0.
               fdwfft(j,k,i) = 0.
            end do
         end do
      end do
!$omp enddo


!$omp single
      call itf_fft_set( F_Gni-Lam_pil_w-Lam_pil_e,'QCOS',pri )
!$omp end single

!$omp do
      do k=1,F_Nkl                ! do forward fft in X direction
         call itf_fft_drv( fdwfft(1+pil_s,k,1+Lam_pil_w), &
         (F_Maxy-F_Miny+1)*(F_Maxz-F_Minz+1),1,       &
         (F_Maxy-F_Miny+1-pil_s-pil_n), -1 )
      end do
!$omp enddo

!$omp single
! do transpose from (j,k,i) to (k,i,j)
      call rpn_comm_transpose                          &
      ( fdwfft, F_Miny,  F_Maxy,  F_Gnj, (F_Maxz-F_Minz+1), &
      F_Minij, F_Maxij, F_Gni, fdg2, 2, 2 )
      call itf_fft_set( F_Gnj-Lam_pil_s-Lam_pil_n,'QCOS',pri )
!$omp end single

!$omp do
      do k=1,F_nij              ! do forward fft in Y direction
         call itf_fft_drv( fdg2(1,k,1+Lam_pil_s),     &
         (F_Maxz-F_Minz+1)*(F_Maxij-F_Minij+1),1, &
         (F_Maxz-F_Minz+1), -1 )
      end do
!$omp enddo

!$omp do
      do jj=1,G_nj+2            ! do filter in X-Y direction
         do ii=1,F_nij
            do kk=F_Minz,F_Maxz
               fdg2(kk,ii,jj)=fdg2(kk,ii,jj)*fxy(ii+F_nij0,jj)
            end do
         end do
      end do
!$omp enddo

!$omp do
      do k=1,F_nij              ! do backward fft in Y direction
         call itf_fft_drv( fdg2(1,k,1+Lam_pil_s),     &
         (F_Maxz-F_Minz+1)*(F_Maxij-F_Minij+1),1, &
         (F_Maxz-F_Minz+1), +1 )
      end do
!$omp enddo

!$omp single
! do backward transpose from (k,i,j) to (j,k,i)
      call rpn_comm_transpose                          &
      ( fdwfft, F_Miny,  F_Maxy,  F_Gnj, (F_Maxz-F_Minz+1), &
      F_Minij, F_Maxij, F_Gni, fdg2, -2, 2 )
      call itf_fft_set( F_Gni-Lam_pil_w-Lam_pil_e,'QCOS',pri )
!$omp end single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! inverse projection ( r = x * w )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$omp do
      do k=1, F_Nkl               ! do backward fft in X direction
         call itf_fft_drv( fdwfft(1+pil_s,k,1+Lam_pil_w), &
         (F_Maxy-F_Miny+1)*(F_Maxz-F_Minz+1),1,       &
         (F_Maxy-F_Miny+1-pil_s-pil_n), +1 )
      end do
!$omp enddo
!$omp end parallel

! do backward transpose from (j,k,i) to (i,j,k)
      call rpn_comm_transpose                         &
           ( Ldiff3D, F_Minx, F_Maxx, F_Gni, (F_Maxy-F_Miny+1), &
                     F_Minz, F_Maxz, F_Nk,  fdwfft, -1, 2 )

      do kk=2,F_Nk
         fld3d(1:l_ni,1:l_nj,kk) = &
         fld3d(1:l_ni,1:l_nj,kk) + &
         prof(kk)*SNGL(Ldiff3D(1:l_ni,1:l_nj,kk))*spn_wt
      end do
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_fld

