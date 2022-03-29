!--------------------------------- LICENCE BEGIN -------------------------------
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
!---------------------------------- LICENCE END ---------------------------------

!**s/r dcmip_vrd_main - Call variables and Apply vertical diffusion over column
!                       on Winds, Potential temperature and Tracers

      subroutine dcmip_vrd_main

      use glb_ld
      use gmm_vt1
      use gmm_pw
      use lun

      implicit none

#include <arch_specific.hf>

      !object
      !============================================================
      !     Call variables and Apply vertical diffusion over column
      !     on Winds, Potential temperature and Tracers
      !============================================================
!
!---------------------------------------------------------------------
!
      if (Lun_out>0) write (Lun_out,1000)

      !Diffuse Winds, Potential temperature and Tracers
      !------------------------------------------------
      call dcmip_vrd_drv (ut1,vt1,zdt1,wt1,tt1,pw_pt_plus,l_minx,l_maxx,l_miny,l_maxy,l_nk)

      return

 1000 format( &
      /,'APPLY VERTICAL DIFFUSION: (S/R DCMIP_VRD_MAIN)',   &
      /,'==============================================',/,/)

      end subroutine dcmip_vrd_main

!--------------------------------------------------------------------------------

!**s/r dcmip_vrd_drv - Apply vertical diffusion over column on Winds, Potential temperature and Tracers
!                      (Based on eqspng_drv)

      subroutine dcmip_vrd_drv (F_u,F_v,F_zd,F_w,F_tv,F_pr_t,Minx,Maxx,Miny,Maxy,Nk)

      use canonical
      use dcmip_vrd_coef
      use glb_ld
      use mem_tracers
      use tdpack, only : cappa_8
      use tr3d

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_u,F_v,F_zd,F_w,F_tv
      real, dimension(Minx:Maxx,Miny:Maxy,NK), intent(in)    :: F_pr_t

      !______________________________________________________________________
      !        |                                             |           |   |
      ! NAME   |             DESCRIPTION                     | DIMENSION |I/O|
      !--------|---------------------------------------------|-----------|---|
      ! F_u    | x component of velocity                     | 3D (Nk)   |i/o|
      ! F_v    | y component of velocity                     | 3D (Nk)   |i/o|
      ! F_zd   | dZeta/dt                                    | 3D (Nk)   |i/o|
      ! F_w    | dz/dt                                       | 3D (Nk)   |i/o|
      ! F_tv   | Virtual temperature                         | 3D (Nk)   |i/o|
      ! F_pr_t | pressure THERMO                             | 3D (Nk)   |i  |
      !________|_____________________________________________|___________|___|
      !
      !    There are nlev levels that will be diffused.
      !    There are nlev-1 coefficient passed in namelist by user.
      !
      !    The boundary conditions are flux=0. This is achieved by
      !    putting coef=0 at boundaries (see drawing below).
      !    Also the index km and kp make the wind derivatives zero
      !    at boundary.
      !
      !    ---- this u equals u(1) (see km in loop)   \
      !                                               |
      !    ==== coef(1)=0, Top boundary.              | this derivative = 0
      !                                               |
      !    ---- u(1) first levels diffused            /
      !
      !    ==== coef(2)=first coef passed by user
      !
      !    ---- u(2)
      !
      !        ...
      !
      !    ---- u(nlev-1)
      !
      !    ==== coef(nlev)=last coef passed by user
      !
      !    ---- u(nlev) last level diffused           \
      !                                               |
      !    ==== coef(nlev+1)=0, Bottom boundary.      | this derivative = 0
      !                                               |
      !    ---- this u equal u(nlev) (see kp in loop) /
      !
      !_____________________________________________________________________
      !

      !object
      !=====================================================================================
      !     Apply vertical diffusion over column on Winds, Potential temperature and Tracers
      !     (Based on eqspng_drv)
      !=====================================================================================

      integer :: i,j,k,n,istat
      real, pointer, dimension(:,:,:) :: tr
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: pp,th_p
!
!---------------------------------------------------------------------
!
      !-----------------
      !Diffuse U,V,W,ZD
      !-----------------
      if (Dcmip_wd_L) then

         call dcmip_vrd_fld (F_u,  uref,Dcmip_ref_wd,Dcmip_cp_wd_m,Dcmip_cm_wd_m,Minx,Maxx,Miny,Maxy,Nk,1,l_niu,1,l_nj )
         call dcmip_vrd_fld (F_v,  vref,Dcmip_ref_wd,Dcmip_cp_wd_m,Dcmip_cm_wd_m,Minx,Maxx,Miny,Maxy,Nk,1,l_ni ,1,l_njv)

         call dcmip_vrd_fld (F_w,  wref,Dcmip_ref_wd,Dcmip_cp_wd_t,Dcmip_cm_wd_t,Minx,Maxx,Miny,Maxy,Nk,1,l_ni,1,l_nj)
         call dcmip_vrd_fld (F_zd,zdref,Dcmip_ref_wd,Dcmip_cp_wd_t,Dcmip_cm_wd_t,Minx,Maxx,Miny,Maxy,Nk,1,l_ni,1,l_nj)

      end if

      !-----------------------------
      !Diffuse Potential temperature
      !-----------------------------
      if (Dcmip_th_L) then

         !Calculate Potential temperature from Virtual temperature
         !--------------------------------------------------------
         do k = 1,Nk
         do j = 1,l_nj
         do i = 1,l_ni

            !Theta=t/pp; pp=(pr/p0)**cappa; p0=1
            !-----------------------------------
            pp  (i,j,k) = F_pr_t(i,j,k)**cappa_8

            th_p(i,j,k) = F_tv(i,j,k)/pp(i,j,k)

         end do
         end do
         end do

         call dcmip_vrd_fld (th_p,thref,Dcmip_ref_th,Dcmip_cp_th_t,Dcmip_cm_th_t,Minx,Maxx,Miny,Maxy,Nk,1,l_ni,1,l_nj)

         do k = 1,Nk
         do j = 1,l_nj
         do i = 1,l_ni
            F_tv(i,j,k) = pp(i,j,k)*th_p(i,j,k)
         end do
         end do
         end do

      end if

      !---------------
      !Diffuse Tracers
      !---------------
      if (Dcmip_tr_L) then

      do n=1,Tr3d_ntr

         if (Tr3d_name_S(n)=="HU") istat = tr_get('HU:P',tr)
         if (Tr3d_name_S(n)=="QC") istat = tr_get('QC:P',tr)
         if (Tr3d_name_S(n)=="RW") istat = tr_get('RW:P',tr)

         if (Tr3d_name_S(n)=="HU") call dcmip_vrd_fld (tr,qvref,Dcmip_ref_tr,Dcmip_cp_tr_t,Dcmip_cm_tr_t,Minx,Maxx,Miny,Maxy,Nk,1,l_ni,1,l_nj)
         if (Tr3d_name_S(n)=="QC") call dcmip_vrd_fld (tr,qcref,Dcmip_ref_tr,Dcmip_cp_tr_t,Dcmip_cm_tr_t,Minx,Maxx,Miny,Maxy,Nk,1,l_ni,1,l_nj)
         if (Tr3d_name_S(n)=="RW") call dcmip_vrd_fld (tr,qrref,Dcmip_ref_tr,Dcmip_cp_tr_t,Dcmip_cm_tr_t,Minx,Maxx,Miny,Maxy,Nk,1,l_ni,1,l_nj)

         do k=1,Nk
            do j=1,l_nj
               do i=1,l_ni
                  tr(i,j,k)=max(0.,tr(i,j,k))
               end do
            end do
         end do

      end do

      end if

      return

      end subroutine dcmip_vrd_drv

!--------------------------------------------------------------------------------

!**s/r dcmip_vrd_fld - Apply vertical diffusion over column on a given fld (Based on eqspng_drv)

      subroutine dcmip_vrd_fld (F_fld,F_ref_fld,F_ref_on,F_cp,F_cm,Minx,Maxx,Miny,Maxy,Nk,F_i0,F_in,F_j0,F_jn)

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in)    :: Minx,Maxx,Miny,Maxy,Nk,F_i0,F_in,F_j0,F_jn
      real,    intent(inout) :: F_fld(Minx:Maxx,Miny:Maxy,Nk)
      real,    intent(in)    :: F_ref_fld(Minx:Maxx,Miny:Maxy,Nk),F_cm(Nk),F_cp(Nk),F_ref_on

      !object
      !==============================================================================
      !     Apply vertical diffusion over column on a given fld (Based on eqspng_drv)
      !==============================================================================

      real, dimension(:,:,:), pointer :: w_fld,p_fld
      integer :: i,j,k,km,kp
!
!---------------------------------------------------------------------
!
      allocate (w_fld(F_i0:F_in,F_j0:F_jn,Nk),p_fld(F_i0:F_in,F_j0:F_jn,Nk))

      p_fld (F_i0:F_in,F_j0:F_jn,1:Nk) = F_fld (F_i0:F_in,F_j0:F_jn,1:Nk) - F_ref_on * F_ref_fld(F_i0:F_in,F_j0:F_jn,1:Nk)

      do k=1,Nk
         kp=min(Nk,k+1)
         km=max(1,k-1)
         do j=F_j0,F_jn
            do i=F_i0,F_in
               w_fld(i,j,k)=p_fld(i,j,k)+(F_cp(k)*(p_fld(i,j,kp)-p_fld(i,j,k )) &
                                         -F_cm(k)*(p_fld(i,j,k )-p_fld(i,j,km)))
            end do
         end do
      end do

      F_fld(F_i0:F_in,F_j0:F_jn,1:Nk) = w_fld(F_i0:F_in,F_j0:F_jn,1:Nk) + F_ref_on * F_ref_fld(F_i0:F_in,F_j0:F_jn,1:Nk)

      deallocate(w_fld,p_fld)
!
!---------------------------------------------------------------------
!
      return

      end subroutine dcmip_vrd_fld
