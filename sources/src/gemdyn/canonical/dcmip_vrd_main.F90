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

      use gmm_vt1
      use glb_ld
      use lun
      use gmm_itf_mod

      implicit none

#include <arch_specific.hf>

      !object
      !============================================================
      !     Call variables and Apply vertical diffusion over column
      !     on Winds, Potential temperature and Tracers
      !============================================================

      integer istat

      if (Lun_out>0) write (Lun_out,1000)

      !Get GMM variables
      !-----------------
      istat = gmm_get(gmmk_ut1_s,  ut1)
      istat = gmm_get(gmmk_vt1_s,  vt1)
      istat = gmm_get(gmmk_wt1_s,  wt1)
      istat = gmm_get(gmmk_zdt1_s, zdt1)
      istat = gmm_get(gmmk_tt1_s,  tt1)
      istat = gmm_get(gmmk_st1_s,  st1)

      !Diffuse Winds, Potential temperature and Tracers
      !------------------------------------------------
      call dcmip_vrd_drv (ut1,vt1,zdt1,wt1,tt1,st1,l_minx,l_maxx,l_miny,l_maxy,l_nk)

      return

 1000 format( &
      /,'APPLY VERTICAL DIFFUSION: (S/R DCMIP_VRD_MAIN)',   &
      /,'==============================================',/,/)

      end subroutine dcmip_vrd_main

!--------------------------------------------------------------------------------

!**s/r dcmip_vrd_drv - Apply vertical diffusion over column on Winds, Potential temperature and Tracers
!                      (Based on eqspng_drv)

      subroutine dcmip_vrd_drv (F_u,F_v,F_zd,F_w,F_tv,F_s,Minx,Maxx,Miny,Maxy,Nk)

      use canonical
      use dcmip_options
      use dcmip_vrd_coef
      use tdpack, only : cappa_8
      use glb_ld
      use tr3d
      use ver
      use gmm_itf_mod

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_u, F_v, F_zd, F_w, F_tv
      real, dimension(Minx:Maxx,Miny:Maxy),    intent(inout) :: F_s

      !object
      !=====================================================================================
      !     Apply vertical diffusion over column on Winds, Potential temperature and Tracers
      !     (Based on eqspng_drv)
      !=====================================================================================

      !arguments
      !______________________________________________________________________
      !        |                                             |           |   |
      ! NAME   |             DESCRIPTION                     | DIMENSION |I/O|
      !--------|---------------------------------------------|-----------|---|
      ! F_u    | x component of velocity                     | 3D (Nk)   |i/o|
      ! F_v    | y component of velocity                     | 3D (Nk)   |i/o|
      ! F_zd   | dZeta/dt                                    | 3D (Nk)   |i/o|
      ! F_w    | dz/dt                                       | 3D (Nk)   |i/o|
      ! F_tv   | Virtual temperature                         | 3D (Nk)   |i/o|
      ! F_s    | log (surface pressure/pref)                 | 2D (1)    |i  |
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

      !---------------------------------------------------------------

      integer :: i,j,k,n,istat
      real, pointer, dimension(:,:,:) :: tr
      real  pp(Minx:Maxx,Miny:Maxy,Nk),th_p(Minx:Maxx,Miny:Maxy,Nk)

      !---------------------------------------------------------------

      !-----------------
      !Diffuse U,V,W,ZD
      !-----------------
      if (Dcmip_wd_L) then

         istat = gmm_get(gmmk_uref_s , uref )
         istat = gmm_get(gmmk_vref_s , vref )
         istat = gmm_get(gmmk_wref_s , wref )
         istat = gmm_get(gmmk_zdref_s, zdref)

         call dcmip_vrd_fld (F_u,  uref,Dcmip_ref_wd,Dcmip_cp_wd_m,Dcmip_cm_wd_m,Minx,Maxx,Miny,Maxy,Nk)
         call dcmip_vrd_fld (F_v,  vref,Dcmip_ref_wd,Dcmip_cp_wd_m,Dcmip_cm_wd_m,Minx,Maxx,Miny,Maxy,Nk)

         call dcmip_vrd_fld (F_w,  wref,Dcmip_ref_wd,Dcmip_cp_wd_t,Dcmip_cm_wd_t,Minx,Maxx,Miny,Maxy,Nk)
         call dcmip_vrd_fld (F_zd,zdref,Dcmip_ref_wd,Dcmip_cp_wd_t,Dcmip_cm_wd_t,Minx,Maxx,Miny,Maxy,Nk)

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

            !Theta=t/pi; pi=(p/p0)**cappa; p=exp(a+b*s); p0=1.
            !-------------------------------------------------
            pp  (i,j,k) = exp(cappa_8*(Ver_z_8%t(k)+Ver_b_8%t(k)*F_s(i,j)))

            th_p(i,j,k) = F_tv(i,j,k)/pp(i,j,k)

         end do
         end do
         end do

         istat = gmm_get(gmmk_thref_s, thref)

         call dcmip_vrd_fld (th_p,thref,Dcmip_ref_th,Dcmip_cp_th_t,Dcmip_cm_th_t,Minx,Maxx,Miny,Maxy,Nk)

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

         nullify (tr)
         istat = gmm_get('TR/'//trim(Tr3d_name_S(n))//':P',tr)

         if (Tr3d_name_S(n)=="HU") istat = gmm_get(gmmk_qvref_s, qvref)
         if (Tr3d_name_S(n)=="QC") istat = gmm_get(gmmk_qcref_s, qcref)
         if (Tr3d_name_S(n)=="RW") istat = gmm_get(gmmk_qrref_s, qrref)

         if (Tr3d_name_S(n)=="HU") call dcmip_vrd_fld (tr,qvref,Dcmip_ref_tr,Dcmip_cp_tr_t,Dcmip_cm_tr_t,Minx,Maxx,Miny,Maxy,Nk)
         if (Tr3d_name_S(n)=="QC") call dcmip_vrd_fld (tr,qcref,Dcmip_ref_tr,Dcmip_cp_tr_t,Dcmip_cm_tr_t,Minx,Maxx,Miny,Maxy,Nk)
         if (Tr3d_name_S(n)=="RW") call dcmip_vrd_fld (tr,qrref,Dcmip_ref_tr,Dcmip_cp_tr_t,Dcmip_cm_tr_t,Minx,Maxx,Miny,Maxy,Nk)

!$omp parallel do private(i,j,k) shared(tr)
         do k=1,Nk
            do j=1,l_nj
               do i=1,l_ni
                  tr(i,j,k)=max(0.,tr(i,j,k))
               end do
            end do
         end do
!$omp end parallel do

      end do

      end if

      return

      end subroutine dcmip_vrd_drv

!--------------------------------------------------------------------------------

!**s/r dcmip_vrd_fld - Apply vertical diffusion over column on a given fld (Based on eqspng_drv)

      subroutine dcmip_vrd_fld (F_fld,F_ref_fld,F_ref_on,F_cp,F_cm,Minx,Maxx,Miny,Maxy,Nk)

      use glb_ld

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer Minx,Maxx,Miny,Maxy,Nk

      real F_fld(Minx:Maxx,Miny:Maxy,Nk),F_ref_fld(Minx:Maxx,Miny:Maxy,Nk), &
           F_cm(Nk),F_cp(Nk),F_ref_on

      !object
      !==============================================================================
      !     Apply vertical diffusion over column on a given fld (Based on eqspng_drv)
      !==============================================================================

      !------------------------------------------------

      real, dimension(:,:,:), pointer :: w_fld,p_fld
      integer :: i,j,k,km,kp

      !------------------------------------------------

      allocate (w_fld(l_ni,l_nj,Nk),p_fld(l_ni,l_nj,Nk))

      p_fld (1:l_ni,1:l_nj,1:Nk) = F_fld (1:l_ni,1:l_nj,1:Nk) - F_ref_on * F_ref_fld(1:l_ni,1:l_nj,1:Nk)

!$omp parallel do private(kp,km,i,j,k) shared(w_fld,p_fld)
      do k=1,Nk
         kp=min(Nk,k+1)
         km=max(1,k-1)
         do j=1,l_nj
            do i=1,l_ni
               w_fld(i,j,k)=p_fld(i,j,k)+(F_cp(k)*(p_fld(i,j,kp)-p_fld(i,j,k )) &
                                         -F_cm(k)*(p_fld(i,j,k )-p_fld(i,j,km)))
            end do
         end do
      end do
!$omp end parallel do

      F_fld(1:l_ni,1:l_nj,1:Nk) = w_fld(1:l_ni,1:l_nj,1:Nk) + F_ref_on * F_ref_fld(1:l_ni,1:l_nj,1:Nk)

      deallocate(w_fld,p_fld)

      return

      end subroutine dcmip_vrd_fld
