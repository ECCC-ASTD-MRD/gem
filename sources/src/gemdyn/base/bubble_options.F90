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
!---------------------------------- LICENCE END ---------------------------------
module bubble_options
   use dynkernel_options
   use HORgrid_options
   use dyn_fisl_options
   use VERgrid_options
   use lun
   use glb_ld
   use cstv
   use dcst
   use tdpack
   implicit none
   public
   save

   !#
   integer :: bubble_ni = 101
   namelist /bubble_cfgs/ bubble_ni
   !#
   integer :: bubble_nj = 1
   namelist /bubble_cfgs/ bubble_nj
   !#
   integer :: bubble_nk = 100
   namelist /bubble_cfgs/ bubble_nk
   !#
   real :: bubble_dx = 10.
   namelist /bubble_cfgs/ bubble_dx
   !#
   real :: bubble_dz = 10.
   namelist /bubble_cfgs/ bubble_dz
   !#
   real :: bubble_theta = 303.16
   namelist /bubble_cfgs/  bubble_theta
   !#
   integer :: bubble_rad = 25
   namelist /bubble_cfgs/ bubble_rad
   !#
   integer :: bubble_ictr = -1
   namelist /bubble_cfgs/ bubble_ictr
   !#
   integer :: bubble_kctr = -1
   namelist /bubble_cfgs/ bubble_kctr
   !#
   logical :: bubble_gaus_L = .false.
   namelist /bubble_cfgs/ bubble_gaus_L

contains

      integer function bubble_nml (F_unf)
      implicit none

      integer F_unf

      logical nml_must
      character(len=64) :: nml_S
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         bubble_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=bubble_cfgs)
         return
      end if

      bubble_nml= -1 ; nml_must= .true. ; nml_S= 'bubble_cfgs'

      rewind(F_unf)
      read (F_unf, nml=bubble_cfgs, end= 1001, err=1003)
      bubble_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         bubble_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (bubble_nml < 0 ) return
      if ((Lun_out>=0).and.(bubble_nml==0)) write (Lun_out, 6004) trim(nml_S)

      ! establish horizontal grid configuration
      ! (must absolutely be done here)
      Dcst_rayt_8 = Dcst_rayt_8*0.1d0 ! an accuracy problem
      Dcst_inv_rayt_8 = Dcst_inv_rayt_8 * 10.d0 ! an accuracy problem
      Grd_typ_S='LU'
      Grd_ni = bubble_ni ; Grd_nj = bubble_nj
      Grd_dx = (bubble_dx/Dcst_rayt_8)*(180./pi_8)
      Grd_dy = Grd_dx
      Grd_latr = 0.
      Grd_lonr = (bubble_ni/2 + 20) * Grd_dx
      Grd_maxcfl = 3

      bubble_nml=0

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      return
      end function bubble_nml
!
!-------------------------------------------------------------------
!
      integer function bubble_cfg()
      implicit none
#include <arch_specific.hf>

      integer k
      real*8 c1_8,Exner_8,height_8,pres_8,pref_8,ptop_8,htop_8
!
!     ---------------------------------------------------------------
!
      bubble_cfg = -1

      ! establish vertical grid configuration
      pref_8 = 1.d5

      G_nk   = bubble_nk
      htop_8 = G_nk*bubble_dz

      if ( hyb(1) < 0 ) then

        !isentropic case
         c1_8=grav_8/(cpd_8*bubble_theta)
         Exner_8=1.d0-c1_8*htop_8
         ptop_8 = Exner_8**(1.d0/cappa_8)*pref_8
!        Uniform distribution of levels in terms of height
         do k=1,G_nk
            height_8=htop_8*(1.d0-(dble(k)-.5d0)/G_nk)
            Exner_8=1.d0-c1_8*height_8
            pres_8=Exner_8**(1.d0/cappa_8)*pref_8
            hyb(k)=(pres_8-ptop_8)/(pref_8-ptop_8)
            hyb(k) = hyb(k) + (1.-hyb(k))*ptop_8/pref_8
         end do

      else

         do k=1024,1,-1
            if(hyb(k) < 0) G_nk=k-1
         end do

      end if

      if (bubble_ictr < 0) bubble_ictr = int(float(Grd_ni-1)*0.5)+1
      if (bubble_kctr < 0) bubble_kctr = G_nk - bubble_rad - 1

      bubble_cfg = 1

      return
      end function bubble_cfg
!
!-------------------------------------------------------------------
!

!**s/r bubble_data - generates initial condition for Robert's bubble
!                    experiment (Robert 1993 JAS)
!
      subroutine bubble_data ( F_u, F_v, F_t, F_s, F_q, F_topo,&
                               Mminx,Mmaxx,Mminy,Mmaxy,nk )
      use gmm_geof
      use geomh
      use gmm_itf_mod
      use glb_pil
      use ptopo
      use type_mod
      use ver
      implicit none
#include <arch_specific.hf>

      integer Mminx,Mmaxx,Mminy,Mmaxy,nk
      real F_u    (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_v    (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_t    (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_s    (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_q    (Mminx:Mmaxx,Mminy:Mmaxy,nk+1)

      integer :: i,j,k,istat,ii,err
      real*8 :: pp, ex, theta, r,rad
!
!     ---------------------------------------------------------------
!
      istat = gmm_get (gmmk_sls_s     ,   sls )

      sls   (:,:) = 0.0
      F_topo(:,:) = 0.0
      F_s   (:,:) = 0.0
      F_u (:,:,:) = 0.0
      F_v (:,:,:) = 0.0
      F_q (:,:,:) = 0.0
!
!---------------------------------------------------------------------
!     Initialize temperature
!---------------------------------------------------------------------
!
      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

         if(.not.bubble_gaus_L) then
         do k=1,g_nk
            do j=1,l_nj
            do i=1,l_ni
               ii=i+l_i0-1
               theta=bubble_theta
               if ( (((ii)-bubble_ictr)**2 +((k)-bubble_kctr)**2) < bubble_rad**2 ) then
                  theta=theta+0.5d0
               end if
               pp = exp(Ver_a_8%t(k)+Ver_b_8%t(k)*F_s(i,j))
               ex = (pp/Cstv_pref_8)**cappa_8
               F_t(i,j,k)=theta*ex
            end do
            end do
         end do
         else
            !Gaussian
            do k=1,g_nk
               do j=1,l_nj
               do i=1,l_ni
                  ii=i+l_i0-1
                  theta=bubble_theta
                  r = sqrt( (dble((ii)-bubble_ictr)**2 + dble((k)-bubble_kctr)**2) )*dble(bubble_dx)
                  rad = bubble_rad*dble(bubble_dx)
                  if ( r <= rad ) then
                     theta = bubble_theta + 0.5d0
                  else
                     theta = bubble_theta + 0.5d0 * exp( -(r - rad)**2 / 100.d0**2 );
                  end if
                  pp = exp(Ver_a_8%t(k)+Ver_b_8%t(k)*F_s(i,j))
                  ex = (pp/Cstv_pref_8)**cappa_8
                  F_t(i,j,k)=theta*ex
               end do
               end do
            end do
         end if

      end if

      err=0
      if ((l_north).and.(l_nj-2*pil_n+1<1)) err=-1
      if ((l_east ).and.(l_ni-2*pil_e+1<1)) err=-1
      if ((l_south).and.(2*pil_s>l_nj)    ) err=-1
      if ((l_west ).and.(2*pil_w>l_ni)    ) err=-1
      call gem_error(err,'ABORT in BUBBLE_DATA',&
          'Partitionning NOT allowed for MIRROR')
!
 9000 format(/,'CREATING INPUT DATA FOR BUBBLE THEORETICAL CASE' &
            /,'================================================')
!
!     -----------------------------------------------------------------
!
      return
      end subroutine bubble_data

end module bubble_options
