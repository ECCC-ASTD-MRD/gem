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

!**s/r theo_hdif_main - applies horizontal diffusion on a given set of fields
!
      subroutine theo_hdif_main
      use gmm_vt1
      use HORgrid_options
      use glb_ld
      use rmn_gmm
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

!author
!
!revision

      integer istat
!     _________________________________________________________________
!
      istat = gmm_get(gmmk_ut1_s,ut1)
      istat = gmm_get(gmmk_zdt1_s,zdt1)
      istat = gmm_get(gmmk_wt1_s,wt1)

!
!*******************************************
!  1. One-dimensional Horizontal diffusion *
!*******************************************
!
      call theo_hdif_ctl ( ut1, l_minx,l_maxx,l_miny,l_maxy, G_nk)
      call theo_hdif_ctl (zdt1, l_minx,l_maxx,l_miny,l_maxy, G_nk)
      call theo_hdif_ctl ( wt1, l_minx,l_maxx,l_miny,l_maxy, G_nk)
!
!     _________________________________________________________________
!
      return
      end

!**s/r theo_hdif_ctl - applies horizontal explicit diffusion
!
      subroutine theo_hdif_ctl(F_f2dif, Minx,Maxx,Miny,Maxy, NK)
      use gem_options
      use glb_ld
      use rmn_gmm
      use theo_dif
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer       Minx,Maxx,Miny,Maxy,Nk
      real          F_f2dif (Minx:Maxx,Miny:Maxy,Nk)
!
!AUTHOR    C. Girard
!

      integer nn, mm
      real wk1(l_minx:l_maxx,l_miny:l_maxy,Nk)
      real(kind=REAL64) :: nu_dif, lnr, pwr

      real(kind=REAL64), parameter :: epsilon = 1.0d-12, pt25=0.25d0
!     __________________________________________________________________
!
      lnr = hdif_lnr
      pwr = hdif_pwr

      nu_dif = pt25*lnr**(2.d0/pwr)
      nu_dif  = min ( nu_dif, pt25-epsilon )
      if (nu_dif < 1.0e-10) return

      nn = pwr/2

      call rpn_comm_xch_halo ( F_f2dif, l_minx,l_maxx,l_miny,l_maxy,&
           l_ni,l_nj, Nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      do mm=1,nn

         call theo_hdif(F_f2dif, wk1, l_minx,l_maxx,l_miny,l_maxy,&
                                                   Nk, nu_dif, mm,nn )
         if (mm /= nn) then
              call rpn_comm_xch_halo( wk1, l_minx,l_maxx,l_miny,l_maxy,&
                                      l_ni,l_nj, Nk, G_halox,G_haloy,  &
                                      G_periodx,G_periody,l_ni,0)
         end if

      end do
!     __________________________________________________________________
!
      return
      end

!**s/r theo_hdif

      subroutine theo_hdif(rfd,sfd,Minx,Maxx,Miny,Maxy,lnk,nu_dif,m,n)
      use gem_options
      use glb_ld
      use rmn_gmm
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      integer Minx,Maxx,Miny,Maxy,lnk,m,n
      real rfd (Minx:Maxx,Miny:Maxy,lnk), &
           sfd (Minx:Maxx,Miny:Maxy,lnk)
      real(kind=REAL64) nu_dif

!AUTHOR    C. Girard
!
!revision


      integer i,j,k,id,jd,iff,jf,i0,in,j0,jn
      real wk(l_minx:l_maxx,l_miny:l_maxy)
      real(kind=REAL64), parameter :: two = 2.d0
!
!----------------------------------------------------------------------
!
      id = 1    + pil_w
      jd = 1    + pil_s
      iff= l_ni - pil_e
      jf = l_nj - pil_n

      do k=1,lnk
         i0= id-1 ; in= iff+1
         j0= jd-1 ; jn= jf +1
         if (m == 1) then
            sfd(i0:in,j0:jn,k) = rfd(i0:in,j0:jn,k)
         else if (m == 2) then
            sfd(i0:in,j0:jn,k) = rfd(i0:in,j0:jn,k) - sfd(i0:in,j0:jn,k)
         else
            i0 = i0+west  ; in= in-east
            j0 = j0+south ; jn= jn-north
            sfd(i0:in,j0:jn,k) = rfd(i0:in,j0:jn,k) - sfd(i0:in,j0:jn,k)
         end if
         if (m == n) then
            do j=jd,jf
               do i=id,iff
                  rfd(i,j,k)= rfd(i,j,k) +  &
                              nu_dif*(sfd(i+1,j,k)+sfd(i-1,j,k)-two*sfd(i,j,k))
               end do
            end do
         else
            do j=jd,jf
               do i=id,iff
                  wk(i,j) = nu_dif*(sfd(i+1,j,k)+sfd(i-1,j,k)-2.d0*sfd(i,j,k))
               end do
            end do
            sfd(id:iff,jd:jf,k) = rfd(id:iff,jd:jf,k) + wk(id:iff,jd:jf)
         end if
      end do
!
!----------------------------------------------------------------------
!
      return
      end
