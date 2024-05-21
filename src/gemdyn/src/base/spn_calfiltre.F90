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

!*s/r spn_calfiltre - compute a filter for spectral nudging

      subroutine spn_calfiltre
      use dcst
      use glb_ld
      use glb_pil
      use HORgrid_options
      use spn_options
      use tdpack
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer i,j,il,jl
      integer ni_trunc, ni_truncx, nj_trunc
      real(kind=REAL64) nix, njx, nkx, nk_cut
      real(kind=REAL64) WXL, WXS, DX, DY
!
!----------------------------------------------------------------------
!
      DX = Grd_dx*pi_8*Dcst_rayt_8/(180.*1000.)
      DY = Grd_dy*pi_8*Dcst_rayt_8/(180.*1000.)
      WXL= Spn_cutoff_scale_large
      WXS= Spn_cutoff_scale_small

      if (Lun_out > 0) write(Lun_out,1000) WXL,WXS

      ni_trunc = int(DX*(G_ni-2*Grd_extension)/WXL)
      nj_trunc = int(DY*(G_nj-2*Grd_extension)/WXL)
      ni_truncx= int(DX*(G_ni-2*Grd_extension)/WXS)

      nk_cut = float(ni_truncx)/float(ni_trunc)

      allocate (Spn_flt(Spn_22n,G_nj))
      Spn_flt= 0.
      do j=1+Grd_extension,G_nj-Grd_extension
         jl= j-Grd_extension-1
         do i=1+Spn_22pil_w,Spn_22n-Spn_22pil_e
            il= i+Spn_22n0-Grd_extension-2
            nix = dble(il)/dble(2.*ni_trunc)
            njx = dble(jl)/dble(2.*nj_trunc)
            nkx = sqrt(nix*nix + njx*njx)
            if ( nkx > nk_cut ) then
               Spn_flt(i,j)= 0.0
            else if ( nkx > 1.0 .and. nkx<=nk_cut) then
               Spn_flt(i,j)=(cos( (pi_8/2.0)* ((nkx-1.)/(nk_cut-1.))))**2
            else
               Spn_flt(i,j)= 1.0
            end if
         end do
      end do
      
 1000 format(/' Spn_calfiltre, Large,Small cutoff_scales= ',2f7.2)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_calfiltre


