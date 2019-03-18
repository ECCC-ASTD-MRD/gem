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

      subroutine spn_calfiltre ( F_nis, F_njs )
      use dcst
      use spn_work_mod
      use HORgrid_options
      use spn_options
      use tdpack
      use lun
      use glb_pil
      implicit none
#include <arch_specific.hf>

      integer F_nis,F_njs
!
!author
!     Minwei Qian (CCRD) & Bernard Dugas, Syed Husain  (MRB)  - summer 2015
!
!revision
! v4_80 - Qian, Dugas, Hussain            - initial version
! v4_80 - Baek - clarification


      integer i,j
      integer ni_trunc, ni_truncx, nj_trunc
      real*8 nix, njx, nkx, nk_cut
      real*8 WXL, WXS, DX, DY
!
!----------------------------------------------------------------------
!
      DX = Grd_dx*pi_8*Dcst_rayt_8/(180.*1000.)
      DY = Grd_dy*pi_8*Dcst_rayt_8/(180.*1000.)
      WXL= Spn_cutoff_scale_large
      WXS= Spn_cutoff_scale_small

      if (Lun_out > 0) write(Lun_out,1000) WXL
      if (Lun_out > 0) write(Lun_out,1001) WXS

      ni_trunc = int(DX*F_nis/WXL)
      nj_trunc = int(DY*F_njs/WXL)
      ni_truncx= int(DX*F_nis/WXS)

      nk_cut = float(ni_truncx)/float(ni_trunc)

      ! BAEK debug
      if (Lun_out > 0) write(Lun_out,*) "ni_trunc: ",ni_trunc
      if (Lun_out > 0) write(Lun_out,*) "ni_truncx: ",ni_truncx
      if (Lun_out > 0) write(Lun_out,*) "nj_trunc: ",nj_trunc
      if (Lun_out > 0) write(Lun_out,*) "nk_cut: ",nk_cut

      fxy = 0.

      do j=0,F_njs-1
         do i=0,F_nis-1
            nix = dble(i)/dble(ni_trunc)
            njx = dble(j)/dble(nj_trunc)

            nkx = sqrt(nix*nix + njx*njx)

            if ( nkx > nk_cut ) then

               fxy(i+Lam_pil_w+1,j+Lam_pil_s+1) = 0.0

            else if ( nkx > 1.0 ) then

               fxy(i+Lam_pil_w+1,j+Lam_pil_s+1) = &
               (cos( (pi_8/2.0) * ((nkx-1.)/(nk_cut-1.)) ))**2

            else

               fxy(i+Lam_pil_w+1,j+Lam_pil_s+1) = 1.0

            end if

         end do
      end do

 1000 format(/' In SPN_CALFILTRE, Spn_cutoff_scale_L = ',f7.2/)
 1001 format(/' In SPN_CALFILTRE, Spn_cutoff_scale_S = ',f7.2/)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_calfiltre


