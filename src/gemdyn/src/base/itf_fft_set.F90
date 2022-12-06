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

!**s/r itf_fft_set

      subroutine itf_fft_set ( F_dim, F_type_S, F_pri_8 )
      use tdpack
      use glb_ld
      use glb_pil
      use gem_fft
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_type_S
      integer, intent(in) ::  F_dim
      real(kind=REAL64), intent(out) ::  F_pri_8
      integer :: istat, npts

      ! With an FFTW backend, this routine needs only to save:
      ! * The transform type
      ! * The transform size
      ! * The normalization constant

      ! In particular, no calculation of trigonometric factors
      ! is necessary.
!
!----------------------------------------------------------------------
!
      istat = 0
      select case (trim(F_type_S))
         case ('PERIODIC')
            npts    = F_dim
            F_pri_8 = dble(npts) / ( 2.0d0 * pi_8 )
            
         case ('SIN')
            npts    = F_dim + 1
            F_pri_8 = dble(npts)/(G_xg_8(G_ni-Lam_pil_e)-G_xg_8(Lam_pil_w-1))
         case ('QSIN' , 'QCOS')
            npts    = F_dim
            F_pri_8 = dble(npts)/(G_xg_8(G_ni-Lam_pil_e)-G_xg_8(Lam_pil_w))
         case DEFAULT
            npts = -1
            istat = -1
            
      end select
      call handle_error(istat,'itf_fft_set', &
           'received invalid F_type_S'//trim(F_type_S))

      Fft_type_S = F_type_S
      Fft_n      = npts
!
!----------------------------------------------------------------------
!
      return
      end subroutine itf_fft_set
