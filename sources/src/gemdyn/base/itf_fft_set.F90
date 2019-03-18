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
!
      subroutine itf_fft_set ( F_dim, F_type_S, F_pri_8 )
      use tdpack
      use glb_ld
      use glb_pil
      use fft
      implicit none
#include <arch_specific.hf>
!
      character(len=*) F_type_S
      integer       F_dim
      real*8        F_pri_8
!
!author
!     Michel Desgagne - spring 2012
!
!revision
! v4_50 - Desgagne M.       - initial version


      integer, save :: same = -1
      integer i, npts, istat
      real*8  del, angle
      real*8, parameter :: half=0.5, one=1.0, two=2.0
!     __________________________________________________________________
!
      istat = 0
      select case (trim(F_type_S))
         case ('PERIODIC')
            npts    = F_dim
            F_pri_8 = dble(npts) / ( two * pi_8 )
         case ('SIN')
            npts    = F_dim + 1
            F_pri_8 = dble(npts)/(G_xg_8(G_ni-Lam_pil_e)-G_xg_8(Lam_pil_w-1))
         case ('COS')
            npts    = F_dim - 1
         case ('QSIN' , 'QCOS')
            npts    = F_dim
            F_pri_8 = dble(npts)/(G_xg_8(G_ni-Lam_pil_e)-G_xg_8(Lam_pil_w))
         case DEFAULT
            istat = -1
      end select
      call handle_error(istat,'itf_fft_set', &
           'received invalid F_type_S'//trim(F_type_S))

      if ((npts==same) .and. (trim(F_type_S)==trim(Fft_type_S))) return

      Fft_type_S = F_type_S
      Fft_n      = npts
      same       = npts

      if (trim(Fft_type_S) /= 'PERIODIC') then
         Fft_m      = Fft_n/2
         Fft_nstore = Fft_n + 2

         if (associated(Fft_ssin)) deallocate(Fft_ssin,stat=istat)
         if (associated(Fft_ccos)) deallocate(Fft_ccos,stat=istat)
         if (associated(Fft_qsin)) deallocate(Fft_qsin,stat=istat)

         allocate (Fft_ssin(Fft_n-Fft_m-1), &
                   Fft_ccos(Fft_n-Fft_m-1), Fft_qsin(0:Fft_m-1))

         del = acos( - one )/Fft_n

         do i=1,Fft_n-Fft_m-1
            angle = i * del
            Fft_ccos( i ) = cos( angle )
            Fft_ssin( i ) = sin( angle )
         end do

         do i=0,Fft_m-1
            Fft_qsin( i ) = sin( ( i + half ) * del )
         end do
      end if

      call setfft8( Fft_n )
!     __________________________________________________________________
!
      return
      end
