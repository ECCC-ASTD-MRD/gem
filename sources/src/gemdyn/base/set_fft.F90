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
!**s/r set_fft - Check for fast fourier transforms

      integer function set_fft ()
      use HORgrid_options
      use dyn_fisl_options
      use glb_ld
      use lun
      use glb_pil
      use fft
      implicit none

      integer npts,onept,next_down
!
!     ---------------------------------------------------------------
!
      set_fft    = 0
      Fft_fast_L = .false.

      if (trim(Sol_type_S) == 'DIRECT' ) then

         if (Lun_out > 0) write(Lun_out,1000)

         onept= 0
         if (Grd_yinyang_L) onept= 1
         npts= G_ni-Lam_pil_w-Lam_pil_e+onept

         call itf_fft_nextfactor2 ( npts, next_down )

         if ( npts /= G_ni-Lam_pil_w-Lam_pil_e+onept ) then
            if (Lun_out > 0) write (Lun_out,3001) &
            G_ni-Lam_pil_w-Lam_pil_e,npts-onept,next_down-onept
            if (sol_fft_L) then
               if (Lun_out > 0) write (Lun_out,3002)
               set_fft    = -1
!               call gem_error ( -1,'SET_FFT', '-- ABORTING --')
            end if
         else
            Fft_fast_L = .true.
         end if

         if (Lun_out > 0) write(Lun_out,*) 'Fft_fast_L = ',Fft_fast_L

      end if

 1000 format( &
      /,'VERIFYING G_NI FOR FFT FACTORIZATION (S/R SET_FFT)', &
      /,'==================================================')
 3001 format (' Fft_fast_L = .false. ====> G_NI = ',i6,' NOT FACTORIZABLE' &
              /' Neighboring factorizable G_NIs are: ',i6,' and',i6)
 3002 format (' USER REQUESTED MANDATORY FFT with sol_fft_L=.t. but'/&
              ' G_NI IS NOT FACTORIZABLE')
!
!     ---------------------------------------------------------------
!
      return
      end function set_fft
