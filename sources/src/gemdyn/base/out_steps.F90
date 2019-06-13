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
!**s/r out_steps -

      subroutine out_steps
      use step_options
      use init_options
      use gem_options
      use out_options
      use cstv
      use out_mod
      use out3
      use outp
      use out_listes
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=16) :: datev
      integer istep,step0,stepf
      integer, save :: marker
      logical, save :: done = .false., dgflt_L = .true.
      real,   parameter :: eps=1.e-12
      real(kind=REAL64), parameter :: OV_day = 1.0d0/86400.0d0
      real(kind=REAL64)  dayfrac
!
!     ---------------------------------------------------------------
!
      if (.not. done) then
         nullify(outd_sorties,outp_sorties,outc_sorties,outp_lasstep)
         marker= Lctl_step - 1
      end if

      if ( (.not.Init_mode_L) .and. dgflt_L) marker= Lctl_step - 1
      dgflt_L = Init_mode_L

      if (marker < Lctl_step) then

         step0 = Lctl_step
         stepf = Lctl_step + 50
         marker= stepf

         if (associated(outd_sorties)) deallocate (outd_sorties)
         if (associated(outp_sorties)) deallocate (outp_sorties)
         if (associated(outc_sorties)) deallocate (outc_sorties)
         if (associated(outp_lasstep)) deallocate (outp_lasstep)
         allocate (outd_sorties(0:MAXSET,step0:stepf))
         allocate (outp_sorties(0:MAXSET,step0:stepf))
         allocate (outc_sorties(0:MAXSET,step0:stepf))
         allocate (outp_lasstep(MAXSET,step0:stepf))

         outd_sorties(0,:)= 0
         outp_sorties(0,:)= 0
         outc_sorties(0,:)= 0

         do istep = step0, stepf
            if (.not.( Init_mode_L .and.            &
            (istep-Step_initial) >= Init_halfspan)) &
            call out_thistep (outd_sorties(0,istep),istep,MAXSET,'DYN')
            if (       Init_mode_L .and.            &
            (istep-Step_initial) >= Init_halfspan+1) cycle
            call out_thistep (outp_sorties(0,istep),istep,MAXSET,'PHY')
            call out_thistep (outc_sorties(0,istep),istep,MAXSET,'CHM')
         end do

      end if

      Out_dateo = Out3_date
      if ( lctl_step < 0 ) then  ! adjust Out_dateo because ip2=npas=0
         dayfrac = dble(lctl_step-Step_delay) * Cstv_dt_8 * OV_day
         call incdatsd (datev,Step_runstrt_S,dayfrac)
         call datp2f   (Out_dateo,datev)
      end if

      Out_ip2  = int (dble(lctl_step) * Out_deet / 3600. + eps)
      Out_ip2  = max (0, Out_ip2  )
      Out_npas = max (0, Lctl_step)

      Out_ip3  = 0
      if (Out3_ip3 == -1) Out_ip3 = max (0, Lctl_step)
      if (Out3_ip3 > 0 ) Out_ip3 = Out3_ip3

      Out_typvar_S = 'P'
      if (Lctl_step < 0) Out_typvar_S = 'I'

      done = .true.
!
!     ---------------------------------------------------------------
!
      return
      end
