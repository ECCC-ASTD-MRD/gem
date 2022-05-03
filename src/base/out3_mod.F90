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

module out3
   use dimout, only: MAXELEM
   implicit none
   public
   save
!
!__________________________________________________________________________
!                                                                           |
!  OUTPUT GLOBAL CONTROL VARIABLES                                          |
!__________________________________________________________________________ |
!                    |                                                      |
! NAME               | DESCRIPTION                                          |
!--------------------|----------------------------------------------------- |
! Out3_xnbits_s      | list of variables that will require a compaction     |
!                    | that is different from Out3_nbitg                    |
! Out3_conv_s        | list of variables that will require a different      |
!                    | conversion from the default set in the model         |
! Out3_filt_s        | list of variables what will require filtering;       |
!                    | for each variable, there is a corresponding          |
!                    | number of filter passes (Out3_filtpass) and the      |
!                    | filter coefficient (Out3_filtcoef)                   |
! Out3_date          | date time stamp of the run                           |
! Out3_sfcdiag_L     ! to output physics surface diagnostic level           |
!                    |        water/ice phase for saturation                |
!                    | .false. for water phase only for saturation          |
! Out3_filtpass      | number of filtering passes (see Out3_filt_s)         |
! Out3_filtcoef      | filter coefficient (see Out3_filt_s)                 |
! Out3_convmult      | muliplication factor in the unit conversion          |
! Out3_convadd       | adding factor in the unit conversion                 |
! Out3_unit_S        | output name in units from Out3_close_interval_S      |
! Out3_accavg_L      | Output control for accumulater,averaging by DYN      |
! --------------------------------------------------------------------------

   character(len=16) :: Out3_xnbits_S(MAXELEM),Out3_filt_S(MAXELEM), Out3_conv_S(MAXELEM)
   character(len=4)  :: Out3_unit_S
   logical :: Out3_accavg_L, Out3_sfcdiag_L, Out3_ezcoll_L
   integer :: Out3_conv_max, Out3_date, Out3_lieb_nk
   integer :: Out3_filtpass(MAXELEM),Out3_xnbits(MAXELEM)
   integer :: Out3_filtpass_max, Out3_xnbits_max
   integer :: Out3_comm_id, Out3_comm_setno, Out3_iome, &
              Out3_comm_io, Out3_iobcast
   real    :: Out3_filtcoef(MAXELEM), Out3_close_interval
   real    :: Out3_convmult(MAXELEM),Out3_convadd(MAXELEM)

end module out3
