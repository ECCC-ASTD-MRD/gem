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

!**s/p set_filt - initialize list of variables for filtering before output
      integer function set_filt (F_argc,F_argv_S,F_cmdtyp_S,F_v1,F_v2)
      use gem_options
      use glb_ld
      use lun
      use out3
      implicit none
#include <arch_specific.hf>
!
      integer F_argc,F_v1,F_v2
      character(len=*) F_argv_S(0:F_argc),F_cmdtyp_S
!
!author Vivian Lee - rpn - July 2004
!
!revision
! v3_20 - Lee V.            - initial MPI version
!
!object
!       initialization of the common blocks OUT3. This function is
!       called when the keyword "filtre" is found in the first word
!       of the directives in the input file given in the statement
!       "process_f_callback". This feature is enabled by the
!       ARMNLIB "rpn_fortran_callback" routine (called in "srequet")
!       which allows a different way of passing user directives than
!       the conventional FORTRAN namelist. This function will process
!       the following example command read from the named input file.
!
! ie:   filtre([UU,VV,TT],coef,.1,pass,3)
!
!       The "rpn_fortran_callback" routine will process the above
!       statement and return 5 arguments to this function. For more
!       information to how this is processed, see "SREQUET".
!
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_argc       I    - number of elements in F_argv_S
! F_argv_S     I    - array of elements received
!                     if F_argv_S(ii) contains "[", the value in this
!                     argument indicates number of elements following it
! F_cmdtyp_S   I    - character command type - not used
! F_v1         I    - integer parameter 1 - not used
! F_v2         I    - integer parameter 2 - not used
!----------------------------------------------------------------


      character(len=5) :: stuff_S
      integer pass, varmax
      real coefficient
      integer i, j, ii, jj
!
!----------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,*)
      if (Lun_out > 0) write(Lun_out,*) F_argv_S
      set_filt=0
      Out3_filtpass_max = Out3_filtpass_max + 1
      if (Out3_filtpass_max > MAXELEM) then
          if (Lun_out > 0) write(Lun_out,*)'SET_FILT WARNING: Too many definitions filter'
          Out3_filtpass_max = Out3_filtpass_max - 1
          set_filt = 1
          return
      end if

      if (index(F_argv_S(1),'[') > 0) then
          stuff_S=F_argv_S(1)
          read(stuff_S(2:4),*) varmax
      else
         if (Lun_out > 0) then
            write(Lun_out,*) 'set_filt WARNING: syntax incorrect'
         end if
         set_filt=1
         Out3_filtpass_max = Out3_filtpass_max - 1
        return
      end if
!
!     Obtain coefficient,npass values for filter
!
      pass=0
      coefficient=0.0
      do i=varmax+2, F_argc
         if (F_argv_S(i) == 'coef') then
            read(F_argv_S(i+1),*) coefficient
         else if (F_argv_S(i) == 'pass') then
            read(F_argv_S(i+1),*) pass
         end if
      end do

      if (pass <= 0) then
          if (Lun_out > 0) write(Lun_out,*) 'set_filt WARNING: Number of passes not chosen'
         set_filt=1
          Out3_filtpass_max = Out3_filtpass_max - 1
         return
      end if

      if (coefficient == 0.0) then
          if (Lun_out > 0) write(Lun_out,*) 'set_filt WARNING: Filtering coefficient=0.0'
         set_filt=1
          Out3_filtpass_max = Out3_filtpass_max - 1
         return
      end if
!
!     Store variables in variable sets
!
      j = Out3_filtpass_max + varmax

      if (j > MAXELEM) then
          if (Lun_out > 0) write(Lun_out,*) 'set_filt WARNING: too many variables to filt'
          set_filt=1
          return
      end if
!
      jj=Out3_filtpass_max
      do ii=1,varmax
         jj = jj+1
         Out3_filt_S(jj)    = F_argv_S(ii+1)
         Out3_filtpass(jj)  = pass
         Out3_filtcoef(jj)  = coefficient
      end do
      Out3_filtpass_max = jj
!
!----------------------------------------------------------------
!
      return
      end
