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

!**s/p set_xnbit - initialize list of variables required different compaction
!

!
      integer function set_xnbit (F_argc,F_argv_S,F_cmdtyp_S,F_v1,F_v2)
!
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
!author Vivian Lee - rpn - JULY 2004
!
!revision
! v3_20 - Lee V.            - initial MPI version
!
!object
!       initialization of the common blocks OUT3. This function is
!       called when the keyword "xnbit" is found in the first word
!       of the directives in the input file given in the statement
!       "process_f_callback". This feature is enabled by the
!       ARMNLIB "rpn_fortran_callback" routine (called in "srequet")
!       which allows a different way of passing user directives than
!       the conventional FORTRAN namelist. This function will process
!       the following example command read from the named input file.
!
! ie:   xnbit([Z0,MX],bits,24)
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
!
!
      character(len=5) :: stuff_S
      integer bits,varmax
      integer i, j, ii, jj
!
!----------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,*)
      if (Lun_out > 0) write(Lun_out,*) F_argv_S
      set_xnbit=0
      Out3_xnbits_max = Out3_xnbits_max + 1
      if (Out3_xnbits_max > MAXELEM) then
         if (Lun_out > 0) then
            write(Lun_out,*)'set_xnbit WARNING: Too many definitions XNBITS'
         end if
         Out3_xnbits_max = Out3_xnbits_max - 1
         set_xnbit = 1
         return
      end if

      if (index(F_argv_S(1),'[') > 0) then
         stuff_S=F_argv_S(1)
         read(stuff_S(2:4),*) varmax
      else
         if (Lun_out > 0) then
            write(Lun_out,*) 'set_xnbit WARNING: syntax incorrect'
         end if
         set_xnbit=1
         Out3_xnbits_max = Out3_xnbits_max - 1
         return
      end if
!
!     Obtain compaction bits request
!
      bits=0
      do i=varmax+2, F_argc
         if (F_argv_S(i) == 'bits') read(F_argv_S(i+1),*) bits
      end do

      if (bits <= 0) then
         if (Lun_out > 0) then
            write(Lun_out,*) 'set_xnbit WARNING: Number of bits not chosen'
         end if
         set_xnbit=1
         Out3_xnbits_max = Out3_xnbits_max - 1
         return
      end if
!
!     Store variables in variable sets
!
      j = Out3_xnbits_max + varmax
      if (j > MAXELEM) then
         if (Lun_out > 0) then
             write(Lun_out,*) 'set_xnbit WARNING: too many variables for xnbits'
         end if
         set_xnbit=1
         return
      end if
!
      jj=Out3_xnbits_max
      do ii=1,varmax
         jj = jj+1
         Out3_xnbits_S(jj)    = F_argv_S(ii+1)
         Out3_xnbits(jj)  = bits
      end do
      Out3_xnbits_max = jj
!
!----------------------------------------------------------------
!
      return
      end
