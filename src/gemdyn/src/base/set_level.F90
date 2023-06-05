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

!**s/r set_level - initialization of common block LEVEL
!

!
      integer function set_level (F_argc,F_argv_S,F_cmdtyp_S,F_v1,F_v2)
!
      use gem_options
      use glb_ld
      use lun
      use dimout
      use out_mod
      use levels
      implicit none
#include <arch_specific.hf>
!
        integer F_argc,F_v1,F_v2
        character(len=*) F_argv_S(0:F_argc),F_cmdtyp_S
!
!author Vivian Lee - RPN - April 1999
!
!revision
! v2_00 - Lee V.         - initial MPI version
! v2_10 - Lee V.         - corrected so that all eta levels are
! v2_10                    outputted when a "-1" is indicated
! v2_21 - Dugas B.       - use convip
! v2_30 - Lee V.         - reduced dimension of Level_typ to 1
! v2_32 - Lee V.         - levset is now an ID defined by user, not the
! v2_32                    actual "set" number forced to be in sequence
! v3_01 - Lee V.         - new ip1 encoding (kind=5 -- unnormalized)
! v3_02 - Lee V.         - eliminate levels repeated in one level set
! v3_21 - Lee V.         - bug correction when levdesc=1
! v4_40 - Lee V.         - modification in order to select the "bot" levels
!
!
!object
!       initialization of the common block LEVEL. This function is
!       called when the keyword "levels" is found in the first word
!       of the directives in the input file given in the statement
!       "process_f_callback". This feature is enabled by the
!       ARMNLIB "rpn_fortran_callback" routine (called in "srequet")
!       which allows a different way of passing user directives than
!       the conventional FORTRAN namelist. This function will process
!       the following example command read from the named input file.
!
! ie:  levels=1,pres,[1000.,925.,850.];
!
!       The "rpn_fortran_callback" routine will process the above
!       statement and return 5 arguments to this function. For more
!       information to how this is processed, see "SREQUET".
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_argc       I    - number of elements in F_argv_S
! F_argv_S     I    - array of elements received
!                     if F_argv_S(ii) contains "[", the value in this
!                     argument indicates number of elements following it.
! F_cmdtyp_S   I    - character command type - not used
! F_v1         I    - integer parameter 1 - not used
! F_v2         I    - integer parameter 2 - not used
!----------------------------------------------------------------
!
!Notes:
!
!    levels=levelset#,pres/eta/arbitrary,{list};
! ie:  levels=2,eta,[1,5,10];
!      levels=3,eta,<1,28,2>;
!      levels=4,eta,-1;
!
!      Should label the levelset# sequentially: 1,2,3,....
!      'eta'       - model levels (eta)
!      'pres'      - pressure (hPa)
!      '-1' with "eta" levels will give all model levels.
!      [a,b,c] means level a,b and c are requested
!      <a,b,c> means levels a to b, incrementing every c are requested
!

!
!*
      logical :: press_L,eta_L
      character(len=5) :: stuff_S
      integer :: i,j,k,ii,idx,levset,num,levdesc
      integer, dimension(size(Level_allpres)) :: ip1_stub
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) then
          write(Lun_out,*)
          write(Lun_out,*) F_argv_S(0),'=',F_argv_S(1),',',F_argv_S(2),',',(F_argv_S(i),i=3,F_argc)
      end if
      Level_momentum = G_nk
      Level_thermo = G_nk
      set_level = 0
      read( F_argv_S(1), * ) levset
      Level_sets = Level_sets + 1

      if (Level_sets > MAXSET) then

          if (Lun_out > 0) then
            write(Lun_out,*)'SET_LEVEL WARNING: Too many sets of LEVELS'
          end if
          Level_sets = Level_sets -1
          set_level=1
          return

      end if

      j=Level_sets

      i=0
      Level_id(j)=levset

!      i is the counter for the number of levels requested

      press_L = .false.
      eta_L   = .false.
      levdesc = -1

      each_lev: do ii=2,F_argc

         if (index(F_argv_S(ii),'[') > 0) then

            stuff_S=F_argv_S(ii)
            read( stuff_S(2:4), * ) num

         else if (F_argv_S(ii) == 'eta') then
            if (press_L) then
                if (Lun_out > 0) then
                  write(Lun_out,*) 'SET_LEVEL WARNING:  Pressure levels are already defined'
                end if
                Level_sets = Level_sets -1
                set_level=1
                return
            end if
            levdesc = 1
            eta_L = .true.

         else if (F_argv_S(ii) == 'pres') then
            if (eta_L) then
               if (Lun_out > 0) then
                  write(Lun_out,*) 'SET_LEVEL WARNING: Model levels are already defined'
               end if
               Level_sets = Level_sets -1
               set_level=1
               return
            end if

            levdesc = 2
            press_L = .true.

         else if (F_argv_S(ii) == 'bot') then
            if (press_L) then
                if (Lun_out > 0) then
                  write(Lun_out,*) 'SET_LEVEL WARNING: Pressure levels are already defined'
                end if
                Level_sets = Level_sets -1
                set_level=1
                return
            end if
            levdesc = 3
            eta_L = .true.

         else if (levdesc == 1) then

            Level_typ_S(j)='M'
            i = i+1
            read( F_argv_S(ii), * ) Level(i,j)

            if (Level(i,j) == -1) then

!              request for all model eta levels
!              will put in the max levels which is the number of thermo

                  i = i-1
                  do idx=1,Level_thermo
                     i = i+1
                     Level(i,j) = float( idx )
                  end do
            end if

         else if (levdesc == 3) then
                  Level_typ_S(j)='M'
                  i = i+1
                  read( F_argv_S(ii), * ) Level(i,j)
                  k=nint(Level(i,j))
!              request for model levels close and include the surface
                  if (k > 0.and. k < Level_thermo) then
                     i=i-1
                     do idx=k,1,-1
                        i = i+1
                        Level(i,j) = -1.0*idx + 1.0
                     end do
                  else
                     if (Lun_out > 0) then
                        write(Lun_out,*) 'SET_LEVEL WARNING: Level index out of range'
                     end if

                     i = i - 1
                  end if
         else if (levdesc == 2) then

                  Level_typ_S(j)='P'
                  i = i+1
                  read( F_argv_S(ii), * ) Level(i,j)
                  Level_allpres(Level_npres+i) = Level(i,j)

         else

                  if (Lun_out > 0) then
                     write(Lun_out,*) 'SET_LEVEL WARNING: Level type not recognizable'
                  end if
                  Level_sets = Level_sets -1
                  set_level=1
                  return

         end if

      end do each_lev

      if (i > MAXLEV) then

         if (Lun_out > 0) then
            write(Lun_out,*)'SET_LEVEL WARNING: Requested levels > MAXLEV'
         end if
         Level_sets = Level_sets -1
         set_level = 1
         return

      end if

      if (i == 0) then

         if (Lun_out > 0) then
            write(Lun_out,*)'SET_LEVEL WARNING: No levels requested'
         end if
         Level_sets = Level_sets -1
         set_level = 1
         return

      end if

!     Eliminate repeated levels in one Level set
      ip1_stub = 0
      call sortlev(Level(:,j),ip1_stub(1:i),i,Level_max(Level_sets))

!     Eliminate repeated levels in the full pressure list
      if (Level_typ_S(j) == 'P') then
         call sortlev(Level_allpres(1:Level_npres+i),ip1_stub(1:Level_npres+i), &
              Level_npres+i,Level_npres)
      end if

      if (Lun_out > 0) then

         write(Lun_out,*) ' Level_set(',j,') : Level_id=',Level_id(j)
         write(Lun_out,*) ' Level_type=',Level_typ_S(j)
         write(Lun_out,*) ' Level=',(Level(i,j),i=1,Level_max(j))

      end if
!
!     ---------------------------------------------------------------
!
      return

    end function set_level


