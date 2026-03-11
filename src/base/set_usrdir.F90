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

!**s/p set_usrdir - initialization of common block

      integer function set_usrdir (F_argc,F_argv_S,F_cmdtyp_S,F_v1,F_v2)
      use lun
      use outusrdir
      implicit none
#include <arch_specific.hf>

      integer F_argc,F_v1,F_v2
      character(len=*) F_argv_S(0:F_argc),F_cmdtyp_S
      

!object
!	initialization of the module outusrdir. This function is
!       called when the keyword "usrdir" is found in the first word
!       of the directives in the input file given in the statement
!       "process_f_callback". This feature is enabled by the
!       ARMNLIB "rpn_fortran_callback" routine (called in "srequet")
!       which allows a different way of passing user directives than
!       the conventional FORTRAN namelist. This function will process
!       the following example command read from the named input file.
!
!   ie: usrdir=1,usr1;
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
! F_cmdtyp_S   I    - character command type - not used
! F_v1         I    - integer parameter 1 - not used
! F_v2         I    - integer parameter 2 - not used
!----------------------------------------------------------------
!
!Notes:
!
! examples:
! usrdir=1,usr1;
! usrdir=2,usr2;
!
! general syntax
! usrdir=usrdirid,usr[1-9];
!
!      usrdirid  - number to identify usrdirset to relate to sortie statement
!
      integer i,j, usrdirset, len_str
      character(len=OutUsrdir_name_lenght) :: temp
      character(len=1), dimension(9) :: OK_list_S
!
!-------------------------------------------------------------------
!
      if (Lun_out > 0) then
          write(Lun_out,*)
          write(Lun_out,*) F_argv_S(0),'=',F_argv_S(1),',',F_argv_S(2),',',(F_argv_S(i),i=3,F_argc)
      end if

      set_usrdir = 0

      OK_list_S = ['1','2','3','4','5','6','7','8','9']

      print*,'F_argv_S(1)=',F_argv_S(1)
      read(F_argv_S(1),*) usrdirset
      
      OutUsrdir_sets = OutUsrdir_sets + 1
      if (OutUsrdir_sets > OUTUSRDIR_MAX) then
          if (Lun_out > 0) then
             write(Lun_out,*)'SET_USRDIR WARNING: Too many usrdir definitions, maximum is ',OUTUSRDIR_MAX
          end if
          OutUsrdir_sets = OutUsrdir_sets - 1
          set_usrdir = 1
          return
      end if

      OutUsrdir_id(OutUsrdir_sets) = usrdirset
      
      ! Supprimer les guillemets (' et ou ") de F_argv_S(2)
      OutUsrdir_name_S(OutUsrdir_sets) = ""
      len_str = len_trim(F_argv_S(2))
      j = 1
      do i = 1, len_str
         if (F_argv_S(2)(i:i) /= '"' .and. F_argv_S(2)(i:i) /= "'") then
            OutUsrdir_name_S(OutUsrdir_sets)(j:j) = F_argv_S(2)(i:i)
            j = j + 1
         end if
      end do

      if(OutUsrdir_name_S(OutUsrdir_sets)(1:3) /= 'usr' .or. &
           .not. any(OutUsrdir_name_S(OutUsrdir_sets)(4:4) == OK_list_S)  .or. &
           len_trim(OutUsrdir_name_S(OutUsrdir_sets)) /= 4 )then
         write(Lun_out,*)'SET_USRDIR WARNING: value for usrdir must be usr[1-9], got ',&
              trim(OutUsrdir_name_S(OutUsrdir_sets))
         OutUsrdir_name_S(OutUsrdir_sets)=''
         OutUsrdir_sets = OutUsrdir_sets - 1
         set_usrdir = 1
      endif
      
!-------------------------------------------------------------------
      return
      end
