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

module cpl_nml_mod

   private
   public :: cpl_nml

contains

      integer function cpl_nml (F_namelistf_S, F_unout)
      use cpl_mod
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_namelistf_S
      integer, intent(in)          :: F_unout

!authors    Francois Roy -- spring 2015
! 
!revision
! v4_80 - Roy, F.  - initial version
!
!Purpose
!  Default configuration and reading namelist coupling

#include <rmn/WhiteBoard.hf>

      namelist /coupling/ cpl_ocn_L, cpl_wav_L, cplocn_ocnf_nsprd, cplocn_debug_L, cplocn_iweight_L

      integer, external :: fnom,wkoffit
      character*60 name_list
      logical print_L,found_namelist
      integer nrec,unf,err,err_open
!
!-------------------------------------------------------------------
!
      cpl_nml = -1
      print_L = (F_unout>0)

      if ((F_namelistf_S.eq.'print').or.(F_namelistf_S.eq.'PRINT')) then
         cpl_nml = 1
         if (print_L) then
            write(F_unout,2000)
            write(F_unout,NML=coupling)
         endif
         return
      endif
!
! Defaults values for coupling namelist variables
!
      cpl_ocn_L = .false.
      cpl_wav_L = .false.
      cplocn_ocnf_nsprd = 0
      cplocn_debug_L = .false.
      cplocn_iweight_L = .false.

      unf            = 0
      found_namelist = .false.
      err            = wkoffit (F_namelistf_S)

      if (err.ge.-1) then
         
         err_open= fnom (unf, F_namelistf_S, 'SEQ+OLD' , nrec)
         if (err_open.eq.0) then
            
            name_list = 'coupling'
            read (unf, nml=coupling, end= 333, err= 90)
            
            found_namelist = .true.
            goto 333
90          if (print_L) then
               write (F_unout, 1500) trim(name_list),trim(F_namelistf_S)
               return
            endif
            
         endif
333      call fclos (unf)
         
      endif
      if ( (err.lt.-1) .or. (err_open.ne.0) ) then
         if (print_L) write (F_unout, 1600) trim(F_namelistf_S)
         cpl_nml= 0
         return
      endif

      if (.not.found_namelist) then
         if (print_L) write (F_unout,1200)
         cpl_nml= 0
         return
      endif

      cpl_nml= 0
      if (cpl_ocn_L) cpl_nml= 1

!      cplwav = cpl_wav_L

1200 format (/3X,61('*')/3x, &
             'NAMELIST &coupling NOT available; running without coupling'&
             /3x,61('*'))
1500 format (/,' NAMELIST ',a,' INVALID IN FILE: ',a/)
1600 format (/,' NAMELIST FILE ',a,' NOT AVAILABLE: RUNNING WITHOUT COUPLING')
 2000 FORMAT (/4x,'COUPLING NAMELIST :',/,4x,21('='))
!
!-------------------------------------------------------------------
!
      return
      end function cpl_nml

end module cpl_nml_mod
