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

!**s/r pe_zero_topo - Initialize processor topology

      subroutine pe_zero_topo (F_npx, F_npy )
      use clib_itf_mod
      use path
      use cpus_options
      use version
      implicit none
#include <arch_specific.hf>

      integer F_npx, F_npy

      include 'gemdyn_version.inc'
#include <rmnlib_basics.hf>

      integer,external :: cpus_nml

      character(len=50) :: DSTP,name_S,arch_S,compil_S,user_S
      character(len=4056) :: fn
      logical :: is_official_L
      integer :: err, unf
!
!-------------------------------------------------------------------
!
      err= clib_mkdir (Path_output_S)

      call  open_status_file3 (trim(Path_output_S)//'/status_MOD.dot')
      call write_status_file3 ('_status=ABORT' )

      call atm_model_getversion2(name_S,Version_number_S,DSTP,arch_S, &
           compil_S,user_S,is_official_L)
      if (is_official_L) then
         Version_title_S = trim(name_S)//' --- Release of: '//trim(DSTP)
      else
         Version_title_S = trim(name_S)//' --- '//trim(user_S)//' Build: '//trim(DSTP)
      end if

      err = exdb(trim(Version_title_S),trim(Version_number_S),'NON')

!
! Read namelist &cpus from file model_settings
!
      unf= 0
      fn = trim(Path_work_S)//'/model_settings.nml'
      if (fnom (unf,trim(fn), 'SEQ+OLD', 0) == 0) then
         write (6, 6000) trim( fn )
         if (cpus_nml (unf) == 1 ) then
            F_npx = Cpus_npex
            F_npy = Cpus_npey
            err = cpus_nml (-1)
         else
            F_npx = 0
            F_npy = 0
            write (6, 8000)
         end if
         err= fclos(unf)
      else
         write (6, 6001) trim( fn)
      end if

      write (6,1001) trim(GEMDYN_NAME_S),trim(GEMDYN_VERSION_S), &
                     trim(GEMDYN_DSTP_S),trim(GEMDYN_EC_ARCH_S)

 1001 format (/3x,60('*')/3x,'Package: ',a,5x,'version: ',a/ &
               3x,'Release of: ',a,5x,'COMPILER: 'a/3x,60('*')/)
 6000 format (' READING &cpus namelists from FILE: '/A)
 6001 format (/,' Namelist FILE: ',A,' NOT AVAILABLE'/)
 8000 format (/,'========= ERROR IN S/R PE_ZERO_TOPO ============='/)
!
!-------------------------------------------------------------------
!
      return
      end
