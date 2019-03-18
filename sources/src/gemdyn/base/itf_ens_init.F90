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

!**s/r itf_ens_init - initialize ensemble prevision system

      subroutine itf_ens_init
      use ens_options
      use glb_ld
      use lun
      use path
      implicit none
#include <arch_specific.hf>

#include <rmnlib_basics.hf>

      integer err_nml, ier, unf
!
!-------------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1000)
!
! Read namelist &ensembles from file model_settings
!
      unf= 0 ; err_nml= -1
      if (fnom (unf,Path_nml_S, 'SEQ+OLD', 0) == 0) then
         if (Lun_out > 0) write (6, 6000) trim( Path_nml_S )
         err_nml= ens_nml(unf)
         ier= fclos(unf)
      else
         if (Lun_out > 0) write (6, 6001) trim( Path_nml_S )
      end if

      if (err_nml==1) then
         if (Lun_out > 0) write(Lun_out,1010)
         return
      end if

      call handle_error (err_nml,'itf_ens_init','Problem with ens_nml')

      call ens_setmem (l_ni, l_nj, l_nk, Lun_out)

 1000 format(/,'INITIALIZATION OF ENSEMBLES (S/R ITF_ENS_INIT)'/(46('=')))
 1010 format(/,'NO ENSEMBLES REQUIRED       (S/R ITF_ENS_INIT)'/(46('=')))
 6000 format (' READING &ensembles namelists from FILE: '/A)
 6001 format (/,' Namelist FILE: ',A,' NOT AVAILABLE'/)
!
!-------------------------------------------------------------------
!
      return
      end
