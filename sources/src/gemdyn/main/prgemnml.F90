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
!TODO: opening unit 6 may be unsafe, better to pass unit to printing functions
#define FILEID 60

subroutine prgemnml()
   use lun
   implicit none
#include <arch_specific.hf>
   print *,'creating gemdict.nml'
   open(FILEID, file="gemdict.nml", access='SEQUENTIAL', form='FORMATTED')
   Lun_out=FILEID
   call nml_gem()
   call nml_theo()
   close(FILEID)
   return
end subroutine prgemnml

subroutine nml_gem()
   use HORgrid_options
   use VERgrid_options
   use dynkernel_options
   use dyn_expo_options
   use gem_options
   use adv_options
   use grdc_options
   use hvdif_options
   use init_options
   use inp_options
   use lam_options
   use out_options
   use spn_options
   use step_options
   implicit none
   integer err

   err = dynKernel_nml(-2)
   err = dyn_fisl_nml (-2)
   err = adv_nml      (-2)
   err = hvdif_nml    (-2)
   err = step_nml     (-2)
   err = gem_nml      (-2)
   err = HORgrid_nml  (-2)
   err = grdc_nml     (-2)
   err = VERgrid_nml  (-2)
   err = init_nml     (-2)
   err = inp_nml      (-2)
   err = lam_nml      (-2)
   err = out_nml      (-2)
   err = spn_nml      (-2)
   return
end subroutine nml_gem

subroutine nml_theo()
   use gem_options
   use HORgrid_options
   use step_options
   use theo_options
   use bubble_options
   use mtn_options
   implicit none
   write(FILEID, nml=theo_cfgs)
   write(FILEID, nml=bubble_cfgs)
   write(FILEID, nml=mtn_cfgs)
   return
end subroutine nml_theo
