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

!**s/r pw_init - Initialization of pw_*_moins variables

      subroutine pw_init
      use gmm_pw
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>


      integer istat

!     ________________________________________________________________
!
      istat = gmm_get ( gmmk_pw_uu_plus_s  , pw_uu_plus  )
      istat = gmm_get ( gmmk_pw_vv_plus_s  , pw_vv_plus  )
      istat = gmm_get ( gmmk_pw_tt_plus_s  , pw_tt_plus  )
      istat = gmm_get ( gmmk_pw_gz_plus_s  , pw_gz_plus  )
      istat = gmm_get ( gmmk_pw_pm_plus_s  , pw_pm_plus  )
      istat = gmm_get ( gmmk_pw_pt_plus_s  , pw_pt_plus  )
      istat = gmm_get ( gmmk_pw_me_plus_s  , pw_me_plus  )
      istat = gmm_get ( gmmk_pw_p0_plus_s  , pw_p0_plus  )
      istat = gmm_get ( gmmk_pw_uu_moins_s , pw_uu_moins )
      istat = gmm_get ( gmmk_pw_vv_moins_s , pw_vv_moins )
      istat = gmm_get ( gmmk_pw_tt_moins_s , pw_tt_moins )
      istat = gmm_get ( gmmk_pw_gz_moins_s , pw_gz_moins )
      istat = gmm_get ( gmmk_pw_pm_moins_s , pw_pm_moins )
      istat = gmm_get ( gmmk_pw_pt_moins_s , pw_pt_moins )
      istat = gmm_get ( gmmk_pw_me_moins_s , pw_me_moins )
      istat = gmm_get ( gmmk_pw_p0_moins_s , pw_p0_moins )

      pw_uu_moins = pw_uu_plus
      pw_vv_moins = pw_vv_plus
      pw_tt_moins = pw_tt_plus
      pw_gz_moins = pw_gz_plus
      pw_pm_moins = pw_pm_plus
      pw_pt_moins = pw_pt_plus
      pw_me_moins = pw_me_plus
      pw_p0_moins = pw_p0_plus
!     ________________________________________________________________
!
      return
      end
