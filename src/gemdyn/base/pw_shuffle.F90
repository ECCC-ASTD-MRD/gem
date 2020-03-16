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

!**s/r pw_shuffle - Rename time level :P -> :M
!
      subroutine pw_shuffle
      use gmm_itf_mod
      use gmm_pw
      implicit none
#include <arch_specific.hf>

      integer gmmstat
      character(len=8) , dimension(2), parameter :: pw_uulist = ['PW_UU:P','PW_UU:M']
      character(len=8) , dimension(2), parameter :: pw_vvlist = ['PW_VV:P','PW_VV:M']
      character(len=8) , dimension(2), parameter :: pw_ttlist = ['PW_TT:P','PW_TT:M']
      character(len=8) , dimension(2), parameter :: pw_pmlist = ['PW_PM:P','PW_PM:M']
      character(len=8) , dimension(2), parameter :: pw_ptlist = ['PW_PT:P','PW_PT:M']
      character(len=8) , dimension(2), parameter :: pw_gzlist = ['PW_GZ:P','PW_GZ:M']
      character(len=8) , dimension(2), parameter :: pw_melist = ['PW_ME:P','PW_ME:M']
      character(len=8) , dimension(2), parameter :: pw_p0list = ['PW_P0:P','PW_P0:M']
!     ________________________________________________________________
!
      gmmstat = gmm_shuffle (pw_uulist)
      gmmstat = gmm_shuffle (pw_vvlist)
      gmmstat = gmm_shuffle (pw_ttlist)
      gmmstat = gmm_shuffle (pw_pmlist)
      gmmstat = gmm_shuffle (pw_ptlist)
      gmmstat = gmm_shuffle (pw_gzlist)
      gmmstat = gmm_shuffle (pw_melist)
      gmmstat = gmm_shuffle (pw_p0list)

      gmmstat = gmm_get(gmmk_pw_uu_plus_s  ,pw_uu_plus)
      gmmstat = gmm_get(gmmk_pw_vv_plus_s  ,pw_vv_plus)
      gmmstat = gmm_get(gmmk_pw_tt_plus_s  ,pw_tt_plus)
      gmmstat = gmm_get(gmmk_pw_pm_plus_s  ,pw_pm_plus)
      gmmstat = gmm_get(gmmk_pw_pt_plus_s  ,pw_pt_plus)
      gmmstat = gmm_get(gmmk_pw_gz_plus_s  ,pw_gz_plus)
      gmmstat = gmm_get(gmmk_pw_me_plus_s  ,pw_me_plus)
      gmmstat = gmm_get(gmmk_pw_p0_plus_s  ,pw_p0_plus)

      gmmstat = gmm_get(gmmk_pw_uu_moins_s  ,pw_uu_moins)
      gmmstat = gmm_get(gmmk_pw_vv_moins_s  ,pw_vv_moins)
      gmmstat = gmm_get(gmmk_pw_tt_moins_s  ,pw_tt_moins)
      gmmstat = gmm_get(gmmk_pw_pt_moins_s  ,pw_pt_moins)
      gmmstat = gmm_get(gmmk_pw_gz_moins_s  ,pw_gz_moins)
      gmmstat = gmm_get(gmmk_pw_pm_moins_s  ,pw_pm_moins)
      gmmstat = gmm_get(gmmk_pw_me_moins_s  ,pw_me_moins)
      gmmstat = gmm_get(gmmk_pw_p0_moins_s  ,pw_p0_moins)

!     ________________________________________________________________
!
      return
      end
