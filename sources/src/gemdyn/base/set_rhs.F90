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

!**s/r set_rhs - initialization of the VMM common for right-hand side of
!               equations
!
#define SPY_VMM_CREATE spy_vmm_create

!
      subroutine set_rhs
      use gmm_rhsc
      use gmm_orh
      use glb_ld
      use lun
      use gmm_itf_mod
      use var_gmm
      implicit none
#include <arch_specific.hf>
!
!author
!     Gravel/Roch   - rpn - august 1993
!
!revision
! v2_00 - Desgagne/Lee   - initial MPI version (from setrhs v1_03)
! v2_21 - J. P. Toviessi - rename some model output variable
! v4_05 - Girard C.         - Open top
!
!object
!     See above id.
!
!arguments
!	none
!

!
      integer :: istat
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1000)
!
! Assign the names of the variables
!
      gmmk_rhsu_s = 'RHSU'
      gmmk_rhsv_s = 'RHSV'
      gmmk_rhst_s = 'RHST'
      gmmk_rhsc_s = 'RHSC'
      gmmk_rhsw_s = 'RHSW'
      gmmk_rhsf_s = 'RHSF'
      gmmk_rhsp_s = 'RHSP'
      gmmk_rhsb_s = 'RHSB'
!
      gmmk_ruw1_s = 'RUW1'
      gmmk_ruw2_s = 'RUW2'
      gmmk_rvw1_s = 'RVW1'
      gmmk_rvw2_s = 'RVW2'
!
      gmmk_orhsu_s = 'ORHU'
      gmmk_orhsv_s = 'ORHV'
      gmmk_orhst_s = 'ORHT'
      gmmk_orhsc_s = 'ORHC'
      gmmk_orhsw_s = 'ORHW'
      gmmk_orhsf_s = 'ORHF'
!
      nullify (rhsu,rhsv,rhst,rhsc,rhsw,rhsf,rhsp,rhsb)
      nullify (ruw1,rvw1,ruw2,rvw2)
      nullify (orhsu,orhsv,orhst,orhsc,orhsw,orhsf)
!
      istat = gmm_create(gmmk_rhsu_s,rhsu,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhsu)'
      istat = gmm_create(gmmk_rhsv_s,rhsv,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhsv)'
      istat = gmm_create(gmmk_rhst_s,rhst,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhst)'
      istat = gmm_create(gmmk_rhsc_s,rhsc,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhsc)'
      istat = gmm_create(gmmk_rhsw_s,rhsw,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhsw)'
      istat = gmm_create(gmmk_rhsf_s,rhsf,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhsf)'
      istat = gmm_create(gmmk_rhsp_s,rhsp,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhsp)'
      istat = gmm_create(gmmk_rhsb_s,rhsb,meta2d   , GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rhsb)'
!
      istat = gmm_create(gmmk_ruw1_s,ruw1,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(ruw1)'
      istat = gmm_create(gmmk_rvw1_s,rvw1,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rvw1)'
      istat = gmm_create(gmmk_ruw2_s,ruw2,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(ruw2)'
      istat = gmm_create(gmmk_rvw2_s,rvw2,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(rvw2)'
!
      istat = gmm_create(gmmk_orhsu_s,orhsu,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(orhsu)'
      istat = gmm_create(gmmk_orhsv_s,orhsv,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(orhsv)'
      istat = gmm_create(gmmk_orhst_s,orhst,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(orhst)'
      istat = gmm_create(gmmk_orhsc_s,orhsc,meta3d_nk, GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(orhsc)'
      istat = gmm_create(gmmk_orhsw_s,orhsw,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(orhsw)'
      istat = gmm_create(gmmk_orhsf_s,orhsf,meta3d_nk,GMM_FLAG_IZER)
      if (GMM_IS_ERROR(istat)) print *,'set_rhs ERROR at gmm_create(orhsf)'
!
      istat = gmm_get (gmmk_rhsw_s, rhsw)
      rhsw  = 0.

 1000 format( &
      /,'INITIALIZATION OF RIGHT-HAND SIDE COMDECK (S/R SET_RHS)', &
      /,'=======================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
