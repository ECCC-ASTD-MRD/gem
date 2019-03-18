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

!**s/r ta2t1tx -  Transfer variables ta into t1
!
      subroutine ta2t1tx
      use gmm_vt1
      use gmm_vta
      use gem_options
      use tr3d
      use gmm_itf_mod
      implicit none

!author
!     Alain Patoine - april 94
!
!revision
! v2_00 - Desgagne M.       - initial MPI version
! v2_30 - Edouard  S.       - remove pi' at the top
! v2_31 - Desgagne M.       - remove treatment of HU and QC and
!                             re-introduce tracers
! v3_00 - Desgagne & Lee    - Lam configuration
! v3_21 - Lee V.            - remove TR2D
! V4    - Girard-Plante     - Staggered version
! v4_05 - Lepine M.         - VMM replacement with GMM
! V4_11 - Plante A.         - Add diag var on vertical scope
!                             recompute fip

#include <arch_specific.hf>

      integer n,istat
      real, pointer, dimension(:,:,:) :: tr,tra
!
!     ---------------------------------------------------------------
!
      istat = gmm_get(gmmk_uta_s,uta)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(uta)'
      istat = gmm_get(gmmk_ut1_s,ut1)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(ut1)'
      istat = gmm_get(gmmk_vta_s,vta)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(vta)'
      istat = gmm_get(gmmk_vt1_s,vt1)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(vt1)'
      istat = gmm_get(gmmk_wta_s,wta)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(wta)'
      istat = gmm_get(gmmk_wt1_s,wt1)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(wt1)'
      istat = gmm_get(gmmk_tta_s,tta)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(tta)'
      istat = gmm_get(gmmk_tt1_s,tt1)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(tt1)'
      istat = gmm_get(gmmk_zdta_s,zdta)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(zdta)'
      istat = gmm_get(gmmk_zdt1_s,zdt1)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(zdt1)'
      istat = gmm_get(gmmk_sta_s,sta)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(sta)'
      istat = gmm_get(gmmk_st1_s,st1)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(st1)'
      istat = gmm_get(gmmk_qt1_s,qt1)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(qt1)'
      istat = gmm_get(gmmk_qta_s,qta)
      if (GMM_IS_ERROR(istat)) print *,'ta2t1tx ERROR at gmm_get(qt1)'

      ut1  = uta ; vt1  = vta ; wt1  = wta
      tt1  = tta ; zdt1 = zdta; qt1  = qta
      st1  = sta

      do n=1,Tr3d_ntr
         istat = gmm_get('DIGF_'//trim(Tr3d_name_S(n))      , tra)
         istat = gmm_get('TR/'  //trim(Tr3d_name_S(n))//':P', tr )
         tr = tra
      end do
!
!     ---------------------------------------------------------------
!
      return
      end
