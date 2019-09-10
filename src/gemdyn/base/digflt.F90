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

!**s/r digflt -  Compute digitally filtered fields
!
     subroutine digflt
      use step_options
      use gmm_vt1
      use gmm_vta
      use init_options
      use glb_ld
      use tr3d
      use gmm_itf_mod
      implicit none

#include <arch_specific.hf>

      integer i, j, k, n, istat
      real dfcoef
      real, pointer, dimension(:,:,:) :: tr,tra
!     __________________________________________________________________
!
      dfcoef = Init_dfco ( abs( (Init_halfspan - Step_kount ) ) )

      istat = gmm_get(gmmk_uta_s,uta)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(uta)'
      istat = gmm_get(gmmk_ut1_s,ut1)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(ut1)'
      istat = gmm_get(gmmk_vta_s,vta)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(vta)'
      istat = gmm_get(gmmk_vt1_s,vt1)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(vt1)'
      istat = gmm_get(gmmk_tta_s,tta)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(tta)'
      istat = gmm_get(gmmk_tt1_s,tt1)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(tt1)'
      istat = gmm_get(gmmk_zdta_s,zdta)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(zdta)'
      istat = gmm_get(gmmk_zdt1_s,zdt1)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(zdt1)'
      istat = gmm_get(gmmk_sta_s,sta)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(sta)'
      istat = gmm_get(gmmk_st1_s,st1)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(st1)'
      istat = gmm_get(gmmk_wta_s,wta)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(wta)'
      istat = gmm_get(gmmk_wt1_s,wt1)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(wt1)'
      istat = gmm_get(gmmk_qt1_s,qt1)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(qt1)'
      istat = gmm_get(gmmk_qta_s,qta)
      if (GMM_IS_ERROR(istat)) print *,'digflt ERROR at gmm_get(qta)'

      do k= 1, l_nk
      do j= 1, l_nj
      do i= 1, l_ni
         tta (i,j,k)   =  tta(i,j,k)   + dfcoef *  tt1(i,j,k)
         uta (i,j,k)   =  uta(i,j,k)   + dfcoef *  ut1(i,j,k)
         vta (i,j,k)   =  vta(i,j,k)   + dfcoef *  vt1(i,j,k)
         zdta(i,j,k)   = zdta(i,j,k)   + dfcoef * zdt1(i,j,k)
         wta (i,j,k)   =  wta(i,j,k)   + dfcoef *  wt1(i,j,k)
         qta (i,j,k)   =  qta(i,j,k)   + dfcoef *  qt1(i,j,k)
      end do
      end do
      end do

      do j= 1, l_nj
      do i= 1, l_ni
         sta (i,j)        = sta(i,j)        + dfcoef * st1(i,j)
         qta (i,j,l_nk+1) = qta(i,j,l_nk+1) + dfcoef * qt1(i,j,l_nk+1)
      end do
      end do

!***************************************************************
!     Passive tracers (no passive tracers in linear model)
!***************************************************************

      do n=1,Tr3d_ntr

         istat = gmm_get('DIGF_'//trim(Tr3d_name_S(n))      , tra)
         istat = gmm_get('TR/'  //trim(Tr3d_name_S(n))//':P', tr )

         if ( Init_dftr_L ) then
            do k=1,G_nk
            do j=1,l_nj
            do i=1,l_ni
               tra(i,j,k) = tra(i,j,k) + dfcoef * tr(i,j,k)
            end do
            end do
            end do
         elseif ( Step_kount == Init_halfspan ) then
            tra(1:l_ni,1:l_nj,1:G_nk) = tr(1:l_ni,1:l_nj,1:G_nk)
         end if

      end do
!     __________________________________________________________________
!
      return
      end

