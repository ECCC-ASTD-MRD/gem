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

!**s/r slabsym - symmetrical boundary conditions for theoretical
!               cases
!
      subroutine slabsym ()
      use gmm_vt1
      use gem_options
      use glb_ld
      use lun
      use tr3d
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>

!author
!     Gravel              - spring 2003 (after MC2 v_4.9.3)
!
!revision
! v3_11 - Gravel             - initial version
! v3_30 - Lee V              - changed t0=> t1,nest is after t02t1
! v4_05 - Lepine M.          - VMM replacement with GMM


      type(gmm_metadata) :: mymeta
      integer err,i,j,k,n
      integer jin, jj
      real, pointer, dimension(:,:,:) :: tr
!
!----------------------------------------------------------------------
!
      err = gmm_get(gmmk_ut1_s,ut1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(ut1)'
      err = gmm_get(gmmk_vt1_s,vt1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(vt1)'
      err = gmm_get(gmmk_wt1_s,wt1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(wt1)'
      err = gmm_get(gmmk_tt1_s,tt1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(tt1)'
      err = gmm_get(gmmk_zdt1_s,zdt1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(zdt1)'
      err = gmm_get(gmmk_st1_s,st1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(st1)'
      err = gmm_get(gmmk_wt1_s,wt1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(wt1)'
      err = gmm_get(gmmk_qt1_s,qt1,mymeta)
      if (GMM_IS_ERROR(err)) print *,'slabsym ERROR at gmm_get(qt1)'
!
      if (l_north) then
         do k=1,G_nk
            jin = l_nj-pil_n-1
            jj  = l_nj-pil_n
            do i=1,l_ni
               vt1  (i,jj,k) = vt1  (i,jin,k)
            end do
            jin = l_nj-pil_n
            do j=1,pil_n
            jj  = l_nj-pil_n+j
            do i=1,l_ni
               vt1 (i,jj,k) = vt1 (i,jin,k)
               tt1 (i,jj,k) = tt1 (i,jin,k)
               wt1 (i,jj,k) = wt1 (i,jin,k)
               zdt1(i,jj,k) = zdt1(i,jin,k)
               qt1 (i,jj,k) = qt1 (i,jin,k)
            end do
            do i=1,l_niu
               ut1 (i,jj,k) = ut1 (i,jin,k)
            end do
            end do
         end do
         jin = l_nj-pil_n
         do j=1,pil_n
         jj  = l_nj-pil_n+j
         do i=1,l_ni
            st1(i,jj)        = st1(i,jin)
            qt1(i,jj,G_nk+1) = qt1(i,jin,G_nk+1)
         end do
         end do
      end if
!
      if (l_south) then
         do k=1,G_nk
            jin = pil_s+1
            do j=1,pil_s
            jj  = pil_s-j+1
            do i=1,l_ni
               vt1 (i,jj,k) = vt1 (i,jin,k)
               wt1 (i,jj,k) = wt1 (i,jin,k)
               tt1 (i,jj,k) = tt1 (i,jin,k)
               zdt1(i,jj,k) = zdt1(i,jin,k)
               qt1 (i,jj,k) = qt1 (i,jin,k)
            end do
            do i=1,l_niu
               ut1 (i,jj,k) = ut1 (i,jin,k)
            end do
            end do
         end do
         jin = pil_s+1
         do j=1,pil_s
         jj  = pil_s-j+1
         do i=1,l_ni
            st1(i,jj)        = st1(i,jin)
            qt1(i,jj,G_nk+1) = qt1(i,jin,G_nk+1)
         end do
         end do
      end if
!
      do n=1,Tr3d_ntr
         nullify(tr)
         err = gmm_get('TR/'//trim(Tr3d_name_S(n))//':P',tr,mymeta)
         if (err == 0) then
         if (l_north)      then
            do k=1,G_nk
               jin = l_nj-pil_n
               do j=1,pil_n
                  jj  = l_nj-pil_n+j
                  do i=1,l_ni
                     tr(i,jj,k) = tr(i,jin,k)
                  end do
               end do
            end do
         end if
         if (l_south) then
            do k=1,G_nk
               jin = pil_s+1
               do j=1,pil_s
                  jj  = pil_s-j+1
                  do i=1,l_ni
                     tr(i,jj,k) = tr(i,jin,k)
                  end do
               end do
            end do
         end if
         end if
      end do
!
!----------------------------------------------------------------------
      return
      end
