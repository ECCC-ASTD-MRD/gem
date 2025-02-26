
module series_geop_mod
   use tdpack_const, only: PI
   use phygridmap, only: phydim_ni, phydim_nj
   use phymem, only: phymem_find, phymem_getdata
   use series_options
   use series_xst_mod, only: series_xst_geo
   implicit none
   private

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   public :: series_geop

contains

   subroutine series_geop()
!!!#include <arch_specific.hf>
      !@object Prepares "first record" output for time series.
      !@author Andre Methot - cmc - june 1994 v0_14
      !@revision
      ! v2_00 - Desgagne M.     - initial MPI version
      ! v2_20 - Lee V.          - extract geophysical fields for time-series 
      ! v2_20                     from physics permanent bus,not VMM variables
      ! v3_11 - A. Plante       - Adjust code for LAM time-series
      ! v3_20 - Winger K.       - correct time series handling in climate mode
      ! v3_30 - Winger K.       - Change serset to serset8 for HEURE
      ! v3_30 - Desgagne M.     - Remove Mem_phyncore_L
      ! 2017-01, S.Chamberland  - Major revision
      !@description
      !      This subroutine is part of time serie's package
      !      initialisation. It extracts and produce output of constant
      !      fields to be used by the unwrapper.
      integer :: j, istat, idxv1(1),k1
      real, pointer, contiguous :: ptr1d(:), ptr2d(:,:)
      real :: prcon, w1(phydim_ni)

      include "surface.cdk"
      !---------------------------------------------------------------
      
      call msg(MSG_INFO, PKGNAME_S//'Extracting Geop fields')

      prcon = 180./pi

      w1 = 1.
      do j= 1, phydim_nj
         call series_xst_geo(w1, 'MA', j)
      end do

      istat = phymem_find(idxv1, 'DLAT', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr1d)
            istat = phymem_getdata(ptr1d, idxv1(1), j)
            if (istat < 0) cycle
            w1(1:phydim_ni) = ptr1d(:) * prcon
            call series_xst_geo(w1, 'LA', j)
         end do
      endif

      istat = phymem_find(idxv1, 'DLON', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr1d)
            istat = phymem_getdata(ptr1d, idxv1(1), j)
            if (istat < 0) cycle
            w1(1:phydim_ni) = ptr1d(:) * prcon
            where(w1 < 0.) w1 = w1 + 360.
            call series_xst_geo(w1, 'LO', j)
         end do
      endif

      !#TODO: loop over var
      istat = phymem_find(idxv1, 'Z0', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr2d)
            istat = phymem_getdata(ptr2d, idxv1(1), j)
            if (istat < 0) cycle
            k1 = min(indx_agrege, size(ptr2d,2))
            w1(1:phydim_ni) = ptr2d(:,k1)
            call series_xst_geo(w1, 'ZP', j)
         end do
      endif

      istat = phymem_find(idxv1, 'MG', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr1d)
            istat = phymem_getdata(ptr1d, idxv1(1), j)
            if (istat < 0) cycle
            w1(1:phydim_ni) = ptr1d(:)
            call series_xst_geo(w1, 'MG', j)
         end do
      endif

      istat = phymem_find(idxv1, 'LHTG', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr1d)
            istat = phymem_getdata(ptr1d, idxv1(1), j)
            if (istat < 0) cycle
            w1(1:phydim_ni) = ptr1d(:)
            call series_xst_geo(w1, 'LH', j)
         end do
      endif

      istat = phymem_find(idxv1, 'ALVIS', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr2d)
            istat = phymem_getdata(ptr2d, idxv1(1), j)
            if (istat < 0) cycle
            k1 = min(indx_agrege, size(ptr2d,2))
            w1(1:phydim_ni) = ptr2d(:,k1)
            call series_xst_geo(w1, 'AL', j)
         end do
      endif

      istat = phymem_find(idxv1, 'SNODP', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr2d)
            istat = phymem_getdata(ptr2d, idxv1(1), j)
            if (istat < 0) cycle
            k1 = min(indx_agrege, size(ptr2d,2))
            w1(1:phydim_ni) = ptr2d(:,k1)
            call series_xst_geo(w1, 'SD', j)
         end do
      endif

      istat = phymem_find(idxv1, 'TWATER', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr2d)
            istat = phymem_getdata(ptr2d, idxv1(1), j)
            if (istat < 0) cycle
            k1 = min(indx_agrege, size(ptr2d,2))
            w1(1:phydim_ni) = ptr2d(:,k1)
            call series_xst_geo(w1, 'TM', j)
         end do
      endif

      istat = phymem_find(idxv1, 'TSOIL', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr2d)
            istat = phymem_getdata(ptr2d, idxv1(1), j)
            if (istat < 0) cycle
            k1 = min(indx_agrege, size(ptr2d,2))
            w1(1:phydim_ni) = ptr2d(:,k1)
            call series_xst_geo(w1, 'TP', j)
         end do
      endif

      istat = phymem_find(idxv1, 'GLSEA', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr2d)
            istat = phymem_getdata(ptr2d, idxv1(1), j)
            if (istat < 0) cycle
            k1 = min(indx_agrege, size(ptr2d,2))
            w1(1:phydim_ni) = ptr2d(:,k1)
            call series_xst_geo(w1, 'GL', j)
         end do
      endif

      istat = phymem_find(idxv1, 'WSOIL', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      if (istat > 0) then
         do j= 1, phydim_nj
            nullify(ptr2d)
            istat = phymem_getdata(ptr2d, idxv1(1), j)
            if (istat < 0) cycle
            k1 = min(indx_agrege, size(ptr2d,2))
            w1(1:phydim_ni) = ptr2d(:,k1)
            call series_xst_geo(w1, 'HS', j)
         end do
      endif
      !---------------------------------------------------------------
      return
   end subroutine series_geop

end module series_geop_mod
