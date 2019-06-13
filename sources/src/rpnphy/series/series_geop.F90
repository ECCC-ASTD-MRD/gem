!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

module series_geop_mod
   use tdpack_const, only: PI
   use phygridmap, only: phydim_ni, phydim_nj
   use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
   use series_options
   use series_xst_mod, only: series_xst_geo
   implicit none
   private

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

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
      integer :: j, istat
      real :: prcon, w1(phydim_ni)
      type(phymetaplus) :: metaplus
      !---------------------------------------------------------------

#define BUSSLICE(META,JJ) META%vptr(1:phydim_ni,JJ)

      call msg(MSG_INFO, PKGNAME_S//'Extracting Geop fields')

      prcon = 180./pi

      w1 = 1.
      do j= 1, phydim_nj
         call series_xst_geo(w1, 'MA', j)
      end do

      istat = phygetmetaplus(metaplus, 'DLAT', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j) * prcon
         call series_xst_geo(w1, 'LA', j)
      end do

      istat = phygetmetaplus(metaplus, 'DLON', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j) * prcon
         where(w1 < 0.) w1 = w1 + 360.
         call series_xst_geo(w1, 'LO', j)
      end do

      !#TODO: loop over var
      istat = phygetmetaplus(metaplus, 'Z0', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'ZP', j)
      end do

      istat = phygetmetaplus(metaplus, 'MG', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'MG', j)
      end do

      istat = phygetmetaplus(metaplus, 'LHTG', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'LH', j)
      end do

      istat = phygetmetaplus(metaplus, 'ALVIS', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'AL', j)
      end do

      istat = phygetmetaplus(metaplus, 'SNODP', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'SD', j)
      end do

      istat = phygetmetaplus(metaplus, 'TWATER', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'TM', j)
      end do

      istat = phygetmetaplus(metaplus, 'TSOIL', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'TP', j)
      end do

      istat = phygetmetaplus(metaplus, 'GLSEA', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'GL', j)
      end do

      istat = phygetmetaplus(metaplus, 'WSOIL', 'V', 'P', F_quiet=.true., &
           F_shortmatch=.false.)
      do j= 1, phydim_nj
         w1(1:phydim_ni) = BUSSLICE(metaplus,j)
         call series_xst_geo(w1, 'HS', j)
      end do
      !---------------------------------------------------------------
      return
   end subroutine series_geop

end module series_geop_mod
