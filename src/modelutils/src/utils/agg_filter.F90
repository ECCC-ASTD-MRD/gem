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

!/@
module agg_filter_mod
   implicit none
   private
   !@objective Aggregation filter
   !@author  Stephane Chamberland, 2012-02
   !@description
   ! Public functions
   public :: aggreg
   !@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   interface aggreg
      module procedure aggreg_r4_2d
      module procedure aggreg_r4_3d
   end interface

contains


   !/@
   subroutine aggreg_r4_2d(F_data,F_i0,F_j0,F_in,F_jn,F_di,F_dj)
      implicit none
      !@objective Aggregate by doing a simple average
      !@arguments
      real,pointer :: F_data(:,:)
      integer,intent(inout) :: F_i0,F_j0,F_in,F_jn
      integer,intent(in) :: F_di,F_dj
      !@author
      !@/
      integer :: l_ijk(2),u_ijk(2),ii,jj,di,dj,&
           iout,jout,ii0,jj0,ii1,jj1,di2,dj2,ni2,nj2
      ! ---------------------------------------------------------------------
      l_ijk = lbound(F_data)
      u_ijk = ubound(F_data)
      F_i0 = min(max(l_ijk(1),F_i0),u_ijk(1))
      F_j0 = min(max(l_ijk(2),F_j0),u_ijk(2))
      F_in = min(max(F_i0,F_in),u_ijk(1))
      F_jn = min(max(F_j0,F_jn),u_ijk(1))
      ni2 = F_in-F_i0+1
      nj2 = F_jn-F_j0+1
      di = min(max(1,F_di),ni2)
      dj = min(max(1,F_dj),nj2)
      if (all((/di,dj/) == 1)) return
      if (di >= ni2 .and. dj >= nj2) then
         F_data(F_i0,F_j0) = sum(F_data(F_i0:F_in,F_j0:F_jn))/float(ni2*nj2)
         return
      endif
      !TODO: ?should aggreg be done centered with +-dij or +-dij/2 points?
      jout = F_j0 - 1
      do jj = F_j0,F_jn,dj
         jout = jout + 1
         jj0 = min(jj,F_jn-dj+1)
         jj1 = min(jj0+dj-1,F_jn)
         dj2 = jj1-jj0+1
         iout = F_i0 - 1
         do ii = F_i0,F_in,di
            iout = iout + 1
            ii0 = min(ii,F_in-di+1)
            ii1 = min(ii0+di-1,F_in)
            di2 = ii1-ii0+1
            F_data(iout,jout) = sum(F_data(ii0:ii1,jj0:jj1))/float(di2*dj2)
         enddo
      enddo
      F_in = iout
      F_jn = jout
      ! ---------------------------------------------------------------------
      return
   end subroutine aggreg_r4_2d


   !/@
   subroutine aggreg_r4_3d(F_data,F_i0,F_j0,F_in,F_jn,F_di,F_dj)
      implicit none
      !@objective Aggregate by doing a simple average
      !@arguments
      real,pointer :: F_data(:,:,:)
      integer,intent(inout) :: F_i0,F_j0,F_in,F_jn
      integer,intent(in) :: F_di,F_dj
      !@author
      !@/
      integer :: l_ijk(3),u_ijk(3),ii,jj,kk,di,dj,&
           iout,jout,ii0,jj0,ii1,jj1,di2,dj2,ni2,nj2
      real :: rtmp
      ! ---------------------------------------------------------------------
      l_ijk = lbound(F_data)
      u_ijk = ubound(F_data)
      F_i0 = min(max(l_ijk(1),F_i0),u_ijk(1))
      F_j0 = min(max(l_ijk(2),F_j0),u_ijk(2))
      F_in = min(max(F_i0,F_in),u_ijk(1))
      F_jn = min(max(F_j0,F_jn),u_ijk(1))
      ni2 = F_in-F_i0+1
      nj2 = F_jn-F_j0+1
      di = min(max(1,F_di),ni2)
      dj = min(max(1,F_dj),nj2)
      if (all((/di,dj/) == 1)) return
      if (di >= ni2 .and. dj >= nj2) then
         do kk=l_ijk(3),u_ijk(3)
            F_data(F_i0,F_j0,kk) = sum(F_data(F_i0:F_in,F_j0:F_jn,kk))/float(ni2*nj2)
         enddo
         return
      endif
      jout = F_j0 - 1
      do jj = F_j0,F_jn,dj
         jout = jout + 1
         jj0 = min(jj,F_jn-dj+1)
         jj1 = min(jj0+dj-1,F_jn)
         dj2 = jj1-jj0+1
         iout = F_i0 - 1
         do ii = F_i0,F_in,di
            iout = iout + 1
            ii0 = min(ii,F_in-di+1)
            ii1 = min(ii0+di-1,F_in)
            di2 = ii1-ii0+1
            rtmp = 1./float(di2*dj2)
            do kk=l_ijk(3),u_ijk(3)
               F_data(iout,jout,kk) = sum(F_data(ii0:ii1,jj0:jj1,kk))*rtmp
            enddo
         enddo
      enddo
      F_in = iout
      F_jn = jout
      ! ---------------------------------------------------------------------
      return
   end subroutine aggreg_r4_3d

end module agg_filter_mod
