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

!**function samevert -

      integer function samevert ( zam_8,zbm_8,lvm, zat_8,zbt_8,lvt, &
                                  analtopo, F_topo, ni,nj,i0,in,j0,jn)
      use glb_ld
      use lun
      use ver
      use type_mod
      implicit none
#include <arch_specific.hf>

      integer ni,nj,i0,in,j0,jn,lvm,lvt
      real analtopo(ni,nj),F_topo (ni,nj)
      real*8 zam_8(lvm),zbm_8(lvm),zat_8(lvt),zbt_8(lvt)

!author
! v3_30 - Desgagne M        - from samelevel
!
!revision
! v4_60 - Lee V             - range to be checked and return an integer
!
!object
!     Compare analysis levels to model levels and topographies between
!     analysis and the model. Returns 3 types of codes:
!     samevert > 0 - means Vertical interpolation is needed.
!     samevert = 0 - means NO vertical interpolation is needed.
!     samevert < 0 - error in this function


      integer :: i,j,k,err,ndiff,sumdiff
      real, parameter :: DIFSIG=1.e-5
!
! ---------------------------------------------------------------------
!
      samevert = 0
      if (Lun_out > 0) write(Lun_out,*) 'SAMEVERT CHECK'
      if (in > ni.or.jn > nj.or.i0 < 1.or.j0 < 1) then
          if (Lun_out > 0) then
             write(lun_out,*) ' ******* ERROR   ********'
             write(Lun_out,*) 'RANGE TO VERIFY IS GREATER THAN DIMENSIONS'
             write(lun_out,*) ' i0,in,j0,jn=',i0,in,j0,jn,' in (',ni,',',nj,')'
          end if
          samevert=-1
          return
      end if

      if (G_nk /= lvm) then
          samevert = 1
          return
      end if

      do k=1,lvm
         if (abs(zam_8(k)-Ver_a_8%m(k)) > DIFSIG .or. &
             abs(zbm_8(k)-Ver_b_8%m(k)) > DIFSIG) samevert=samevert+1
      end do


      do k=1,lvt
         if (abs(zat_8(k)-Ver_a_8%t(k)) > DIFSIG .or. &
             abs(zbt_8(k)-Ver_b_8%t(k)) > DIFSIG) samevert=samevert+1
      end do
      if (samevert > 0) then
          return
      end if

!     Check if analysis and model have the same topography

      ndiff   = 0
      sumdiff = 0
      do j= j0, jn
      do i= i0, in
         if (abs(F_topo(i,j)-analtopo(i,j)) > 25.) then
            ndiff = ndiff+1
         end if
      end do
      end do

      call rpn_comm_allreduce (ndiff,sumdiff,1,"MPI_INTEGER", &
                                            "MPI_SUM","grid",err)
      samevert = sumdiff
!
! ---------------------------------------------------------------------
!
      return
      end

