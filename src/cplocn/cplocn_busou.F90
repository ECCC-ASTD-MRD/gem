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

!**s/r cplocn_busou

      subroutine cplocn_busou
      use phy_itf, only: phy_get
      use cpl_mod
      use cplocn_mod
      use rmn_gmm
      implicit none
#include <arch_specific.hf>
#include <rmn/msg.h>

!
!authors    Michel Desgagne - Spring 2008
! 
!revision
! v3_31 - Desgagne M.          - initial MPI version
! v4_06 - Lepine M.            - VMM replacement with GMM
! v???? - Roy F. Belanger J-M  - coupling with nemo option
!                              - default C_coupling_L (MoGSL)
! v4_7  - Roy F. Desgagne M.   - rewritten from itf_cpl_fillatm
!                              - now executed in physics world

! Purpose: Pack information to send to ocean model

      integer istat
      real, dimension(:,:  ), pointer :: ptr2d
      real, dimension(:,:,:), pointer :: ptr3d
      real, dimension(cpl_drv_lni,cpl_drv_lnj,cplocn_n_fldou), target :: busou_tmp
      character(len=GMM_MAXNAMELENGTH) bus_name
      character(len=1)  bpath

      integer ivar,kstr,kend
!     ________________________________________________________________

      do ivar=1,cplocn_n_fldou

         bus_name = cplocn_cvou_N(ivar)
         bpath    = cplocn_cvou_G(ivar)
         kstr     = cplocn_cvou_K(1,ivar)
         kend     = cplocn_cvou_K(2,ivar)

         if ( cplocn_cvou_G(ivar) == 'P' ) then ! physical grid

            ptr3d => busou_tmp(cpl_drv_i0:cpl_drv_in, &
                               cpl_drv_j0:cpl_drv_jn, ivar:ivar)

 
            istat = phy_get(ptr3d, bus_name,                          &
                         F_npath='V',F_bpath=bpath,                &
                         F_start=(/  1,  1,kstr/),                  &
                         F_end=  (/ -1, -1,kend/) )

         else ! dynamic grid

          if ( cplocn_cvou_S(ivar) == 'P0A' ) then ! special case for 2D surface pressure

            istat = gmm_get(bus_name,ptr2d)
            busou_tmp(cpl_drv_i0:cpl_drv_in, cpl_drv_j0:cpl_drv_jn, ivar) = &
                ptr2d(cpl_drv_i0:cpl_drv_in,cpl_drv_j0:cpl_drv_jn)

          else

            istat = gmm_get(bus_name,ptr3d)
            busou_tmp(cpl_drv_i0:cpl_drv_in, cpl_drv_j0:cpl_drv_jn, ivar:ivar) = &
                ptr3d(cpl_drv_i0:cpl_drv_in,cpl_drv_j0:cpl_drv_jn,kstr:kend)

          endif

         endif ! test grid type

         if (istat < 0) goto 998

         ! clean up the values before accumulating them
         where (isnan(busou_tmp(:,:,ivar))) busou_tmp(:,:,ivar) = 0.
         where (abs(busou_tmp(:,:,ivar))>1e10) busou_tmp(:,:,ivar) = 0.

      enddo

      ! accumulate
      cplocn_it = cplocn_it + 1
      ocn_busou(1:cpl_drv_lni,1:cpl_drv_lnj,1:cplocn_n_fldou) = &
      ocn_busou(1:cpl_drv_lni,1:cpl_drv_lnj,1:cplocn_n_fldou) + &
      busou_tmp(1:cpl_drv_lni,1:cpl_drv_lnj,1:cplocn_n_fldou)

      ! divide when ready
      if ( cplocn_rap_dt == 1 ) then
         cplocn_it = 0 ! no need to divide by 1
      elseif (cplocn_it == cplocn_rap_dt ) then
         ocn_busou(:,:,:) = ocn_busou(:,:,:) / real(cplocn_rap_dt)
         cplocn_it = 0
      endif

      goto 999

 998  call handle_error(istat,'cplocn_busou','Problems filling ocn_busou')
 999  continue

!     ________________________________________________________________
!
      return
      end

