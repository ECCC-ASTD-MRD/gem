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
      use rmn_gmm
      implicit none
#include <arch_specific.hf>
#include <rmn/msg.h>

      include "cpl.cdk"
      include "cplocn.cdk"
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

      integer unout, istat
      logical print_L
      real, dimension(:,:,:), pointer :: ptr3d
      character(len=20) bus_name
      character(len=1)  bpath

      integer ivar,ni,nj,nk,kstr,kend,gni,gnj,gnk
      integer, external :: msg_getUnit
!     ________________________________________________________________

      ni=cpl_drv_lni
      nj=cpl_drv_lnj
      nk=cpl_drv_gnk-1

      gni=cpl_drv_gni
      gnj=cpl_drv_gnj
      gnk=1

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)

      !* FLUX COUPLING OCEAN BUS OUT
      !* Ice-Ocean independant
      !* 1- FB  - SW down                          p
      !* 2- FI  - LW down                          p
      !* 3- RT  - Precipitation                    p
      !* Ocean model flux calculation
      !* 4- TT  - Air temperature                  d
      !* 5- UU  - Wind x component                 d
      !* 6- VV  - Wind y component                 d
      !* 7- QA  - Specific humidity                d
      !* 8- PX  - First momentum level pressure    d
      !* 9- PX  - First thermo   level pressure    d
      !*10- P0  - Ground level pressure            d
      !* Atmospheric model flux calculated (DE-ACTIVATED)
      !*11- SHO - Sensible heat flux over water    ?
      !*12- SHI - Sensible heat flux over ice      ?
      !*13- LHO - Latent heat flux over water      ?
      !*14- LHI - Latent heat flux over ice        ?
      !*15- TXO - Wind stress x component (water)  ?
      !*16- TYO - Wind stress y component (water)  ?
      !*17- TXI - Wind stress x component (ice)    ?
      !*18- TYI - Wind stress y component (ice)    ?

      do ivar=1,cplocn_n_fldou

         ptr3d => ocn_busou(cpl_drv_i0:cpl_drv_in, &
                            cpl_drv_j0:cpl_drv_jn, ivar:ivar)

         SELECT CASE ( cplocn_cvou_S(ivar) )

         CASE ( 'FBA' )

           bus_name='flusolis'
           bpath='P'
           kstr=1
           kend=1

         CASE ( 'FIA' )

           bus_name='fdsi'
           bpath='P'
           kstr=1
           kend=1

         CASE ( 'RTA' )

           bus_name='rt'
           bpath='P'
           kstr=1
           kend=1

         CASE ( 'TTA' )

           bus_name='PW_TT:P'
           bpath='D'
           kstr=nk
           kend=nk

         CASE ( 'UUA' )

           bus_name='PW_UU:P'
           bpath='D'
           kstr=nk
           kend=nk

         CASE ( 'VVA' )

           bus_name='PW_VV:P'
           bpath='D'
           kstr=nk
           kend=nk

         CASE ( 'QQA' )

           bus_name='TR/HU:P'
           bpath='D'
           kstr=nk
           kend=nk

         CASE ( 'PMA' )

           bus_name='PW_PM:P'
           bpath='D'
           kstr=nk
           kend=nk

         CASE ( 'PTA' )

           bus_name='PW_PT:P'
           bpath='D'
           kstr=nk
           kend=nk

         CASE ( 'P0A' )

           bus_name='PW_P0:P'
           bpath='D'
           kstr=1
           kend=1

         CASE DEFAULT

           istat=-1
           call handle_error(istat,'cplocn_busou','WRONG NAME TAG')

         END SELECT

         istat = phy_get(ptr3d, bus_name,                          &
                         F_npath='V',F_bpath=bpath,                &
                         F_start=(/  1,  1,kstr/),                  &
                         F_end=  (/ -1, -1,kend/) )

         if (istat < 0) goto 998

      enddo

      goto 999

 998  call handle_error(istat,'cplocn_busou','Problems filling ocn_busou')
 999  continue

!     ________________________________________________________________
!
      return
      end

