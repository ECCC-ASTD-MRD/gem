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


      subroutine cplocn_update (F_f2u, F_name_S, F_pos, F_ni, cplu, &
                                F_rho, F_u, F_v, F_vmod, F_cmu)
      use rmn_gmm
      use cpl_mod
      use cplocn_mod
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_name_S
      integer, intent(in)          :: F_ni
      integer, dimension(2,F_ni),           intent(in)    :: F_pos
      real   , dimension  (F_ni),           intent(inout) :: F_f2u
      real   , dimension  (F_ni),           intent(in)    :: F_rho, F_u, F_v, F_vmod
      real   , dimension  (F_ni),           intent(inout) :: F_cmu
      logical,                              intent(out)   :: cplu

!authors    Francois Roy -- spring 2015
! 
!revision
! v4_80 - Roy, F.  - initial version

!Purpose
! Update field F_f2u

      integer i,idrv,jdrv,iv,ivU,ivV,ier
      real, parameter :: eps_cplocn_blnd = 1.e-4
      real  mask, var, tau, minval, maxval, maskmax
      real, dimension(:,:), pointer :: surf0

      integer ii

!     ---------------------------------------------------------------

      maskmax=0.
      cplu = .false.

!$omp critical
      ier = gmm_get(gmmk_ocn_busin_s , ocn_busin)
!$omp end critical

      select case(F_name_S)
         case ('GLI','TMO','I8I','SDI')
!$omp critical
            if ( F_name_S == 'GLI' ) then
               ier = gmm_get(gmmk_gli_0_s,gli_0)
               iv  = icvin_GLI
               minval=0.
               maxval=1.
               surf0 => gli_0
            endif
            if ( F_name_S == 'TMO' ) then
               ier = gmm_get(gmmk_tmo_0_s,tmo_0)
               iv  = icvin_TMO
               minval=-9999.
               maxval=9999.
               surf0 => tmo_0
            endif
            if ( F_name_S == 'I8I' ) then
               ier = gmm_get(gmmk_i8i_0_s,i8i_0)
               iv  = icvin_I8I
               minval=0.
               maxval=9999.
               surf0 => i8i_0
            endif
            if ( F_name_S == 'SDI' ) then
               ier = gmm_get(gmmk_sdi_0_s,sdi_0)
               iv  = icvin_SDI
               minval=0.
               maxval=9999.
               surf0 => sdi_0
            endif
!$omp end critical
            do i=1,F_ni
               idrv = cpl_drv_i0 + F_pos(1,i) - 1 
               jdrv = cpl_drv_j0 + F_pos(2,i) - 1
               mask = ocn_busin(idrv,jdrv,icvin_MCP)
               maskmax=max(mask,maskmax)
               var  = ocn_busin(idrv,jdrv,iv)
               if ( mask > eps_cplocn_blnd ) then
                  if (surf0(idrv,jdrv) == cplocn_missval)  &
                      surf0(idrv,jdrv) = F_f2u(i)
                  F_f2u(i) = surf0(idrv,jdrv) * (1.-mask) &
                           + var              *     mask
                  F_f2u(i) = min(maxval,max(minval,F_f2u(i)))
               endif
            enddo
         case ('UVI','UVO')
            if ( F_name_S == 'UVI' ) then
               ivU = icvin_UUI
               ivV = icvin_VVI
            endif
            if ( F_name_S == 'UVO' ) then
               ivU = icvin_UUO
               ivV = icvin_VVO
            endif
            do i=1,F_ni
               idrv = cpl_drv_i0 + F_pos(1,i) - 1 
               jdrv = cpl_drv_j0 + F_pos(2,i) - 1
               mask = ocn_busin(idrv,jdrv,icvin_MCP)
               maskmax=max(mask,maskmax)
               F_f2u(i)= F_f2u(i) * (1.-mask) +    &
                 SQRT((F_u(i)-ocn_busin(idrv,jdrv,ivU))**2  &
                     +(F_v(i)-ocn_busin(idrv,jdrv,ivV))**2) &
                                  *     mask
            enddo
         case ('FRI','FRO')
            if ( F_name_S == 'FRI' ) then
               ivU = icvin_TXI
               ivV = icvin_TYI
            endif
            if ( F_name_S == 'FRO' ) then
               ivU = icvin_TXO
               ivV = icvin_TYO
            endif
            do i=1,F_ni
               idrv = cpl_drv_i0 + F_pos(1,i) - 1 
               jdrv = cpl_drv_j0 + F_pos(2,i) - 1
               mask = ocn_busin(idrv,jdrv,icvin_MCP)
               maskmax=max(mask,maskmax)
               tau = SQRT(ocn_busin(idrv,jdrv,ivU)**2 + &
                          ocn_busin(idrv,jdrv,ivV)**2)
               F_f2u(i)= F_f2u(i) * (1.-mask) +    &
                         SQRT(tau/F_rho(i))        &
                                  *     mask
               F_cmu(i)= F_cmu(i) * (1.-mask) +    &
                         F_f2u(i)*F_f2u(i)/max(F_vmod(i),1.e-6) &
                                  *     mask
            enddo
         case ('T4I','T4O')
            if ( F_name_S == 'T4I' ) iv = icvin_T4I
            if ( F_name_S == 'T4O' ) iv = icvin_T4O
            do i=1,F_ni
               idrv = cpl_drv_i0 + F_pos(1,i) - 1 
               jdrv = cpl_drv_j0 + F_pos(2,i) - 1
               mask = ocn_busin(idrv,jdrv,icvin_MCP)
               maskmax=max(mask,maskmax)
               F_f2u(i)= F_f2u(i) * (1.-mask) +          &
                         (ocn_busin(idrv,jdrv,iv)**0.25) &
                                  *     mask
            enddo
         case ('ZMO', 'ZMI', 'ZHO', 'ZHI')
            ! The log value was transfered
            if ( F_name_S == 'ZMO' ) iv = icvin_ZMO
            if ( F_name_S == 'ZMI' ) iv = icvin_ZMI
            if ( F_name_S == 'ZHO' ) iv = icvin_ZHO
            if ( F_name_S == 'ZHI' ) iv = icvin_ZHI
            do i=1,F_ni
               idrv = cpl_drv_i0 + F_pos(1,i) - 1 
               jdrv = cpl_drv_j0 + F_pos(2,i) - 1
               mask = ocn_busin(idrv,jdrv,icvin_MCP)
               maskmax=max(mask,maskmax)
               F_f2u(i)= F_f2u(i) * (1.-mask) +          &
                         exp(ocn_busin(idrv,jdrv,iv))    &
                                  *     mask
            end do
         case DEFAULT
            select case(F_name_S)
               case ('SHO')
                  iv = icvin_SHO
               case ('LHO')
                  iv = icvin_LHO
               case ('ZUO')
                  iv = icvin_ZUO
               case ('ZVO')
                  iv = icvin_ZVO
               case ('ZTO')
                  iv = icvin_ZTO
               case ('ZQO')
                  iv = icvin_ZQO
               case ('SHI')
                  iv = icvin_SHI
               case ('LHI')
                  iv = icvin_LHI
               case ('ZUI')
                  iv = icvin_ZUI
               case ('ZVI')
                  iv = icvin_ZVI
               case ('ZTI')
                  iv = icvin_ZTI
               case ('ZQI')
                  iv = icvin_ZQI
               case ('I7I')
                  iv = icvin_I7I
               case ('QSO')
                  iv = icvin_QSO
               case ('QSI')
                  iv = icvin_QSI
               case ('ILO')
                  iv = icvin_ILO
               case ('ILI')
                  iv = icvin_ILI
               case ('TMW')         !special name to update TSURF in water routine
                  iv = icvin_TMO    !rather than TWATER in itf_sfc_main routine (impact?)
               case DEFAULT
                  write(6,*) 'cplocn_update problem '
                  call flush(6)
                  stop 'SHOULD NEVER COME HERE'
            end select
            do i=1,F_ni
               idrv = cpl_drv_i0 + F_pos(1,i) - 1 
               jdrv = cpl_drv_j0 + F_pos(2,i) - 1
               mask = ocn_busin(idrv,jdrv,icvin_MCP)
               maskmax=max(mask,maskmax)
               F_f2u(i)= F_f2u(i) * (1.-mask) +    &
                         ocn_busin(idrv,jdrv,iv)  &
                                  *     mask
            enddo
      end select

      if (maskmax > 0.) cplu=.true.
!
!     ---------------------------------------------------------------
!
      return
      end subroutine cplocn_update
