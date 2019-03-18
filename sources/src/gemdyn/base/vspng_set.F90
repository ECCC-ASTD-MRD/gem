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

!**s/r vspng_set - Vertical sponge setup

      subroutine vspng_set
      use dcst
      use HORgrid_options
      use hvdif_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use ver
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer i,j,k
      real*8, dimension(:), allocatable :: weigh
      real*8 pis2_8,pbot_8,delp_8,c_8,nutop
!
!     ---------------------------------------------------------------
!
      Vspng_niter = 0
      if (Vspng_coeftop < 0.) Vspng_nk = 0
      if (Vspng_nk <= 0) return

      if (Lun_out > 0) write (Lun_out,1001)

      Vspng_nk = min(G_nk,Vspng_nk)
      pis2_8   = pi_8/2.0d0

      pbot_8 = exp(Ver_z_8%m(Vspng_nk+1))
      delp_8 = pbot_8 - exp(Ver_z_8%m(1))

      allocate ( weigh(Vspng_nk), Vspng_coef_8(Vspng_nk) )
      do k=1,Vspng_nk
         weigh(k) = (sin(pis2_8*(pbot_8-exp(Ver_z_8%m(k)))/(delp_8)))**2
      end do

      i=G_ni/2
      j=G_nj/2
      c_8 = min ( G_xg_8(i+1) - G_xg_8(i), G_yg_8(j+1) - G_yg_8(j) )

      nutop = Vspng_coeftop*Cstv_dt_8/(Dcst_rayt_8*c_8)**2
      Vspng_niter = int(8.d0*nutop+0.9999999)

      if (Lun_out > 0) then
         write (Lun_out,2002) Vspng_coeftop,Vspng_nk,nutop,Vspng_niter
         write (Lun_out,3001)
      end if

      nutop = dble(Vspng_coeftop)/max(1.d0,dble(Vspng_niter))

      do k=1,Vspng_nk
         Vspng_coef_8(k) = weigh(k) * nutop
         if (Lun_out > 0) write (Lun_out,2005) &
                       Vspng_coef_8(k)      ,&
                       Vspng_coef_8(k)*Cstv_dt_8/(Dcst_rayt_8*c_8)**2,&
                       exp(Ver_z_8%m(k)),k
         Vspng_coef_8(k)= Vspng_coef_8(k) * Cstv_dt_8 * Dcst_inv_rayt_8**2
      end do

 1001 format(/,'INITIALIZATING SPONGE LAYER PROFILE ',  &
               '(S/R VSPNG_SET)',/,51('='))
 2002 format('  SPONGE LAYER PROFILE BASED ON: Vspng_coeftop=',1pe10.2, &
             '  m**2 AND Vspng_nk=',i3/'  Nu_top=',1pe14.6, &
             '  Vspng_niter=',i8)
 2005 format(1pe14.6,1pe14.6,f11.2,i8)
 3001 format('     Coef           Nu            Pres      Level')
!
!     ---------------------------------------------------------------
!
      return
      end
