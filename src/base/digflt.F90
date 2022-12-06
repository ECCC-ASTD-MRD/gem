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

     subroutine digflt
      use step_options
      use gmm_vt1
      use gmm_vta
      use init_options
      use glb_ld
      use mem_tracers
      implicit none

#include <arch_specific.hf>

      integer i, j, k
      real dfcoef
!     __________________________________________________________________
!
      dfcoef = Init_dfco ( abs( (Init_halfspan - Step_kount ) ) )

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

      if ( Init_dftr_L ) then
         trdf = trdf + dfcoef * trt1
      else
         if ( Step_kount == Init_halfspan ) trdf = trt1
      endif
!     __________________________________________________________________
!
      return
      end

