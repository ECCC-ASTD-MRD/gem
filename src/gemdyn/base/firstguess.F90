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

!**s/r firstguess - Copy data from time level t1 that will be used as a
!                   first guess at time level t0

      subroutine firstguess()
      use gmm_itf_mod
      use gmm_vt0
      use gmm_vt1
      use tr3d
      implicit none
#include <arch_specific.hf>

      integer :: istat, k
      real, pointer, contiguous, dimension(:,:,:) :: plus, minus
!
!     ---------------------------------------------------------------
!
      do k=1,Tr3d_ntr
         nullify (plus, minus)
         istat = gmm_get('TR/'//trim(Tr3d_name_S(k))//':M', minus)
         istat = gmm_get('TR/'//trim(Tr3d_name_S(k))//':P', plus )
         minus = plus
      end do

      tt0  = tt1
      zdt0 = zdt1
      wt0  = wt1
      ut0  = ut1
      vt0  = vt1
      qt0  = qt1
      st0  = st1
!
!     ---------------------------------------------------------------
!
      return
      end
