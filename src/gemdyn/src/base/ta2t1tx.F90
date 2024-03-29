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

!**s/r ta2t1tx -  Transfer variables ta into t1

      subroutine ta2t1tx ()
      use glb_ld
      use gmm_vt1
      use gmm_vta
      use gem_options
      use mem_tracers
      use outp
      implicit none
!
!     ---------------------------------------------------------------
!
       ut1(1:l_ni,1:l_nj,:)=  uta(1:l_ni,1:l_nj,:)
       vt1(1:l_ni,1:l_nj,:)=  vta(1:l_ni,1:l_nj,:)
       wt1(1:l_ni,1:l_nj,:)=  wta(1:l_ni,1:l_nj,:)
       tt1(1:l_ni,1:l_nj,:)=  tta(1:l_ni,1:l_nj,:)
      zdt1(1:l_ni,1:l_nj,:)= zdta(1:l_ni,1:l_nj,:)
       qt1(1:l_ni,1:l_nj,:)=  qta(1:l_ni,1:l_nj,:)
       st1(1:l_ni,1:l_nj  )=  sta(1:l_ni,1:l_nj  )
      trt1 = trdf
      udiag(:,:) = diag_dgf(:,:,1)
      vdiag(:,:) = diag_dgf(:,:,2)
      tdiag(:,:) = diag_dgf(:,:,3)
      qdiag(:,:) = diag_dgf(:,:,4)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine ta2t1tx
