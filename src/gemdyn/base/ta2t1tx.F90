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

      subroutine ta2t1tx
      use gmm_vt1
      use gmm_vta
      use gem_options
      use mem_tracers
      implicit none

#include <arch_specific.hf>
!
!     ---------------------------------------------------------------
!
      ut1  = uta ; vt1  = vta ; wt1  = wta
      tt1  = tta ; zdt1 = zdta; qt1  = qta
      st1  = sta ; trt1 = trdf
!
!     ---------------------------------------------------------------
!
      return
      end
