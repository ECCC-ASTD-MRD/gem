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

!**s/r  sol_mat_vec - call 2D or 3D matrix_vector subroutine
!
      subroutine  sol_matvec ( wk22, wk11, Minx, Maxx, Miny, Maxy,  &
                               nil,njl, nk, minx1,maxx1,minx2,maxx2 )
      use dyn_fisl_options
      use sol
      implicit none
#include <arch_specific.hf>
!
      integer Minx,Maxx,Miny,Maxy, njl,nil,nk, minx1,maxx1,minx2,maxx2
      real*8  wk11(*),wk22(*)
!
!author
!     Abdessamad Qaddouri - decembre 2013
!
!revision
! v4_6 - Qaddouri A.      - initial version
!
!
      integer i0,in,j0,jn,minx3,maxx3
      real*8 fdg1   ( minx1:maxx1, minx2:maxx2,nk ), &
             wint_8 (  Minx:Maxx ,  Miny:Maxy ,nk ), &
             wint_81(  Minx:Maxx ,  Miny:Maxy ,nk )
!
!     ---------------------------------------------------------------
!
      i0 = 1   + sol_pil_w
      in = nil - sol_pil_e
      j0 = 1   + sol_pil_s
      jn = njl - sol_pil_n

      minx3= 0 ; maxx3= nk+1

      call  tab_vec ( wint_81, Minx,Maxx,Miny,Maxy,nk, &
                      wk11   , i0,in,j0,jn, -1 )

      if (Sol_type_S(11:12) == '2D') then

         call mat_vecs2D ( wint_81, wint_8, Minx, Maxx, Miny, Maxy,nil, &
                           njl,minx1, maxx1, minx2, maxx2,nk, fdg1 )
      else

         call mat_vecs3D ( wint_81, wint_8, Minx, Maxx, Miny, Maxy,nil, &
                           njl,minx1, maxx1, minx2, maxx2,nk,minx3,maxx3)
      end if

      call  tab_vec ( wint_8 , Minx,Maxx,Miny,Maxy,nk, &
                      wk22   , i0,in,j0,jn, +1 )
!
!     ---------------------------------------------------------------
!
      return
      end

