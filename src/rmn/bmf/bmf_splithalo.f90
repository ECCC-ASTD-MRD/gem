!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2005  Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
subroutine bmf_splithalo(il,ir,jl,jr,gil,gir,gjl,gjr)
  use bmf_modsplit
!
! BMF_MODSPLIT, L. Corbeil, 31/08/2001
!ARGUMENTS
!  il,ir      size of right and left halo of the local tile
!  jl,jr      size of upper and lower halo of the local tile
!  gil,gir      size of right and left halo of the local tile
!  gjl,gjr      size of upper and lower halo of the local tile
! 
! To specify halo values like MC2 does, i.e. when a call to
! a bmf_splitwr* will be executed, we will consider that for
! ni points, gil+gir points are in the halo and will not be
! considered for the calculation of the size of the splitted
! domain. Also, the size of the local domain will be modified
! according to global halos (if the tile is on the border of
! the global domain) and local  halos (elsewhere)
!

  implicit none
  integer il,ir,jl,jr,gil,gir,gjl,gjr
  bmf_haloileft=il
  bmf_haloiright=ir
  bmf_halojleft=jl
  bmf_haloiright=jr
  bmf_ghaloileft=gil
  bmf_ghaloiright=gir
  bmf_ghalojleft=gjl
  bmf_ghaloiright=gjr

  return 
end subroutine

