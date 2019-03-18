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
module bmf_modsplit
   save
   character(len=1024), allocatable, dimension(:) :: split_files
   integer, allocatable, dimension(:) :: split_unit
   logical split_started
   integer bmf_npex,bmf_npey
   integer bmf_haloileft, bmf_haloiright
   integer bmf_halojleft, bmf_halojright
   integer bmf_ghaloileft, bmf_ghaloiright
   integer bmf_ghalojleft, bmf_ghalojright
   integer bmf_nig, bmf_njg

   type bmf_hole
   character(len=4) bmf_holename 
   integer bmf_holet1, bmf_holet2
   integer bmf_holei0, bmf_holeisize
   integer bmf_holej0, bmf_holejsize
   end type bmf_hole

   type bmf_holelist
          type(bmf_hole) bmf_trou
          type(bmf_holelist), pointer :: trou_suivant
   end type bmf_holelist

   type(bmf_holelist), pointer :: holelist

end module bmf_modsplit
