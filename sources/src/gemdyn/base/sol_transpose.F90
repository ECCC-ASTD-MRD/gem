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

!**s/r sol_transpose - Establish memory layout for the elliptic solver
!                      and setup for the 2 level data transposition

      integer function sol_transpose ( F_npex, F_npey, F_checkparti_L )
      use dyn_fisl_options
      use glb_ld
      use ldnh
      use lun
      use trp
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_checkparti_L
      integer, intent(in) :: F_npex, F_npey

      logical, external :: decomp
      integer, parameter :: lowest = 2
      integer :: minx, maxx, n, npartiel, n0, n1
!
!     ---------------------------------------------------------------
!
      sol_transpose = 0
      if (Lun_out > 0) write (Lun_out,1000)

! Establishing local dimensions and computing arena (data topology) for:
!          G_ni distributed on Ptopo_npex PEs and
!          G_nj distributed on Ptopo_npey PEs both without halo

      if (Lun_out > 0) then
         if (Sol_type_S == 'DIRECT') then
            write(Lun_out,1002) ' Transpose 1===>2 for SOLVER (no halo):', &
                                ' G_ni distributed on F_npex PEs', G_ni,F_npex
         else
            write(Lun_out,1002) ' Memory layout for SOLVER (no halo):', &
                                ' G_ni distributed on F_npex PEs', G_ni,F_npex
         end if
      end if

      if (.not. decomp (G_ni, ldnh_minx, ldnh_maxx, ldnh_ni, npartiel, 0, n0, &
                        .true., .true., F_npex, lowest, F_checkparti_L, 0 ))       &
      sol_transpose = -1

      if (Lun_out > 0) then
         if (Sol_type_S == 'DIRECT') then
            write(Lun_out,1002) ' Transpose 1===>2 for SOLVER (no halo):', &
                                ' G_nj distributed on F_npey PEs', G_nj,F_npey
         else
            write(Lun_out,1002) ' Memory layout for SOLVER (no halo):', &
                                ' G_nj distributed on F_npey PEs', G_nj,F_npey
         end if
      end if

      if (.not. decomp (G_nj, ldnh_miny, ldnh_maxy, ldnh_nj, npartiel, 0, n1, &
                        .false.,.true., F_npey, lowest, F_checkparti_L, 0 ))       &
      sol_transpose = -1

      trp_22n = -1

      if (Sol_type_S == 'DIRECT') then

! Transpose 1===>2:Schm_nith distributed on F_npex PEs (no halo)
!                  G_nj distributed on F_npey PEs (original layout)
!                  G_ni NOT distributed
!         initial layout  : (l_minx:l_maxx,    l_miny:l_maxy    ,G_nk)
!         transpose layout: (l_miny:l_maxy,trp_12smin:trp_12smax,G_ni)

         if (Lun_out > 0) write(Lun_out,1002) ' Transpose 1===>2 for SOLVER (no halo):', &
                  ' Schm_nith distributed on F_npex PEs', Schm_nith,F_npex

         if (.not. decomp (Schm_nith, minx, maxx, n, npartiel, 0, n0, &
                           .true., .true., F_npex, -1, F_checkparti_L, 3 ))     &
         sol_transpose = -1

         trp_12smin = minx ! most likely = 1 since no halo
         trp_12smax = maxx
         trp_12sn   = n
         trp_12sn0  = n0

! Transpose 2===>2:G_nk distributed on F_npex PEs (no halo)
!                  G_nj NOT distributed
!                  G_ni distributed on F_npey PEs (no Halo)
!  initial layout  : (    l_miny:l_maxy    ,trp_12smin:trp_12smax,G_ni)
!  transpose layout: (trp_12smin:trp_12smax,trp_22min :trp_22max ,G_nj)

         if (Lun_out > 0) write(Lun_out,1002) ' Transpose 2===>2 for SOLVER (no halo):', &
                                              ' G_ni distributed on F_npey PEs', G_ni,F_npey

         if (.not. decomp (G_ni, minx, maxx, n, npartiel, 0, n0, &
                           .false., .true., F_npey, lowest, F_checkparti_L, 0 )) &
         sol_transpose = -1

         trp_22min = minx ! most likely = 1 since no halo
         trp_22max = maxx
         trp_22n   = n
         trp_22n0  = n0

      end if

      if (sol_transpose < 0)  then
         if  (Lun_out > 0) then
            write(lun_out,*) 'SOL_TRANSPOSE: ILLEGAL DOMAIN PARTITIONING'
         end if
      end if

 1000 format (/' SOL_TRANSPOSE: checking SOLVER dimension partitionning:')
 1002 format (a/a45,i6,' /',i5)
!
!     ---------------------------------------------------------------
!
      return
      end

