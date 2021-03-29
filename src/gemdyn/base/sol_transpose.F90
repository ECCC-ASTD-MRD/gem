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
      use sol
      use trp
      use numa
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_checkparti_L
      integer, intent(in) :: F_npex, F_npey

      logical, external :: decomp, sol_numa_space
      integer, parameter :: lowest = 2
      integer, dimension(F_npex) :: lnis, nis
      integer, dimension(F_npey) :: lnjs
      integer :: i, minx, maxx, n, npartiel, n0, err1, err2, dim
      real    :: sumk
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
            write(Lun_out,1002) ' Horizontal partitionnong (no halo) for DIRECT Solver:', &
                                ' G_ni distributed on F_npex PEs', G_ni,F_npex
         else
            write(Lun_out,1002) ' Memory layout for SOLVER (no halo):', &
                                ' G_ni distributed on F_npex PEs', G_ni,F_npex
         end if
      end if

      if (.not. decomp (G_ni, ldnh_minx, ldnh_maxx, lnis, npartiel, 0, ldnh_i0, &
               .true., .true., F_npex, lowest, F_checkparti_L, 0 )) sol_transpose = -1
      ldnh_ni= lnis(1)

      if (Lun_out > 0) then
         if (Sol_type_S == 'DIRECT') then
            write(Lun_out,1002) ' Horizontal partitionnong (no halo) for DIRECT Solver:', &
                                ' G_nj distributed on F_npey PEs', G_nj,F_npey
         else
            write(Lun_out,1002) ' Memory layout for SOLVER (no halo):', &
                                ' G_nj distributed on F_npey PEs', G_nj,F_npey
         end if
      end if

      if (.not. decomp (G_nj, ldnh_miny, ldnh_maxy, lnjs, npartiel, 0, ldnh_j0,&
               .false.,.true., F_npey, lowest, F_checkparti_L, 0 )) sol_transpose = -1
      ldnh_nj= lnjs(1)

      trp_22n = -1

      if (Sol_type_S == 'DIRECT') then

! Transpose 1===>2:Schm_nith distributed on F_npex PEs (no halo)
!                  G_nj distributed on F_npey PEs (original layout)
!                  G_ni NOT distributed
!         initial layout  : (l_minx:l_maxx,    l_miny:l_maxy    ,G_nk)
!         transpose layout: (l_miny:l_maxy,trp_12smin:trp_12smax,G_ni)

         if (Lun_out > 0) write(Lun_out,1002) ' Transpose 1===>2 for SOLVER (no halo):', &
                  ' Schm_nith distributed on F_npex PEs', Schm_nith,F_npex

         err1= 0
         if (decomp (Schm_nith, minx, maxx, lnis, npartiel, 0, n0, &
                     .true., .true., F_npex, -1, F_checkparti_L, 2 )) then
            if (F_checkparti_L) then
               n= lnis(1)
               sumk= 0
               do i=1,F_npex
                  if (lnis(i)<1) sumk= sumk+1
               end do
               if (sumk>Schm_nith*.1) then
                  err1= -1
                  if (Lun_out > 0) write(Lun_out,1003) 'NPEX too large for current G_NK transpose'
               endif
            else
               n= lnis(1) ; lnis= 0
               lnis(Ptopo_mycol+1)= n
               call rpn_comm_Allreduce (lnis,nis,F_npex,"MPI_INTEGER", &
                                        "MPI_SUM","EW",err2)
               sumk= 0
               do i=1,F_npex
                  if (nis(i)<1) sumk= sumk+1
               end do
               if (sumk>Schm_nith*.1) then
                  err1= -1
                  if (Lun_out > 0) write(Lun_out,1003) 'NPEX too large for current G_NK transpose'
               endif
            endif
         else
            err1= -1
         endif

         call rpn_comm_Allreduce (err1,sol_transpose,1,"MPI_INTEGER", &
                                  "MPI_SUM","grid",err2)

         if (sol_transpose < 0) goto 9988

         trp_12smin = minx ! most likely = 1 since no halo
         trp_12smax = maxx
         trp_12sn   = n
         trp_12sn0  = n0

         if (sol_one_transpose_L) sol_one_transpose_L= sol_numa_space (F_checkparti_L)
         
         if (.not.sol_one_transpose_L) then
! Transpose 2===>2:G_nk distributed on F_npex PEs (no halo)
!                  G_nj NOT distributed
!                  G_ni distributed on F_npey PEs (no Halo)
!  initial layout  : (    l_miny:l_maxy    ,trp_12smin:trp_12smax,G_ni)
!  transpose layout: (trp_12smin:trp_12smax,trp_22min :trp_22max ,G_nj)

            if (Lun_out > 0) write(Lun_out,1002) ' Transpose 2===>2 for SOLVER (no halo):', &
                                              ' G_ni distributed on F_npey PEs', G_ni,F_npey

            if (.not. decomp (G_ni, minx, maxx, lnis, npartiel, 0, n0, &
                .false., .true., F_npey, lowest, F_checkparti_L, 0 )) sol_transpose = -1
            n= lnis(1)
            trp_22min = minx ! most likely = 1 since no halo
            trp_22max = maxx
            trp_22n   = n
            trp_22n0  = n0
            dim = (trp_12smax-trp_12smin+1)*(trp_22max-trp_22min+1)*G_nj
            if (.not.allocated(Sol_ai_8)) &
            allocate (Sol_ai_8(dim),Sol_bi_8(dim),Sol_ci_8(dim))
         endif

      end if

 9988 if (sol_transpose < 0)  then
         if  (Lun_out > 0) then
            write(lun_out,*) 'SOL_TRANSPOSE: ILLEGAL DOMAIN PARTITIONING'
         end if
      else
         if (sol_one_transpose_L) then
            if (Lun_out > 0) write(Lun_out,'(" SOL_TRANSPOSE: Using one-transpose direct solver")')
         else
            if (Lun_out > 0) write(Lun_out,'(" SOL_TRANSPOSE: Using two-transpose direct solver")')
         end if
      end if

 1000 format (/' SOL_TRANSPOSE: checking SOLVER dimension partitionning:')
 1002 format (a/a45,i6,' /',i5)
 1003    format (/a/)
!
!     ---------------------------------------------------------------
!
      return
      end

