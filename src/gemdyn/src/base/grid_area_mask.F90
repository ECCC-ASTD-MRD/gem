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

!**s/r grid_area_mask - Evaluate area and mask

      subroutine grid_area_mask (F_area_8,F_mask_8,F_area_mask_8,F_ni,F_nj)
      use adz_mem
      use adz_options
      use geomh
      use glb_ld
      use glb_pil
      use HORgrid_options
      use lun
      use ptopo
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,                       intent(in)  :: F_ni,F_nj
      real(kind=REAL64) , dimension(F_ni,F_nj), intent(out) :: F_area_8
      real(kind=REAL64) , dimension(F_ni,F_nj), intent(out) :: F_mask_8
      real(kind=REAL64) , dimension(F_ni,F_nj), intent(out) :: F_area_mask_8

      !object
      !===========================
      !     Evaluate area and mask
      !===========================

      include 'mpif.h'
      include 'rpn_comm.inc'

      integer :: i,j,n,r_j1,r_j2,r_i1,r_i2,r_ish,r_jsh,i_north,j_north, &
                 i_middle,j_middle,err,i1,j1,i2,j2,i_sh,j_sh,comm

      integer, parameter :: MAXB=100000

      real(kind=REAL64) :: blon_8(MAXB+1),blat_8(MAXB+1),iu_8(0:G_ni),jv_8(0:G_nj), &
                g_mask_8(G_ni,G_nj),l_mask_8(l_minx:l_maxx,l_miny:l_maxy), &
                r_ic1_8,r_jc1_8,r_ic2_8,r_jc2_8,h1_8,h2_8,h3_8,ib_8,jb_8, &
                ic1_8,jc1_8,ic2_8,jc2_8,slope_8,sf_8(2),sp_8(2), &
                rlat_8,rlon_8,s_8(2,2),x_8,y_8
      real(kind=REAL64) spa_8(l_ni,l_nj,2)

      real(kind=REAL64) :: gathV(2,Ptopo_numproc)

      logical :: almost_zero,nj_even_L,ni_even_L,between_L,line_L
      logical :: glbsum_L
!
!---------------------------------------------------------------------
!
      glbsum_L=.true.

      do j = 1,F_nj
         do i = 1,F_ni
            F_area_8(i,j) = geomh_hx_8 * cos(geomh_y_8(j)) * geomh_hy_8
         end do
      end do

      comm = RPN_COMM_comm ('GRID')

      !-------
      !GEM LAM
      !-------
      if (.not.Grd_yinyang_L) then

         F_mask_8 = 1.d0

         F_area_mask_8 = F_area_8

         return

      end if

      !------------
      !GEM Yin-Yang
      !------------

      !Prepare list of Mask Boundary points
      !------------------------------------
      do n = 1,MAXB+1

         rlat_8 = -0.25*pi_8

         rlon_8 = (dble(n-1)/dble(MAXB))*0.5d0*pi_8 + pi_8 !Range [0,2 PI]

         x_8 = rlon_8 - pi_8
         y_8 = rlat_8

         call smat(s_8,rlon_8,rlat_8,x_8,y_8)

         rlon_8 = rlon_8 + pi_8

         blon_8(n) = rlon_8
         blat_8(n) = rlat_8

      end do

      !Initialize Global staggered grid indices
      !----------------------------------------
      do i = 0,G_ni
         iu_8(i) = .5d0 + i
      end do

      do j = 0,G_nj
         jv_8(j) = .5d0 + j
      end do

      ni_even_L = .false.

      if ((G_ni/2)*2==G_ni) ni_even_L = .true.

      i_middle = int((G_ni+1)/2) !West Range

      !Find NORTHERN cell containing (-PI/2,PI/4)
      !------------------------------------------
      ic1_8 = 1+(0.50*pi_8-G_xg_8(1))/geomh_hx_8 !Range [0,2 PI]
      jc1_8 = 1+(0.25*pi_8-G_yg_8(1))/geomh_hy_8

   L1:do j = (G_nj+1)/2,G_nj
         do i = 1,i_middle

            if (iu_8(i-1)<=ic1_8.and.ic1_8<=iu_8(i).and. &
                jv_8(j-1)<=jc1_8.and.jc1_8<=jv_8(j)) exit L1

         end do
      end do L1

      i_north = i
      j_north = j

      !Check if (-3 PI/4,0) is at cell_corner
      !--------------------------------------
      nj_even_L = .false.

      if ((G_nj/2)*2==G_nj) nj_even_L = .true.

      j_middle = int((G_nj+1)/2)

      if (nj_even_L) j_middle = j_middle + 1 !North Range

      g_mask_8 = 0.

      !---------------------------------------------------------
      !Summation of trapezoids to estimate area of Mask Boundary
      !---------------------------------------------------------
      do n = 1,MAXB

         !Convert lon,lat to grid indices for 2 consecutive Mask Boundary points
         !----------------------------------------------------------------------
         ic1_8 = 1+(blon_8(n  )-G_xg_8(1))/geomh_hx_8
         jc1_8 = 1+(blat_8(n  )-G_yg_8(1))/geomh_hy_8

         ic2_8 = 1+(blon_8(n+1)-G_xg_8(1))/geomh_hx_8
         jc2_8 = 1+(blat_8(n+1)-G_yg_8(1))/geomh_hy_8

         !Find cell locations
         !-------------------
         j1 = nint(jc1_8)
         j2 = nint(jc2_8)

         i1 = nint(ic1_8)
         i2 = nint(ic2_8)

         !Set P1=Mask Boundary point(n) and P2=Mask Boundary point(n+1)
         !-------------------------------------------------------------
         r_i1 = i1
         r_j1 = j1

         r_ic1_8 = ic1_8
         r_jc1_8 = jc1_8

         r_i2 = i2
         r_j2 = j2

         r_ic2_8 = ic2_8
         r_jc2_8 = jc2_8

         r_ish = 0
         r_jsh = 0

         !Estimate slope between 2 consecutives Mask Boundary points
         !----------------------------------------------------------
         slope_8 = (jc2_8-jc1_8)/(ic2_8-ic1_8)

         between_L = .true.

         !----------------------------------------------------------
         !Calculate area between 2 consecutives Mask Boundary points
         !----------------------------------------------------------
         do while (between_L)

            !Calculate underlying area between P1 and P2 in the same cell location
            !---------------------------------------------------------------------
            if (r_i1==r_i2.and.r_j1==r_j2) then

               if (r_j1/=j_middle.or.nj_even_L) then
                  h1_8 = r_jc1_8    - jv_8(r_j1-1)
                  h2_8 = r_jc2_8    - jv_8(r_j1-1)
                  h3_8 = jv_8(r_j1) - jv_8(r_j1-1)
               else
                  h1_8 = r_jc1_8    - j_middle
                  h2_8 = r_jc2_8    - j_middle
                  h3_8 = jv_8(r_j1) - j_middle
               end if

               g_mask_8(r_i1,r_j1) = g_mask_8(r_i1,r_j1) + abs(0.5*(h1_8 + h2_8)*(r_ic2_8 - r_ic1_8))

               if (almost_zero(r_jc2_8-jv_8(r_j1))) &
               g_mask_8(r_i1,r_j1) = g_mask_8(r_i1,r_j1) + abs(0.5*(h2_8 + h3_8)*(iu_8(r_i1) - r_ic2_8))

               !We go to next Mask Boundary point
               !---------------------------------
               if (r_i2==i2.and.r_j2==j2) then

                  between_L = .false.

                  cycle

               end if

               !Set P1=P2 and P2=Mask Boundary point(n+1)
               !-----------------------------------------
               r_i1 = r_i2 + r_ish
               r_j1 = r_j2 + r_jsh

               r_ic1_8 = r_ic2_8
               r_jc1_8 = r_jc2_8

               r_i2 = i2
               r_j2 = j2

               r_ic2_8 = ic2_8
               r_jc2_8 = jc2_8

               r_ish = 0
               r_jsh = 0

            !If points P1 and P2 are not at the same cell location,
            !find which of (ib,jv(r_j1)) or (iu(r_i1),jb) are in the cell location of P1
            !---------------------------------------------------------------------------
            else

               !Calculate i intersection
               !------------------------
               jb_8 = slope_8 * (iu_8(r_i1) - r_ic1_8) + r_jc1_8

               !Calculate j intersection
               !------------------------
               ib_8 = (jv_8(r_j1) - r_jc1_8)/slope_8 + r_ic1_8

               !Point (ib,jb) is in the cell location of P1
               !-------------------------------------------
               if (((iu_8(r_i1-1)<=ib_8.and.ib_8<=iu_8(r_i1)).or.almost_zero(abs(iu_8(r_i1)-ib_8))).and. &
                   ((jv_8(r_j1-1)<=jb_8.and.jb_8<=jv_8(r_j1)).or.almost_zero(abs(jv_8(r_j1)-jb_8)))) then

                  r_i2 = r_i1
                  r_j2 = r_j1

                  r_ic2_8 = ib_8
                  r_jc2_8 = jb_8

                  r_ish   = 1
                  r_jsh   = 1

               !Point (ib,jv(r_j1)) is in the cell location of P1
               !-------------------------------------------------
               else if ((iu_8(r_i1-1)<=ib_8.and.ib_8<=iu_8(r_i1)).or.almost_zero(abs(iu_8(r_i1)-ib_8))) then

                  r_i2 = r_i1
                  r_j2 = r_j1

                  r_ic2_8 = ib_8
                  r_jc2_8 = jv_8(r_j1)

                  r_ish   = 0
                  r_jsh   = 1

               !Point (iu(r_i1),jb) is in the cell location of P1
               !-------------------------------------------------
               else if ((jv_8(r_j1-1)<=jb_8.and.jb_8<=jv_8(r_j1)).or.almost_zero(abs(jv_8(r_j1)-jb_8))) then

                  r_i2 = r_i1
                  r_j2 = r_j1

                  r_ic2_8 = iu_8(r_i1)
                  r_jc2_8 = jb_8

                  r_ish   = 1
                  r_jsh   = 0

               else

                  call handle_error (-1,'grid_area_mask','BOUNDARY: We are lost')

               end if

            end if

         end do

      end do

      !(-3 PI/4,0) is not a cell_corner: Complete area of those cells by symmetry
      !--------------------------------------------------------------------------
      if (.not.nj_even_L) then
         do i = 1,i_middle
            g_mask_8(i,j_middle) = 2.*g_mask_8(i,j_middle)
         end do
      end if

      !Complete area of NORTHERN cell containing (-PI/2,PI/4)
      !------------------------------------------------------
      ic2_8 = 1+(blon_8(MAXB+1)-G_xg_8(1))/geomh_hx_8
      jc2_8 = 1+(blat_8(MAXB+1)-G_yg_8(1))/geomh_hy_8

   L2:do j = j_north,j_north
         do i = i_middle,1,-1
            if (.not.almost_zero(g_mask_8(i,j))) then
                g_mask_8(i,j) = g_mask_8(i,j) + abs((iu_8(i_north)-ic2_8)*(jc2_8-jv_8(j_north-1)))
                exit L2
            end if
         end do
      end do L2

      !Fill up NORTH-WEST quarter of Mask
      !----------------------------------
      do j = j_middle,j_north

         i = i_middle

         line_L = .true.

         !Fill up the line
         !----------------
         do while (line_L)

            if (.not.almost_zero(g_mask_8(i,j))) then

               line_L = .false.

               cycle

            end if

            if (j==j_north) g_mask_8(i,j) = abs((iu_8(i)-iu_8(i-1))*(jc2_8-jv_8(j_north-1)))
            if (j/=j_north) g_mask_8(i,j) = 1.d0

            i = i-1

            if (i<1) then

               line_L = .false.

               cycle

            end if

         end do

      end do

      j_sh = 0
      if (nj_even_L) j_sh = 1

      !Fill up SOUTH-WEST quarter of Mask
      !----------------------------------
      do j = j_middle,j_north
         do i = 1,i_middle
            g_mask_8(i,2*j_middle-j-j_sh) = g_mask_8(i,j)
         end do
      end do

      i_sh = 0
      if (ni_even_L) i_sh = 1

      !Fill up EAST half of Mask
      !-------------------------
      do j = 1,j_north
         do i = 1,i_middle
            g_mask_8(2*i_middle-i+i_sh,j) = g_mask_8(i,j)
         end do
      end do

      !Distribute Global Mask
      !----------------------
      call RPN_COMM_dist (g_mask_8,1,G_ni,1,G_nj,G_ni,G_nj,1,0,0,2, &
                          l_mask_8,l_minx,l_maxx,l_miny,l_maxy,0,0,G_periodx,G_periody,err)

      !Preparation to set area of Mask closer to 2*PI
      !----------------------------------------------
      sp_8 = 0.d0
      spa_8 = 0.d0

      do j = 1+pil_s,l_nj-pil_n

         do i = 1+pil_w,l_ni-pil_e

            l_mask_8(i,j) = max(0.d0,min(1.d0,l_mask_8(i,j)))

            if (.NOT.almost_zero(l_mask_8(i,j)-1.d0)) then
               sp_8(1) = sp_8(1) + l_mask_8(i,j)*geomh_area_8(i,j)
               spa_8(i,j,1) =      l_mask_8(i,j)*geomh_area_8(i,j)
            else
               sp_8(2) = sp_8(2) + l_mask_8(i,j)*geomh_area_8(i,j)
               spa_8(i,j,2) =      l_mask_8(i,j)*geomh_area_8(i,j)
            end if

         end do

      end do

      sf_8 = 0.d0


      !MPI reduce to get the global value for sp into sf
      !-------------------------------------------------
      if ( glbsum_L ) then
           call glbsum8 (sf_8,spa_8,1,l_ni,1,l_nj,2,          &
                   1+Lam_pil_w, G_ni-Lam_pil_e, 1+Lam_pil_s, G_nj-Lam_pil_n)
      else
           call MPI_Allgather(sp_8,2,MPI_DOUBLE_PRECISION,gathV,2,MPI_DOUBLE_PRECISION,comm,err)
           do i=1,2
              sf_8(i) = sum(gathV(i,:))
           end do

      end if

      if (Lun_out>0) then
         write(Lun_out,1001)
         write(Lun_out,1002) 'MASK AREA Iteration#0  = ',sf_8(1) + sf_8(2)
         write(Lun_out,1002) 'MASK AREA SF(1)        = ',sf_8(1)
         write(Lun_out,1002) 'MASK AREA SF(2)        = ',sf_8(2)
         write(Lun_out,1002) 'MASK AREA target       = ',2.0*pi_8
         write(Lun_out,1002) 'MASK AREA error        = ',sf_8(1) + sf_8(2) - 2.0*pi_8
         write(Lun_out,*) ''
      end if

      !Do Mask correction
      !------------------
      sp_8 = 0.d0
      spa_8 = 0.d0

      do j = 1+pil_s,l_nj-pil_n

         do i = 1+pil_w,l_ni-pil_e

            x_8 = l_mask_8(i,j)*(2.d0*pi_8 - sf_8(2))/sf_8(1)

            if (.NOT.almost_zero(l_mask_8(i,j)-1.d0)) l_mask_8(i,j) = min( 1.d0, x_8 )

            if (.NOT.almost_zero(l_mask_8(i,j)-1.d0)) then
               sp_8(1) = sp_8(1) + l_mask_8(i,j)*geomh_area_8(i,j)
               spa_8(i,j,1) =      l_mask_8(i,j)*geomh_area_8(i,j)
            else
               sp_8(2) = sp_8(2) + l_mask_8(i,j)*geomh_area_8(i,j)
               spa_8(i,j,2) =      l_mask_8(i,j)*geomh_area_8(i,j)
            end if

         end do

      end do

      sf_8 = 0.d0

      !MPI reduce to get the global value for sp into sf
      !-------------------------------------------------
      if ( glbsum_L ) then
           call glbsum8 (sf_8,spa_8,1,l_ni,1,l_nj,2,          &
                   1+Lam_pil_w, G_ni-Lam_pil_e, 1+Lam_pil_s, G_nj-Lam_pil_n)
      else
           call MPI_Allgather(sp_8,2,MPI_DOUBLE_PRECISION,gathV,2,MPI_DOUBLE_PRECISION,comm,err)
           do i=1,2
              sf_8(i) = sum(gathV(i,:))
           end do
      end if

      if (Lun_out>0) then
         write(Lun_out,1002) 'MASK AREA Iteration#1  = ',sf_8(1) + sf_8(2)
         write(Lun_out,1002) 'MASK AREA SF(1)        = ',sf_8(1)
         write(Lun_out,1002) 'MASK AREA SF(2)        = ',sf_8(2)
         write(Lun_out,1002) 'MASK AREA target       = ',2.0*pi_8
         write(Lun_out,1002) 'MASK AREA error        = ',sf_8(1) + sf_8(2) - 2.0*pi_8
      end if

      do j = 1,F_nj
         do i = 1,F_ni
            F_area_mask_8(i,j) = F_area_8(i,j) * l_mask_8(i,j)
         end do
      end do

      do j = 1,F_nj
         do i = 1,F_ni
            F_mask_8(i,j) = l_mask_8(i,j)
         end do
      end do

 1001 format(/,'INITIALIZATING MASK BOUNDARY (S/R GRID_AREA_MASS)', &
             /,'=================================================')
 1002 format(2X,A25,e22.15)
!
!---------------------------------------------------------------------
!
      return
      end
