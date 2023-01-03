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

!**s/r mass_tr_PM - Evaluate Mass of Tracer TIME P/M and Mass of Flux_out/Flux_in
!                   for all tracers using Bermejo-Conde (assuming in Mixing Ratio)

      subroutine mass_tr_PM ( F_mass_p_8,F_mass_m_8,F_mass_fo_8,F_mass_fi_8,F_bc,        &
                              F_air_mass_p,F_air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,&
                              F_i0,F_in,F_j0,F_jn,F_k0,F_ntr_bc )

      use adz_mem
      use adz_options
      use ctrl
      use dynkernel_options
      use omp_timing
      use geomh
      use HORgrid_options
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,           intent(in)                       :: F_ntr_bc                    !Number of tracers using B-C
      type(bc),          intent(in),  dimension(F_ntr_bc) :: F_bc                        !Pointers for tracers using B-C
      real(kind=REAL64), intent(out), dimension(F_ntr_bc) :: F_mass_p_8,F_mass_m_8       !Mass of Tracer TIME P/M
      real(kind=REAL64), intent(out), dimension(F_ntr_bc) :: F_mass_fo_8,F_mass_fi_8     !Mass of Flux_out/_in
      integer,           intent(in)                       :: F_minx,F_maxx,F_miny,F_maxy !Dimension H
      integer,           intent(in)                       :: F_nk                        !Number of vertical levels
      integer,           intent(in)                       :: F_i0,F_in,F_j0,F_jn,F_k0    !Scope of operator Tracer
      real, intent(in),dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: F_air_mass_p,&     !Air mass TIME P
                                                                      F_air_mass_m       !Air mass TIME M

      !object
      !===================================================================
      !     Evaluate Mass of Tracer TIME P/M and Mass of Flux_out/Flux_in
      !     for all tracers using Bermejo-Conde (assuming in Mixing Ratio)
      !===================================================================

      include 'mpif.h'
      include 'rpn_comm.inc'
      integer :: i,j,k,err,n,comm
      real(kind=REAL64), dimension(F_ntr_bc,4) :: c_mass_8,gc_mass_8
      real(kind=REAL64) :: gathV1(2*F_ntr_bc,Ptopo_numproc*Ptopo_ncolors), &
                           gathV2(4*F_ntr_bc,Ptopo_numproc*Ptopo_ncolors)
      logical :: LAM_L,BC_LAM_Aranami_L
!
!---------------------------------------------------------------------
!
!      call gtmg_start (15, 'MASS__', 74)

      LAM_L = .not.Grd_yinyang_L

      BC_LAM_Aranami_L = LAM_L.and.Adz_BC_LAM_flux==1

      c_mass_8 = 0.0d0

!      call gtmg_start (18, 'SOMME_', 15)

      do n=1,F_ntr_bc

         !Evaluate Local Mass of Tracer TIME P/M
         !--------------------------------------
         if (Schm_autobar_L) then

            do j=F_j0,F_jn
               do i=F_i0,F_in
                  c_mass_8(n,1) = c_mass_8(n,1) + F_bc(n)%p(i,j,1) * geomh_area_mask_8(i,j)
                  c_mass_8(n,2) = c_mass_8(n,2) + F_bc(n)%m(i,j,1) * geomh_area_mask_8(i,j)
               end do
            end do

         else

            do k=F_k0,F_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in
                     c_mass_8(n,1) = c_mass_8(n,1) + F_bc(n)%p(i,j,k) * F_air_mass_p(i,j,k) * geomh_mask_8(i,j)
                     c_mass_8(n,2) = c_mass_8(n,2) + F_bc(n)%m(i,j,k) * F_air_mass_m(i,j,k) * geomh_mask_8(i,j)
                  end do
               end do
            end do

         end if

         !Evaluate Local Mass of Flux_out/Flux_in
         !---------------------------------------
         if (LAM_L.and.BC_LAM_Aranami_L.and..not.Ctrl_theoc_L) then

            if (Schm_autobar_L) then

               do j=Adz_j0b,Adz_jnb
                  do i=Adz_i0b,Adz_inb
                     c_mass_8(n,3) = c_mass_8(n,3) + F_bc(n)%fo(i,j,1) * geomh_area_mask_8(i,j)
                     c_mass_8(n,4) = c_mass_8(n,4) + F_bc(n)%fi(i,j,1) * geomh_area_mask_8(i,j)
                  end do
               end do

            else

               do k=1,F_nk
                  do j=Adz_j0b,Adz_jnb
                     do i=Adz_i0b,Adz_inb
                        c_mass_8(n,3) = c_mass_8(n,3) + F_bc(n)%fo(i,j,k) * F_air_mass_m(i,j,k) * geomh_mask_8(i,j)
                        c_mass_8(n,4) = c_mass_8(n,4) + F_bc(n)%fi(i,j,k) * F_air_mass_m(i,j,k) * geomh_mask_8(i,j)
                     end do
                  end do
               end do

            end if

         end if

      end do

!      call gtmg_stop  (18)

!      call gtmg_start (19, 'REDUCE', 15)

      comm = RPN_COMM_comm ('MULTIGRID')

      !Evaluate Global Mass
      !----------------------------------------
      if (LAM_L.and.BC_LAM_Aranami_L.and..not.Ctrl_theoc_L) then

         call MPI_Allgather(c_mass_8,4*F_ntr_bc,MPI_DOUBLE_PRECISION,gathV2,4*F_ntr_bc,MPI_DOUBLE_PRECISION,comm,err)

         do j=1,4
         do i=1,F_ntr_bc

            n = (j-1)*F_ntr_bc + i
            gc_mass_8(i,j) = sum(gathV2(n,:))

         end do
         end do

         F_mass_p_8 (1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,1)
         F_mass_m_8 (1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,2)
         F_mass_fo_8(1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,3)
         F_mass_fi_8(1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,4)

      else

         call MPI_Allgather(c_mass_8,2*F_ntr_bc,MPI_DOUBLE_PRECISION,gathV1,2*F_ntr_bc,MPI_DOUBLE_PRECISION,comm,err)

         do j=1,2
         do i=1,F_ntr_bc

            n = (j-1)*F_ntr_bc + i
            gc_mass_8(i,j) = sum(gathV1(n,:))

         end do
         end do

         F_mass_p_8 (1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,1)
         F_mass_m_8 (1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,2)
         F_mass_fo_8(1:F_ntr_bc) = 0.0
         F_mass_fi_8(1:F_ntr_bc) = 0.0

      end if

!      call gtmg_stop (19)

!      call gtmg_stop (15)
!
!---------------------------------------------------------------------
!
      return
      end subroutine mass_tr_PM
