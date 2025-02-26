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

      subroutine mass_tr_PM_hlt ( F_mass_p_8,F_mass_m_8,F_mass_fo_8,F_mass_fi_8,F_bc,        &
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
      use omp_lib
      use masshlt
      use mem_tstp

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
      integer :: i,j,k,err,n,comm,dim
      real(kind=REAL64), dimension(F_ntr_bc,4) :: c_mass_8
      real(kind=REAL64) :: gathV1(2*F_ntr_bc,Ptopo_numproc*Ptopo_ncolors), &
                           gathV2(4*F_ntr_bc,Ptopo_numproc*Ptopo_ncolors)
      real(kind=REAL64), dimension(:,:), pointer :: gc_mass_8
      logical :: LAM_L,BC_LAM_Aranami_L

      OMP_max_threads=OMP_get_max_threads()
      thread_sum2(1:4,1:F_ntr_bc,0:OMP_max_threads-1) => WS1_8(1:) ; dim= 4*OMP_max_threads*F_ntr_bc
      gc_mass_8(1:F_ntr_bc,1:4) => WS1_8(dim+1:) ; dim= dim+4*F_ntr_bc

!
!---------------------------------------------------------------------
!
!      call gtmg_start (15, 'MASS__', 74)

      comm = COMM_multigrid

      LAM_L = .not.Grd_yinyang_L

      BC_LAM_Aranami_L = LAM_L.and.Adz_BC_LAM_flux==1

      c_mass_8 = 0.0d0

!      call gtmg_start (18, 'SOMME_', 15)

      do n=1,F_ntr_bc

         !Evaluate Local Mass of Tracer TIME P/M
         !--------------------------------------
         if (Dynamics_autobar_L) then

!$omp do
            do j=F_j0,F_jn
               do i=F_i0,F_in
                  c_mass_8(n,1) = c_mass_8(n,1) + F_bc(n)%p(i,j,1) * geomh_area_mask_8(i,j)
                  c_mass_8(n,2) = c_mass_8(n,2) + F_bc(n)%m(i,j,1) * geomh_area_mask_8(i,j)
               end do
            end do
!$omp end do nowait 
            thread_sum2(1,n,OMP_get_thread_num())=c_mass_8(n,1)
            thread_sum2(2,n,OMP_get_thread_num())=c_mass_8(n,2)

         else

!$omp do collapse(2)
            do k=F_k0,F_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in
                     c_mass_8(n,1) = c_mass_8(n,1) + F_bc(n)%p(i,j,k) * F_air_mass_p(i,j,k) * geomh_mask_8(i,j)
                     c_mass_8(n,2) = c_mass_8(n,2) + F_bc(n)%m(i,j,k) * F_air_mass_m(i,j,k) * geomh_mask_8(i,j)
                  end do
               end do
            end do
!$omp end do nowait
            thread_sum2(1,n,OMP_get_thread_num())=c_mass_8(n,1)
            thread_sum2(2,n,OMP_get_thread_num())=c_mass_8(n,2)

         end if

         !Evaluate Local Mass of Flux_out/Flux_in
         !---------------------------------------
         if (LAM_L.and.BC_LAM_Aranami_L.and..not.Ctrl_theoc_L) then

            if (Dynamics_autobar_L) then

!$omp do
               do j=Adz_j0b,Adz_jnb
                  do i=Adz_i0b,Adz_inb
                     c_mass_8(n,3) = c_mass_8(n,3) + F_bc(n)%fo(i,j,1) * geomh_area_mask_8(i,j)
                     c_mass_8(n,4) = c_mass_8(n,4) + F_bc(n)%fi(i,j,1) * geomh_area_mask_8(i,j)
                  end do
               end do
!$omp end do nowait 
            thread_sum2(3,n,OMP_get_thread_num())=c_mass_8(n,3)
            thread_sum2(4,n,OMP_get_thread_num())=c_mass_8(n,4)

            else

!$omp do collapse(2)
               do k=1,F_nk
                  do j=Adz_j0b,Adz_jnb
                     do i=Adz_i0b,Adz_inb
                        c_mass_8(n,3) = c_mass_8(n,3) + F_bc(n)%fo(i,j,k) * F_air_mass_m(i,j,k) * geomh_mask_8(i,j)
                        c_mass_8(n,4) = c_mass_8(n,4) + F_bc(n)%fi(i,j,k) * F_air_mass_m(i,j,k) * geomh_mask_8(i,j)
                     end do
                  end do
               end do
!$omp end do nowait 
            thread_sum2(3,n,OMP_get_thread_num())=c_mass_8(n,3)
            thread_sum2(4,n,OMP_get_thread_num())=c_mass_8(n,4)

            end if

         end if

      end do

!      call gtmg_stop  (18)

!      call gtmg_start (19, 'REDUCE', 15)


      !Evaluate Global Mass
      !----------------------------------------
      if (LAM_L.and.BC_LAM_Aranami_L.and..not.Ctrl_theoc_L) then
!$OMP BARRIER
!$omp single
         do n=1,F_ntr_bc
         c_mass_8(n,1)=sum(thread_sum2(1,n,:))
         c_mass_8(n,2)=sum(thread_sum2(2,n,:))
         c_mass_8(n,3)=sum(thread_sum2(3,n,:))
         c_mass_8(n,4)=sum(thread_sum2(4,n,:))
         enddo
         call MPI_Allgather(c_mass_8,4*F_ntr_bc,MPI_DOUBLE_PRECISION,gathV2,4*F_ntr_bc,MPI_DOUBLE_PRECISION,comm,err)

         do j=1,4
         do i=1,F_ntr_bc

            n = (j-1)*F_ntr_bc + i
             gc_mass_8(i,j) = sum(gathV2(n,:))

         end do
         end do
!$omp end single

         F_mass_p_8 (1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,1)
         F_mass_m_8 (1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,2)
         F_mass_fo_8(1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,3)
         F_mass_fi_8(1:F_ntr_bc) = gc_mass_8(1:F_ntr_bc,4)

      else
!$OMP BARRIER
!$omp single
         do n=1,F_ntr_bc
         c_mass_8(n,1)=sum(thread_sum2(1,n,:))
         c_mass_8(n,2)=sum(thread_sum2(2,n,:))
         enddo
         call MPI_Allgather(c_mass_8,2*F_ntr_bc,MPI_DOUBLE_PRECISION,gathV1,2*F_ntr_bc,MPI_DOUBLE_PRECISION,comm,err)

         do j=1,2
         do i=1,F_ntr_bc

            n = (j-1)*F_ntr_bc + i
            gc_mass_8(i,j) = sum(gathV1(n,:))

         end do
         end do
!$omp end single

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
      end subroutine mass_tr_PM_hlt
