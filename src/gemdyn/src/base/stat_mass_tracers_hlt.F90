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

!**s/r stat_mass_tracers_hlt - Calculate and print the mass of each tracer scaled by area

      subroutine stat_mass_tracers_hlt (F_time,F_comment_S)

      use adz_mem
      use adz_options
      use dyn_fisl_options
      use HORgrid_options
      use lun
      use ptopo
      use tdpack, only: grav_8
      use tr3d
      use omp_lib
      use mem_tstp
      use masshlt
      use mem_tracers

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,          intent(in) :: F_time       !Time 0 or Time 1
      character(len=*), intent(in) :: F_comment_S  !Comment

      !object
      !===============================================================
      !     Calculate and print the mass of each tracer scaled by area
      !===============================================================

      integer :: i,j,k,n,k0,count,i0_c,in_c,j0_c,jn_c,i0_sb,in_sb,j0_sb,jn_sb,k0_sb,dim
      real(kind=REAL64) :: tracer_8
      logical :: do_subset_GY_L
      real, pointer, dimension (:,:,:) :: fld_tr
      character(len=21) :: type_S
      character(len=7)  :: time_S
      character(len=8)  :: in_S
!
!---------------------------------------------------------------------
!
      OMP_max_threads=OMP_get_max_threads()
      thread_sum(1:2,0:OMP_max_threads-1) => WS1_8(1:) ; dim= 2*OMP_max_threads
      g_avg_8(1:2) => WS1_8(dim+1:) 

      !Set CORE limits
      !---------------
      i0_c = 1+pil_w ; j0_c = 1+pil_s ; in_c = l_ni-pil_e ; jn_c = l_nj-pil_n ; k0 = Adz_k0t 

      !Terminator: Initialization
      !--------------------------
!$omp single
      call canonical_terminator_0 (count)
!$omp end single

      if (F_time==1) time_S = "TIME T1"
      if (F_time==0) time_S = "TIME T0"

      !Reset Air Mass at appropriate time
      !----------------------------------
      call get_air_mass_hlt (air_mass,F_time,l_minx,l_maxx,l_miny,l_maxy,l_nk,k0)

      type_S = "Mass of Mixing  (WET)"
      if (Schm_dry_mixing_ratio_L) type_S = "Mass of Mixing  (DRY)"

!$omp single
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ''
!$omp end single

      !-------------
      !Treat tracers
      !-------------
      do n=1,Tr3d_ntr

         if (F_time==1) then
            fld_tr=>tracers_P(n)%pntr
            in_S = 'TR/'//trim(Tr3d_name_S(n))//':P'
         end if
         if (F_time==0) then
            fld_tr=>tracers_M(n)%pntr
            in_S = 'TR/'//trim(Tr3d_name_S(n))//':M'
         end if

         !Terminator: Prepare CLY
         !-----------------------
!$omp single
         call canonical_terminator_1 (fld_tr,in_S,count,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk)
!$omp end single

         !Print the mass of each tracer scaled by area (CORE)
         !---------------------------------------------------
         call mass_tr_hlt (tracer_8,fld_tr,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0_c,in_c,j0_c,jn_c,k0)

!$omp single
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,'  C= ', &
                                                                  tracer_8/Adz_gc_area_8,Tr3d_name_S(n)(1:4),F_comment_S
!$omp end single

         do_subset_GY_L = Grd_yinyang_L.and.adz_pil_sub_s_g /= -1

         !Print the mass of each tracer scaled by area (SUBSET of YIN)
         !------------------------------------------------------------
         if (do_subset_GY_L) then

            i0_sb = 1+Adz_pil_sub_w ; j0_sb = 1+Adz_pil_sub_s ; in_sb = l_ni-Adz_pil_sub_e; jn_sb = l_nj-Adz_pil_sub_n ; k0_sb = Adz_k0t_sub

            if (Ptopo_couleur == 0) then 
!$omp do collapse(2)
             do k=k0,l_nk
               do j=j0_sb,jn_sb
                 do i=i0_sb,in_sb
                  w_tr(i,j,k) = fld_tr(i,j,k)
                 end do
               end do
             end do
!$omp enddo 
            else
!$omp do collapse(2)
             do k=k0,l_nk
               do j=j0_sb,jn_sb
                 do i=i0_sb,in_sb
                  w_tr(i,j,k) = 0.
                 end do
               end do
             end do
!$omp enddo 
            endif

            call mass_tr_hlt (tracer_8,w_tr,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0_sb,in_sb,j0_sb,jn_sb,k0_sb)

!$omp single
            if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,' SB= ', &
                                                                     tracer_8/Adz_gs_area_8,Tr3d_name_S(n)(1:4),F_comment_S
!$omp end single

         end if

         !Terminator: Print the mass of CLY scaled by area (CORE)
         !-------------------------------------------------------
!$omp single
         call canonical_terminator_2 (air_mass,tracer_8,count,l_minx,l_maxx,l_miny,l_maxy,l_nk,k0,Lun_out, &
                                      Ptopo_couleur,type_S,time_S,F_comment_S)
!$omp end single

      end do

      !-----------------
      !Treat Air density
      !-----------------
!$omp do collapse(2)
      do k=k0,l_nk
        do j=j0_c,jn_c
          do i=i0_c,in_c
            w_tr(i,j,k) = 1.
          end do
        end do
      end do
!$omp enddo 

      type_S = "Mass of Density (WET)"
      if (Schm_dry_mixing_ratio_L) type_S = "Mass of Density (DRY)"

      !Print Total Air mass scaled by area (CORE)
      !------------------------------------------
      call mass_tr_hlt (tracer_8,w_tr,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0_c,in_c,j0_c,jn_c,k0)

!$omp single
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,'  C= ', &
                                                               tracer_8/Adz_gc_area_8,'RHO ',F_comment_S

      !Same scaling as in STAT_PSADJ
      !-----------------------------
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,'  C= ', &
                                                               tracer_8/Adz_gc_area_8*grav_8,'RHOg',F_comment_S
!$omp end single

      !Print Total Air mass scaled by area (SUBSET of YIN)
      !---------------------------------------------------
      if (do_subset_GY_L) then


         if (Ptopo_couleur == 0) then 
!$omp do collapse(2)
          do k=k0_sb,l_nk
            do j=j0_sb,jn_sb
              do i=i0_sb,in_sb
               w_tr(i,j,k) = 1.
              end do
            end do
          end do
!$omp enddo 
         else
!$omp do collapse(2)
          do k=k0,l_nk
            do j=j0_sb,jn_sb
              do i=i0_sb,in_sb
               w_tr(i,j,k) = 0.
              end do
            end do
          end do
!$omp enddo 
         endif

         call mass_tr_hlt (tracer_8,w_tr,air_mass,l_minx,l_maxx,l_miny,l_maxy,l_nk,i0_sb,in_sb,j0_sb,jn_sb,k0_sb)

!$omp single
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,' SB= ', &
                                                                  tracer_8/Adz_gs_area_8,'RHO ',F_comment_S

         !Same scaling as in STAT_PSADJ
         !-----------------------------
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1002) 'TRACERS: ',type_S,time_S,' SB= ', &
                                                                  tracer_8/Adz_gs_area_8*grav_8,'RHOg',F_comment_S
!$omp end single

      end if

!$omp single
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ''
!$omp end single
!
!---------------------------------------------------------------------
!
      return

1002  format(1X,A9,A21,1X,A7,A5,E19.12,1X,A4,1X,A16)

      end subroutine stat_mass_tracers_hlt
