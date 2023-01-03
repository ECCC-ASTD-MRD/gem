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

!**s/r psadj_hlt - Adjust surface pressure for conservation in gem-H

      subroutine psadj_hlt (F_kount)
      use adz_mem
      use cstv
      use dyn_fisl_options
      use dynkernel_options
      use geomh
      use gmm_geof
      use gmm_vt0
      use gmm_pw
      use HORgrid_options
      use init_options
      use metric
      use psadjust
      use masshlt
      use ptopo
      use tdpack
      use tr3d
      use mem_tstp
      use omp_timing
      use omp_lib
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_kount

      include 'mpif.h'
      include 'rpn_comm.inc'
      logical :: LAM_L
      integer :: err,i,j,k,n,MAX_iteration,empty_i,comm,dim,dimH,dimV,dimV1
      real :: empty
      real, dimension(:,:  ), pointer :: delq
      real, dimension(:,:,:), pointer :: qt0i
      real(kind=REAL64), dimension(:,:  ), pointer :: fl_0_8
      real(kind=REAL64) :: l_avg_8(2),g_avg_ps_1_8,g_avg_ps_0_8,g_avg_fl_0_8,substract_8,sum_lcl(2)
      real(kind=REAL64) :: gathV(2,Ptopo_numproc*Ptopo_ncolors),gathS(Ptopo_numproc*Ptopo_ncolors)

!---------------------------------------------------------------------
!
      OMP_max_threads=OMP_get_max_threads()
      dimH= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      dimV= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*l_nk
      dimV1= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*(l_nk+1)
      delq  (l_minx:l_maxx,l_miny:l_maxy) => WS1(1:)           ; dim=dimH
      qt0i  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1) => WS1(dim+1:); dim=dim+dimV1
      p0_0_8(l_minx:l_maxx,l_miny:l_maxy) => WS1_8(    1:) ; dim=dimH
      p0_1_8(l_minx:l_maxx,l_miny:l_maxy) => WS1_8(dim+1:) ; dim=dim+dimH
      p0_dry_0_8(l_minx:l_maxx,l_miny:l_maxy) => WS1_8(dim+1:) ; dim=dim+dimH
      p0_dry_1_8(l_minx:l_maxx,l_miny:l_maxy) => WS1_8(dim+1:) ; dim=dim+dimH
      fl_0_8(l_minx:l_maxx,l_miny:l_maxy) => WS1_8(dim+1:) ; dim= 2*dim+dimH
      thread_sum(1:2,0:OMP_max_threads-1) => WS1_8(dim+1:) ; dim= dim+2*OMP_max_threads
      g_avg_8(1:2) => WS1_8(dim+1:) ; dim= dim+2

      
      !Update Pressure PW_MOINS in accordance with TIME M
      !--------------------------------------------------
      call pressure_hlt ( pw_pm_moins,pw_pt_moins,pw_p0_moins,pw_log_pm,pw_log_pt, &
                          pw_pm_moins_8,pw_p0_moins_8, &
                          l_minx,l_maxx,l_miny,l_maxy,l_nk,0 )

      LAM_L = .not.Grd_yinyang_L

      if ( Schm_psadj_print_L ) call stat_psadj_hlt(0,"AFTER DYNSTEP")

      if ( Schm_psadj == 0 ) return

      if ( Cstv_dt_8*F_kount <= Iau_period ) return

      if ( LAM_L .and. Init_mode_L ) return

      substract_8 = 1.0d0
      if (LAM_L) substract_8 = 0.0d0

      comm = RPN_COMM_comm ('MULTIGRID')

      !LAM: Estimate FLUX_out/FLUX_in and Evaluate water tracers at TIME M
      !-------------------------------------------------------------------
      if (LAM_L) then
         Adz_flux => Adz_flux_3CWP_PS

         !Estimate FLUX_out/FLUX_in using Tracer=1 based on Aranami et al. (2015)
         !-----------------------------------------------------------------------
         call gtmg_start (32, 'C_BCFLUX_PS', 10)
         call adz_BC_LAM_Aranami_hlt (empty,Adz_pb,Adz_num_b,1,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,empty_i,MAXTR3D+1)
         call gtmg_stop (32)
!$omp do
         do j=1,l_nj
            do i=1,l_ni
               fl_0_8(i,j) = 0.0d0
            end do
         end do
!$omp end do nowait


         !Obtain Wet Surface pressure minus Cstv_pref_8 at TIME P
         !-------------------------------------------------------
         if (Schm_psadj==1) then

           if (Adz_k0t==1) then
!$omp do
            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e
                  p0_1_8(i,j) = pw_p0_plus_8(i,j)
               end do
            end do
!$omp end do
           else
!$omp do
            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e
                  p0_1_8(i,j) = pw_pm_plus_8(i,j,l_nk+1) - pw_pm_plus_8(i,j,Adz_k0t)
               end do
            end do
!$omp end do
           endif

         endif !Schm_psadj==1

         !Compute Dry surface pressure at TIME P (Schm_psadj==2)
         !------------------------------------------------------
         if (Schm_psadj==2) then
            call dry_sfc_pressure_hlt_8(p0_dry_1_8,pw_pm_plus_8,pw_p0_plus_8,Adz_k0t,'P')
         !Obtain Wet Surface pressure minus Cstv_pref_8
         !---------------------------------------------
!$omp do
            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e
                  p0_1_8(i,j) = p0_dry_1_8(i,j) - substract_8*Cstv_pref_8
               end do
            end do
!$omp end do
         endif !Schm_psadj==2

         !Estimate air mass on CORE at TIME P
         !-----------------------------------
         sum_lcl= 0.d0
!$omp do
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               sum_lcl(1)= sum_lcl(1) +  p0_1_8(i,j) * geomh_area_mask_8(i,j)
            end do
         end do
!$omp end do nowait
         thread_sum(1,OMP_get_thread_num()) = sum_lcl(1)
!$OMP BARRIER

!$omp single
         l_avg_8(1) = sum(thread_sum(1,:))
         call MPI_Allgather(l_avg_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)
         g_avg_8(1) = sum(gathS)
!$omp end single

         g_avg_ps_1_8 = g_avg_8(1) * PSADJ_scale_8
         
      else

         g_avg_ps_1_8 = PSADJ_g_avg_ps_initial_8

      end if

      MAX_iteration = 10
      if (Schm_psadj==1) MAX_iteration = 1
     
      do n = 1, MAX_iteration  

        if (n==1) then
!$omp do collapse(2)
          do k=1,l_nk+1
            do j=1+pil_s,l_nj-pil_n
              do i=1+pil_w,l_ni-pil_e
                 qt0i(i,j,k) = qt0(i,j,k) 
              end do
            end do
          end do
!$omp end do
        endif



        !Obtain Wet Surface pressure minus Cstv_pref_8 at time M (Schm_psadj==1)
        !-----------------------------------------------------------------------
        if (Schm_psadj==1) then

          if (Adz_k0t==1) then
!$omp do
           do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               p0_0_8(i,j) = pw_p0_moins_8(i,j) - substract_8*Cstv_pref_8 
            end do
           end do
!$omp end do
          else
!$omp do
           do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               p0_0_8(i,j) = pw_pm_moins_8(i,j,l_nk+1) - pw_pm_moins_8(i,j,Adz_k0t) - substract_8*Cstv_pref_8 
            end do
           end do
!$omp end do
          endif

        endif !Schm_psadj==1 

        !Compute Dry surface pressure at TIME M (Schm_psadj==2)
        !------------------------------------------------------
        if (Schm_psadj==2) then
        call dry_sfc_pressure_hlt_8(p0_dry_0_8,pw_pm_moins_8,pw_p0_moins_8,Adz_k0t,'M')
!$omp do
          !Obtain Dry Surface pressure minus Cstv_pref_8
          !---------------------------------------------
          do j=1+pil_s,l_nj-pil_n
             do i=1+pil_w,l_ni-pil_e
               p0_0_8(i,j) = p0_dry_0_8(i,j) - substract_8*Cstv_pref_8
             end do
          end do
!$omp end do
        endif

        !Obtain FLUX on NEST+CORE at TIME M
        !----------------------------------
        if (LAM_L) then
!$omp do
         do j=Adz_j0b,Adz_jnb
            do k=1,l_nk
               do i=Adz_i0b,Adz_inb
                  fl_0_8(i,j) = fl_0_8(i,j) + (pw_pm_moins_8(i,j,k+1) - pw_pm_moins_8(i,j,k)) &
                                   * (Adz_flux(1)%fi(i,j,k) - Adz_flux(1)%fo(i,j,k))
               end do
            end do
         end do
!$omp end do
        end if

        sum_lcl(1)= 0.d0
        sum_lcl(2)= 0.d0
        g_avg_8(1)= 0.d0
        g_avg_8(2)= 0.d0
!$omp do
        do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            sum_lcl(1)= sum_lcl(1) +  p0_0_8(i,j) * geomh_area_mask_8(i,j)
         end do
        end do
!$omp end do nowait
        thread_sum(1,OMP_get_thread_num()) = sum_lcl(1)

        if (LAM_L) then
!$omp do
         do j=Adz_j0b,Adz_jnb
            do i=Adz_i0b,Adz_inb
               sum_lcl(2) = sum_lcl(2) + fl_0_8(i,j) * geomh_area_mask_8(i,j)
            end do
         end do
!$omp end do nowait
         thread_sum(2,OMP_get_thread_num()) = sum_lcl(2)
!$OMP BARRIER
!$omp single
         l_avg_8(1) = sum(thread_sum(1,:))
         l_avg_8(2) = sum(thread_sum(2,:))
         call MPI_Allgather(l_avg_8,2,MPI_DOUBLE_PRECISION,gathV,2,MPI_DOUBLE_PRECISION,comm,err)
         do i=1,2
            g_avg_8(i) = sum(gathV(i,:))
         end do
!$omp end single
        else
!$OMP BARRIER
!$omp single
         l_avg_8(1) = sum(thread_sum(1,:))
         call MPI_Allgather(l_avg_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)
         g_avg_8(1) = sum(gathS)
!$omp end single
        endif

        g_avg_ps_0_8 = g_avg_8(1) * PSADJ_scale_8
        g_avg_fl_0_8 = g_avg_8(2) * PSADJ_scale_8

        !Correct surface pressure in order to preserve air mass on CORE (taking FLUX mass into account in LAM)
        !-----------------------------------------------------------------------------------------------------
!$omp do
        do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            if (fis0(i,j) > 1.) then ! bad coding ???
               pw_p0_moins_8(i,j) = pw_p0_moins_8(i,j) + &
               ((g_avg_ps_1_8 - g_avg_ps_0_8) + g_avg_fl_0_8)*PSADJ_fact_8
            else
               pw_p0_moins_8(i,j) = pw_p0_moins_8(i,j) + &
               ((g_avg_ps_1_8 - g_avg_ps_0_8) + g_avg_fl_0_8)*Cstv_psadj_8
            end if
         end do
        end do
!$omp end do

        if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then
!$omp do
        do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            st0(i,j) = log(pw_p0_moins_8(i,j)/Cstv_pref_8)
         end do
        end do
!$omp end do
        else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then
!$omp do
        do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            qt0(i,j,l_nk+1) = rgasd_8*Cstv_Tstr_8 * &
                    (log(pw_p0_moins_8(i,j))-GVM%lg_pstar_8(i,j,l_nk+1))
         end do
        end do
!$omp end do

        !Adjust q at the other vertical levels at TIME M
        !-----------------------------------------------
!$omp do
        do j=1+pil_s,l_nj-pil_n
          do i=1+pil_w,l_ni-pil_e
            delq(i,j) = qt0(i,j,l_nk+1) - qt0i(i,j,l_nk+1)
          end do
        end do
!$omp end do

!$omp do collapse(2)
        do k=1,l_nk
          do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
                qt0(i,j,k) =  qt0i(i,j,k) + delq(i,j)
            end do
          end do
        end do
!$omp end do
        endif

        !Update Pressure PW_MOINS in accordance with TIME M
        !--------------------------------------------------
        call pressure_hlt ( pw_pm_moins,pw_pt_moins,pw_p0_moins,pw_log_pm,pw_log_pt, &
                            pw_pm_moins_8,pw_p0_moins_8, &
                            l_minx,l_maxx,l_miny,l_maxy,l_nk,0 )
      end do !END MAX_iteration

!$omp single
      if ( Schm_psadj_print_L.and.Lun_out>0.and.Ptopo_couleur==0 ) then
         write(Lun_out,*)    ''
         if (Schm_psadj==1) write(Lun_out,*) 'PSADJ: Conservation WET air'
         write(Lun_out,*)    ''
         write(Lun_out,*)    '------------------------------------------------------------------------------'
         write(Lun_out,1004) 'PSADJ STATS: N=',g_avg_ps_0_8,' O=',g_avg_ps_1_8,' F=',g_avg_fl_0_8
         write(Lun_out,*)    '------------------------------------------------------------------------------'

      end if
!$omp end single

      if ( Schm_psadj_print_L ) call stat_psadj_hlt(0,"AFTER PSADJ")

 1004 format(1X,A15,E19.12,A3,E19.12,A3,E19.12)
!
!---------------------------------------------------------------------
!
      return
      end subroutine psadj_hlt
