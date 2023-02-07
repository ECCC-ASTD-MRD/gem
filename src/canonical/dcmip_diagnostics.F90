!--------------------------------- LICENCE BEGIN -------------------------------
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

!**s/r dcmip_diagnostics - Evaluate Error norms L1,L2,LMASS,L_inf for DCMIP's cases

      subroutine dcmip_diagnostics (F_my_step)

      use adz_options
      use canonical
      use cstv
      use dcmip_options
      use geomh
      use glb_ld
      use gmm_vt1
      use HORgrid_options
      use lun
      use mem_tracers
      use ptopo
      use step_options
      use tr3d

      implicit none

      !arguments
      !---------
      integer, intent(in) :: F_my_step

      !object
      !==============================================================================
      !   CASE  11: Evaluate Error norms L1,L2,LMASS,L_inf for Q1,Q2,Q3,Q4 at final T
      !   CASE  12: Evaluate Error norms L1,L2,LMASS,L_inf for Q1          at final T
      !   CASE  13: Evaluate Error norms L1,L2,LMASS,L_inf for Q1,Q2,Q3,Q4 at final T
      !   CASE 161: Evaluate Error norms L1,L2,LMASS,L_inf for CLY         at each  T
      !   CASE 163: Maximum vertical velocity                              at each  T
      !==============================================================================

      integer :: istat,i,j,k,ierr,n
      real, pointer, dimension(:,:,:) :: cl,cl2,tr,tr_r

      real, parameter :: CLY_REF = 4.*10.**(-6)

      real(kind=REAL64) :: norm_1_8,norm_2_8,norm_m_8,norm_inf_8,abs_8,w1_8,w2_8, &
                           s_err_m_8,s_ref_m_8,s_err_inf_8,s_ref_inf_8,s_max_8,   &
                           g_err_m_8,g_ref_m_8,g_err_inf_8,g_ref_inf_8,g_max_8,   &
                           s_err_1_8,s_ref_1_8,s_err_2_8,s_ref_2_8,               &
                           g_err_1_8,g_ref_1_8,g_err_2_8,g_ref_2_8

      character(len= 12) :: name_S
      character(len= 9)  :: communicate_S

!      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk) :: air_mass

      logical :: almost_zero
!
!---------------------------------------------------------------------
!
      if (Adz_verbose == 0) return

      if (.NOT.(Dcmip_case>=11.and.Dcmip_case<=13).and.Dcmip_case/=161.and.Dcmip_case/=163) return

      if ((Dcmip_case>=11.and.Dcmip_case<=13).and.F_my_step<Step_total) return

      !Reset Air Mass at TIME P
      !------------------------
      call get_air_mass (air_mass,1,l_minx,l_maxx,l_miny,l_maxy,l_nk,1)

      !--------------------------------------------------------------------------
      !CASE 11: Pure advection - 3D deformational flow
      !CASE 12: Pure advection - 3D Hadley-like meridional circulation
      !CASE 13: Pure advection - 2D solid-body rotation of thin cloud-like tracer
      !--------------------------------------------------------------------------
      !Evaluate Error norms L1,L2,L_inf for ALL TRACERS
      !------------------------------------------------
      if (Dcmip_case>=11.and.Dcmip_case<=13) then

         if (Dcmip_case==11.and.Lun_out>0.and.Ptopo_couleur==0) write (Lun_out,1003)
         if (Dcmip_case==12.and.Lun_out>0.and.Ptopo_couleur==0) write (Lun_out,1004)
         if (Dcmip_case==13.and.Lun_out>0.and.Ptopo_couleur==0) write (Lun_out,1005)

         !Error norms: L1,L2,L_inf for ALL TRACERS
         !----------------------------------------
         do n=1,Tr3d_ntr

            if (.NOT.((Tr3d_name_S(n)(1:2)=='Q1').or. &
                      (Tr3d_name_S(n)(1:2)=='Q2').or. &
                      (Tr3d_name_S(n)(1:2)=='Q3').or. &
                      (Tr3d_name_S(n)(1:2)=='Q4'))) cycle

            tr => tracers_P(n)%pntr

            if (Tr3d_name_S(n)(1:2)=='Q1') tr_r => q1ref
            if (Tr3d_name_S(n)(1:2)=='Q2') tr_r => q2ref
            if (Tr3d_name_S(n)(1:2)=='Q3') tr_r => q3ref
            if (Tr3d_name_S(n)(1:2)=='Q4') tr_r => q4ref

            s_err_1_8 = 0.; s_err_2_8 = 0.; s_err_m_8 = 0.; s_err_inf_8 = 0.
            s_ref_1_8 = 0.; s_ref_2_8 = 0.; s_ref_m_8 = 0.; s_ref_inf_8 = 0.

            do k = 1,G_nk

               do j = 1+pil_s,l_nj-pil_n

                  do i = 1+pil_w,l_ni-pil_e

                     s_err_1_8 = s_err_1_8 + abs(tr(i,j,k) - tr_r(i,j,k))    * air_mass(i,j,k) * geomh_mask_8(i,j)
                     s_ref_1_8 = s_ref_1_8 + abs(            tr_r(i,j,k))    * air_mass(i,j,k) * geomh_mask_8(i,j)

                     s_err_2_8 = s_err_2_8 +    (tr(i,j,k) - tr_r(i,j,k))**2 * air_mass(i,j,k) * geomh_mask_8(i,j)
                     s_ref_2_8 = s_ref_2_8 +    (            tr_r(i,j,k))**2 * air_mass(i,j,k) * geomh_mask_8(i,j)

                     s_err_m_8 = s_err_m_8 +    (tr(i,j,k) - tr_r(i,j,k))    * air_mass(i,j,k) * geomh_mask_8(i,j)
                     s_ref_m_8 = s_ref_m_8 +    (            tr_r(i,j,k))    * air_mass(i,j,k) * geomh_mask_8(i,j)

                     w1_8 = abs( tr(i,j,k) - tr_r(i,j,k) )
                     w2_8 = abs(             tr_r(i,j,k) )

                     s_err_inf_8 = max( s_err_inf_8, w1_8 )
                     s_ref_inf_8 = max( s_ref_inf_8, w2_8 )

                  end do

               end do

            end do

            communicate_S = "GRID"
            if (Grd_yinyang_L) communicate_S = "MULTIGRID"

            call RPN_COMM_allreduce(s_err_1_8,g_err_1_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
            call RPN_COMM_allreduce(s_ref_1_8,g_ref_1_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
            call RPN_COMM_allreduce(s_err_2_8,g_err_2_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
            call RPN_COMM_allreduce(s_ref_2_8,g_ref_2_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
            call RPN_COMM_allreduce(s_err_m_8,g_err_m_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
            call RPN_COMM_allreduce(s_ref_m_8,g_ref_m_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)

            call RPN_COMM_allreduce(s_err_inf_8,g_err_inf_8,1,"MPI_double_precision","MPI_MAX",communicate_S,ierr)
            call RPN_COMM_allreduce(s_ref_inf_8,g_ref_inf_8,1,"MPI_double_precision","MPI_MAX",communicate_S,ierr)

            !Evaluate Norms
            !--------------
            norm_1_8 = 10.**8; norm_2_8 = 10.**8; norm_m_8 = 10.**8; norm_inf_8 = 10.**8

            if (.not.almost_zero(g_ref_1_8)  ) norm_1_8   =     (g_err_1_8/g_ref_1_8)
            if (.not.almost_zero(g_ref_2_8)  ) norm_2_8   = sqrt(g_err_2_8/g_ref_2_8)
            if (.not.almost_zero(g_ref_m_8)  ) norm_m_8   =     (g_err_m_8/g_ref_m_8)
            if (.not.almost_zero(g_ref_inf_8)) norm_inf_8 = g_err_inf_8/g_ref_inf_8

            !Print Norms
            !-----------
            if (Lun_out>0.and.Ptopo_couleur==0) then
                if (Tr3d_name_S(n)(1:2)=='Q1') name_S = 'TRACER Q1  '
                if (Tr3d_name_S(n)(1:2)=='Q2') name_S = 'TRACER Q2  '
                if (Tr3d_name_S(n)(1:2)=='Q3') name_S = 'TRACER Q3  '
                if (Tr3d_name_S(n)(1:2)=='Q4') name_S = 'TRACER Q4  '
                write (Lun_out,*) ' +++++++++++++++++++++++++++++++++++++++++++'
                write (Lun_out,1001) name_S,'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM L1    = ',norm_1_8
                write (Lun_out,1001) name_S,'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM L2    = ',norm_2_8
                write (Lun_out,1001) name_S,'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM LMASS = ',norm_m_8
                write (Lun_out,1001) name_S,'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM L_inf = ',norm_inf_8
                write (Lun_out,*) ' +++++++++++++++++++++++++++++++++++++++++++'
                write (Lun_out,*) ' '
            end if

         end do

      end if

      !----------------------------------------------------------
      !CASE 161: Baroclinic wave with Toy Terminal Chemistry
      !----------------------------------------------------------
      !Evaluate Error norms L2,L_inf,relative Mass change for CLY
      !----------------------------------------------------------
      if (Dcmip_case==161) then

         if (Lun_out>0) write (Lun_out,1000)

         istat = tr_get('CL:P',cl)
         istat = tr_get('CL2:P',cl2)

         cly(1:l_ni,1:l_nj,1:G_nk) = cl(1:l_ni,1:l_nj,1:G_nk) + 2.0d0 * cl2(1:l_ni,1:l_nj,1:G_nk)

         s_err_1_8 = 0.; s_err_2_8 = 0.; s_err_m_8 = 0.; s_err_inf_8 = 0.
         s_ref_1_8 = 0.; s_ref_2_8 = 0.; s_ref_m_8 = 0.; s_ref_inf_8 = 0.

         do k=1,G_nk
            do j = 1+pil_s,l_nj-pil_n
            do i = 1+pil_w,l_ni-pil_e

               s_err_1_8 = s_err_1_8 + abs(cly(i,j,k) - cly_REF)     * air_mass(i,j,k) * geomh_mask_8(i,j)
               s_ref_1_8 = s_ref_1_8 + abs(             cly_REF)     * air_mass(i,j,k) * geomh_mask_8(i,j)

               s_err_2_8 = s_err_2_8 +     (cly(i,j,k) - cly_REF)**2 * air_mass(i,j,k) * geomh_mask_8(i,j)
               s_ref_2_8 = s_ref_2_8 +     (             cly_REF)**2 * air_mass(i,j,k) * geomh_mask_8(i,j)

               s_err_m_8 = s_err_m_8 +     (cly(i,j,k) - cly_REF)    * air_mass(i,j,k) * geomh_mask_8(i,j)
               s_ref_m_8 = s_ref_m_8 +     (             cly_REF)    * air_mass(i,j,k) * geomh_mask_8(i,j)

               w1_8 = abs(cly(i,j,k) - cly_REF)
               w2_8 = abs(             cly_REF)

               s_err_inf_8 = max(s_err_inf_8, w1_8)
               s_ref_inf_8 = max(s_ref_inf_8, w2_8)

            end do
            end do
         end do

         communicate_S = "GRID"
         if (Grd_yinyang_L) communicate_S = "MULTIGRID"

         call RPN_COMM_allreduce(s_err_1_8,g_err_1_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
         call RPN_COMM_allreduce(s_ref_1_8,g_ref_1_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
         call RPN_COMM_allreduce(s_err_2_8,g_err_2_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
         call RPN_COMM_allreduce(s_ref_2_8,g_ref_2_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
         call RPN_COMM_allreduce(s_err_m_8,g_err_m_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)
         call RPN_COMM_allreduce(s_ref_m_8,g_ref_m_8,1,"MPI_double_precision","MPI_SUM",communicate_S,ierr)

         call RPN_COMM_allreduce(s_err_inf_8,g_err_inf_8,1,"MPI_double_precision","MPI_MAX",communicate_S,ierr)
         call RPN_COMM_allreduce(s_ref_inf_8,g_ref_inf_8,1,"MPI_double_precision","MPI_MAX",communicate_S,ierr)

         !Evaluate Norms
         !--------------
         norm_1_8 = 10.**8; norm_2_8 = 10.**8; norm_m_8 = 10.**8; norm_inf_8 = 10.**8

         if (.not.almost_zero(g_ref_1_8)  ) norm_1_8   =     (g_err_1_8  /g_ref_1_8)
         if (.not.almost_zero(g_ref_2_8)  ) norm_2_8   = sqrt(g_err_2_8  /g_ref_2_8)
         if (.not.almost_zero(g_ref_m_8)  ) norm_m_8   =     (g_err_m_8  /g_ref_m_8)
         if (.not.almost_zero(g_ref_inf_8)) norm_inf_8 =      g_err_inf_8/g_ref_inf_8

         !Print Norms
         !-----------
         if (Lun_out>0.and.Ptopo_couleur==0) then
            write (Lun_out,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            write (Lun_out,1001) "TRACER CLY ",'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM L1    = ',norm_1_8
            write (Lun_out,1001) "TRACER CLY ",'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM L2    = ',norm_2_8
            write (Lun_out,1001) "TRACER CLY ",'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM LMASS = ',norm_m_8
            write (Lun_out,1001) "TRACER CLY ",'TIME (days) = ',(F_my_step*Cstv_dt_8)/3600./24.,' ERROR NORM L_inf = ',norm_inf_8
            write (Lun_out,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            write (Lun_out,*) ' '
         end if

      end if

      !----------------------------------
      !CASE 163: Supercell (Small Planet)
      !----------------------------------
      !Evaluate Maximum vertical velocity
      !----------------------------------
      if (Dcmip_case==163) then

         if (Lun_out>0) write (Lun_out,1002)

         s_max_8 = 0.0D0

         do k = 1,G_nk
            do j = 1+pil_s,l_nj-pil_n
               do i = 1+pil_w,l_ni-pil_e

                  abs_8 = abs(wt1(i,j,k))

                  s_max_8 = max(s_max_8, abs_8)

               end do
            end do
         end do

         communicate_S = "GRID"
         if (Grd_yinyang_L) communicate_S = "MULTIGRID"

         call RPN_COMM_allreduce(s_max_8,g_max_8,1,"MPI_double_precision","MPI_MAX",communicate_S,ierr)

         !Print Norms
         !-----------
         if (Lun_out>0.and.Ptopo_couleur==0) then
            write (Lun_out,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            write (Lun_out,1001) "FIELD WT1  ",'TIME (days) = ',(F_my_step*Cstv_dt_8*Dcmip_X)/3600./24., &
                                 ' MAX.VERT.VELOCITY= ',g_max_8
            write (Lun_out,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            write (Lun_out,*) ' '
         end if

      end if
!
!---------------------------------------------------------------------
!
      return

 1000 format( &
      /,'EVALUATE ERROR NORMS CLY (Moist Baroclinic wave with Toy Chemistry): (S/R DCMIP_DIAGNOSTICS)',   &
      /,'============================================================================================',/,/)
 1001 format(A12,1X,A14,F10.4,A20,E19.12)
 1002 format( &
      /,'EVALUATE MAXIMUM VELOCITY WT1 (Supercell): (S/R DCMIP_DIAGNOSTICS)',   &
      /,'==================================================================',/,/)
 1003 format( &
      /,'EVALUATE ERROR NORMS Q1-Q2-Q3-Q4 (Pure advection - 3D deformational flow): (S/R DCMIP_DIAGNOSTICS)',   &
      /,'==================================================================================================',/,/)
 1004 format( &
      /,'EVALUATE ERROR NORMS Q1 (Pure advection - 3D Hadley-like meridional circulation): (S/R DCMIP_DIAGNOSTICS)',   &
      /,'=========================================================================================================',/,/)
 1005 format( &
      /,'EVALUATE ERROR NORMS Q1-Q2-Q3-Q4 (Pure advection - 2D solid-body rot. of thin cloud-like tracer): ',   &
        '(S/R DCMIP_DIAGNOSTICS)',   &
      /,'==================================================================================================',   &
        '=======================',   &
      /,/)

      end subroutine dcmip_diagnostics
