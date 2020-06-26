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

!**s/r wil_init - Prepare initial conditions (u,v,gz,s,topo) for Williamson cases

      subroutine wil_init (F_u,F_v,F_gz,F_s,F_topo,F_q,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

      use canonical
      use cstv
      use dynkernel_options
      use dyn_fisl_options
      use glb_ld
      use mem_tracers
      use step_options
      use tdpack
      use ver
      use wil_options

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in)  :: Mminx,Mmaxx,Mminy,Mmaxy,Nk
      logical, intent(in)  :: F_stag_L ! Staggered uv if .T. / Scalar uv if .F.
      real,    intent(out) :: F_u   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar if F_stag_L=.F.
                              F_v   (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar if F_stag_L=.F.
                              F_gz  (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !HEIGHT (m)
                              F_s   (Mminx:Mmaxx,Mminy:Mmaxy),    &
                              F_topo(Mminx:Mmaxx,Mminy:Mmaxy),    &
                              F_q   (Mminx:Mmaxx,Mminy:Mmaxy,Nk+1)

      !object
      !--------------------|-------------------------------------------------|
      ! NAME               | DESCRIPTION                                     |
      !--------------------|-------------------------------------------------|
      ! Williamson_case    |1=Advection 2D (see Williamson_NAIR/Terminator_L)|
      !                    |2=Steady-state nonlinear zonal geostrophic flow  |
      !                    |  Williamson et al.,1992,JCP,102,211-224         |
      !                    |5=Zonal flow over an isolated mountain           |
      !                    |  Williamson et al.,1992,JCP,102,211-224         |
      !                    |6=Rossby-Haurwitz wave                           |
      !                    |  Williamson et al.,1992,JCP,102,211-224         |
      !                    |7=The 21 December 1978 case                      |
      !                    |  Williamson et al.,1992,JCP,102,211-224         |
      !                    |8=Galewsky's barotropic wave                     |
      !                    |  Galewsky et al.,2004,Tellus,56A,429-440        |
      !                    |9=The Matsuno baroclinic wave test case          |
      !                    |  Shamir et al.,2019,GMD,12,2181-2193            |
      !--------------------|-------------------------------------------------|
      ! Williamson_NAIR    |Use when Williamson_case=1                       |
      !--------------------|-------------------------------------------------|
      !                    |=0=Solid body rotation of a cosine bell          |
      !                    |   Williamson et al.,1992,JCP,102,211-224        |
      !                    |=1=Deformational Non-divergent winds             |
      !                    |   Lauritzen et al.,2012,GMD,5,887-901           |
      !                    |=2=Deformational divergent winds                 |
      !                    |   Lauritzen et al.,2012,GMD,5,887-901           |
      !                    |=3=Deformational Flow for Circular vortex        |
      !                    |   Nair and Machenhauer,2002,MWR,130,649-667     |
      !----------------------------------------------------------------------|
      ! Wil._Terminator_L  |Do Terminator chemistry if T                     |
      !                    |Lauritzen et al.,2015,GMD,8,1299-1313            |
      !----------------------------------------------------------------------|

      integer :: istat,istat1,istat2,istat3,istat4,i,j,k
      real, dimension (:,:,:), pointer :: cl,cl2,q1,q2,q3,q4
      real, dimension (Mminx:Mmaxx,Mminy:Mmaxy) :: topo_case5
      real, parameter :: CLY_REF = 4.*10.**(-6)
!
!---------------------------------------------------------------------
!
      if (Williamson_case==7) return

      !Setup 2D Advection runs
      !-----------------------
      if (Williamson_case==1) then

         if (Williamson_Nair==0.or.Williamson_Nair==3) then

            !Initialize Q1
            !-------------
            istat = tr_get('Q1:P',q1)

            if (istat<=0.and.Williamson_Nair==0) call handle_error(-1,'WIL_INIT','Tracer Q1 missing')
            if (istat<=0.and.Williamson_Nair==3) call handle_error(-1,'WIL_INIT','Tracer Q1 missing')

            if (Williamson_Nair==0) call wil_case1(q1,Mminx,Mmaxx,Mminy,Mmaxy,Nk,0,Lctl_step)
            if (Williamson_Nair==3) call wil_case1(q1,Mminx,Mmaxx,Mminy,Mmaxy,Nk,5,Lctl_step)

            !Initialize Q1 REFERENCE
            !-----------------------
            q1ref(1:l_ni,1:l_nj,1:Nk) = q1(1:l_ni,1:l_nj,1:Nk)

            else if (Williamson_Nair==1.or.Williamson_Nair==2) then

            !Initialize Q1,Q2,Q3,Q4
            !----------------------
            istat1 = tr_get('Q1:P',q1)
            istat2 = tr_get('Q2:P',q2)
            istat3 = tr_get('Q3:P',q3)
            istat4 = tr_get('Q4:P',q4)

            if (istat1<=0.or.istat2<=0.or.istat3<=0.or.istat4<=0) call handle_error(-1,'WIL_INIT','Tracers Q1/Q2/Q3/Q4 missing')

            call wil_case1(q1,Mminx,Mmaxx,Mminy,Mmaxy,Nk,1,Lctl_step)
            call wil_case1(q2,Mminx,Mmaxx,Mminy,Mmaxy,Nk,2,Lctl_step)
            call wil_case1(q3,Mminx,Mmaxx,Mminy,Mmaxy,Nk,3,Lctl_step)
            call wil_case1(q4,Mminx,Mmaxx,Mminy,Mmaxy,Nk,4,Lctl_step)

            !Initialize Q1,Q2,Q3,Q4 REFERENCE
            !--------------------------------
            q1ref(1:l_ni,1:l_nj,1:Nk) = q1(1:l_ni,1:l_nj,1:Nk)
            q2ref(1:l_ni,1:l_nj,1:Nk) = q2(1:l_ni,1:l_nj,1:Nk)
            q3ref(1:l_ni,1:l_nj,1:Nk) = q3(1:l_ni,1:l_nj,1:Nk)
            q4ref(1:l_ni,1:l_nj,1:Nk) = q4(1:l_ni,1:l_nj,1:Nk)

         else

            call handle_error (-1,'INIT_BAR','Williamson_Nair not valid')

         end if

         !Setup Terminator's case: Initialize CL/CL2 tracers
         !--------------------------------------------------
         istat1 = tr_get('CL:P',cl)
         istat2 = tr_get('CL2:P',cl2)

         if (Williamson_Terminator_L.and.(istat1<=0.or.istat2<=0)) call handle_error (-1,'INIT_BAR','Tracers CL/CL2 missing')

         if (istat1>0.and.istat2>0) then

            call wil_Terminator_0 (cl,cl2,Mminx,Mmaxx,Mminy,Mmaxy,Nk)

            !Initialize CLY
            !--------------
            cly(1:l_ni,1:l_nj,1:Nk) = cl(1:l_ni,1:l_nj,1:Nk) + 2.0d0 * cl2(1:l_ni,1:l_nj,1:Nk)

            !Initialize CLY REFERENCE
            !------------------------
            clyref(1:l_ni,1:l_nj,1:Nk) =  CLY_REF

         end if

         call wil_uvcase1 (F_u,F_v,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L,0)

      end if

      !Setup Williamson Case 2: Steady state nonlinear geostrophic flow
      !----------------------------------------------------------------
      if (Williamson_case==2) then
          call wil_case2   (F_gz,   Mminx,Mmaxx,Mminy,Mmaxy,Nk)
          call wil_uvcase2 (F_u,F_v,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)
      end if

      !Setup Williamson Case 5: Zonal Flow over an isolated mountain
      !-------------------------------------------------------------
      if (Williamson_case==5) then
          call wil_case5   (F_gz,topo_case5,Mminx,Mmaxx,Mminy,Mmaxy,Nk)
          call wil_uvcase5 (F_u,F_v,        Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)
      end if

      !Setup Williamson Case 6: Rossby-Haurwitz wave
      !---------------------------------------------
      if (Williamson_case==6) then
          call wil_case6   (F_gz,   Mminx,Mmaxx,Mminy,Mmaxy,Nk)
          call wil_uvcase6 (F_u,F_v,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)
      end if

      !Setup Williamson Case 8 == Galewsky's Case: Barotropic wave
      !-----------------------------------------------------------
      if (Williamson_case==8) then
          call wil_case8   (F_gz,   Mminx,Mmaxx,Mminy,Mmaxy,Nk)
          call wil_uvcase8 (F_u,F_v,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)
      end if

      !Setup Williamson Case 9 == The Matsuno baroclinic wave test case
      !----------------------------------------------------------------
      if (Williamson_case==9) then
          call wil_case9   (F_gz,   Mminx,Mmaxx,Mminy,Mmaxy,Nk,Lctl_step)
          call wil_uvcase9 (F_u,F_v,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L,Lctl_step)
      end if

      !Initialize log(surface pressure) and Topo
      !-----------------------------------------
      F_s    = 0.
      F_topo = 0.

      if (Williamson_case==1) return

      if (Williamson_case==5) F_topo(1:l_ni,1:l_nj) = topo_case5(1:l_ni,1:l_nj)*grav_8

      !Initialize log(surface pressure)
      !--------------------------------
      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

         do j=1,l_nj
           do i=1,l_ni
              F_s(i,j) = (grav_8*F_gz(i,j,1)-F_topo(i,j)) &
                         /(Rgasd_8*Cstv_Tstr_8) &
                         +Ver_z_8%m(1)-Cstv_Zsrf_8
           end do
         end do

      !Initialize PHI perturbation in q
      !--------------------------------
      else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

        do k=1,G_nk+1
           do j=1,l_nj
              do i=1,l_ni
                 F_q(i,j,k) = grav_8*F_gz(i,j,1) - 1.0d0/Cstv_invFI_8
              end do
            end do
        end do

      end if
!
!---------------------------------------------------------------------
!
      return
      end
