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

!**s/r dcmip_init - Prepare initial conditions (u,v,w,t,zd,s,q,topo) for DCMIP 2012/2016 runs

      subroutine dcmip_init (F_u, F_v, F_w, F_t, F_zd, F_s, F_q, F_topo, &
                             Mminx,Mmaxx,Mminy,Mmaxy,Nk, &
                             F_stag_L,F_trprefix_S,F_trsuffix_S )
      use canonical
      use dcmip_options
      use gem_options
      use inp_mod
      use gmm_geof
      use gmm_pw
      use dyn_fisl_options
      use glb_ld
      use tr3d
      use gmm_itf_mod
      use tdpack, only : rgasd_8, rgasv_8

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      character(len=*) F_trprefix_S, F_trsuffix_S
      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk
      real F_u (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar if F_stag_L=.F.
           F_v (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Scalar if F_stag_L=.F.
           F_w (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_t (Mminx:Mmaxx,Mminy:Mmaxy,Nk), & !Real T if F_stag_L=.F.
           F_zd(Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_q (Mminx:Mmaxx,Mminy:Mmaxy,Nk+1), &
           F_topo(Mminx:Mmaxx,Mminy:Mmaxy)

      logical F_stag_L ! Staggered uv Virtual Temperature if .T.
                       !    Scalar uv    Real Temperature if .F.

      !object
      !======================================================================|
      !Prepare initial conditions for DCMIP 2012/2016 runs                   |
      !----------------------------------------------------------------------|
      !  case DCMIP 2012   | Pure advection                                  |
      !                    | ------------------------------------------------|
      !                    | 11: 3D deformational flow                       |
      !                    | 12: 3D Hadley-like meridional circulation       |
      !                    | 13: 2D solid-body rotation of thin cloud-like   |
      !                    |     tracer in the presence of orography         |
      !                    | ------------------------------------------------|
      !                    | 20: Steady-state at rest in presence of oro.    |
      !                    | ------------------------------------------------|
      !                    | Gravity waves, Non-rotating small-planet        |
      !                    | ------------------------------------------------|
      !                    | 21: Mountain waves over a Schaer-type mountain  |
      !                    | 22: As 21 but with wind shear                   |
      !                    | 31: Gravity wave along the equator              |
      !                    | ------------------------------------------------|
      !                    | Rotating planet: Hydro. to non-hydro. scales (X)|
      !                    | ------------------------------------------------|
      !                    | 41X: Dry Baroclinic Instability Small Planet    |
      !                    | ------------------------------------------------|
      !                    | 43 : Moist Baroclinic Instability Simple physics|
      !                    | ------------------------------------------------|
      ! case DCMIP 2016    | 161: Baroclinic wave with Toy Terminal Chemistry|
      !                    | 162: Tropical cyclone                           |
      !                    | 163: Supercell (Small Planet)                   |
      !--------------------|-------------------------------------------------|
      !DCMIP_2012: https://www.earthsystemcog.org/projects/dcmip-2012/       |
      !DCMIP_2016: https://www.earthsystemcog.org/projects/dcmip-2016/       |
      !======================================================================|

      !--------------------------------------------------------------------------

      integer istat,istat1,istat2,istat3,istat4, &
              Deep,Pertt,Pert,Moist,Shear,Tracers,i,j,k
      real, pointer, dimension(:,:,:) :: cl,cl2,qv,qc,qr,q1,q2,q3,q4
      real, dimension(1,1,1) :: empty
      real  s_u(Mminx:Mmaxx,Mminy:Mmaxy,Nk),s_v(Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
             tv(Mminx:Mmaxx,Mminy:Mmaxy,Nk), tt(Mminx:Mmaxx,Mminy:Mmaxy,Nk)

      real(8), parameter :: Rd   = Rgasd_8, & ! cte gaz - air sec   [J kg-1 K-1]
                            Rv   = Rgasv_8    ! cte gaz - vap eau   [J kg-1 K-1]

      real(8)  zvir

      !--------------------------------------------------------------------------

      if (Schm_sleve_L ) call handle_error (-1,'DCMIP_init','  SLEVE not available YET  ')

      !Prescribed d(Zeta)dot and dz/dt
      !-------------------------------
      Inp_zd_L = .true.
      Inp_w_L  = .true.

      !Obtain specific humidity
      !------------------------
      istat = gmm_get (trim(F_trprefix_S)//'HU'//trim(F_trsuffix_S),qv)

      !Initialization QC/QR for Precipitation
      !--------------------------------------
      if (Dcmip_prec_type/=-1) then

         istat1 = gmm_get (trim(F_trprefix_S)//'QC'//trim(F_trsuffix_S),qc)
         istat2 = gmm_get (trim(F_trprefix_S)//'RW'//trim(F_trsuffix_S),qr)

         if (istat1/=0.or.istat2/=0) call handle_error(-1,'DCMIP_INIT','Tracers QC/RW required when Precipitation')

         qc = 0. !ZERO Cloud water mixing ratio
         qr = 0. !ZERO Rain  water mixing ratio

         istat = gmm_get(gmmk_art_s, art)
         istat = gmm_get(gmmk_wrt_s, wrt)

         art = 0. !ZERO Averaged precipitation rate
         wrt = 0. !ZERO Averaged precipitation rate (WORK FIELD)

      end if

      !DCMIP 2016: Baroclinic wave with Toy Terminal Chemistry
      !-------------------------------------------------------
      if (Dcmip_case==161) then

          istat1= gmm_get (trim(F_trprefix_S)//'CL'//trim(F_trsuffix_S), cl )
          istat2= gmm_get (trim(F_trprefix_S)//'CL2'//trim(F_trsuffix_S),cl2)

          if (istat1/=0.or.istat2/=0) call handle_error(-1,'DCMIP_INIT','Tracers CL/CL2 required when Chemistry')

          !--------------------------------------------------------------------------
          Deep  = 0           !Deep atmosphere (no=0)
          Pertt = 0           !Type of perturbation (exponential=0/stream function=1)
          Moist = Dcmip_moist !Moist=1/Dry=0 Initial conditions
          !--------------------------------------------------------------------------

          call dcmip_baroclinic_wave_2016 (F_u,F_v,F_w,tv,F_zd,F_s,F_topo,qv,cl,cl2, &
                                           Mminx,Mmaxx,Mminy,Mmaxy,Nk,Deep,Moist,Pertt,Dcmip_X,F_stag_L)

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure

      !DCMIP 2016: Tropical cyclone
      !----------------------------
      else if (Dcmip_case==162) then

          call dcmip_tropical_cyclone (F_u,F_v,F_w,tv,F_zd,F_s,F_topo,qv,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure

      !DCMIP 2016: Supercell (Small planet)
      !------------------------------------
      else if (Dcmip_case==163) then

          istat = gmm_get(gmmk_thbase_s,thbase)

          !---------------------------------------------------------
          Pert = 1 !Thermal perturbation included (0 = no / 1 = yes)
          !---------------------------------------------------------

          call dcmip_supercell (F_u,F_v,F_w,tv,F_zd,F_s,F_topo,qv,Pert,thbase,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure

          !Initialize u,v,zd,w,tv,qv,qc,rw,theta REFERENCE for Vertical diffusion
          !----------------------------------------------------------------------
          istat = gmm_get(gmmk_uref_s , uref )
          istat = gmm_get(gmmk_vref_s , vref )
          istat = gmm_get(gmmk_wref_s , wref )
          istat = gmm_get(gmmk_zdref_s, zdref)
          istat = gmm_get(gmmk_qvref_s, qvref)
          istat = gmm_get(gmmk_qcref_s, qcref)
          istat = gmm_get(gmmk_qrref_s, qrref)
          istat = gmm_get(gmmk_thref_s, thref)

          !Prepare UREF/VREF (on staggered grids) for DCMIP_VRD
          !----------------------------------------------------
          if (.not.F_stag_L) then

          call hwnd_stag (s_u,s_v,F_u,F_v,Mminx,Mmaxx,Mminy,Mmaxy,Nk,.true.)

           uref(1:l_ni-1,1:l_nj,  1:Nk) =    s_u(1:l_ni-1,1:l_nj,  1:Nk)
           vref(1:l_ni,  1:l_nj-1,1:Nk) =    s_v(1:l_ni,  1:l_nj-1,1:Nk)

          else

             uref(1:l_ni-1,1:l_nj,  1:Nk) = F_u(1:l_ni-1,1:l_nj,  1:Nk)
             vref(1:l_ni,  1:l_nj-1,1:Nk) = F_v(1:l_ni,  1:l_nj-1,1:Nk)

          end if

           wref(1:l_ni,  1:l_nj,  1:Nk) =   F_w (1:l_ni,  1:l_nj,  1:Nk)
          zdref(1:l_ni,  1:l_nj,  1:Nk) =   F_zd(1:l_ni,  1:l_nj,  1:Nk)
          qvref(1:l_ni,  1:l_nj,  1:Nk) =     qv(1:l_ni,  1:l_nj,  1:Nk)
          qcref(1:l_ni,  1:l_nj,  1:Nk) =     qc(1:l_ni,  1:l_nj,  1:Nk)
          qrref(1:l_ni,  1:l_nj,  1:Nk) =     qr(1:l_ni,  1:l_nj,  1:Nk)
          thref(1:l_ni,  1:l_nj,  1:Nk) = thbase(1:l_ni,  1:l_nj,  1:Nk)

      !DCMIP 2012: Steady-State Atmosphere at Rest in the Presence of Orography
      !------------------------------------------------------------------------
      else if (Dcmip_case==20) then

          !Set initial conditions according to prescribed mountain
          !-------------------------------------------------------
          call dcmip_steady_state_mountain (F_u,F_v,F_zd,tv,qv,F_topo,F_s,Mminx,Mmaxx,Mminy,Mmaxy,Nk,.true.,F_stag_L)

          if (Vtopo_L) then

              istat = gmm_get(gmmk_topo_low_s , topo_low )
              istat = gmm_get(gmmk_topo_high_s, topo_high)

              topo_low (1:l_ni,1:l_nj) = 0.
              topo_high(1:l_ni,1:l_nj) = F_topo (1:l_ni,1:l_nj)
              F_topo   (1:l_ni,1:l_nj) = 0.

              !Reset initial conditions according to topo_low
              !----------------------------------------------
              call dcmip_steady_state_mountain (F_u,F_v,F_zd,tv,qv,F_topo,F_s,Mminx,Mmaxx,Mminy,Mmaxy,Nk,.false.,F_stag_L)

          end if

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure
          F_w = 0. !ZERO Dz/Dt

      !DCMIP 2012: Pure 3D Advection
      !-----------------------------
      else if (Dcmip_case>=11.and.Dcmip_case<=13) then

          !Get tracers Q1,Q2,Q3,Q4
          !-----------------------
          istat1 = gmm_get(trim(F_trprefix_S)//'Q1'//trim(F_trsuffix_S),q1)
          istat2 = gmm_get(trim(F_trprefix_S)//'Q2'//trim(F_trsuffix_S),q2)
          istat3 = gmm_get(trim(F_trprefix_S)//'Q3'//trim(F_trsuffix_S),q3)
          istat4 = gmm_get(trim(F_trprefix_S)//'Q4'//trim(F_trsuffix_S),q4)

          if ((istat1/=0.or.istat2/=0.or.istat3/=0.or.istat4/=0).and.Dcmip_case==11) goto 999
          if ((istat1/=0)                                       .and.Dcmip_case==12) goto 999
          if ((istat1/=0.or.istat2/=0.or.istat3/=0.or.istat4/=0).and.Dcmip_case==13) goto 999

          !3D deformational flow
          !---------------------
          if (Dcmip_case==11) call dcmip_tracers11_transport (F_u,F_v,F_zd,tv,qv,F_topo,F_s,q1,q2,q3,q4, &
                                                              Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

          !3D Hadley-like meridional circulation
          !-------------------------------------
          if (Dcmip_case==12) call dcmip_tracers12_transport (F_u,F_v,F_zd,tv,qv,F_topo,F_s,q1,&
                                                              Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

          !2D solid-body rotation of thin cloud-like tracer in the presence of orography
          !-----------------------------------------------------------------------------
          if (Dcmip_case==13) call dcmip_tracers13_transport (F_u,F_v,F_zd,tv,qv,F_topo,F_s,q1,q2,q3,q4, &
                                                              Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

          !Store REFERENCE at initial time
          !-------------------------------
          istat = gmm_get(gmmk_q1ref_s,q1ref)
          istat = gmm_get(gmmk_q2ref_s,q2ref)
          istat = gmm_get(gmmk_q3ref_s,q3ref)
          istat = gmm_get(gmmk_q4ref_s,q4ref)

          if (Dcmip_case>  0) q1ref(1:l_ni,1:l_nj,1:Nk) = q1(1:l_ni,1:l_nj,1:Nk)
          if (Dcmip_case/=12) q2ref(1:l_ni,1:l_nj,1:Nk) = q2(1:l_ni,1:l_nj,1:Nk)
          if (Dcmip_case/=12) q3ref(1:l_ni,1:l_nj,1:Nk) = q3(1:l_ni,1:l_nj,1:Nk)
          if (Dcmip_case/=12) q4ref(1:l_ni,1:l_nj,1:Nk) = q4(1:l_ni,1:l_nj,1:Nk)

      !DCMIP 2012: Mountain waves over a Schaer-type mountain on a small planet
      !------------------------------------------------------------------------
      else if (Dcmip_case==21.or.Dcmip_case==22) then

          if (Dcmip_case==21) Shear = 0 !Without wind shear
          if (Dcmip_case==22) Shear = 1 !With    wind shear

          call dcmip_Schaer_mountain (F_u,F_v,F_zd,tv,qv,F_topo,F_s,Mminx,Mmaxx,Mminy,Mmaxy,Nk,Shear,.true.,F_stag_L)

          if (Vtopo_L) then

              istat = gmm_get(gmmk_topo_low_s , topo_low )
              istat = gmm_get(gmmk_topo_high_s, topo_high)

              topo_low (1:l_ni,1:l_nj) = 0.
              topo_high(1:l_ni,1:l_nj) = F_topo (1:l_ni,1:l_nj)
              F_topo   (1:l_ni,1:l_nj) = 0.

              !Reset initial conditions according to topo_low
              !----------------------------------------------
              call dcmip_Schaer_mountain (F_u,F_v,F_zd,tv,qv,F_topo,F_s,Mminx,Mmaxx,Mminy,Mmaxy,Nk,Shear,.false.,F_stag_L)

          end if

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure
          F_w = 0. !ZERO Dz/Dt

      !DCMIP 2012: Gravity wave on a small planet along the equator
      !------------------------------------------------------------
      else if (Dcmip_case==31) then

          istat = gmm_get(gmmk_thbase_s,thbase)

          call dcmip_gravity_wave (F_u,F_v,F_zd,tv,qv,F_topo,F_s,thbase,Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_stag_L)

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure
          F_w = 0. !ZERO Dz/Dt

      !DCMIP 2012: Dry Baroclinic Instability on a Small Planet with dynamic tracers
      !------------------------------------------------------------------------------------
      !            Each Dcmip_case=41X has his own Dcmip_X
      !            Dynamical Tracers: Potential temperature and Ertel's potential vorticity
      !------------------------------------------------------------------------------------
      else if (Dcmip_case==410.or. &
              Dcmip_case==411.or. &
              Dcmip_case==412.or. &
              Dcmip_case==413 ) then

          !-------------------------------------------------------
          Moist   = Dcmip_moist ! Moist=1/Dry=0 Initial conditions
          Tracers = 1           ! Tracers=1/No Tracers=0
          !-------------------------------------------------------

          !Dynamical Tracers: Potential temperature and Ertel's potential vorticity
          !------------------------------------------------------------------------
          istat1 = gmm_get(trim(F_trprefix_S)//'Q1'//trim(F_trsuffix_S),q1)
          istat2 = gmm_get(trim(F_trprefix_S)//'Q2'//trim(F_trsuffix_S),q2)

          if ((istat1/=0.or.istat2/=0).and.Tracers==1) call handle_error(-1,'DCMIP_INIT','Tracers Q1/Q2 required when Dcmip_case=41X')

          call dcmip_baroclinic_wave_2012 (F_u,F_v,F_zd,tv,qv,F_topo,F_s,q1,q2, &
                                           Mminx,Mmaxx,Mminy,Mmaxy,Nk,Moist,Dcmip_X,Tracers,F_stag_L)

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure
          F_w = 0. !ZERO Dz/Dt

      !DCMIP 2012: Moist Baroclinic Instability driven by Simple Physics
      !-----------------------------------------------------------------
      else if (Dcmip_case==43) then

          !-------------------------------------------------------
          Moist   = Dcmip_moist ! Moist=1/Dry=0 Initial conditions
          Tracers = 0           ! Tracers=1/No Tracers=0
          !-------------------------------------------------------

          call dcmip_baroclinic_wave_2012 (F_u,F_v,F_zd,tv,qv,F_topo,F_s,empty,empty, &
                                           Mminx,Mmaxx,Mminy,Mmaxy,Nk,Moist,Dcmip_X,Tracers,F_stag_L)

          F_q = 0. !ZERO log of non-hydrostatic perturbation pressure
          F_w = 0. !ZERO Dz/Dt

      else

          call handle_error(-1,'dcmip_init','DCMIP_CASE 2012/2016 not available')

      end if

      call rpn_comm_xch_halo (F_topo,Mminx,Mmaxx,Mminy,Mmaxy,l_ni,l_nj,1, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !Recover Real Temperature
      !------------------------
      zvir = (Rv/Rd) - 1 ! Constant for virtual temp. calc. is approx. 0.608

      do k=1,Nk
      do j=1,l_nj
      do i=1,l_ni
         tt(i,j,k) = tv(i,j,k)/(1.d0 + zvir * qv(i,j,k))
      end do
      end do
      end do

      !Estimate U-V on scalar grids and Real Temperature
      !-------------------------------------------------
      if (F_stag_L) then

         istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
         istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
         istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt_plus)

         call hwnd_stag ( pw_uu_plus,pw_vv_plus,F_u,F_v, &
                          Mminx,Mmaxx,Mminy,Mmaxy,Nk,.false. )

         pw_tt_plus(1:l_ni,1:l_nj,1:Nk) = tt(1:l_ni,1:l_nj,1:Nk)

         F_t(1:l_ni,1:l_nj,1:Nk) = tv(1:l_ni,1:l_nj,1:Nk)

      else

         F_t(1:l_ni,1:l_nj,1:Nk) = tt(1:l_ni,1:l_nj,1:Nk)

      end if

      !---------------------------------------------------------------

      return

  999 call handle_error(-1,'DCMIP_INIT','Inappropriate list of tracers')

      end
