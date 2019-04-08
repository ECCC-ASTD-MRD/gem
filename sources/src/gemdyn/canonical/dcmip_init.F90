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

      subroutine dcmip_init
      use canonical
      use dcmip_options
      use gem_options
      use inp_mod
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use dyn_fisl_options
      use glb_ld
      use tr3d
      use gmm_itf_mod
      use tdpack, only : rgasd_8, rgasv_8
      implicit none
#include <arch_specific.hf>

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

      integer istat,istat1,istat2,istat3,istat4,err(2), &
              Deep,Pertt,Pert,Moist,Shear,Tracers,i,j,k
      real, pointer, dimension(:,:,:) :: cl,cl2,qv,qc,qr,q1,q2,q3,q4
      real, dimension(1,1,1) :: empty
      real  s_u(l_minx:l_maxx,l_miny:l_maxy,G_nk),s_v(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
             tv(l_minx:l_maxx,l_miny:l_maxy,G_nk), tt(l_minx:l_maxx,l_miny:l_maxy,G_nk)

      real(8), parameter :: Rd   = Rgasd_8, & ! cte gaz - air sec   [J kg-1 K-1]
                            Rv   = Rgasv_8    ! cte gaz - vap eau   [J kg-1 K-1]

      real(8)  zvir

      !--------------------------------------------------------------------------

      if (Schm_sleve_L ) call handle_error (-1,'DCMIP_init','  SLEVE not available YET  ')

      istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
      istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt_plus)
      istat = gmm_get (gmmk_ut1_s ,ut1 )
      istat = gmm_get (gmmk_vt1_s ,vt1 )
      istat = gmm_get (gmmk_wt1_s ,wt1 )
      istat = gmm_get (gmmk_tt1_s ,tt1 )
      istat = gmm_get (gmmk_zdt1_s,zdt1)
      istat = gmm_get (gmmk_st1_s ,st1 )
      istat = gmm_get (gmmk_sls_s ,sls )
      istat = gmm_get (gmmk_fis0_s,fis0)
      istat = gmm_get (gmmk_qt1_s ,qt1 )

      !Prescribed d(Zeta)dot and dz/dt
      !-------------------------------
      Inp_zd_L = .true.
      Inp_w_L  = .true.

      !Obtain specific humidity
      !------------------------
      istat = gmm_get ('TR/HU:P',qv)

      !Initialization QC/QR for Precipitation
      !--------------------------------------
      if (Dcmip_prec_type/=-1) then
         err= 0
         err(1) = gmm_get ('TR/QC:P',qc)
         err(2) = gmm_get ('TR/RW:P',qr)

         call gem_error(minval(err),'DCMIP_INIT','Tracers QC/RW required when Precipitation')

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
         err= 0
         err(1)= gmm_get ('TR/CL:P', cl )
         err(2)= gmm_get ('TR/CL2:P',cl2)

         call gem_error(minval(err),'DCMIP_INIT','Tracers CL/CL2 required when Chemistry')

          !--------------------------------------------------------------------------
          Deep  = 0           !Deep atmosphere (no=0)
          Pertt = 0           !Type of perturbation (exponential=0/stream function=1)
          Moist = Dcmip_moist !Moist=1/Dry=0 Initial conditions
          !--------------------------------------------------------------------------

         call dcmip_baroclinic_wave_2016 (ut1,vt1,wt1,tv,zdt1,st1,fis0,qv,cl,cl2, &
                  l_minx,l_maxx,l_miny,l_maxy,G_nk,Deep,Moist,Pertt,Dcmip_X,.true.)

         qt1 = 0. !ZERO log of non-hydrostatic perturbation pressure

      !DCMIP 2016: Tropical cyclone
      !----------------------------
      else if (Dcmip_case==162) then

         call dcmip_tropical_cyclone (ut1,vt1,wt1,tv,zdt1,st1,fis0,qv,&
                               l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         qt1 = 0. !ZERO log of non-hydrostatic perturbation pressure

      !DCMIP 2016: Supercell (Small planet)
      !------------------------------------
      else if (Dcmip_case==163) then

         istat = gmm_get(gmmk_thbase_s,thbase)

          !---------------------------------------------------------
          Pert = 1 !Thermal perturbation included (0 = no / 1 = yes)
          !---------------------------------------------------------

         call dcmip_supercell (ut1,vt1,wt1,tv,zdt1,st1,fis0,qv,Pert,thbase,&
                                    l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         qt1 = 0. !ZERO log of non-hydrostatic perturbation pressure

         !     Initialize u,v,zd,w,tv,qv,qc,rw,theta REFERENCE for Vertical diffusion
         !----------------------------------------------------------------------
         istat = gmm_get(gmmk_uref_s , uref )
         istat = gmm_get(gmmk_vref_s , vref )
         istat = gmm_get(gmmk_wref_s , wref )
         istat = gmm_get(gmmk_zdref_s, zdref)
         istat = gmm_get(gmmk_qvref_s, qvref)
         istat = gmm_get(gmmk_qcref_s, qcref)
         istat = gmm_get(gmmk_qrref_s, qrref)
         istat = gmm_get(gmmk_thref_s, thref)

         !     Prepare UREF/VREF (on staggered grids) for DCMIP_VRD
         !----------------------------------------------------
         if (.not..true.) then

            call hwnd_stag (s_u,s_v,ut1,vt1,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

            uref(1:l_ni-1,1:l_nj,  1:G_nk) =    s_u(1:l_ni-1,1:l_nj,  1:G_nk)
            vref(1:l_ni,  1:l_nj-1,1:G_nk) =    s_v(1:l_ni,  1:l_nj-1,1:G_nk)

         else

            uref(1:l_ni-1,1:l_nj,  1:G_nk) = ut1(1:l_ni-1,1:l_nj,  1:G_nk)
            vref(1:l_ni,  1:l_nj-1,1:G_nk) = vt1(1:l_ni,  1:l_nj-1,1:G_nk)

         end if

          wref(1:l_ni,  1:l_nj,  1:G_nk) =   wt1 (1:l_ni,  1:l_nj,  1:G_nk)
         zdref(1:l_ni,  1:l_nj,  1:G_nk) =   zdt1(1:l_ni,  1:l_nj,  1:G_nk)
         qvref(1:l_ni,  1:l_nj,  1:G_nk) =     qv(1:l_ni,  1:l_nj,  1:G_nk)
         qcref(1:l_ni,  1:l_nj,  1:G_nk) =     qc(1:l_ni,  1:l_nj,  1:G_nk)
         qrref(1:l_ni,  1:l_nj,  1:G_nk) =     qr(1:l_ni,  1:l_nj,  1:G_nk)
         thref(1:l_ni,  1:l_nj,  1:G_nk) = thbase(1:l_ni,  1:l_nj,  1:G_nk)

      !DCMIP 2012: Steady-State Atmosphere at Rest in the Presence of Orography
      !------------------------------------------------------------------------
      else if (Dcmip_case==20) then

          !Set initial conditions according to prescribed mountain
          !-------------------------------------------------------
         call dcmip_steady_state_mountain (ut1,vt1,zdt1,tv,qv,fis0,st1,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.,.true.)

         if (Vtopo_L) then

            istat = gmm_get(gmmk_topo_low_s , topo_low )
            istat = gmm_get(gmmk_topo_high_s, topo_high)

            topo_low (1:l_ni,1:l_nj) = 0.
            topo_high(1:l_ni,1:l_nj) = fis0 (1:l_ni,1:l_nj)
            fis0   (1:l_ni,1:l_nj) = 0.

              !Reset initial conditions according to topo_low
              !----------------------------------------------
            call dcmip_steady_state_mountain (ut1,vt1,zdt1,tv,qv,fis0,st1,&
                           l_minx,l_maxx,l_miny,l_maxy,G_nk,.false.,.true.)

         end if

         qt1 = 0. !ZERO log of non-hydrostatic perturbation pressure
         wt1 = 0. !ZERO Dz/Dt

      !DCMIP 2012: Pure 3D Advection
      !-----------------------------
      else if (Dcmip_case>=11.and.Dcmip_case<=13) then

          !Get tracers Q1,Q2,Q3,Q4
          !-----------------------
         istat1 = gmm_get('TR/Q1:P',q1)
         istat2 = gmm_get('TR/Q2:P',q2)
         istat3 = gmm_get('TR/Q3:P',q3)
         istat4 = gmm_get('TR/Q4:P',q4)

         if ((istat1/=0.or.istat2/=0.or.istat3/=0.or.istat4/=0).and.Dcmip_case==11) goto 999
         if ((istat1/=0)                                       .and.Dcmip_case==12) goto 999
         if ((istat1/=0.or.istat2/=0.or.istat3/=0.or.istat4/=0).and.Dcmip_case==13) goto 999

          !3D deformational flow
          !---------------------
         if (Dcmip_case==11) call dcmip_tracers11_transport (ut1,vt1,zdt1,tv,qv,fis0,st1,q1,q2,q3,q4, &
                                                             l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

          !3D Hadley-like meridional circulation
          !-------------------------------------
         if (Dcmip_case==12) call dcmip_tracers12_transport (ut1,vt1,zdt1,tv,qv,fis0,st1,q1,&
                                                             l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

          !2D solid-body rotation of thin cloud-like tracer in the presence of orography
          !-----------------------------------------------------------------------------
         if (Dcmip_case==13) call dcmip_tracers13_transport (ut1,vt1,zdt1,tv,qv,fis0,st1,q1,q2,q3,q4, &
                                                             l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

          !Store REFERENCE at initial time
          !-------------------------------
         istat = gmm_get(gmmk_q1ref_s,q1ref)
         istat = gmm_get(gmmk_q2ref_s,q2ref)
         istat = gmm_get(gmmk_q3ref_s,q3ref)
         istat = gmm_get(gmmk_q4ref_s,q4ref)

         if (Dcmip_case>  0) q1ref(1:l_ni,1:l_nj,1:G_nk) = q1(1:l_ni,1:l_nj,1:G_nk)
         if (Dcmip_case/=12) q2ref(1:l_ni,1:l_nj,1:G_nk) = q2(1:l_ni,1:l_nj,1:G_nk)
         if (Dcmip_case/=12) q3ref(1:l_ni,1:l_nj,1:G_nk) = q3(1:l_ni,1:l_nj,1:G_nk)
         if (Dcmip_case/=12) q4ref(1:l_ni,1:l_nj,1:G_nk) = q4(1:l_ni,1:l_nj,1:G_nk)

      !DCMIP 2012: Mountain waves over a Schaer-type mountain on a small planet
      !------------------------------------------------------------------------
         else if (Dcmip_case==21.or.Dcmip_case==22) then

            if (Dcmip_case==21) Shear = 0 !Without wind shear
            if (Dcmip_case==22) Shear = 1 !With    wind shear

            call dcmip_Schaer_mountain (ut1,vt1,zdt1,tv,qv,fis0,st1,&
                 l_minx,l_maxx,l_miny,l_maxy,G_nk,Shear,.true.,.true.)

         if (Vtopo_L) then

            istat = gmm_get(gmmk_topo_low_s , topo_low )
            istat = gmm_get(gmmk_topo_high_s, topo_high)

            topo_low (1:l_ni,1:l_nj) = 0.
            topo_high(1:l_ni,1:l_nj) = fis0 (1:l_ni,1:l_nj)
            fis0   (1:l_ni,1:l_nj) = 0.

              !Reset initial conditions according to topo_low
              !----------------------------------------------
            call dcmip_Schaer_mountain (ut1,vt1,zdt1,tv,qv,fis0,st1,&
               l_minx,l_maxx,l_miny,l_maxy,G_nk,Shear,.false.,.true.)

         end if

         qt1 = 0.               !ZERO log of non-hydrostatic perturbation pressure
         wt1 = 0.               !ZERO Dz/Dt

      !DCMIP 2012: Gravity wave on a small planet along the equator
      !------------------------------------------------------------
      else if (Dcmip_case==31) then

         istat = gmm_get(gmmk_thbase_s,thbase)

         call dcmip_gravity_wave (ut1,vt1,zdt1,tv,qv,fis0,st1,thbase,&
                              l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         qt1 = 0.               !ZERO log of non-hydrostatic perturbation pressure
         wt1 = 0.               !ZERO Dz/Dt

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
         Moist   = Dcmip_moist  ! Moist=1/Dry=0 Initial conditions
         Tracers = 1            ! Tracers=1/No Tracers=0
          !-------------------------------------------------------

          !Dynamical Tracers: Potential temperature and Ertel's potential vorticity
          !------------------------------------------------------------------------
         err= 0
         err(1) = gmm_get('TR/Q1:P',q1)
         err(2) = gmm_get('TR/Q2:P',q2)

         call gem_error(minval(err),'DCMIP_INIT','Tracers Q1/Q2 required when Dcmip_case=41X')

         call dcmip_baroclinic_wave_2012 (ut1,vt1,zdt1,tv,qv,fis0,st1,q1,q2, &
               l_minx,l_maxx,l_miny,l_maxy,G_nk,Moist,Dcmip_X,Tracers,.true.)

         qt1 = 0.         !ZERO log of non-hydrostatic perturbation pressure
         wt1 = 0.               !ZERO Dz/Dt

      !DCMIP 2012: Moist Baroclinic Instability driven by Simple Physics
      !-----------------------------------------------------------------
      else if (Dcmip_case==43) then

          !-------------------------------------------------------
         Moist   = Dcmip_moist  ! Moist=1/Dry=0 Initial conditions
         Tracers = 0            ! Tracers=1/No Tracers=0
          !-------------------------------------------------------

         call dcmip_baroclinic_wave_2012 (ut1,vt1,zdt1,tv,qv,fis0,st1,empty,empty, &
                       l_minx,l_maxx,l_miny,l_maxy,G_nk,Moist,Dcmip_X,Tracers,.true.)

         qt1 = 0. !ZERO log of non-hydrostatic perturbation pressure
         wt1 = 0.               !ZERO Dz/Dt

      else

         call gem_error(-1,'dcmip_init','DCMIP_CASE 2012/2016 not available')

      end if

      call rpn_comm_xch_halo (fis0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,1, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !Recover Real Temperature
      !------------------------
      zvir = (Rv/Rd) - 1 ! Constant for virtual temp. calc. is approx. 0.608

      do k=1,G_nk
      do j=1,l_nj
      do i=1,l_ni
         tt(i,j,k) = tv(i,j,k)/(1.d0 + zvir * qv(i,j,k))
      end do
      end do
      end do

      !Estimate U-V on scalar grids and Real Temperature
      !-------------------------------------------------
!!$      if (.true.) then

         call hwnd_stag ( pw_uu_plus,pw_vv_plus,ut1,vt1, &
                          l_minx,l_maxx,l_miny,l_maxy,G_nk,.false. )

         pw_tt_plus(1:l_ni,1:l_nj,1:G_nk) = tt(1:l_ni,1:l_nj,1:G_nk)

         tt1(1:l_ni,1:l_nj,1:G_nk) = tv(1:l_ni,1:l_nj,1:G_nk)

!!$      else
!!$
!!$         tt1(1:l_ni,1:l_nj,1:G_nk) = tt(1:l_ni,1:l_nj,1:G_nk)
!!$
!!$      end if

      !---------------------------------------------------------------

      return

  999 call gem_error(-1,'DCMIP_INIT','Inappropriate list of tracers')

      end
