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

!**s/r dcmip_init - Prepare initial conditions (u,v,w,tv,zd,s,q,topo) for DCMIP 2012/2016 runs

      subroutine dcmip_init ()

      use canonical
      use cstv
      use dcmip_options
      use dyn_fisl_options
      use dynkernel_options
      use gem_options
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use inp_mod
      use mem_tracers
      use tdpack, only : rgasd_8, rgasv_8, grav_8, cpd_8
      use tr3d
      use ver

      implicit none

#include <arch_specific.hf>

      !object
      !======================================================================|
      !Prepare initial conditions for DCMIP 2012/2016 runs                   |
      !----------------------------------------------------------------------|
      !NOTE: U,V output on Staggered grids                                   |
      !----------------------------------------------------------------------|
      ! case DCMIP 2012    | Pure advection                                  |
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

      integer :: istat1,istat2,istat3,istat4,err(2), &
                 Deep,Pertt,Pert,Moist,Shear,Tracers,i,j,k,iter,niter

      real, pointer, dimension(:,:,:) :: cl,cl2,qv,qc,qr,q1,q2,q3,q4
      real, dimension(1,1,1) :: empty

      real, dimension(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) :: zmom_8,lg_pstar,log_pt,log_pm,pt_plus,pm_plus
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk)     :: tt
      real, dimension(l_minx:l_maxx,l_miny:l_maxy)          :: ps

      real(kind=REAL64), parameter :: Rd = Rgasd_8, & ! cte gaz - air sec   [J kg-1 K-1]
                                      Rv = Rgasv_8    ! cte gaz - vap eau   [J kg-1 K-1]

      real(kind=REAL64) :: zvir,aaa_8

      logical :: GEM_P_L
!
!---------------------------------------------------------------------
!
      if (Schm_sleve_L ) call gem_error(-1,'DCMIP_INIT','  SLEVE not available YET  ')

      GEM_P_L = trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P'

      !Prescribed d(Zeta)dot and dz/dt
      !-------------------------------
      Inp_zd_L = .true.
      Inp_w_L  = .true.

      !Obtain specific humidity
      !------------------------
      qv => tracers_P(Tr3d_hu)%pntr

      !Initialization QC/QR for Precipitation
      !--------------------------------------
      if (Dcmip_prec_type/=-1) then

         err(1) = tr_get('QC:P',qc)
         err(2) = tr_get('RW:P',qr)

         if (err(1)<=0.or.err(2)<=0) call handle_error(-1,'DCMIP_INIT','Tracers QC/RW required when Precipitation')

         qc = 0. !ZERO Cloud water mixing ratio
         qr = 0. !ZERO Rain  water mixing ratio

         art = 0. !ZERO Averaged precipitation rate
         wrt = 0. !ZERO Averaged precipitation rate (WORK FIELD)

      end if

      !DCMIP 2016: Baroclinic wave with Toy Terminal Chemistry
      !-------------------------------------------------------
      if (Dcmip_case==161) then

         err(1) = tr_get('CL:P',cl)
         err(2) = tr_get('CL2:P',cl2)

         if (err(1)<=0.or.err(2)<=0) call handle_error (-1,'DCMIP_INIT','Tracers CL/CL2 required when Chemistry')

         !--------------------------------------------------------------------------
         Deep  = 0           !Deep atmosphere (no=0)
         Pertt = 0           !Type of perturbation (exponential=0/stream function=1)
         Moist = Dcmip_moist !Moist=1/Dry=0 Initial conditions
         !--------------------------------------------------------------------------

         call dcmip_baroclinic_wave_2016 (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,cl,cl2, &
                                          l_minx,l_maxx,l_miny,l_maxy,G_nk,Deep,Moist,Pertt,Dcmip_X,.true.)

      !DCMIP 2016: Tropical cyclone
      !----------------------------
      else if (Dcmip_case==162) then

         call dcmip_tropical_cyclone (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

      !DCMIP 2016: Supercell (Small planet)
      !------------------------------------
      else if (Dcmip_case==163) then

         !---------------------------------------------------------
         Pert = 1 !Thermal perturbation included (0 = no / 1 = yes)
         !---------------------------------------------------------

         call dcmip_supercell (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,Pert,thbase,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         !----------------------------------------------------------------------
         !Initialize u,v,zd,w,tv,qv,qc,rw,theta REFERENCE for Vertical diffusion
         !----------------------------------------------------------------------

         !Prepare UREF/VREF (on staggered grids) for DCMIP_VRD
         !----------------------------------------------------
         uref(1:l_ni-1,1:l_nj,  1:G_nk) =    ut1(1:l_ni-1,1:l_nj,  1:G_nk)
         vref(1:l_ni,  1:l_nj-1,1:G_nk) =    vt1(1:l_ni,  1:l_nj-1,1:G_nk)

          wref(1:l_ni, 1:l_nj,  1:G_nk) =    wt1(1:l_ni,  1:l_nj,  1:G_nk)
         zdref(1:l_ni, 1:l_nj,  1:G_nk) =   zdt1(1:l_ni,  1:l_nj,  1:G_nk)
         qvref(1:l_ni, 1:l_nj,  1:G_nk) =     qv(1:l_ni,  1:l_nj,  1:G_nk)
         qcref(1:l_ni, 1:l_nj,  1:G_nk) =     qc(1:l_ni,  1:l_nj,  1:G_nk)
         qrref(1:l_ni, 1:l_nj,  1:G_nk) =     qr(1:l_ni,  1:l_nj,  1:G_nk)
         thref(1:l_ni, 1:l_nj,  1:G_nk) = thbase(1:l_ni,  1:l_nj,  1:G_nk)

      !DCMIP 2012: Steady-State Atmosphere at Rest in the Presence of Orography
      !------------------------------------------------------------------------
      else if (Dcmip_case==20) then

         !Set initial conditions according to prescribed mountain
         !-------------------------------------------------------
         call dcmip_steady_state_mountain (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.,.true.)

         if (Vtopo_L) then

            topo_low (1:l_ni,1:l_nj) = 0.
            topo_high(1:l_ni,1:l_nj) = fis0(1:l_ni,1:l_nj)
            fis0     (1:l_ni,1:l_nj) = 0.

            !Reset initial conditions according to topo_low
            !----------------------------------------------
            call dcmip_steady_state_mountain (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,l_minx,l_maxx,l_miny,l_maxy,G_nk,.false.,.true.)

         end if

      !DCMIP 2012: Pure 3D Advection
      !-----------------------------
      else if (Dcmip_case>=11.and.Dcmip_case<=13) then

         !Get tracers Q1,Q2,Q3,Q4
         !-----------------------
         istat1 = tr_get('Q1:P',q1)
         istat2 = tr_get('Q2:P',q2)
         istat3 = tr_get('Q3:P',q3)
         istat4 = tr_get('Q4:P',q4)

         if (((istat1<=0.or.istat2<=0.or.istat3<=0.or.istat4<=0).and.Dcmip_case==11).or. &
             ((istat1<=0)                                       .and.Dcmip_case==12).or. &
             ((istat1<=0.or.istat2<=0.or.istat3<=0.or.istat4<=0).and.Dcmip_case==13)) then

             call handle_error(-1,'DCMIP_INIT','Inappropriate list of tracers')

         end if

         !3D deformational flow
         !---------------------
         if (Dcmip_case==11) call dcmip_tracers11_transport (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,q1,q2,q3,q4, &
                                                             l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         !3D Hadley-like meridional circulation
         !-------------------------------------
         if (Dcmip_case==12) call dcmip_tracers12_transport (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,q1,&
                                                             l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         !2D solid-body rotation of thin cloud-like tracer in the presence of orography
         !-----------------------------------------------------------------------------
         if (Dcmip_case==13) call dcmip_tracers13_transport (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,q1,q2,q3,q4, &
                                                             l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

         !Store REFERENCE at initial time
         !-------------------------------
         if (Dcmip_case>  0) q1ref(1:l_ni,1:l_nj,1:G_nk) = q1(1:l_ni,1:l_nj,1:G_nk)
         if (Dcmip_case/=12) q2ref(1:l_ni,1:l_nj,1:G_nk) = q2(1:l_ni,1:l_nj,1:G_nk)
         if (Dcmip_case/=12) q3ref(1:l_ni,1:l_nj,1:G_nk) = q3(1:l_ni,1:l_nj,1:G_nk)
         if (Dcmip_case/=12) q4ref(1:l_ni,1:l_nj,1:G_nk) = q4(1:l_ni,1:l_nj,1:G_nk)

      !DCMIP 2012: Mountain waves over a Schaer-type mountain on a small planet
      !------------------------------------------------------------------------
      else if (Dcmip_case==21.or.Dcmip_case==22) then

         if (Dcmip_case==21) Shear = 0 !Without wind shear
         if (Dcmip_case==22) Shear = 1 !With    wind shear

         call dcmip_Schaer_mountain (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,l_minx,l_maxx,l_miny,l_maxy,G_nk,Shear,.true.,.true.)

         if (Vtopo_L) then

            topo_low (1:l_ni,1:l_nj) = 0.
            topo_high(1:l_ni,1:l_nj) = fis0(1:l_ni,1:l_nj)
            fis0     (1:l_ni,1:l_nj) = 0.

            !Reset initial conditions according to topo_low
            !----------------------------------------------
            call dcmip_Schaer_mountain (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,l_minx,l_maxx,l_miny,l_maxy,G_nk,Shear,.false.,.true.)

         end if

      !DCMIP 2012: Gravity wave on a small planet along the equator
      !------------------------------------------------------------
      else if (Dcmip_case==31) then

         call dcmip_gravity_wave (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,thbase,thfull,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

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
         err(1) = tr_get('Q1:P',q1)
         err(2) = tr_get('Q2:P',q2)

         if (err(1)<=0.or.err(2)<=0) call handle_error(-1,'DCMIP_INIT','Tracers Q1/Q2 required when Dcmip_case=41X')

         call dcmip_baroclinic_wave_2012 (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,q1,q2, &
                                          l_minx,l_maxx,l_miny,l_maxy,G_nk,Moist,Dcmip_X,Tracers,.true.)

      !DCMIP 2012: Moist Baroclinic Instability driven by Simple Physics
      !-----------------------------------------------------------------
      else if (Dcmip_case==43) then

         !-------------------------------------------------------
         Moist   = Dcmip_moist ! Moist=1/Dry=0 Initial conditions
         Tracers = 0           ! Tracers=1/No Tracers=0
         !-------------------------------------------------------

         call dcmip_baroclinic_wave_2012 (ut1,vt1,wt1,zdt1,tt1,qv,fis0,st1,ps,empty,empty, &
                                          l_minx,l_maxx,l_miny,l_maxy,G_nk,Moist,Dcmip_X,Tracers,.true.)

      else

         call handle_error(-1,'dcmip_init','DCMIP_CASE 2012/2016 not available')

      end if

      call rpn_comm_xch_halo (fis0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,1, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !GEM-H: Estimate Q (Pressure deviation) as in INP_BASE_H and Virtual temperature (T31 only)
      !------------------------------------------------------------------------------------------
      if (.not.GEM_P_L) then

         zmom_8(:,:,0) = Ver_z_8%m(0)

         do k=1,G_nk
            zmom_8(:,:,k) = Ver_z_8%m(k) + Ver_b_8%m(k)*fis0(:,:)/grav_8
         end do

         zmom_8(:,:,G_nk+1) = fis0(:,:)/grav_8

         lg_pstar(:,:,G_nk+1) = log(1.d5) - grav_8*zmom_8(:,:,G_nk+1)/(rgasd_8*Cstv_Tstr_8)

         do k=G_nk,1,-1
            lg_pstar(:,:,k) = lg_pstar(:,:,k+1) + grav_8*(zmom_8(:,:,k+1)-zmom_8(:,:,k))/(rgasd_8*Cstv_Tstr_8)
         end do

         niter = 1
         if (Dcmip_case==31) niter = 5

         do iter = 1,niter

            !Obtain pressure deviation Q from Surface Pressure and Virtual Temperature
            !-------------------------------------------------------------------------
            aaa_8 = rgasd_8 * Cstv_Tstr_8

            do j=1,l_nj
               do i=1,l_ni
                  qt1(i,j,G_nk+1) = aaa_8 * log(ps(i,j)/1.e5)
               end do
            end do

            aaa_8 = grav_8 * Cstv_Tstr_8

            do k=G_nk,1,-1
               do j=1,l_nj
                  do i=1,l_ni
                     qt1(i,j,k) = qt1(i,j,k+1) + aaa_8 * (zmom_8(i,j,k+1) - zmom_8(i,j,k))/tt1(i,j,k)
                  end do
               end do
            end do

            do k=1,G_nk+1
               do j=1,l_nj
                  do i=1,l_ni
                     qt1(i,j,k) = qt1(i,j,k) + grav_8*zmom_8(i,j,k)
                  end do
               end do
            end do

            if (niter>1) then

               !Obtain Momentum and Thermo pressures
               !------------------------------------
               do k=1,G_nk

                  if(k == 1) then
                     log_pm(1:l_ni,1:l_nj,k) = (qt1(1:l_ni,1:l_nj,k)/(rgasd_8*Cstv_Tstr_8)+lg_pstar(1:l_ni,1:l_nj,k))
                  end if

                  pm_plus(1:l_ni,1:l_nj,k) = exp(log_pm(1:l_ni,1:l_nj,k))

                  log_pm(1:l_ni,1:l_nj,k+1) = (qt1(1:l_ni,1:l_nj,k+1)/(rgasd_8*Cstv_Tstr_8)+lg_pstar(1:l_ni,1:l_nj,k+1))

                  if(k==G_nk) &
                  ps(1:l_ni,1:l_nj) = exp(log_pm(1:l_ni,1:l_nj,G_nk+1))

               end do

               do k=1,G_nk
                  log_pt(1:l_ni,1:l_nj,k) = 0.5*(log_pm(1:l_ni,1:l_nj,k+1)+log_pm(1:l_ni,1:l_nj,k))
               end do

               log_pt(1:l_ni,1:l_nj,G_nk+1) = log_pm(1:l_ni,1:l_nj,G_nk+1)

               do k=1,G_nk
                  pt_plus(1:l_ni,1:l_nj,k) = exp(log_pt(1:l_ni,1:l_nj,k))
               end do

               !Obtain revised Virtual Temperature !qv==0
               !-----------------------------------------
               do k=1,G_nk
                  do j=1,l_nj
                     do i=1,l_ni
                        tt1(i,j,k) = thfull(i,j,k) * (pt_plus(i,j,k)/Cstv_pref_8) ** (rgasd_8/cpd_8)
                     end do
                  end do
               end do

            end if

         end do

      !GEM-P: ZERO log of non-hydrostatic perturbation pressure
      !--------------------------------------------------------
      else

         qt1 = 0.

      end if

      !Recover Real Temperature
      !------------------------
      zvir = (Rv/Rd) - 1 ! Constant for virtual temp. calc. is approx. 0.608

      do k=1,G_nk
      do j=1,l_nj
      do i=1,l_ni
         tt(i,j,k) = tt1(i,j,k)/(1.d0 + zvir * qv(i,j,k))
      end do
      end do
      end do

      !Estimate U-V on scalar grids and Real Temperature in PW comdeck
      !---------------------------------------------------------------
      call hwnd_stag ( pw_uu_plus,pw_vv_plus,ut1,vt1, &
                       l_minx,l_maxx,l_miny,l_maxy,G_nk,.false. )

      pw_tt_plus(1:l_ni,1:l_nj,1:G_nk) = tt(1:l_ni,1:l_nj,1:G_nk)
!
!---------------------------------------------------------------------
!
      return
      end
