!-------------------------------------- LICENCE BEGIN --------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module linoz
   implicit none
   private
   public :: linoz3

contains

   !/@*
   subroutine linoz3(d,v,f,dsiz,vsiz,fsiz,dt,kount,trnch,ni,nkm1,nk)
      use iso_c_binding
      use debug_mod, only: init2nan
      use tdpack_const, only: CAPPA, CONSOL2, GRAV, PI, STEFAN, RGASD
      use phy_options
      use phybus
      use linoz_mod
      implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      !@object Stratospheric Ozone chemistry
      ! Produces the ozone tendency due to stratospheric 
      ! photochemistry (based on LINOZ scheme)
      !@arguments
      ! n        horizonal index
      ! nkm1     vertical  index nk-1
      ! nk       vertical  index
      ! dt       timestep
      ! trnch    index of vertical slice n*nk
      ! kount    number of timesteps
      ! v        volatile bus
      ! f        permanent bus
      ! d        dynamics bus

      integer,    intent(in)    :: ni,nkm1,nk,dsiz,vsiz,fsiz,kount,trnch
      real,       intent(in)    :: dt
      real,target,intent(inout) :: d(dsiz), v(vsiz), f(fsiz)

      !@author J. de Grandpre (ARQI): February 2013
      !@revisions
      ! * I. Ivanova (2015) - Complete remake
      !*@/
#include <msg.h>

      include "phyinput.inc"
      include "ccc_tracegases.cdk"

      external ccc_tracedata

      ! declaration of local arguments
      !
      !          - input -
      ! o3_vmr   linoz tracer (t+)        3D  o3    vmr (mole /mole)
      ! ch4_vmr  linoz tracer (t+)        3D  ch4   vmr (mole /mole)
      ! n2o_vmr  linoz tracer (t+)        3D  n2o   vmr (mole /mole)
      ! f11_vmr  linoz tracer (t+)        3D  cfc11 vmr (mole /mole)
      ! f12_vmr  linoz tracer (t+)        3D  cfc12 vmr (mole /mole)
      ! o3fk_vmr climatology F-K  zonavg  2D  o3    vmr (mole /mole)
      ! tt       temperature (K) mid-thermo layer
      ! qq       specific humidity (kg /kg)
      ! ps       surface pressure
      ! shtj     sigma levels on momentum/flux levels (nk)
      !
      !          - input/output -
      ! zo3lplus  linoz tracer (t*)        3D  o3    mmr (micro g /kg air)
      ! ziarplus age of air               3D  air   seconds

      real, dimension(ni)       :: ps
      real, dimension(ni, nkm1) :: tt, qq, aird                 ! nk-1 layers
      real, dimension(ni, nk)   :: shtj,s_qrt                   ! nk "flux" levels

      real, dimension(ni, nkm1) :: o3_vmr,o3fk_vmr,o3c_vmr,ch4_vmr,n2o_vmr,f11_vmr,f12_vmr ! nk-1 layers
      real, dimension(ni, nkm1) :: o1_tend, o4_tend, o6_tend, o7_tend
      real, dimension(ni, nkm1) :: o3_tend, ch4_tend, n2o_tend, f11_tend, f12_tend
      real, dimension(ni, nkm1) :: o3_over
      real, dimension(ni, nkm1) :: o3new, ch4new, f11new, f12new, n2onew
      real, dimension(ni, nkm1) :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14 !nk-1 layers

      integer :: i, k, ni2
      real :: o3fk_over,o3c_over,ch4_over,n2o_over,f11_over,f12_over,ptop,pbot
      real :: delage, delrlx
      real :: dp

#define PHYPTRDCL
#include "linoz_ptr.hf"

      !----------------------------------------------------------------
      if (.not.llinoz) return

      call msg_toall(MSG_DEBUG, 'linoz [BEGIN]')

#undef PHYPTRDCL
#include "linoz_ptr.hf"

      !   call init2nan(o3_vmr,ch4_vmr,n2o_vmr,f11_vmr,f12_vmr,o3fk_vmr)
      !   call init2nan(tt,qq,aird)
      !   call init2nan(ps)
      !   call init2nan(shtj)

      call init2nan(o3new,o3_over,o3_tend,o1_tend,o4_tend,o6_tend,o7_tend)

      !   call init2nan(ch4new,n2onew,f11new,f12new)
      !   call init2nan(ch4_tend,n2o_tend,f11_tend,f12_tend)

      !   call init2nan(zch4chmtd,zn2ochmtd,zf11chmtd,zf12chmtd)
      !   call init2nan(zch4col,zn2ocol,zf11col,zf12col)
      !   call init2nan(zch4tc,zn2otc,zf11tc,zf12tc)

      c1 = 0.
      c2 = 0.
      c3 = 0.
      c4 = 0.
      c5 = 0.
      c6 = 0.
      c7 = 0.

      c8 = 0.
      c9 = 0.
      c10= 0.
      c11= 0.
      c12= 0.
      c13= 0.
      c14= 0.

      ! Mass mixing ratio kg / kg air  mmr
      ch4_ppm = qch4 * 1.e-6                 ! ppmv 
      rmch4   = ch4_ppm * mwt_ch4 / mwt_air  ! kg / kg air

      n2o_ppm = qn2o * 1.e-6
      rmn2o   = n2o_ppm * mwt_n2o / mwt_air

      f11_ppm = qcfc11 * 1.e-9
      rmf11   = f11_ppm * mwt_f11 / mwt_air

      f12_ppm = qcfc12 * 1.e-9
      rmf12   = f12_ppm * mwt_f12 / mwt_air

      ! Calculate age of air
      !
      ! ---  change in age of air (age in units of years)
      !   --- assuming points with a water vapour mixing ratio
      !         greater than 20 ppmv are in the troposphere (ozone smaller than 150 ppbv)
      !   --- tropospheric points are relaxed back to zero age (here
      !         defined with an offset of 2 years) with a time constant
      !         of 0.1 days
      delage = dt/(365.0*86400.0)
      delrlx = exp(-1.0*dt/8640.0)

      ni2 = int(float(ni)/2.)

      IF_KOUNT0: if (kount == 0) then

         if (radghg_L) then

            ! ERA-5 ozone
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'o3ce')) zo3ce = -1.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'ch4c')) zch4c = 1.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'n2oc')) zn2oc = 1.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'cf1c')) zcf1c = 1.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'cf2c')) zcf2c = 1.
            zch4  = rmch4*zch4c ! kg /kg air mmr
            zn2o  = rmn2o*zn2oc ! kg /kg air mmr
            zcf11 = rmf11*zcf1c ! kg /kg air mmr
            zcf12 = rmf12*zcf2c ! kg /kg air mmr

            zch4  =  zch4 / mwt_ch4 *mwt_air  ! mole /mole vmr 
            zn2o  =  zn2o / mwt_n2o *mwt_air  ! mole /mole vmr
            zcf11 = zcf11 / mwt_f11 *mwt_air  ! mole /mole vmr
            zcf12 = zcf12 / mwt_f12 *mwt_air  ! mole /mole vmr

            zch4  = max( zch4, 1.E-20)      ! mole /mole vmr
            zn2o  = max( zn2o, 1.E-20)      ! mole /mole vmr
            zcf11 = max(zcf11, 1.E-20)      ! mole /mole vmr
            zcf12 = max(zcf12, 1.E-20)      ! mole /mole vmr
            
         endif

         !      if (.not.any(phyinread_list_s(1:phyinread_n) == 'tts' )) ztts  = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'ttce')) zttce = 0.
         !      if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin1')) zlin1 = 0.
         !      if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin2')) zlin2 = 0.
         !      if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin3')) zlin3 = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin4')) zlin4 = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin5')) zlin5 = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin6')) zlin6 = 0.
         if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin7')) zlin7 = 0.

         ! Initialize LINOZ ozone tracer with climatology, if not read from analysis
         if (.not.any(dyninread_list_s == 'o3l') .or.   &
              .not.any(phyinread_list_s(1:phyinread_n) == 'tr/o3l:p')) &
              !                                              zo3lplus = zo3fk * 1E+9
              !micro g /kg air <--      kg /kg air
              ! climato+ghg+linoz+era5+ (Irena I)
              zo3lplus = zo3ce                             !micro g /kg air
         ! climato_phase2 (Paul V)
         !                                              zo3lplus = zo3ce * 1E+9
         !micro g /kg air <--      kg /kg air

         IF_LINGHG1: if (llingh) then

            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin8' )) zlin8  = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin9' )) zlin9  = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin10')) zlin10 = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin11')) zlin11 = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin12')) zlin12 = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin13')) zlin13 = 0.
            if (.not.any(phyinread_list_s(1:phyinread_n) == 'lin14')) zlin14 = 0.

            if (.not.any(dyninread_list_s == 'ch4l') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/ch4l:p')) zch4lplus = zch4         !mole /mole vmr
            if (.not.any(dyninread_list_s == 'n2ol') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/n2ol:p')) zn2olplus = zn2o         !mole /mole vmr
            if (.not.any(dyninread_list_s == 'f11l') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/f11l:p')) zf11lplus = zcf11        !mole /mole vmr
            if (.not.any(dyninread_list_s == 'f12l') .or. &
                 .not.any(phyinread_list_s(1:phyinread_n) == 'tr/f12l:p')) zf12lplus = zcf12        !mole /mole vmr

         end if IF_LINGHG1

      endif IF_KOUNT0


      ! Calculate shtj, sigma at flux levels
      
      do i = 1, ni
         s_qrt(i,1) = zsigw(i,1) / zsigw(i,2)
         s_qrt(i,nk) = 1.0
      enddo
      do k = 2, nkm1
         do i = 1, ni
            s_qrt(i,k) = zsigw(i,k-1) * zsigw(i,k)
         enddo
      enddo

      call vssqrt(shtj,s_qrt,ni*nk)

      do i = 1, ni
         shtj(i,1) = zsigw(i,1) * shtj(i,1)
      enddo


      ! --------------------
      !  Loop on longitudes
      ! -------------------
      DO_I1: do i = 1, ni

         o3fk_over = 0.
         o3c_over = 0.

         ! ----------------------------
         ! Loop on all vertical levels  
         ! ----------------------------

         DO_K1: do k = 1, nkm1 !nk-1 layers, 1=TOA-layer, nk=SFC-layer

            !        if(local_dbg .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
            !           write(lun_out,*) 'linoz: ', i,trnch,k, &
            !              'o3lplus, o3fk_vmr=',                  &
            !               zo3lplus(i,k), o3fk_vmr(i,k)   !,        &
            !              'pplus, shtj, tplus, huplus=',      &

            !               zpplus (i), shtj(i,k), ztplus(i,k), zhuplus(i,k)

            ps  (i)   = zpplus(i)
            tt  (i,k) = ztplus(i,k)
            qq  (i,k) = zhuplus(i,k)
            aird(i,k) = zsigw(i,k) * ps(i) / (tt(i,k) * rgasd) !kg /m3

            !        if(mod(i, ni2)==0 .and. mod(trnch, 32)==0)   &
            !           write(lun_out,*) 'linoz: ', i,trnch, k, &
            !             'tt qq aird ps=', tt(i,k), qq(i,k), aird(i,k),ps(i)


            ! output air density (molec. /cm3)
            zrhoa(i,k) = aird(i,k) 
            !
            !        aird_vmr  = aird_mmr    !mole /mole <-- kg /kg
            !        rho_air   = aird_vmr * avno !molecules (x) mole^-1a <-- mole /mole
            !        rho_air   = aird(i,k) * ???

            ! avno      = 0.6022000000000e+24 ! avogadro's num        atoms mol-1
            ! grav      = .980616e+1          ! m s-2; gravitational acceleration
            ! mwt_air   = 0.2897000000000e+02 ! mol. wgt. of air  g mol-1
            !
            ! Column X
            !
            ! Col X (molec/cm2) = Y [mole/mole] /mwt_air * avno * 1.E+3 (g/kg) * 1E-4 (m2/cm2)
            !
            ! 1E-4 m2/cm2; 1E+3 g/kg => 1E-1
            !       cst_x = grav /mwt_air * avno * 1.0E+1    ! molec/cm2
            ! 
            !      x_nmr = x_vmr * avno                     ! molecules (x) mole^-1a <-- mole /mole
            !      dpgm = 1.0E-1*(p_low-p_up)/grav/mwt_air    ! mole (air) cm-2
            !      x_pcol = 1.0E-15* x_nmr * dpgm           ! 1E+15 (Peta) molecules cm-2
            !      x_col = x_col + x_pcol                   ! 1E+15 (Peta) molecules cm-2

            dp = shtj(i,k+1) - shtj(i,k)
            dp = max(dp*ps(i),0.)

            ptop = shtj(i,k)  *ps(i)       ! pressure (Pa) at the upper interface
            pbot = shtj(i,k+1)*ps(i)       ! pressure (Pa) at the bottom interface

            ! Calculate age of air
            !
            ! ---  change in age of air (age in units of years)
            !   --- assuming points with a water vapour mixing ratio
            !         greater than 20 ppmv are in the troposphere (ozone smaller than 150 ppbv)
            !   --- tropospheric points are relaxed back to zero age (here
            !         defined with an offset of 2 years) with a time constant
            !         of 0.1 days

            ! Troposphere
            !       if (zo3lplus (i,k).lt.150.0E-9) then  ! ppb <-- mole /mole vmr
            ! PBL pressure 900 hPa (p_linoz = 0.1000000000000e+05 !pressure level for stratospheric linoz ozone Pa  (100 hPa)
            if (ptop >= 0.900E+05) then !900 hPa 
               zairplus(i,k) = 2.0 + (zairplus(i,k) - 2.0)*delrlx   !years
            else
               ! Stratosphere
               zairplus(i,k) = zairplus(i,k) + delage               !years
            endif

            ! Set lower limit on species as in BIRA (units mole /mole)
            o3_vmr(i,k) = zo3lplus(i,k) *1E-9 * mwt_air/ mwt_o3     !mole /mole vmr <-- micro g/kg air
            o3_vmr(i,k) = max(o3_vmr(i,k), QepsO3)

            ! F-K ozone climatology
            o3fk_vmr(i,k) = zo3fk(i,k) * mwt_air/ mwt_o3             !mole /mole vmr  <-- kg /kg air

            if (ptop <= 1.E+2) then 
               !F-K+HALOE ozone climatology above 1 hPa
               o3c_vmr(i,k) = zo3fk(i,k) * mwt_air/ mwt_o3          !mole /mole vmr  <-- kg /kg
            else
               ! ERA-5 ozone climatology below 1 hPa
               if (minval(zo3ce) >= 0.) &
                    
                    ! 'climato+ghg+linoz+era5+' (Irena I)
                    
                    o3c_vmr(i,k) = zo3ce(i,k)*1E-9 * mwt_air/ mwt_o3      !mole /mole vmr <-- micro g /kg air 

               ! 'climato_phase2' (Paul V)
               !           o3c_vmr(i,k) = zo3ce(i,k) * mwt_air/ mwt_o3         !mole /mole vmr <-- kg /kg air

            end if

            ! Total column F-K ozone climatology F-K (D.U.)
            !
            call linoz_xcol(o3fk_vmr(i,k), pbot, ptop, o3fk_over)
            zo3fkcol(i,k) = o3fk_over      !total column DU

            ! Total column ERA5 ozone climatology (D.U.)
            !
            call linoz_xcol(o3c_vmr(i,k), pbot, ptop, o3c_over)
            zo3ccol(i,k) = o3c_over      !total column DU
            !
            ! Ozone Linoz table of coefficients
            !
            !        c1(i, k) = zlin1(i,k)
            !        c1(i, k) = o3fk_vmr(i,k)        ! F-K Ozone instead (mole /mole)
            c1(i, k) = o3c_vmr(i,k)        ! ERA5 Ozone (mole /mole)

            !        c2(i, k) = zlin2(i,k)
            !        c2(i, k) = ztts(i,k)         ! Temp. climatology GEM (degK)
            c2(i, k) = zttce(i,k)         ! Temp. climatology ERA5 11-years (degK)

            !        c3(i, k) = zlin3(i,k)
            !        c3(i, k) = zo3fkcol(i,k)        ! F-K Ozone Total Column (DU)
            c3(i, k) = zo3ccol(i,k)        ! ERA5 Ozone Total Column (DU)

            c4(i, k) = zlin4(i,k) 
            c5(i, k) = zlin5(i,k)
            c6(i, k) = zlin6(i,k)
            c7(i, k) = zlin7(i,k)

            !        if(mod(i, ni2)==0 .and. mod(trnch, 32)==0)   &
            !        if(o3_vmr(i,k) < QepsO3)                                      &
            !           write(lun_out,*) 'linoz: ', i,trnch, k, &
            !             'o3lplus o3fk o3fkcol c1-c7 =', o3_vmr(i,k), o3fk_vmr(i,k), zo3fkcol(i,k), &
            !             c1(i,k), c2(i,k),c3(i,k),c4(i,k),c5(i,k),c6(i,k),c7(i,k)

            !
            ! GHG Linoz table of coefficients
            !
            IF_LINGHG2: if (llingh) then 
               c8(i, k) = zlin8(i,k)
               c9(i, k) = zlin9(i,k)
               c10(i, k) = zlin10(i,k)
               c11(i, k) = zlin11(i,k)
               c12(i, k) = zlin12(i,k)
               c13(i, k) = zlin13(i,k)
               c14(i, k) = zlin14(i,k)

               ch4_vmr(i,k) = zch4lplus(i,k)                    !mole /mole vmr
               n2o_vmr(i,k) = zn2olplus(i,k)                    !mole /mole vmr
               f11_vmr(i,k) = zf11lplus(i,k)                    !mole /mole vmr
               f12_vmr(i,k) = zf12lplus(i,k)                    !mole /mole vmr

               !        if(ch4_vmr(i,k) < Qeps)                                      &
               !           write(lun_out,*) 'linoz zero check: ', i,trnch, k, &
               !             'ch4lplus=', ch4_vmr(i,k)

               !        if(n2o_vmr(i,k) < Qeps)                                      &
               !           write(lun_out,*) 'linoz zero check: ', i,trnch, k, &
               !             'n2olplus=', n2o_vmr(i,k)

               !        if(f11_vmr(i,k) < Qeps)                                      &
               !           write(lun_out,*) 'linoz zero check: ', i,trnch, k, &
               !             'f11lplus=', f11_vmr(i,k)

               !        if(f12_vmr(i,k) < Qeps)                                      &
               !           write(lun_out,*) 'linoz zero check: ', i,trnch, k, &
               !             'f12lplus=', f12_vmr(i,k)

               ! Set lower limit on species as in BIRA (units mole /mole)

               ch4_vmr(i,k) = max(ch4_vmr(i,k) ,Qeps)
               n2o_vmr(i,k) = max(n2o_vmr(i,k) ,Qeps)
               f11_vmr(i,k) = max(f11_vmr(i,k) ,Qeps)
               f12_vmr(i,k) = max(f12_vmr(i,k) ,Qeps)

               !        if(mod(i, ni2)==0 .and. mod(trnch, 32)==0)              &
               !           write(lun_out,*)                                     &
               !              'linoz: ', i,trnch, k,                            &
               !              'ch4, n2o, f11, f12 c8-c14='                  , &
               !               ch4_vmr(i,k), n2o_vmr(i,k)                   , &
               !               f11_vmr(i,k), f12_vmr(i,k)                   , &
               !               c8(i,k), c9(i,k),c10(i,k),c11(i,k),c12(i,k),c13(i,k),c14(i,k)

            end if IF_LINGHG2

         enddo DO_K1

         ! 2D Total Column 
         zo3fktc(i) = zo3fkcol(i,nkm1) !FK
         zo3ctc(i) = zo3ccol(i,nkm1) !ERA5+FK merge

      enddo DO_I1

      !    Linoz tendencies 

      call linoz_tend( &
           o3_vmr ,ch4_vmr,n2o_vmr            , & !input, mole/mole vmr
!!!!                 f11_vmr,f12_vmr,o3fk_vmr            , & !input, mole/mole vmr F-K ozone climato in troposphere
           f11_vmr,f12_vmr,o3c_vmr            , & !input, mole/mole vmr ERA-3 ozone climato in troposphere
           tt, ps, shtj, qq                   , & !input
           c1 ,c2, c3 , c4 , c5 , c6 ,c7      , & !input
           c8 ,c9 ,c10, c11, c12, c13,c14     , & !input
           o3new  ,ch4new ,n2onew             , & !output, mole /mole 
           f11new ,f12new                     , & !output, mole /mole 
           !                     zo1chmtd,  zo4chmtd,  zo6chmtd, zo7chmtd, & !output, mole /mole /sec 
           !                     zo3chmtd, zch4chmtd, zn2ochmtd          , & !output, mole /mole /sec 
           !                     zf11chmtd,zf12chmtd                     , & !output, mole /mole /sec 
           !                     zo3col                                  , &
           o1_tend, o4_tend, o6_tend, o7_tend, &
           o3_tend, ch4_tend, n2o_tend, f11_tend, f12_tend, &
           o3_over, &
           dt, kount, trnch, ni, nkm1, nk)          !input

      ! Ozone Prognostic LINOZ

      DO_I2: do i =1,ni

         do k = 1, nkm1    !1=TOA, nk-1=SFC

            !        if(local_dbg  .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
            !           write(lun_out,*) 'linoz: ', i,trnch,k, &
            !              ' o3vmr,  o3new,  o3chmtd=',  &
            !                o3_vmr(i,k), o3new(i,k), v(o3chmtd + (k-1)*n + i-1)
            !                o3_vmr(i,k),  o3new,  o3_tend  !, &

            ! Update 'plus' values here instead in 'apply_tendencies1'
            o3new(i,k) = max(o3new(i,k), QepsO3)             ! mole /mole vmr
            zo3lplus(i,k) = o3new(i,k) *1E+9 * mwt_o3 /mwt_air  ! micro g /kg air <-- mole /mole vmr

            !        zo1chmtd(i,k) = o1_tend(i,k)
            !        zo4chmtd(i,k) = o4_tend(i,k)
            !        zo6chmtd(i,k) = o6_tend(i,k)
            !        zo7chmtd(i,k) = o7_tend(i,k)

            zo3chmtd(i,k) = o3_tend(i,k)
            zo3col(i,k) = o3_over(i,k)

         end do ! k-loop

         ! 2D Total Column
         zo3tc (i) =zo3col (i,nkm1) !LINOZ

      end do DO_I2

      ! GHG Prognostic LINOZ

      IF_LINGHG3: if (llingh) then

         do i = 1,ni

            ch4_over = 0.
            n2o_over = 0.
            f11_over = 0.
            f12_over = 0.

            do k = 1, nkm1    !1=TOA, nk-1=SFC

               ptop = shtj(i,k)  *ps(i)       ! pressure (Pa) at the upper interface
               pbot = shtj(i,k+1)*ps(i)       ! pressure (Pa) at the bottom interface

               !    Total column prognostic ch4 tracer linoz
               !
               call ghg_xcol(ch4_vmr(i,k), pbot, ptop, ch4_over)
               zch4col(i,k) = ch4_over       !molec cm-2

               !    Total column prognostic n2o tracer linoz
               !
               call ghg_xcol(n2o_vmr(i,k), pbot, ptop, n2o_over)
               zn2ocol(i,k) = n2o_over       !molec cm-2

               !    Total column prognostic f11 tracer linoz
               !
               call ghg_xcol(f11_vmr(i,k), pbot, ptop, f11_over)
               zf11col(i,k) = f11_over       !molec cm-2

               !    Total column prognostic f12 tracer linoz
               !
               call ghg_xcol(f12_vmr(i,k), pbot, ptop, f12_over)
               zf12col(i,k) = f12_over       !molec cm-2

               !        if(local_dbg  .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
               !           write(lun_out,*) 'linoz: ', i,trnch,k, &
               !              'ch4vmr, ch4new, ch4chmtd=', ch4_vmr(i,k), ch4new(i,k), ch4_tend, &
               !              'n2ovmr, n2onew, n2ochmtd=', n2o_vmr(i,k), n2onew(i,k), n2o_tend, &
               !              'f11vmr, f11new, f11tend=', f11_vmr(i,k), f11new(i,k), f11_tend, &
               !              'f12vmr, f12new, f12tend=', f12_vmr(i,k), f12new(i,k), f12_tend

               ! Update 'plus' values here instead in 'apply_tendencies1'
               !        zch4lplus(i,k) = max(zch4lplus(i,k) ,Qeps)
               !        zf11lplus(i,k) = max(zf11lplus(i,k) ,Qeps)
               !        zf12lplus(i,k) = max(zf12lplus(i,k) ,Qeps)
               !        zn2olplus(i,k) = max(zn2olplus(i,k) ,Qeps)


               ! Update 'plus' values here instead in 'apply_tendencies1'
               zch4lplus(i,k) = ch4new(i,k)
               zf11lplus(i,k) = f11new(i,k)
               zf12lplus(i,k) = f12new(i,k)
               zn2olplus(i,k) = n2onew(i,k)

               zch4chmtd(i,k) = ch4_tend(i,k)
               zn2ochmtd(i,k) = n2o_tend(i,k)
               zf11chmtd(i,k) = f11_tend(i,k)
               zf12chmtd(i,k) = f12_tend(i,k)

            end do !k-loop 

            ! 2D Total Column
            zch4tc(i) =zch4col(i,nkm1)
            zn2otc(i) =zn2ocol(i,nkm1)
            zf11tc(i) =zf11col(i,nkm1)
            zf12tc(i) =zf12col(i,nkm1)

         end do    !i-loop

      end if IF_LINGHG3

      ! Update 'plus' values above instead in 'apply_tendencies1'
      !   call apply_tendencies1(d,dsiz,v,vsiz,f,fsiz,o3lplus ,o3chmtd ,n,nkm1)
      !   call apply_tendencies1(d,dsiz,v,vsiz,f,fsiz,ch4lplus,ch4chmtd,n,nkm1)
      !   call apply_tendencies1(d,dsiz,v,vsiz,f,fsiz,n2olplus,n2ochmtd,n,nkm1)
      !   call apply_tendencies1(d,dsiz,v,vsiz,f,fsiz,f11lplus,f11chmtd,n,nkm1)
      !   call apply_tendencies1(d,dsiz,v,vsiz,f,fsiz,f12lplus,f12chmtd,n,nkm1)


      call msg_toall(MSG_DEBUG, 'linoz [END]')
      !----------------------------------------------------------------
      return
   end subroutine linoz3

end module linoz
