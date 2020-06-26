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
   
!/@*
subroutine linoz_tend(o3, ch4, n2o, f11, f12, o3c, &
     temp, ps, shtj, qq,                &
     c1, c2, c3,  c4,  c5,  c6,  c7,    & 
     c8, c9, c10, c11, c12, c13, c14,   &
     o3_new, ch4_new, n2o_new,          &
     f11_new, f12_new,                  &
     do1dt,  do4dt,  do6dt, do7dt,      &
     do3dt,  dch4dt, dn2odt,            &
     df11dt, df12dt,                    &
     colo3,                             &
     timestep, kount, trnch, ni, nkm1, nk)
   use phy_options
   use linoz_mod
   implicit none
#include <arch_specific.hf>
   !@object Stratospheric Ozone chemistry
   ! Produces the ozone tendency due to stratospheric
   ! photochemistry (based on LINOZ scheme)
   !@arguments
   ! ni       horizonal index
   ! nkm1     vertical  index nk-1
   ! nk       vertical  index
   ! timestep timestep
   ! trnch    index of vertical slice n*nk
   ! kount    number of timesteps
   ! v        volatile bus
   ! f        permanent bus
   ! d        dynamics bus
   !----------------------------------------------------
   !   Input:
   !   -----
   !      o3 : o3 mole /mole vmr
   !      ch4 : ch4 mole /mole vmr
   !      n2o : n2o mole /mole vmr
   !      f11 : f11 mole /mole vmr
   !      f12 : f12 mole /mole vmr
   !      o3c: o3 mole /mole vmr F-K climatology
   !      c1 : o3 linoz climatology
   !      c2 : temp
   !      c3 : colo3
   !      c4 : P-L
   !      c5 : d(P-L)/dO3
   !      c6 : d(P-L)/dT
   !      c7 : d(P-L)/dcolo3
   !      c8  : n2o loss frequency (s-1)
   !      c9  : ch4 loss frequency (s-1)
   !      c10 : F11 loss frequency (s-1)
   !      c11 : F12 loss frequency (s-1)
   !      c12 : noy loss frequency (s-1)
   !      c13 : noy prod ratio
   !      c14 : noy loss ratio
   !
   !   Output:
   !   ------
   !      o3_new : o3 (t+dt) mole /mole vmr
   !      ch4_new:
   !      n2o_new:
   !      f11_new:
   !      f12_new:
   !      do3dt  : o3 linoz tendency
   !      dch4dt :
   !      dn2odt :
   !      df11dt :
   !      df12dt :
   !      colo3  : total column DU
   !------------------------------------------------------

   integer, intent(in)             :: kount,trnch,ni,nkm1,nk
   real,    intent(in)             :: timestep
   real, dimension(ni, nkm1), intent(in)    :: o3, n2o,ch4,f11,f12
   real, dimension(ni, nkm1), intent(in)    :: o3c,temp,qq
   real, dimension(ni),       intent(in)    :: ps
   real, dimension(ni, nk),   intent(in)    :: shtj
   real, dimension(ni, nkm1), intent(inout) :: c1,c4
   real, dimension(ni, nkm1), intent(in)    :: c2,c3,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14
   real, dimension(ni, nkm1), intent(out)   :: o3_new, n2o_new, ch4_new, f11_new, f12_new
   real, dimension(ni, nkm1), intent(out)   :: do1dt,do4dt,do6dt,do7dt
   real, dimension(ni, nkm1), intent(out)   :: do3dt,dn2odt,dch4dt,df11dt,df12dt
   real, dimension(ni, nkm1), intent(out)   :: colo3 

   !@author J. de Grandpre (ARQI): February 2013
   !@revisions
   ! * I. Ivanova (2015) - Complete remake
   !*@/

   !  Declaration of local variables.

   integer :: i,k, ni2
   real :: o3_over, ptop, pbot, hu_ppm
   real :: tau,fss, c50, p_linoz_c4

   ni2 = int(float(ni)/2.)

   ! --------------------
   !  Loop on longitudes
   ! -------------------
   DO_I1: do i = 1, ni

      o3_over = 0.

      ! ----------------------------
      ! Loop on all vertical levels
      ! ----------------------------

      DO_K1: do k=1,nkm1

         hu_ppm = consth * qq(i,k)      ! units ppmv <- kg kg-1
         ptop = shtj(i,k)  *ps(i)       ! pressure (Pa) at the upper interface
         pbot = shtj(i,k+1)*ps(i)       ! pressure (Pa) at the bottom interface

         !     if(local_dbg .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
         !      write(lun_out,*) 'linoz_tend: ', i, trnch, k, timestep, kount         , &
         !              'o3vmr,  o3cvmr='                                             , &
         !               o3(i,k), o3c(i,k)      !                                      , &
         !              'c1, c2, c3, c4, c5, c6, c7='                                 , &
         !               c1(i,k), c2(i,k), c3(i,k), c4(i,k), c5(i,k), c6(i,k), c7(i,k), &
         !              'c8, c9, c10, c11, c12, c13, c14 ='                           , &
         !               c8(i,k), c9(i,k), c10(i,k), c11(i,k), c12(i,k), c13(i,k), c14(i,k), &
         !              'temp, pbot, ptop =', &
         !               temp(i,k), pbot,ptop , & 
         !               'ch4, n2o, f11, f12, =', &
         !               ch4(i,k), n2o(i,k), f11(i,k), f12(i,k)

         ! Total column ozone tracer linoz (D.U.)

         call linoz_xcol(o3(i,k), pbot, ptop, o3_over)
         colo3(i,k) = o3_over

         !     if(local_dbg .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
         !        write(lun_out,*) 'linoz_tend: ', i, trnch, k             , &


         !  Fix for preventing division by 0 originating from horizontal
         !  interpolation (JDG,march 2010)

         if (c5(i,k) <= 0.) then 
            c50 = min(c5(i,k), -1.E-19)
         else  
            c50 = max(c5(i,k), 1.E-19)
         end if

         ! Replace at all levels linoz O3 climatology with F-K O3 climatology
         c1(i,k) = o3c(i,k)
         tau = -1./c50
         !     p_linoz_c4 = p_linoz*1.E-2     ! 1hPa (formerly 10hPa)
         p_linoz_c4 = p_linoz*1.E-1     ! 10hPa (formerly 1hPa)


         !  Tendency


         ! Apply LINOZ tendencies only above tropopause, defined as specific humidity 'hu_linoz=10ppmv' or pressure 'p_linoz=100mb'
         !
         ! 1) Apply LINOZ tendencies only above pressure p_linoz=100mb, or 
         ! 2) Apply LINOZ tendencies only above tropopause, defined as specific humidity hu_linoz=10ppmv

         IF_PTOPLIN: if (ptop < p_linoz .or. hu_ppm < hu_linoz) then    ! (pressure less than 100 Pa .or. humid less than 10 ppmv)

            ! Ignore P-L term on RHS above 10hPa
            if (ptop < p_linoz_c4) then !10 hPa
               c4(i,k) = 0.
            end if

            ! Apply LINOZ tendencies below the stratopause 1hPa
            if (ptop < p_linoz*1.E-2 ) then                  !1hPa

            ! Reset O3 to FK-HALOE climatology above stratopause 1 hPa (tau_linoz=2 days)

               o3_new (i,k) = o3c(i,k)                                            !mole /mole vmr

            else

               do4dt(i,k) = c4(i,k)                     /o3(i,k)   !c4 : P-L            cm-3 sec-1 
               do1dt(i,k) = c5(i,k)*(o3(i,k) - c1(i,k) )/o3(i,k)   !c5 : d(P-L)/dO3          sec-1 
               do6dt(i,k) = c6(i,k)*(temp(i,k)-c2(i,k)) /o3(i,k)   !c6 : d(P-L)/dT      cm-3 sec-1 K-1
               do7dt(i,k) = c7(i,k)*(colo3(i,k)-c3(i,k))/o3(i,k)   !c7 : d(P-L)/dcolo3  cm-3 sec-1 DU-1 
 
               fss = c1(i,k) + (c4(i,k) + c6(i,k)*(temp(i,k)-c2(i,k)) + c7(i,k)*(colo3(i,k)-c3(i,k)))*tau                !mole/mole vmr

            ! Formulation explicite (not used)
            !
            !      o3l= (fss - o3l)/tau*timestep + o3l
            !
            ! Utilisation de la formulation analytique  (used here)
            !
            !     (aout 2012; jdg)
            !     df/dt = - c(f - fss)
            !     df/(f-fss) = -cdt
            !     ln{(fnew - fss)/(fold-fss)} = -cdt
            !     fnew = fss + (fold-fss)exp(-cdt)
            !     df = fnew - fold = (fss-fold)(1-exp(-cdt))
            !     fnew = fold + df = fold + (fss-fold)(1-exp(-cdt))
            !
               o3_new(i,k) = o3(i,k) + (fss-o3(i,k))*(1.-exp(-timestep/tau))          !mole /mole vmr

            end if

         else ! IF_PTOPLIN -- below tropopause p_linoz & hu_linoz

            o3_new(i,k) = o3(i,k)                                                  !mole /mole vmr

            if (ptop >= ptop_clim) then     ! (below 400mb)

               !    Ozone relaxation toward the Fortuin & Kelder climatology in the troposphere
               !    Relaxation time set to 7 days below 100 hPa which represent the characteristic
               !    vertical mixing time in the troposphere. It is long enough to allow 
               !    stratospheric intrusion of ozone (specific to LINOZ-2)
               !              
               !    Relaxation O3  to FK climatology      below 400 hPa (tau_linoz=2 days)  

               o3_new(i,k) = o3(i,k) - (o3(i,k) - o3c(i,k)) * timestep/tau_linoz   !mole /mole vmr

            end if ! ptop_clim

         end if IF_PTOPLIN

         ! Ozone tendency LINOZ (mole/mole/sec)

         do3dt(i,k) = (o3_new(i,k) - o3(i,k)) / timestep               !mole/mole/sec


      end do DO_K1
   end do DO_I1
   

   IF_LINGHG: if (llingh) then

      ! --------------------
      !  Loop on longitudes
      ! -------------------
      DO_I2: do i = 1, ni

         ! ----------------------------
         ! Loop on all vertical levels
         ! ----------------------------

         DO_K2: do k=1,nkm1

            ! Apply LINOZ tendencies only above tropopause, defined as specific humidity 'hu_linoz=10ppmv' or pressure 'p_linoz=100mb'
            !
            ! 1) Apply LINOZ tendencies only above pressure p_linoz=100mb, or 
            ! 2) Apply LINOZ tendencies only above tropopause, defined as specific humidity hu_linoz=10ppmv
            !
            hu_ppm = consth * qq(i,k)      ! units ppmv <- kg kg-1
            ptop = shtj(i,k)  *ps(i)       ! pressure (Pa) at the upper interface
            pbot = shtj(i,k+1)*ps(i)       ! pressure (Pa) at the bottom interface

            IF_PTOPLIN2: if (ptop < p_linoz .or. hu_ppm < hu_linoz) then    ! (pressure less than 100 Pa .or. humid less than 10 ppmv)

               ! Autres constituents : fss=0

               ! Formulation explicite (not used)
               !
               !      o3l= (fss - o3l)/tau*timestep + o3l
               !      n2oL = n2oL * ( 1 - timestep*c8)
               !      ch4L = ch4L * ( 1 - timestep*c9)
               !      f11L = f11L * ( 1 - timestep*c10)
               !      f12L = f12L * ( 1 - timestep*c11)

               ch4_new(i,k) = ch4(i,k) * exp( -c9(i,k)*timestep)                  !mole /mole vmr
               f11_new(i,k) = f11(i,k) * exp(-c10(i,k)*timestep)                  !mole /mole vmr
               f12_new(i,k) = f12(i,k) * exp(-c11(i,k)*timestep)                  !mole /mole vmr
               n2o_new(i,k) = n2o(i,k) * exp( -c8(i,k)*timestep)                  !mole /mole vmr

               ! N2O-NOY
               ! dn2o: puit de n2o (valeurs positives)
               ! commente temporairement: division par 0
               ! quand le champ est initialise constant a 1ppbv
               !
               !     dn2o=n2o(i,k)*(1.-exp(-c8(i,k)*timestep))                     !mole /mole vmr
               !     n2o_new(i,k) = n2o(i,k)-dn2o                                  !mole /mole vmr
               !     n2o_new(i,k) = n2o(i,k)                                       !mole /mole vmr

               ! set limits to selected NOY parameters
               ! c12 is NOY photodissociation parameters and c14 is NOY loss ratio
               !      c12(i,k) = max( c12(i,k),1.E-15)
               !      c14(i,k) = max( c14(i,k),1.)

               !      pnoy=dn2o*c13(i,k)
               !      lnoy=2.*c12(i,k)*c14(i,k)*noy(i,k)/(1.+c14(i,k)*noy(i,k))
               !      noyss=pnoy/lnoy
               !
               !      dnoy=(noyss-noy)*(1.-exp(-lnoy*timestep))
               !
               !      noy_new(i,k) = noy(i,k) + dnoy                               !mole /mole vmr

            else ! IF_PTOPLIN2 -- below tropopause (p_linoz & hu_linoz) GHG

               ch4_new(i,k) = ch4(i,k)
               n2o_new(i,k) = n2o(i,k)
               f11_new(i,k) = f11(i,k)
               f12_new(i,k) = f12(i,k)

               !   Relaxation to background ghg not tested!!!
               !        ch4_new(i,k) = ch4(i,k) - (ch4(i,k) - qch4  *1.E-6) * timestep/tau_linoz      !mole /mole vmr
               !        n2o_new(i,k) = n2o(i,k) - (n2o(i,k) - qn2o  *1.E-6) * timestep/tau_linoz      !mole /mole vmr
               !        f11_new(i,k) = f11(i,k) - (f11(i,k) - qcfc11*1.E-9) * timestep/tau_linoz      !mole /mole vmr
               !        f12_new(i,k) = f12(i,k) - (f12(i,k) - qcfc12*1.E-9) * timestep/tau_linoz      !mole /mole vmr

            end if IF_PTOPLIN2


            !     Set lower boundary values as in BIRA-Tropo (units mole /mole)
 
            if (k >= nkm1-2) then  !nk-1 =SFC
               ch4_new(i,k) = 1.76E-06  
               n2o_new(i,k) = 3.22E-07
               f11_new(i,k) = 2.60E-10
               f12_new(i,k) = 5.44E-10
            end if

            ! GHG tendency LINOZ (mole/mole/sec)

            dch4dt(i,k) = (ch4_new(i,k) - ch4(i,k)) / timestep               !mole/mole/sec
            df11dt(i,k) = (f11_new(i,k) - f11(i,k)) / timestep               !mole/mole/sec
            df12dt(i,k) = (f12_new(i,k) - f12(i,k)) / timestep               !mole/mole/sec 
            dn2odt(i,k) = (n2o_new(i,k) - n2o(i,k)) / timestep               !mole/mole/sec 

            !     if(local_dbg  .and. mod(i, ni2)==0 .and. mod(trnch, 32)==0) &
            !       write(lun_out,*) 'linoz_tend: ', i,trnch, k, &
            !        ' o3vmr,  o3new,  o3tend= ',   o3(i,k),  o3_new(i,k),  do3dt(i,k), &
            !        ' ch4vmr, ch4new, ch4tend=',  ch4(i,k), ch4_new(i,k), dch4dt(i,k) !, &
            !        ' n2ovmr, n2onew, n2otend=',  n2o(i,k), n2o_new(i,k), dn2odt(i,k), &
            !        ' f11vmr, f11new, f11tend=',  f11(i,k), f11_new(i,k), df11dt(i,k), &
            !        ' f12vmr, f12new, f12tend=',  f12(i,k), f12_new(i,k), df12dt(i,k)

         end do DO_K2
      end do DO_I2

   end if IF_LINGHG

   !----------------------------------------------------------------
   return
end subroutine linoz_tend
