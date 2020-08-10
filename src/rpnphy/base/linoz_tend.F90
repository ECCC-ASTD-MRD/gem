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
subroutine linoz_tend(o3, ch4, n2o, f11, f12, & 
     colo3,                             &
     temp, ps, shtj, qq,                &
     o3c, c2, c3,  c4,  c5,  c6,  c7,    & 
     c8, c9, c10, c11,                  &
     o3_new, &
     ch4_new, n2o_new, f11_new, f12_new,&
     do1dt,  do4dt,  do6dt, do7dt,      &
     do3dt, &
     dch4dt, dn2odt, df11dt, df12dt,    &
     timestep, ni, nkm1, nk)
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
   !      o3c: o3 mole /mole vmr FK-ERA5 climatology
   !      c2 : temp climatology ERA5
   !      c3 : colo3 climatology
   !      c4 : P-L
   !      c5 : d(P-L)/dO3
   !      c6 : d(P-L)/dT
   !      c7 : d(P-L)/dcolo3
   !      c8  : n2o loss frequency (s-1)
   !      c9  : ch4 loss frequency (s-1)
   !      c10 : F11 loss frequency (s-1)
   !      c11 : F12 loss frequency (s-1)
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

   integer, intent(in)             :: ni,nkm1,nk
   real,    intent(in)             :: timestep
   real, dimension(ni, nkm1), intent(in)    :: o3, n2o,ch4,f11,f12
   real, dimension(ni, nkm1), intent(in)    :: o3c,qq
   real, dimension(ni, nkm1), intent(in)    :: colo3
   real, dimension(ni),       intent(in)    :: ps
   real, dimension(ni, nk),   intent(in)    :: shtj,temp
   real, dimension(ni, nkm1), intent(in)    :: c4,c2,c3,c5,c6,c7
   real, dimension(ni, nkm1), intent(in)    :: c8,c9,c10,c11
   real, dimension(ni, nkm1), intent(out)   :: o3_new, n2o_new, ch4_new, f11_new, f12_new
   real, dimension(ni, nkm1), intent(out)   :: do1dt,do4dt,do6dt,do7dt
   real, dimension(ni, nkm1), intent(out)   :: do3dt,dn2odt,dch4dt,df11dt,df12dt

   !@author J. de Grandpre (ARQI): February 2013
   !@revisions
   ! * I. Ivanova (2015) - Complete remake
   !*@/

   !  Declaration of local variables.

   integer :: i,k, ni2
   real :: ptop, hu_ppm
   real :: tau,fss, c50

   ni2 = int(float(ni)/2.)

   ! --------------------
   !  Loop on longitudes
   ! -------------------
   DO_I1: do i = 1, ni

      ! ----------------------------
      ! Loop on all vertical levels
      ! ----------------------------

      DO_K1: do k=1,nkm1

         hu_ppm = consth * qq(i,k)      ! units ppmv <- kg kg-1
         ptop = shtj(i,k)  *ps(i)       ! pressure (Pa) at the upper interface

         !     if(local_dbg .and. mod(i, ni2)==0 ) &
         !      write(lun_out,*) 'linoz_tend: ', i, k, timestep,          , &
         !              'o3vmr,  o3cvmr='                                             , &
         !               o3(i,k), o3c(i,k)      !                                      , &
         !              'c1, c2, c3, c4, c5, c6, c7='                                 , &
         !               o3c(i,k),c2(i,k), c3(i,k), c4(i,k), c5(i,k), c6(i,k), c7(i,k), &
         !              'c8, c9, c10, c11 ='                           , &
         !               c8(i,k), c9(i,k), c10(i,k), c11(i,k),  &
         !              'temp, ptop =', &
         !               temp(i,k), ptop , & 
         !               'ch4, n2o, f11, f12, =', &
         !               ch4(i,k), n2o(i,k), f11(i,k), f12(i,k)

         !     if(local_dbg .and. mod(i, ni2)==0 ) &
         !        write(lun_out,*) 'linoz_tend: ', i, k             , &


         !  Fix for preventing division by 0 originating from horizontal
         !  interpolation (JDG,march 2010)

         if (c5(i,k) <= 0.) then 
            c50 = min(c5(i,k), -1.E-19)
         else  
            c50 = max(c5(i,k), 1.E-19)
         end if

         tau = -1./c50

         !  Tendency


         ! Apply LINOZ tendencies above tropopause 100hPa&10ppmv and below stratopause 1hPa
         IF_PTOPLIN: if (ptop < p_linoz_tropo .or. hu_ppm < hu_linoz_tropo) then

            ! Relaxation O3  to FK climatology above 1 hPa, tau_linoz_meso=6 hours 
            if (ptop < p_linoz_meso ) then                  !1hPa

               ! Relaxation O3  to FK climatology above 1 hPa, tau_linoz_meso=6 hours 
               !
               !    Ozone relaxation toward FK-HALOE climatology in the mesosphere
               !    Relaxation time set to 6 days above 1 hPa

               o3_new(i,k) = o3(i,k) - (o3(i,k) - o3c(i,k)) * timestep/tau_linoz_meso   !mole /mole vmr

            else

               ! Apply LINOZ tendencies only above tropopause , defined as specific humidity 'hu_linoz_tropo=10ppmv' or pressure 'p_linoz_tropo=100mb'
               do4dt(i,k) = c4(i,k)                      /o3(i,k)   !c4 : P-L            cm-3 sec-1 
               do1dt(i,k) = c5(i,k)*(o3(i,k)   -o3c(i,k))/o3(i,k)   !c5 : d(P-L)/dO3          sec-1 
               do6dt(i,k) = c6(i,k)*(temp(i,k) -c2(i,k)) /o3(i,k)   !c6 : d(P-L)/dT      cm-3 sec-1 K-1
               do7dt(i,k) = c7(i,k)*(colo3(i,k)-c3(i,k)) /o3(i,k)   !c7 : d(P-L)/dcolo3  cm-3 sec-1 DU-1 
 
               fss = o3c(i,k) + (c4(i,k)                     + &
                                c6(i,k)*(temp(i,k)-c2(i,k)) + &
                                c7(i,k)*(colo3(i,k)-c3(i,k)))*tau                !mole/mole vmr

            ! Utilisation de la formulation analytique 
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

         else ! IF_PTOPLIN -- below tropopause p_linoz_tropo & hu_linoz_tropo

            o3_new(i,k) = o3(i,k)                                                  !mole /mole vmr

            if (ptop >= p_linoz_clim) then     ! (below 400mb)

               !    Ozone relaxation toward the Fortuin & Kelder climatology in the troposphere
               !    Relaxation time set to 7 days below 100 hPa which represent the characteristic
               !    vertical mixing time in the troposphere. It is long enough to allow 
               !    stratospheric intrusion of ozone (specific to LINOZ-2)
               !              
               !    Relaxation O3  to FK climatology      below 400 hPa (tau_linoz_tropo=2 days)  

               o3_new(i,k) = o3(i,k) - (o3(i,k) - o3c(i,k)) * timestep/tau_linoz_tropo  !mole /mole vmr

            end if ! p_linoz_clim

         end if IF_PTOPLIN

         ! Set lower limit on species as in BIRA (mole /mole)
         o3_new(i,k) = max(o3_new(i,k), QepsO3)                ! mole /mole vmr

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

            ! Apply LINOZ tendencies above tropopause, defined as specific humidity 'hu_linoz=10ppmv' or pressure 'p_linoz_tropo=100mb'
            hu_ppm = consth * qq(i,k)      ! units ppmv <- kg kg-1
            ptop = shtj(i,k)  *ps(i)       ! pressure (Pa) at the upper interface

            ! 100 Pa 10 ppmv
            IF_PTOPLIN2: if (ptop < p_linoz_tropo .or. hu_ppm < hu_linoz_tropo) then

               ch4_new(i,k) = ch4(i,k) * exp( -c9(i,k)*timestep)                  !mole /mole vmr
               f11_new(i,k) = f11(i,k) * exp(-c10(i,k)*timestep)                  !mole /mole vmr
               f12_new(i,k) = f12(i,k) * exp(-c11(i,k)*timestep)                  !mole /mole vmr
               n2o_new(i,k) = n2o(i,k) * exp( -c8(i,k)*timestep)                  !mole /mole vmr

            else ! IF_PTOPLIN2 -- below tropopause (p_linoz_tropo & hu_linoz_tropo) GHG

               ch4_new(i,k) = ch4(i,k)
               n2o_new(i,k) = n2o(i,k)
               f11_new(i,k) = f11(i,k)
               f12_new(i,k) = f12(i,k)

            end if IF_PTOPLIN2

            ! Set lower limit on species as in BIRA (units mole /mole)
            ch4_new(i,k)= max(ch4_new(i,k),Qeps)
            n2o_new(i,k)= max(n2o_new(i,k),Qeps)
            f11_new(i,k)= max(f11_new(i,k),Qeps)
            f12_new(i,k)= max(f12_new(i,k),Qeps)


            ! Set lower boundary values as in BIRA-Tropo (units mole /mole)
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

            !     if(local_dbg  .and. mod(i, ni2)==0) &
            !       write(lun_out,*) 'linoz_tend: ', i,k, &
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
