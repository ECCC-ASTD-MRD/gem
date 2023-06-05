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
!-------------------------------------- LICENCE END ----------------------------
 
module pbl_height
   implicit none
   private
   public :: pbl_height1

contains

   !/@*
   function pbl_height1(hpbl,temp,hum,uwind,vwind,hghtm,hghtt,sigt,temp_s,hum_s,psfc, &
        z0m,z0t,lat,fcor,n,nk) result(status)
      use tdpack_const, only: CAPPA, CPD, EPS1, GRAV
      use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      !@Arguments
      integer, intent(in) :: n                          !horizontal dimension
      integer, intent(in) :: nk                         !vertical dimension
      real, dimension(n), intent(in) :: temp_s          !surface temperature (K)
      real, dimension(n), intent(in) :: hum_s           !surface specific humidity (kg/kg)
      real, dimension(n), intent(in) :: psfc            !surface pressure (Pa)
      real, dimension(n), intent(in) :: z0m             !momentum mixing length (m)
      real, dimension(n), intent(in) :: z0t             !thermodynamic mixing length (m)
      real, dimension(n), intent(in) :: lat             !latitude (rad)
      real, dimension(n), intent(in) :: fcor            !Coriolis factor (/s)
      real, dimension(n,nk), intent(in) :: temp         !dry bulb temperature (K)
      real, dimension(n,nk), intent(in) :: hum          !specific humidity (kg/kg)
      real, dimension(n,nk), intent(in) :: uwind        !x-component wind (m/s)
      real, dimension(n,nk), intent(in) :: vwind        !y-component wind (m/s)
      real, dimension(n,nk), intent(in) :: hghtm        !momentum-level level heights (m)
      real, dimension(n,nk), intent(in) :: hghtt        !thermodynamic-level level heights (m)
      real, dimension(n,nk), intent(in) :: sigt         !sigma of thermodynamic levels
      real, dimension(n), intent(out) :: hpbl           !PBL depth (m)
      integer :: status                                 !completion status (RMN_OK or RMN_ERR)

      !@Author R. McTaggart-Cowan (Fall 2016)
      !@Object
      !          Calculate the depth of the PBL using the approach proposed
      !          by Seidel et al. (2012; JGR-Atm) based on the work of 
      !          Vogelezang and Holtslag (1996; BLM).  Variable names follow
      !          page 50 of IFS documentation (Chapter 3: Turbulent transport
      !          and interactions with the surface) for Cy40r1.
      !*@/

      ! Local parameter definitions
      real, parameter :: HGHT_NEARSFC=2.,RIB_CRIT=0.25,BETA=100.

      ! Local variable declarations
      integer :: i,k
      real :: svhbl,rib,rib_below,shrsq,hght,hght_below,zp
      real, dimension(n) :: theta_air,spd_air,dir_air,temp_n,hum_n,svn,uwind_n,vwind_n,ustar

      ! Initialize status
      status = RMN_ERR

      ! Compute basic near-surface thermodynamic properties
      if (sl_prelim(temp(:,nk),hum(:,nk),uwind(:,nk),vwind(:,nk),psfc,hghtm(:,nk), &
           spd_air=spd_air,dir_air=dir_air) /= SL_OK) then
         print*, 'Error in pbl_height for surface layer preparation'
         return
      endif
      theta_air = temp(:,nk)*sigt(:,nk)**(-CAPPA)
      if (sl_sfclayer(theta_air,hum(:,nk),spd_air,dir_air,hghtm(:,nk),hghtt(:,nk),temp_s,hum_s, &
           z0m,z0t,lat,fcor,hghtt_diag=HGHT_NEARSFC,hghtm_diag=HGHT_NEARSFC, &
           t_diag=temp_n,q_diag=hum_n,u_diag=uwind_n,v_diag=vwind_n,ue=ustar) /= SL_OK) then
         print*, 'Error in pbl_height for surface layer calculations'
         return
      endif

      ! Compute near-surface dry static energy
      svn = CPD*temp_n*(1+EPS1*hum_n) + GRAV*HGHT_NEARSFC

      ! Scan bulk Richardson numbers until the critical value is reached
      do i=1,n
         rib = 0.
         hght = HGHT_NEARSFC
         k = nk
         do while (rib < RIB_CRIT .and. k > 0)
            if (hghtt(i,k) <= HGHT_NEARSFC) then
               k = k-1
               cycle
            endif
            rib_below = rib
            hght_below = hght
            svhbl = CPD*temp(i,k)*(1+EPS1*hum(i,k)) + GRAV*hghtt(i,k)
            shrsq = (uwind(i,k)-uwind_n(i))**2+(vwind(i,k)-vwind_n(i))**2
            hght = 0.5*(hghtm(i,k)+hghtt(i,k))
            zp = (hghtm(i,k)-HGHT_NEARSFC)**2/(hghtt(i,k)-HGHT_NEARSFC)
            rib = zp*2*GRAV*(svhbl-svn(i)) / &
                 ((svhbl+svn(i)-GRAV*hghtt(i,k)-GRAV*HGHT_NEARSFC)*(shrsq + BETA*ustar(i)**2))
            k = k-1
         enddo
         if (k == 0) then
            hpbl(i) = hght
         else
            hpbl(i) = hght - (rib-RIB_CRIT)/(rib-rib_below)*(hght-hght_below)
         endif
      enddo

      ! Successful completion
      status = RMN_OK
      return
   end function pbl_height1

end module pbl_height
