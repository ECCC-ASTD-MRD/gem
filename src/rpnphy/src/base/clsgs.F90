!-------------------------------------- LICENCE BEGIN -------------------------
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

subroutine clsgs8(thl,tve,qw,qc,frac,fnn,fngauss,fnnonloc,c1,zn,ze, &
     hpbl,hflux,s,ps,gztherm,mg,mrk2,acoef,bcoef,ccoef,vcoef,pblsigs,pblq1,n,nk)
   use tdpack_const
   use phy_options
   use ens_perturb, only: ens_nc2d, ens_spp_get
   implicit none
!!!#include <arch_specific.hf>

   !@Arguments
   integer, intent(in) :: n                             !horizontal dimension
   integer, intent(in) :: nk                            !vertical dimension
   real, dimension(n,nk), intent(in) :: thl             !liquid water potential temperature (K; theta_l)
   real, dimension(n,nk), intent(in) :: tve             !virtual temperature on e-levs (K)
   real, dimension(n,nk), intent(in) :: qw              !total water mixing ratio (kg/kg; q_tot)
   real, dimension(n,nk), intent(in) :: c1              !coefficient C1 in second-order moment closure
   real, dimension(n,nk), intent(in) :: zn              !mixing length (m)
   real, dimension(n,nk), intent(in) :: ze              !dissipation length (m)
   real, dimension(n), intent(in) :: hpbl               !boundary layer height (m)
   real, dimension(n), intent(in) :: hflux              !surface heat flux (W/m2)
   real, dimension(n,nk), intent(in) :: s               !sigma for full levels
   real, dimension(n),    intent(in) :: ps              !surface pressure (Pa)
   real, dimension(n,nk), intent(in) :: gztherm         !height of thermodynamic levels (m)
   real, dimension(n),    intent(in) :: mg              !land-sea mask
   real, dimension(n,ens_nc2d), intent(in) :: mrk2      !Markov chains for stochastic parameters
   real, dimension(n,nk), intent(in) :: acoef           !thermodynamic coefficient A
   real, dimension(n,nk), intent(in) :: bcoef           !thermodynamic coefficient B
   real, dimension(n,nk), intent(in) :: ccoef           !thermodynamic coefficient C
   real, dimension(*), intent(in) :: vcoef              !coefficients for vertical interpolation
   real, dimension(n,nk), intent(inout) :: fnnonloc     !nonlocal cloud fraction
   real, dimension(n,nk), intent(out) :: qc             !PBL cloud water content (kg/kg)
   real, dimension(n,nk), intent(out) :: frac           !cloud fraction
   real, dimension(n,nk), intent(out) :: fnn            !flux enhancement * cloud fraction
   real, dimension(n,nk), intent(out) :: fngauss        !Gaussian cloud fraction
   real, dimension(n,nk), intent(out) :: pblsigs        !Subgrid moisture variance
   real, dimension(n,nk), intent(out) :: pblq1          !Normalized saturation deficit

   !@Author
   !          J. Mailhot (Jun 2002)
   !@Revision
   ! 001      J. Mailhot (Feb 2003) Clipping at upper levels
   ! 002      S. Belair  (Apr 2003) Minimum values of 50 m for ZE and ZN
   !                                in calculation of sigmase.
   ! 003      A-M. Leduc (Jun 2003) Pass ps to blweight ---> blweight2
   ! 004      J. P. Toviessi ( Oct. 2003) - IBM conversion
   !               - calls to vslog routine (from massvp4 library)
   !               - unnecessary calculations removed
   !               - etc.
   ! 005      L. Spacek (Dec 2007) - add "vertical staggering" option
   !                                 change the name to clsgs3
   ! 006      A. Zadra (Oct 2015) -- add land-water mask (MG) to input and
   !               add a new user-defined reduction parameter
   !               (which may or may not depend on the land-water fraction)
   !               to control the flux enhancement factor (fnn)
   ! 007      A. Zadra / R. McT-C (Sep 2016) - Revisions from Adrian Lock
   !               and Jocelyn Mailhot (2012) (non-local scalings,
   !               clips for sigmas and Q1, remove min of 50m for ZE and ZN)
   !@Object Calculate the boundary layer sub-grid-scale cloud properties
   !@Notes
   !          Implicit (i.e. subgrid-scale) cloudiness scheme for unified
   !             description of stratiform and shallow, nonprecipitating
   !             cumulus convection appropriate for a low-order turbulence
   !             model based on Bechtold et al.:
   !            - Bechtold and Siebesma 1998, JAS 55, 888-895 (BS98)
   !            - Cuijpers and Bechtold 1995, JAS 52, 2486-2490
   !            - Bechtold et al. 1995, JAS 52, 455-463 (BCMT95)
   !            - Bechtold et al. 1992, JAS 49, 1723-1744
   !          Revisions for non-local scalings based on:
   !            - Lock and Mailhot 2006, BLM 121, 313-338 (LM06)
   !
   !    The parameters fnn_mask and fnn_reduc are the two parameters
   !    that should be read from the gem_settings
   !       fnn_mask = T/F means "use (do not use) the land-water fraction
   !            to modulate the fnn reduction
   !       fnn_reduc = is the reduction parameter that should vary within
   !            the range 0. to 1.; 1 means "keep the original estimate";
   !            any value smaller than 1 means " multiply the original fnn
   !            by a factor fnn_reduc"
   !
   !            The default values should be
   !              fnn_mask = .false.
   !              fnn_reduc = 1
   !
   !            The values tested and chosen for the RDPS-10km are
   !              fnn_mask = .true.
   !              fnn_reduc = 0.8

   ! Local parameters
   real, parameter :: EPS=1e-10,QCMIN=1e-6,QCMAX=1e-3,WHMIN=0.5,WHMAX=1.2

   ! Local variables
   integer :: j,k
   real :: fnn_weight,sigmas_cu,hfc,lnsig,c1coef
   real(kind=8) :: gravinv
   real, dimension(n) :: wf,fnnreduc
   real, dimension(n,nk) :: dz,dqwdz,dthldz,sigmas,q1,weight,thlm,qwm,wh

   ! Retrieve stochastic information for parameter perturbations
   fnnreduc(:) = ens_spp_get('fnnreduc', mrk2, fnn_reduc)
   
   ! Pre-compute constants and vertical derivatives
   gravinv = 1.0/dble(GRAV)
   call vint_thermo2mom(thlm, thl, vcoef, n, nk)
   call vint_thermo2mom(qwm,  qw,  vcoef, n, nk)
   do k=1,nk-1
      do j=1,n
         lnsig   = log(s(j,k+1)/s(j,k))
         dz(j,k) = -RGASD*tve(j,k)*lnsig*gravinv
      end do
   end do
   dz(:,nk) = 0.0
   call dvrtdf( dthldz, thlm, dz, n, n, n, nk)
   call dvrtdf( dqwdz , qwm , dz, n, n, n, nk)

   ! Compute subgrid standard deviation of supersaturation (s) on e-levels following Eq. 10 of BCMT95
   sigmas = 0.
   do k=1,nk-1
      do j=1,n
         c1coef = sqrt(c1(j,k)*max(zn(j,k),50.)*max(ze(j,k),50.))
         sigmas(j,k) = c1coef * abs(acoef(j,k)*dqwdz(j,k) - bcoef(j,k)*dthldz(j,k))
      end do
   end do

   ! Compute the normalized saturation defecit (Q1)
   q1(:,1:nk-1) = ccoef(:,1:nk-1) / ( sigmas(:,1:nk-1) + EPS )
   q1(:,1) = 0. ; q1(:,nk) = 0. ;
   q1=max(-6.,min(4.,q1))
   sigmas(:,1) = 0.; sigmas(:,nk) = 0.
   pblq1 = q1
   pblsigs = sigmas

   ! Compute cloud properties for local or nonlocal scalings
   do k=2,nk-1
      do j=1,n

         NONLOCAL: if( pbl_nonloc == 'LOCK06' ) then

            ! Compute cloud fractions (LM06)
            if( q1(j,k) > 6.0 ) then
               fngauss(j,k) = 1.0
               fnn(j,k)  = 1.0
            elseif ( q1(j,k) .gt. -6.0 ) then
               ! Represent Bechtold FN function as 0.5(1+tanh(0.8*Q1))
               ! (gives much the same shape as atan but without need for max/min)
               fngauss(j,k) = 0.5*( 1.0 + tanh(0.8*q1(j,k)) )
               ! Use shifted approximate Gaussian for flux enhancement (Eq. 5)
               fnn(j,k)  = 0.5*( 1.0 + tanh(0.8*(q1(j,k)+0.5)) )
            else
               fngauss(j,k) = 0.0
               fnn(j,k)  = 0.0
            endif
            ! Combine Gaussian and non-local cumulus cloud fractions
            frac(j,k) = max(fngauss(j,k),fnnonloc(j,k))

            ! Now increase sigmas to include cumulus contribution in the low cloudiness regime (C<0).
            ! This then gives a less negative Q1 and so more realistic QC
            if ( ccoef(j,k) < 0. .and. fnnonloc(j,k) > 0.0001 .and. fnnonloc(j,k) < 0.4999) then
               !         ! Take inverse of tanh function to find sigmas_cu given FNNONLOC
               sigmas_cu = 1.6 * ccoef(j,k) / alog( fnnonloc(j,k)/(1.-fnnonloc(j,k)) )
               sigmas(j,k) = max( sigmas(j,k), sigmas_cu )
               q1(j,k)     = ccoef(j,k)/sigmas(j,k)
            endif

         else

            ! Compute cloud fraction (BS98, Eq. B1)
            if( q1(j,k) > -1.2) then
               frac(j,k) = max(0.,min(1.,0.5 + 0.36*atan(1.55*q1(j,k))))
            elseif( q1(j,k) >= -6.0) then
               frac(j,k) = exp(q1(j,k)-1.0)
            else
               frac(j,k) = 0.
            endif
            fnnonloc(j,k) = frac(j,k)

         endif NONLOCAL

         ! Compute liquid water specific humidity (BS98, Eq. B2)
         !JM start (modification to LM06 Eq. 6) - no impact
!!$        if( q1(j,k) >= 2.0 ) then
!!$          qc(j,k) = 2.032 + 0.9*(q1(j,k)-2.)
!!$        elseif( q1(j,k) >= 0.0 ) then
         if (q1(j,k) >= 0.0) then
            !JM end
            qc(j,k) = exp(-1.) + 0.66*q1(j,k) + 0.086*q1(j,k)**2
         elseif( q1(j,k) >= -6.0 ) then
            qc(j,k) = exp(1.2*q1(j,k)-1.)
         else
            qc(j,k) = 0.
         endif
         !JM start - no impact
!!$        qc(j,k) = min ( qc(j,k)*sigmas(j,k) , qcmax )
         qc(j,k) = min(qc(j,k)*(sigmas(j,k)+EPS),qcmax)
         !JM end

         ! Compute cloud-induced flux enhancement factor * cloud fraction
         if (pbl_nonloc == 'NIL') then
            TRADE_WIND_CU: if (pbl_cucloud) then
               ! Compute flux enhancement factor (BS98, approx. final Eq. of Appendix B, Fig. 4a)
            fnn(j,k) = 1.0
            if (q1(j,k) < 1.0 .and. q1(j,k) >= -1.68) then
               fnn(j,k) = exp(-0.3*(q1(j,k)-1.0))
            elseif (q1(j,k) < -1.68 .and. q1(j,k) >= -2.5) then
               fnn(j,k) = exp(-2.9*(q1(j,k)+1.4))
            elseif (q1(j,k) < -2.5) then
               fnn(j,k) = 23.9 + exp(-1.6*(q1(j,k)+2.5))
            endif
               ! Adjust for cloud fraction (BS98, Fig. 4b)
            fnn(j,k) = fnn(j,k)*frac(j,k)
            if( q1(j,k).le.-2.39 .and. q1(j,k).ge.-4.0 ) then
               fnn(j,k) = 0.60
            elseif( q1(j,k).lt.-4.0 .and. q1(j,k).ge.-6.0 ) then
               fnn(j,k) = 0.30*( q1(j,k)+6.0 )
            elseif( q1(j,k).lt.-6. ) then
               fnn(j,k) = 0.0
            endif
            else
               ! Without a trade wind cumulus regime, FnN = N as in BS98 Appendix B
               fnn(j,k) = frac(j,k)
            endif TRADE_WIND_CU
         endif

         ! Permit PBL moist processes only under specified conditions
         if (.not.pbl_moistke_legacy_cloud) then
            ! Linear ramp-down of cloud effects from 5 W/m2 to -5 W/m2
            hfc   = 5.
            wf(j) = 0.5*( 1. + hflux(j)/hfc)
            wf(j) = max(0.,min(wf(j),1.))
            frac(j,k)     = wf(j)*frac(j,k)
            fnnonloc(j,k) = wf(j)*fnnonloc(j,k)
            fnn(j,k)      = wf(j)*fnn(j,k)
            qc(j,k)       = wf(j)*qc(j,k)
         endif

         ! Adjust FNN by user-controlled scaling factor potentially land-sea masked
         fnn_weight = 1.
         if (fnn_mask) then
            fnn_weight = fnnreduc(j) + (1.-fnnreduc(j))*max(0.,min(1.,mg(j)))
         else
            fnn_weight = fnnreduc(j)
         endif
         fnn(j,k) = fnn_weight*fnn(j,k)

      end do
   end do

   ! Fill array extrema
   frac(1:n,1) = 0. ; frac(1:n,nk) = 0.
   fnn (1:n,1) = 0. ; fnn (1:n,nk) = 0.
   qc  (1:n,1) = 0. ; qc  (1:n,nk) = 0.
   fnnonloc(1:n,1) = 0.; fnnonloc(1:n,nk) = 0.
   fngauss(1:n,1) = 0.; fngauss(1:n,nk) = 0.

   ! Create a vertical function to eliminate PBL cloud effects at upper levels
   call blweight2(weight,s,ps,n,nk)

   ! Do not allow PBL cloud effects to reach beyond 1.5 times the PBL depth
   if (.not.pbl_moistke_legacy_cloud) then
      do k=1,nk
         ! Linear ramp-down of cloud effects from 1.5*hpbl to 2*hpbl
         wh(:,k) = (2.*hpbl(:) - gztherm(:,k))/(0.5*hpbl(:))
         wh(:,k) = max(0., min(1., wh(:,k)))
         weight(:,k) = weight(:,k)*wh(:,k)
      enddo
   endif

   ! Scale all PBL cloud quantities
   frac = frac*weight
   fnnonloc = fnnonloc*weight
   fnn = fnn*weight
   qc = qc*weight
   fngauss = fngauss*weight

   return
end subroutine clsgs8
