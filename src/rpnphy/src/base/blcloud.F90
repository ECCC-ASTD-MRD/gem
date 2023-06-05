!-------------------------------------- LICENCE BEGIN ------------------------
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

subroutine blcloud6(u,v,t,tve,qv,qc,fnn,frac,fngauss,fnnonloc,w_cld,&
     wb_ng,wthl_ng,wqw_ng,uw_ng,vw_ng,f_cs,dudz,dvdz,&
     hpar,frv,z0m,fb_surf,gzmom,ze,s,sw,ps,dudz2,ri,&
     dthv,tau,vcoef,n,nk,ncld)
   use tdpack_const, only: DELTA, GRAV, RGASD
   use phy_options

   implicit none
!!!#include <arch_specific.hf>
   !
   ! Arguments
   integer, intent(in) :: n                             !horizontal dimension
   integer, intent(in) :: nk                            !vertical dimension
   integer, intent(in) :: ncld                          !number of cloud profiles for length scales
   real, intent(in) :: tau                              !time step length (s)
   real, dimension(n,nk), intent(in) :: u               !u-component wind (m/s)
   real, dimension(n,nk), intent(in) :: v               !v-component wind (m/s)
   real, dimension(n,nk), intent(in) :: t               !dry air temperature (K)
   real, dimension(n,nk), intent(in) :: tve             !virt. temperature on e-lev (K)
   real, dimension(n,nk), intent(in) :: qv              !specific humidity (kg/kg)
   real, dimension(n,nk), intent(in) :: qc              !PBL cloud water content (kg/kg)
   real, dimension(n,nk), intent(in) :: fngauss         !Gaussian (local) cloud fraction
   real, dimension(n), intent(in) :: frv                !friction velocity (m/s)
   real, dimension(n), intent(in) :: z0m                !roughness length (m)
   real, dimension(n), intent(in) :: fb_surf            !surface buoyancy flux
   real, dimension(n,nk), intent(in) :: gzmom           !height of momentum levels (m)
   real, dimension(n,nk), intent(in) :: ze              !height of e-levs (m)
   real, dimension(n,nk), intent(in) :: s               !sigma for full levels
   real, dimension(n,nk), intent(in) :: sw              !sigma for working levels
   real, dimension(n), intent(in) :: ps                 !surface pressure (Pa)
   real, dimension(*), intent(in) :: vcoef              !coefficients for vertical interpolation
   real, dimension(n), intent(inout) :: hpar            !height of parcel ascent (m)
   real, dimension(n,nk), intent(inout) :: fnn          !flux enhancement * cloud fraction
   real, dimension(n,nk), intent(out) :: frac           !cloud fraction (nonlocal)
   real, dimension(n,nk,ncld), intent(out) :: w_cld     !cloud-layer velocity scales
   real, dimension(n,nk), intent(out) :: wb_ng          !non-gradient buoyancy flux (nonlocal)
   real, dimension(n,nk), intent(out) :: wthl_ng        !non-gradient theta_l flux (nonlocal)
   real, dimension(n,nk), intent(out) :: wqw_ng         !non-gradient total water flux (nonlocal)
   real, dimension(n,nk), intent(out) :: uw_ng          !non-gradient x-component stress (nonlocal)
   real, dimension(n,nk), intent(out) :: vw_ng          !non-gradient y-component stress (nonlocal)
   real, dimension(n,nk), intent(out) :: dudz           !x-compontent vertical wind shear (s^-1)
   real, dimension(n,nk), intent(out) :: dvdz           !y-compontent vertical wind shear (s^-1)
   real, dimension(n,nk), intent(out) :: dudz2          !square of vertical wind shear (s^-2)
   real, dimension(n,nk), intent(out) :: ri             !gradient Richardson number
   real, dimension(n,nk), intent(out) :: dthv           !buoyancy flux (m2/s2)
   real, dimension(n,nk), intent(out) :: fnnonloc       !nonlocal cloud fraction
   real, dimension(n,nk), intent(out) :: f_cs           !factor for l_s coeff c_s (nonlocal)

   !Author
   !          J. Mailhot (Nov 2000)

   !Revision
   ! 001      A.-M. Leduc (Oct 2001) Automatic arrays
   ! 002      J. Mailhot (Jun 2002) Change calling sequence and rename BLCLOUD1
   ! 003      J. Mailhot (Feb 2003) Change calling sequence and rename BLCLOUD2
   ! 004      L. Spacek (Dec 2007) - add "vertical staggering" option
   ! 005                             change the name to blcloud3
   ! 006      J. Mailhot/A. Lock (Aug 2012) Revisions for non-local scaling option
   !                      Change calling sequence and rename BLCLOUD4

   !Object
   !          Calculate the boundary layer buoyancy parameters (virtual potential
   !          temperature, buoyancy flux) and the vertical shear squared.

   !Notes
   !          Implicit (i.e. subgrid-scale) cloudiness scheme for unified
   !             description of stratiform and shallow, nonprecipitating
   !             cumulus convection appropriate for a low-order turbulence
   !             model based on Bechtold et al.:
   !            - Bechtold and Siebesma 1998, JAS 55, 888-895
   !            - Cuijpers and Bechtold 1995, JAS 52, 2486-2490
   !            - Bechtold et al. 1995, JAS 52, 455-463
   !            - Bechtold et al. 1992, JAS 49, 1723-1744
   !          The boundary layer cloud properties (cloud fraction, cloud water
   !            content) are computed in the companion S/R CLSGS.

   ! Local parameter definitions
   integer, parameter :: IMPLICIT_CLOUD=0
   logical, parameter :: COMPUTE_FROM_CONSERVED=.true.

   ! Local variable declarations
   integer :: j, k
   real, dimension(n,nk) :: thl,qw,alpha,beta,a,b,c,dz,dqwdz,dthldz,coefthl,coefqw,ficelocal, &
        unused,thv,tmom,qvmom,qcmom

   ! Compute ice fraction from temperature profile
   call ficemxp2(ficelocal,unused,unused,t,n,nk)

   ! Compute thermodynamic coefficients following Bechtold and Siebsma (JAS 1998)
   call thermco3(t,qv,qc,sw,ps,t,ficelocal,fnn,thl,qw,a,b,c,alpha,beta, &
        IMPLICIT_CLOUD,COMPUTE_FROM_CONSERVED,n,nk)
   
   ! Compute layer thicknesses
   do k=1,nk-1
      dz(:,k) = -RGASD*tve(:,k)*alog(s(:,k+1)/s(:,k))/GRAV
   end do
   dz(:,nk) = 0.

   ! Compute non-gradient fluxes and PBL cloud fraction in nonlocal formulation
   NONLOCAL_CLOUD: IF ( pbl_nonloc == 'LOCK06' ) THEN
      ! Convert to non-staggered state for nonlocal cloud estimates (FIXME)
      call vint_thermo2mom(tmom,t,vcoef,n,nk)
      call vint_thermo2mom(qvmom,qv,vcoef,n,nk)
      call vint_thermo2mom(qcmom,qc,vcoef,n,nk)
      ! Compute non-local cloud properties
      call nlcalc(fnn,frac,fngauss,fnnonloc,wb_ng,wthl_ng,wqw_ng,w_cld,u,v,uw_ng,vw_ng, &
           f_cs,tmom,qvmom,qcmom,frv,z0m,fb_surf,hpar,tau,gzmom,ze,ps,s,n,n,nk,size(w_cld,dim=3))
   else
      w_cld = 0.
      f_cs = 1.
   endif NONLOCAL_CLOUD

   ! Compute terms of buoyancy flux equation (Eq. 4 of Bechtold and Siebsma (JAS 1998))
   thv = thl + alpha*qw + beta*qc
   coefthl = 1.0 + DELTA*qw  - beta*b*fnn
   coefqw = alpha + beta*a*fnn

   ! Add to WB_NG (non-gradient flux of buoyancy) the contributions from
   ! WTHL_NG and WQW_NG (non-gradient fluxes of theta_l and q_w)
   NONLOCAL_FLUX: IF (pbl_nonloc == 'LOCK06') THEN
      do k=1,nk-1
         do j=1,n
            wb_ng(j,k) = wb_ng(j,k) + ( coefthl(j,k)*wthl_ng(j,k) &
                 + coefqw(j,k)*wqw_ng(j,k) ) &
                 * ( grav / thv(j,k) )
         end do
      end do
      wb_ng(:,nk) = 0.
   endif NONLOCAL_FLUX

   ! Compute vertical derivatives of conserved variables
   call vint_thermo2mom(thl, thl, vcoef, n, nk)
   call vint_thermo2mom(qw,  qw,  vcoef, n, nk)
   call dvrtdf(dthldz,thl,dz,n,n,n,nk)
   call dvrtdf(dqwdz,qw,dz,n,n,n,nk)

   ! Compute buoyancy flux
   do k=1,nk-1
      do j=1,n
         dthv(j,k) = (coefthl(j,k)*dthldz(j,k) + coefqw(j,k)*dqwdz(j,k)) * (GRAV/thv(j,k))
      enddo
   enddo
   dthv(:,nk) = 0.
 
   ! Compute vertical wind shear
   call dvrtdf(dudz,u,dz,n,n,n,nk)
   call dvrtdf(dvdz,v,dz,n,n,n,nk)
   dudz2 = dudz**2 + dvdz**2

   ! Compute gradient Richardson number
   ri = dthv / (dudz2+1e-6)  !FIXME: contains large background shear
   ri(:,nk) = 0.

   return
end subroutine blcloud6
