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
!-------------------------------------- LICENCE END --------------------------

subroutine baktotq8(dt,dqv,dqc,thl,qw,dthl,dqw,qc,s,sw,ps,gztherm,tif,fice,&
     tve,hpbl,hflux,qcbl,fnn,fn,fngauss,fnnonloc,c1,zn,ze,mg,mrk2, &
     vcoef,pblsigs,pblq1,tau,n,nk)
   use tdpack_const, only: CPD, CAPPA, CHLC, CHLF
   use ens_perturb, only: ens_nc2d
   implicit none
!!!#include <arch_specific.hf>

   !@Arguments
   integer, intent(in) :: n                             !horizontal dimension
   integer, intent(in) :: nk                            !vertical dimension
   real, intent(in) :: tau                              !time step length (s)
   real, dimension(n,nk), intent(in) :: thl             !liquid water potential temperature (K; theta_l)
   real, dimension(n,nk), intent(in) :: qw              !total water mixing ratio (kg/kg; q_tot)
   real, dimension(n,nk), intent(in) :: qc              !PBL cloud water content (kg/kg)
   real, dimension(n,nk), intent(in) :: s               !sigma for full levels
   real, dimension(n,nk), intent(in) :: sw              !sigma for working levels
   real, dimension(n),    intent(in) :: ps              !surface pressure (Pa)
   real, dimension(n,nk), intent(in) :: gztherm         !height of thermodynamic levels (m)
   real, dimension(n,nk), intent(in) :: tif             !temperature used for ice fraction (K)
   real, dimension(n,nk), intent(in) :: fice            !ice fraction
   real, dimension(n,nk), intent(in) :: dthl            !tendency of theta_l (K/s)
   real, dimension(n,nk), intent(in) :: dqw             !tendency of q_tot (kg/kg/s)
   real, dimension(n,nk), intent(in) :: tve             !virtual temperature on e-levs (K)
   real, dimension(n), intent(in) :: hpbl               !boundary layer height (m)
   real, dimension(n), intent(in) :: hflux              !surface heat flux (W/m2)
   real, dimension(n,nk), intent(in) :: c1              !constant C1 in second-order moment closure
   real, dimension(n,nk), intent(in) :: zn              !mixing length (m)
   real, dimension(n,nk), intent(in) :: ze              !dissipation length (m)
   real, dimension(n),    intent(in) :: mg              !land-sea mask
   real, dimension(n,ens_nc2d), intent(in) :: mrk2      !Markov chains for stochastic parameters
   real, dimension(*), intent(in) :: vcoef              !coefficients for vertical interpolation   
   real, dimension(n,nk), intent(inout) :: fnn          !flux enhancement * cloud fraction
   real, dimension(n,nk), intent(inout) :: fnnonloc     !nonlocal cloud fraction
   real, dimension(n,nk), intent(out) :: fn             !cloud fraction
   real, dimension(n,nk), intent(out) :: qcbl           !water content of PBL clouds (kg/kg)
   real, dimension(n,nk), intent(out) :: dt             !tendency of dry air temperature (K/s)
   real, dimension(n,nk), intent(out) :: dqv            !tendency of specific humidity (kg/kg/s)
   real, dimension(n,nk), intent(out) :: dqc            !tendency of cloud water content (kg/kg/s)
   real, dimension(n,nk), intent(out) :: fngauss        !Gaussian cloud fraction
   real, dimension(n,nk), intent(out) :: pblsigs        !Subgrid moisture variance
   real, dimension(n,nk), intent(out) :: pblq1          !Normalized saturation deficit

   !@Author  J. Mailhot (Nov 2000)
   !@Revision
   ! 001      A.-M. Leduc (Oct 2001) Automatic arrays
   ! 002      B. Bilodeau and J. Mailhot (Dec 2001) Add a test to
   !                      check the presence of advected explicit cloud water.
   ! 003      J. Mailhot (Nov 2000) Cleanup of routine
   ! 004      J. Mailhot (Feb 2003) - MOISTKE option based on implicit clouds only
   ! 005      A-M. Leduc (Jun 2003) - pass ps to clsgs---> clsgs2
   ! 006      J. P. Toviessi ( Oct. 2003) - IBM conversion
   !               - calls to exponen4 (to calculate power function '**')
   !               - etc.
   ! 007      B. Bilodeau (Dec 2003)   More optimizations for IBM
   !                                   - Call to vspown1
   !                                   - Replace divisions by multiplications
   ! 008      L. Spacek (Dec 2007) - add "vertical staggering" option
   !                                 change the name to baktotq3
   ! 009      A. Zadra (Oct 2015) -- add land-water mask (MG) to input, which
   !                                 is then passed on to CLSGS4
   ! 010      A. Zadra, R. McT-C (Sep 2016) - implement non-local scaling
   !                      deveoped by J. Mailhot/A. Lock (Aug 2012)
   !@Object
   !          Transform conservative variables and their tendencies
   !          back to non-conservative variables and tendencies.
   !          Calculate the boundary layer cloud properties (cloud fraction, cloud
   !          water content, flux enhancement factor).
 
   ! Local parameter definitions
   integer, parameter :: IMPLICIT_CLOUD=0
   logical, parameter :: COMPUTE_FROM_STATE=.false.

   ! Local variables
   real :: cpdinv,tauinv
   real, dimension(n,nk) :: exner,thl_star,qw_star,acoef,bcoef,ccoef,alpha,beta,qcp,unused,qv

   ! Precompute inverses
   CPDINV = 1./CPD
   TAUINV = 1./TAU

   ! Update conserved variables following diffusion
   exner = sw**CAPPA
   thl_star = thl + tau*dthl
   qw_star = qw + tau*dqw

   ! Compute thermodynamic coefficients from conserved variables
   call thermco3(unused,unused,unused,sw,ps,tif,fice,fnn,thl_star,qw_star,acoef,bcoef,ccoef,alpha,beta,&
        IMPLICIT_CLOUD,COMPUTE_FROM_STATE,n,nk)

   ! Retrive updated cloud water content (qcp) from conserved variables
   call clsgs8(thl_star,tve,qw_star,qcp,fn,fnn,fngauss,fnnonloc,c1,zn,ze,hpbl,hflux,s,ps,gztherm,&
        mg,mrk2,acoef,bcoef,ccoef,vcoef,pblsigs,pblq1,n,nk)

   ! Convert back to state variables and tendencies
   qv = qw - max(0.,qc)
   dqc = (max(0.,qcp) - max(0.,qc)) * tauinv
   dqc = max(dqc,-max(0.0,qc)*tauinv) !prevent negative values of qc
   dqc = min(dqc,dqw + max(0.0,qv)*tauinv) !prevent over-depletion of water vapour by PBL cloud
   qcbl = max(0.,qc) + dqc*tau
   dt = exner*dthl + ((CHLC+fice*CHLF)*CPDINV)*dqc
   dqv = max((dqw-dqc),-max(0.0,qv)*tauinv)

   return
end subroutine baktotq8

