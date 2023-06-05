!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------

subroutine thermco3(t,qv,qc,sw,ps,tif,fice,fnn,thl,qw,acoef,bcoef,ccoef,alpha,beta, &
     type,inmode,n,nk)
   use tdpack   
   implicit none
!!!#include <arch_specific.hf>
   
   ! Arguments
   integer, intent(in) :: n                             !horizontal dimension
   integer, intent(in) :: nk                            !vertical dimension
   integer, intent(in) :: type                          !cloud type switch (0=implicit; 1=explicit; 2=both)
   logical, intent(in) :: inmode                        !input mode (T=compute from t,qv,qc; F=compute from thl,qw)
   real, dimension(n,nk), intent(in) :: t               !dry air temperature (K)
   real, dimension(n,nk), intent(in) :: qv              !specific humidity (kg/kg)
   real, dimension(n,nk), intent(in) :: qc              !PBL cloud water content (kg/kg)
   real, dimension(n,nk), intent(in) :: sw              !sigma for working levels
   real, dimension(n), intent(in) :: ps                 !surface pressure (Pa)
   real, dimension(n,nk), intent(in) :: tif             !temperature used for ice fraction (K)
   real, dimension(n,nk), intent(in) :: fice            !ice fraction
   real, dimension(n,nk), intent(in) :: fnn             !flux enhancement * cloud fraction
   real, dimension(n,nk), intent(inout) :: thl          !liquid water potential temperature (K; theta_l)
   real, dimension(n,nk), intent(inout) :: qw           !total water content (kg/kg; q_tot)
   real, dimension(n,nk), intent(out) :: acoef          !thermodynamic coefficient A
   real, dimension(n,nk), intent(out) :: bcoef          !thermodynamic coefficient B
   real, dimension(n,nk), intent(out) :: ccoef          !thermodynamic coefficient C
   real, dimension(n,nk), intent(out) :: alpha          !thermodynamic coefficient alpha
   real, dimension(n,nk), intent(out) :: beta           !thermodynamic coefficient beta

   !Author
   !          J. Mailhot (Nov 1999)
   
   !Revision
   ! 001      J. Mailhot  (Jan 2000) - Changes to add mixed-phase mode
   ! 002      A.-M. Leduc (Oct 2001) Automatic arrays
   ! 003      J. Mailhot  (Jun 2002) - Add cloud type and input mode
   !                       Change calling sequence and rename THERMCO2
   ! 004      A. Plante   (May 2003) - IBM conversion
   !                         - calls to exponen4 (to calculate power function '**')
   !                         - divisions replaced by reciprocals
   !                         - calls to optimized routine mfdlesmx
   ! 005      B. Bilodeau (Aug 2003) - exponen4 replaced by vspown1
   
   !Object
   !          Calculate the thermodynamic coefficients used in the presence of clouds
   !          and the conservative variables.
   
   !Notes
   !          See definitions in:
   !          - Bechtold and Siebesma 1998, JAS 55, 888-895

   ! Local variable declarations
   integer :: j,k
   real, dimension(n,nk) :: pres,exner,qsat,dqsat,th,tl,ffice,tfice,dfice,work

   ! Precompute pressure and Exner function
   do k=1,nk
      pres(:,k) = sw(:,k)*ps
      exner(:,k) = sw(:,k)**CAPPA
   enddo
   ffice = fice
   
   ! Compute conserved variables from state inputs
   COMPUTE_CONSERVED: if (inmode) then
      if ( type .eq. 0 ) call ficemxp2(ffice,tfice,dfice,tif,n,nk)
      th = dble(t) / dble(exner)
      thl = th * (1.0-((CHLC+ffice*CHLF)/CPD) * (qc/t))
      qw = qv + qc
   endif COMPUTE_CONSERVED
   
   ! Compute liquid water temperature
   tl = exner*thl

   ! Compute saturation specific humidity for selected cloud type
   CLOUD_TYPE: select case (type)
      
   case(0) !Implicit clouds
      call ficemxp2(work,tfice,dfice,tif,n,nk)
      do k=1,nk
         do j=1,n
            qsat(j,k) = fqsmx(tl(j,k),pres(j,k),tfice(j,k))
         enddo
      enddo
      call mfdlesmx(work,tl,tfice,dfice,n,nk)
      do k=1,nk
         do j=1,n
            dqsat(j,k) = fdqsmx(qsat(j,k),work(j,k) )
         enddo
      enddo
      
   case (1) !Explicit clouds
      do k=1,nk
         do j=1,n
            qsat(j,k) = foqsa(tl(j,k),pres(j,k))
            dqsat(j,k) = fodqa(qsat(j,k),tl(j,k))
         enddo
      enddo

   case (2) !Combined implicit and explicit clouds
      call ficemxp2(work,tfice,dfice,tif,n,nk)
      call mfdlesmx(work,tl,tfice,dfice,n,nk)
      do k=1,nk
         do j=1,n
            if (fnn(j,k) < 1.0) then
               qsat(j,k) = fqsmx(tl(j,k),pres(j,k),tfice(j,k))
               dqsat(j,k) = fdqsmx(qsat(j,k),work(j,k) )
            else
               qsat(j,k) = foqsa(tl(j,k),pres(j,k))
               dqsat(j,k) = fodqa(qsat(j,k),tl(j,k))
            endif
         enddo
      enddo

   end select CLOUD_TYPE

   ! Compute thermodynamic coefficients following Bechtold and Siebesma (JAS 1998) Appendix A
   acoef = 1.0/( 1.0 + ((CHLC+ffice*CHLF)/CPD)*dqsat)
   bcoef = acoef*exner*dqsat
   ccoef = acoef*(qw-qsat)
   if (inmode) then
      alpha = DELTA*th
      beta = ((CHLC+ffice*CHLF)/CPD)/exner - (1.0+DELTA)*th
   else
      alpha = 0.
      beta = 0.
   endif

   return
end subroutine thermco3
