
module ccc2_lwtragh
   private
   public :: ccc2_lwtragh1
   
contains

subroutine ccc2_lwtragh1(fu, fd, slwf, tauci, omci, &
     taual, taug, bf, urbf, cldfrac, em0, bs, cut, ni, lay, lev)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include "phymkptr.hf"

   integer, intent(in) :: ni, lay, lev
   real, intent(in) :: cut
   real, intent(in) :: slwf(ni), tauci(ni,lay), omci(ni,lay), &
        taual(ni,lay), taug(ni,lay), bf(ni,lev), urbf(ni,lay), &
        cldfrac(ni,lay), em0(ni), bs(ni)
   real, intent(inout) :: fu(ni,2,lev), fd(ni,2,lev)

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)
   !@Revisions
   ! 001    P. Vaillancourt, K. Winger (Sep 2006)- remplace 1-cld-eps par max(1-cld,eps)
   ! 002    M. Lazare (Nov 22/2006) -previous version lwtragh2 for gcm15f:
   !                               - bound of 1.e-10 used for 1-cld
   !                                 instead of subtracting an epsilon.
   !                               - work arays for this routine now
   !                                 are automatic arrays and their
   !                                 space is not passed in.
   ! 003    M.Lazare (Dec 05,2007 - previous version lwtragh3 for gcm15g:
   !                              - support added for emissivity<1.
   ! 004    M. Desgagne, P. Vaillancourt, M. Lazarre (Dec 2008) :
   !        Give appropriate dimensions to xu,xd,dtr,fy and fx;
   !        Initialize fx(I,1,1) and fx(I,2,1) correctly
   ! 005    J. Cole (Feb 11,2009) -  new version for gcm15h:
   !                              - correct bug in specification of abse0.
   !                              - initialize clear-sky fx in top layer.
   !                              - change any work arrays having a middle
   !                                 dimension of "4" to "2", as this is
   !                                 what is used.
   ! 006  P. vaillancourt, A-M. Leduc (Apr 2010) - Integrating the Li/Lazare 2009 version.(lwtragh4)
   !                              - modification of the fu and fy in loop 300.
   !                              - add variables EMBS and ABSEO for  efficiency.
   !@Object OPTICAL DEPTHS CALCULATION
   !        In the g space with interval close 1 (very large optical depth)
   !        or in the case with cloud absorption is very small or the weight
   !        of flux and cooling rate are very small. The cloud radiative
   !        process can be highly simplified. The absorption approximation
   !        method is used and cloud random and maximum overlap is
   !        considered, but cloud scattering and inhomogeneity are ignored.
   !        The exponential source planck function is used which is more
   !        accurate in the region above 200 mb in comparison with linear
   !        source function.
   !@Arguments
   !          - Output -
   ! fu       upward infrared flux
   ! fd       downward infrared flux
   !          - Input -
   ! slwf     input solar flux at model top level for each band
   ! tauci    cloud optical depth for the infrared
   ! omci     cloud single scattering albedo times optical depth
   ! taual    aerosol optical depth for the infrared
   ! taug     gaseous optical depth for the infrared
   ! bf       blackbody intensity integrated over each band at each
   !          level in units w / m^2 / sr.
   ! urbf     u times the difference of log(bf) for two neighbor levels
   !          used for exponential source function (li, 2002 jas p3302)
   ! cldfrac  cloud fraction
   ! em0      surface emission
   ! bs       the blackbody intensity at the surface.
   ! cut                 TO DEFINE    AML    small value (0.001) defined in raddriv
   ! ni          horizontal dimension
   ! lay          number of model levels
   ! lev          number of flux levels (lay+1)     TO CHECK    AML
   !----------------------------------------------------------------------
   ! xu       the emission part in the upward flux transmission
   !          (Li, 2002 JAS p3302)
   ! xd       the emission part in the downward flux transmission
   ! dtr      direct transmission
   ! fy       upward flux for pure clear portion (1) and pure cloud  portion (2)
   ! fx       the same as fy but for the downward flux

   !----------------------------------------------------------------------

   integer  ::i, k, km1, km2, kp1
   real  ::ubeta, epsd, epsu, taul2, cow
   real  ::ctaul2, crtaul2
   real  ::taul1(ni,lay), rtaul1(ni,lay)
   real  ::xu(ni,2,lay), xd(ni,2,lay), dtr(ni,2,lay), fy(ni,2,lev), &
        fx(ni,2,lev)  !# TODO: split these in 2 vars instead of adding a kank... better alignment?
   real(REAL64) :: dtr_vs(ni,lay)
   real  ::embs, abse0
   real :: w1, w2, v1, v2

   real, parameter :: RU = 1.6487213

   !----------------------------------------------------------------------
   !     initialization for first layer. calculate the downward flux in
   !     the second layer
   !     combine the optical properties for the infrared,
   !     1, aerosol + gas; 2, cloud + aerosol + gas.
   !     fd (fu) is down (upward) flux
   !     the overlap between solar and infrared in 4 - 10 um is
   !     considered, slwf is the incoming solar flux
   !     singularity for xd and xu has been considered as li jas 2002
   !----------------------------------------------------------------------

   do k = 2, lev
      km1 = k - 1
      do i = 1, ni
         taul1(i,km1)    =  taual(i,km1) + taug(i,km1)
         rtaul1(i,km1)   =  taul1(i,km1) * ru
         dtr_vs(i,km1)   =  exp(dble(-rtaul1(i,km1)))
      enddo
   enddo

   DO100: do i = 1, ni
      fd(i,1,1)         =  slwf(i)
      fd(i,2,1)         =  slwf(i)
      fx(i,1,1)         =  slwf(i)
      fx(i,2,1)         =  slwf(i)

      dtr(i,1,1)        =  dtr_vs(i,1)
      ubeta             =  urbf(i,1) / (taul1(i,1) + 1.e-20)
      epsd              =  ubeta + 1.0
      epsu              =  ubeta - 1.0

!!$      if (abs(epsd) .gt. 0.001)                                   then
!!$         xd(i,1,1)       = (bf(i,2) - bf(i,1) * dtr(i,1,1)) / epsd
!!$      else
!!$         xd(i,1,1)       =  rtaul1(i,1) * bf(i,1) * dtr(i,1,1)
!!$      endif
!!$      if (abs(epsu) .gt. 0.001)                                   then
!!$         xu(i,1,1)       = (bf(i,2) * dtr(i,1,1) - bf(i,1)) / epsu
!!$      else
!!$         xu(i,1,1)       =  rtaul1(i,1) * bf(i,2) * dtr(i,1,1)
!!$      endif
      
      w1 = MASK_LE(abs(epsd), 0.001)   !# abs(epsd) <= 0.001
      epsd = sign(max(0.001, abs(epsd)), epsd)
      v1 = (bf(i,2) - bf(i,1) * dtr(i,1,1)) / epsd
      v2 = rtaul1(i,1) * bf(i,1) * dtr(i,1,1)
      xd(i,1,1) = (1. - w1) * v1 + w1 * v2 

      w1 = MASK_LE(abs(epsu), 0.001)   !# abs(epsu) <= 0.001
      epsu = sign(max(0.001, abs(epsu)), epsu)
      v1 = (bf(i,2) * dtr(i,1,1) - bf(i,1)) / epsu
      v2 = rtaul1(i,1) * bf(i,2) * dtr(i,1,1)
      xu(i,1,1) = (1. - w1) * v1 + w1 * v2

      fd(i,1,2)         =  fd(i,1,1) * dtr(i,1,1) + xd(i,1,1)

      !#TODO: transform this one too into a weight/mask fn
      if (cldfrac(i,1) .lt. cut)                                  then
         fx(i,1,2)       =  fd(i,1,2)
         fx(i,2,2)       =  fd(i,1,2)
         fd(i,2,2)       =  fd(i,1,2)
         dtr(i,2,1)      = 0.
         xd(i,2,1)       = 0.
         xu(i,2,1)       = 0.
      else
         taul2           =  tauci(i,1) + taul1(i,1)
         cow             =  1.0 - omci(i,1) / taul2
         ctaul2          =  cow * taul2
         crtaul2         =  ctaul2 * ru
         dtr(i,2,1)      =  exp (- crtaul2)
         ubeta           =  urbf(i,1) / (ctaul2)
         epsd            =  ubeta + 1.0
         epsu            =  ubeta - 1.0

!!$         if (abs(epsd) .gt. 0.001)                                 then
!!$            xd(i,2,1)     = (bf(i,2) - bf(i,1) * dtr(i,2,1)) / epsd
!!$         else
!!$            xd(i,2,1)     =  crtaul2 * bf(i,1) * dtr(i,2,1)
!!$         endif
!!$         if (abs(epsu) .gt. 0.001)                                 then
!!$            xu(i,2,1)     = (bf(i,2) * dtr(i,2,1) - bf(i,1)) / epsu
!!$         else
!!$            xu(i,2,1)     =  crtaul2 * bf(i,2) * dtr(i,2,1)
!!$         endif
         
         w1 = MASK_LE(abs(epsd), 0.001)   !# abs(epsd) <= 0.001
         epsd = sign(max(0.001, abs(epsd)), epsd)
         v1 = (bf(i,2) - bf(i,1) * dtr(i,2,1)) / epsd
         v2 = crtaul2 * bf(i,1) * dtr(i,2,1)
         xd(i,2,1) = (1. - w1) * v1 + w1 * v2 

         w1 = MASK_LE(abs(epsu), 0.001)   !# abs(epsu) <= 0.001
         epsu = sign(max(0.001, abs(epsu)), epsu)
         v1 = (bf(i,2) * dtr(i,2,1) - bf(i,1)) / epsu
         v2 = crtaul2 * bf(i,2) * dtr(i,2,1)
         xu(i,2,1) = (1. - w1) * v1 + w1 * v2

         fx(i,1,2)       =  fx(i,1,1) * dtr(i,1,1) + xd(i,1,1)
         fx(i,2,2)       =  fx(i,2,1) * dtr(i,2,1) + xd(i,2,1)
         fd(i,2,2)       =  fx(i,1,2) + &
              cldfrac(i,1) * (fx(i,2,2) - fx(i,1,2))
      endif
   enddo DO100

   DO250: do k = 3, lev
      km1 = k - 1
      km2 = km1 - 1
      do i = 1, ni
         dtr(i,1,km1)    =  dtr_vs(i,km1)
         ubeta           =  urbf(i,km1) / (taul1(i,km1) + 1.e-20)
         epsd            =  ubeta + 1.0
         epsu            =  ubeta - 1.0

!!$         if (abs(epsd) .gt. 0.001)                                 then
!!$            xd(i,1,km1)   = (bf(i,k) - bf(i,km1) * dtr(i,1,km1)) / epsd
!!$         else
!!$            xd(i,1,km1)   =  rtaul1(i,km1) * bf(i,km1) * dtr(i,1,km1)
!!$         endif
!!$         if (abs(epsu) .gt. 0.001)                                 then
!!$            xu(i,1,km1)   = (bf(i,k) * dtr(i,1,km1) - bf(i,km1)) / epsu
!!$         else
!!$            xu(i,1,km1)   =  rtaul1(i,km1) * bf(i,k) * dtr(i,1,km1)
!!$         endif

         w1 = MASK_LE(abs(epsd), 0.001)   !# abs(epsd) <= 0.001
         epsd = sign(max(0.001, abs(epsd)), epsd)
         v1 = (bf(i,k) - bf(i,km1) * dtr(i,1,km1)) / epsd
         v2 = rtaul1(i,km1) * bf(i,km1) * dtr(i,1,km1)
         xd(i,1,km1) = (1. - w1) * v1 + w1 * v2 

         w1 = MASK_LE(abs(epsu), 0.001)   !# abs(epsu) <= 0.001
         epsu = sign(max(0.001, abs(epsu)), epsu)
         v1 = (bf(i,k) * dtr(i,1,km1) - bf(i,km1)) / epsu
         v2 = rtaul1(i,km1) * bf(i,k) * dtr(i,1,km1)
         xu(i,1,km1) = (1. - w1) * v1 + w1 * v2
          
         fd(i,1,k)       =  fd(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)

!!$         w2 = MASK_GE(cldfrac(i,km1), cut)   !# cldfrac(i,km1) >= cut
         
         !#TODO: transform this one too into a weight/mask fn
         if (cldfrac(i,km1) .lt. cut)                              then
            fd(i,2,k)     =  fd(i,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
            fx(i,1,k)     =  fd(i,2,k)
            fx(i,2,k)     =  fd(i,2,k)
            dtr(i,2,km1)  =  0.
            xd(i,2,km1)   =  0.
            xu(i,2,km1)   =  0.
         else
            taul2         =  tauci(i,km1) + taul1(i,km1)
            cow           =  1.0 - omci(i,km1) / taul2
            ctaul2        =  cow * taul2
            crtaul2       =  ctaul2 * ru
            dtr(i,2,km1)  =  exp (- crtaul2)
            ubeta         =  urbf(i,km1) / (ctaul2)
            epsd          =  ubeta + 1.0
            epsu          =  ubeta - 1.0

!!$            if (abs(epsd) .gt. 0.001)                               then
!!$               xd(i,2,km1) = (bf(i,k) - bf(i,km1) * dtr(i,2,km1)) / epsd
!!$            else
!!$               xd(i,2,km1) =  crtaul2 * bf(i,km1) * dtr(i,2,km1)
!!$            endif
!!$            if (abs(epsu) .gt. 0.001)                               then
!!$               xu(i,2,km1) = (bf(i,k) * dtr(i,2,km1) - bf(i,km1)) / epsu
!!$            else
!!$               xu(i,2,km1) =  crtaul2 * bf(i,k) * dtr(i,2,km1)
!!$            endif

            w1 = MASK_LE(abs(epsd), 0.001)   !# abs(epsd) <= 0.001
            epsd = sign(max(0.001, abs(epsd)), epsd)
            v1 = (bf(i,k) - bf(i,km1) * dtr(i,2,km1)) / epsd
            v2 = crtaul2 * bf(i,km1) * dtr(i,2,km1)
            xd(i,2,km1) = (1. - w1) * v1 + w1 * v2 

            w1 = MASK_LE(abs(epsu), 0.001)   !# abs(epsu) <= 0.001
            epsu = sign(max(0.001, abs(epsu)), epsu)
            v1 = (bf(i,k) * dtr(i,2,km1) - bf(i,km1)) / epsu
            v2 = crtaul2 * bf(i,k) * dtr(i,2,km1)
            xu(i,2,km1) = (1. - w1) * v1 + w1 * v2
            
!!$            if (cldfrac(i,km1) .le. cldfrac(i,km2))                 then
!!$               fx(i,1,k)   = ( fx(i,2,km1) + (1.0 - cldfrac(i,km2)) / &
!!$                    max(1.0 - cldfrac(i,km1),1.e-10) * &
!!$                    (fx(i,1,km1) - fx(i,2,km1)) ) * &
!!$                    dtr(i,1,km1) + xd(i,1,km1)
!!$               fx(i,2,k)   =  fx(i,2,km1) * dtr(i,2,km1) + xd(i,2,km1)
!!$            else if (cldfrac(i,km1) .gt. cldfrac(i,km2))            then
!!$               fx(i,1,k)   =  fx(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
!!$               fx(i,2,k)   = (fx(i,1,km1)+cldfrac(i,km2)/cldfrac(i,km1) * &
!!$                    (fx(i,2,km1) - fx(i,1,km1))) * &
!!$                    dtr(i,2,km1) + xd(i,2,km1)
!!$            endif
            
            w1 = MASK_LE(cldfrac(i,km1), cldfrac(i,km2))   !# cldfrac(i,km1) <= cldfrac(i,km2)
            
            v1 = (fx(i,2,km1) + (1.0 - cldfrac(i,km2)) / &
                 max(1.e-10, 1.0 - cldfrac(i,km1)) * &
                 (fx(i,1,km1) - fx(i,2,km1)) ) * &
                 dtr(i,1,km1) + xd(i,1,km1)
            v2 = fx(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
            fx(i,1,k) = w1 * v1 + (1. - w1) * v2

            v1 = fx(i,2,km1) * dtr(i,2,km1) + xd(i,2,km1)
            v2 = (fx(i,1,km1) + cldfrac(i,km2) / max(1.e-10, cldfrac(i,km1)) * &
                 (fx(i,2,km1) - fx(i,1,km1))) * &
                 dtr(i,2,km1) + xd(i,2,km1)
            fx(i,2,k) = w1 * v1 + (1. - w1) * v2
             
            fd(i,2,k) = fx(i,1,k) + cldfrac(i,km1) * (fx(i,2,k) - fx(i,1,k))
         endif
      enddo
   enddo DO250

   do i = 1, ni
      embs             =  em0(i) * bs (i)
      abse0             = 1.0 - em0(i)
      fu(i,1,lev)      =  embs + abse0 * fd(i,1,lev)
      fy(i,1,lev)      =  embs + abse0 * fx(i,1,lev)
      fy(i,2,lev)      =  embs + abse0 * fx(i,2,lev)
      fu(i,2,lev)      =  fy(i,1,lev) + &
           cldfrac(i,lay) * (fy(i,2,lev) - fy(i,1,lev))

      fu(i,1,lay)      =  fu(i,1,lev) * dtr(i,1,lay) + xu(i,1,lay)
      fy(i,1,lay)      =  fy(i,1,lev) * dtr(i,1,lay) + xu(i,1,lay)

!!$      if (cldfrac(i,lay) .lt. cut)                                then
!!$         fy(i,2,lay)    =  fy(i,2,lev) * dtr(i,1,lay) + xu(i,1,lay)
!!$         fu(i,2,lay)    =  fy(i,1,lay)
!!$      else
!!$         fy(i,2,lay)    =  fy(i,2,lev) * dtr(i,2,lay) + xu(i,2,lay)
!!$         fu(i,2,lay)    =  fy(i,1,lay) + &
!!$              cldfrac(i,lay) * (fy(i,2,lay) - fy(i,1,lay))
!!$      endif

      w1 = MASK_GE(cldfrac(i,lay), cut)   !# cldfrac(i,lay) >= cut

!!$      print *,'ccc2a:',i,fy(i,2,lev),dtr(i,1,lay),xu(i,1,lay)
!!$      call flush(6)
      v1 = fy(i,2,lev) * dtr(i,1,lay) + xu(i,1,lay)
!!$      print *,'ccc2b:',i,fy(i,2,lev),dtr(i,2,lay),xu(i,2,lay)
!!$      call flush(6)
      v2 = fy(i,2,lev) * dtr(i,2,lay) + xu(i,2,lay)
!!$      print *,'ccc2c:',i,w1,v1,v2
!!$      call flush(6)
      fy(i,2,lay) = (1. - w1) * v1 + w1 * v2  !# crash here
      
      v1 = fy(i,1,lay)
      v2 = fy(i,1,lay) + cldfrac(i,lay) * (fy(i,2,lay) - fy(i,1,lay))
      fu(i,2,lay) = (1. - w1) * v1 + w1 * v2
      
   enddo

   do k = lev - 2, 1, - 1
      kp1 = k + 1
      do i = 1, ni
         fu(i,1,k)      =  fu(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)

!!$         w2 = MASK_GE(cldfrac(i,k), cut)   !# cldfrac(i,k) >= cut

         if (cldfrac(i,k) .lt. cut)                                then
            fu(i,2,k)    =  fu(i,2,kp1) * dtr(i,1,k) + xu(i,1,k)
            fy(i,1,k)    =  fu(i,2,k)
            fy(i,2,k)    =  fu(i,2,k)
         else

!!$            if (cldfrac(i,k) .lt. cldfrac(i,kp1))                   then
!!$               fy(i,1,k)  = ( fy(i,2,kp1) + (1.0 - cldfrac(i,kp1)) / &
!!$                    (1.0 - cldfrac(i,k)) * (fy(i,1,kp1) - &
!!$                    fy(i,2,kp1)) ) * dtr(i,1,k) + xu(i,1,k)
!!$               fy(i,2,k)  =  fy(i,2,kp1) * dtr(i,2,k) + xu(i,2,k)
!!$            else
!!$               fy(i,1,k)  =  fy(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)
!!$               fy(i,2,k)  = ( fy(i,1,kp1) + cldfrac(i,kp1)/cldfrac(i,k) * &
!!$                    (fy(i,2,kp1) - fy(i,1,kp1)) ) * dtr(i,2,k) + &
!!$                    xu(i,2,k)
!!$            endif

            w1 = MASK_GE(cldfrac(i,k), cldfrac(i,kp1))   !# cldfrac(i,k) >= cldfrac(i,kp1)
            
            v1 = ( fy(i,2,kp1) + (1.0 - cldfrac(i,kp1)) / &
                   max(1.e-10, (1.0 - cldfrac(i,k))) * (fy(i,1,kp1) - &
                    fy(i,2,kp1)) ) * dtr(i,1,k) + xu(i,1,k)
            v2 = fy(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)
            fy(i,1,k) = (1. - w1) * v1 + w1 * v2
            
            v1 = fy(i,2,kp1) * dtr(i,2,k) + xu(i,2,k)
            v2 = ( fy(i,1,kp1) + cldfrac(i,kp1) / max(1.e-10, cldfrac(i,k)) * &
                 (fy(i,2,kp1) - fy(i,1,kp1)) ) * dtr(i,2,k) + xu(i,2,k)
            fy(i,2,k) = (1. - w1) * v1 + w1 * v2

            fu(i,2,k)    =  fy(i,1,k) + &
                 cldfrac(i,k) * (fy(i,2,k) - fy(i,1,k))
         endif
         
!!$         fu(i,2,k) = (1. - w2) * va1 + w2 * va2
!!$         fy(i,1,k) = (1. - w2) * vb1 + w2 * vb2
!!$         fy(i,2,k) = (1. - w2) * vc1 + w2 * vc2
         
      enddo
   enddo

   return
end subroutine ccc2_lwtragh1

end module ccc2_lwtragh
