module condload_mod
   implicit none
   public

contains

!** S/P CONDLOAD
      SUBROUTINE CONDLOAD(QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ, &
                        QNEWIC,QLQOUT,QICOUT,GRAV)
!
!
      implicit none
!!!#include <arch_specific.hf>
!
!
      REAL BOTERM,CONV,DQ,DZ
      REAL ENTERM,GRAV,G1,OLDQ,PPTDRG
      REAL QEST,QICE,QICOUT,QLIQ,QLQOUT,QNEW
      REAL QNEWIC,QNEWLQ,QTOT,RATE
      REAL RATIO3,RATIO4,WAVG,WTW
!
!AUTHOR
!          Kain and Fritsch   (1988)
!
!REVISION
! 001      Stephane Belair (1994)
! 002      A-M Leduc       (2001)  Adaptation to physics 3.7
!
!OBJECT
!          to perform the precipitation fallout
!
!ARGUMENTS
!
!          - INPUT/OUTPUT -
! QLIQ     liquid water in the cloud layer as input
!          liquid water after the precipitation fallout as output
! QICE     cloud ice in the cloud layer as input
!          cloud ice after the precipitation fallout as output
! WTW      W*W of the updraft as input
!          new vertical velocity with the condensation load
!          effect as output
!
!          - INPUT -
! DZ       row number
! BOTERM   buoyancy term for the updraft
! ENTERM   entrainment term for the updraft
! RATE     rate of precipitation fallout production in the
!          updraft (= 0.01)
! QNEWLQ   fresh liquid water condensate in the updraft
! QNEWIC   fresh cloud ice in the updraft
!
!          - Output -
! QLQOUT   precipitation fallout (liquid water)
! QICOUT   precipitation fallout (cloud ice)
!
!
!Notes
!  9/18/88...This precipitation fallout scheme is based on the scheme
!  used by Ogura and Cho (1973). Liquid water fallout from a parcel
!  is calculated using the equation dq=-rate*q*dt, but to simulate a
!  quasi-continuous process, and to eliminate a dependency on vertical
!  resolution, this is expressed as q=q*exp(-rate*dz).
!*
!
      QTOT=QLIQ+QICE
      QNEW=QNEWLQ+QNEWIC
!
!
!
!                              Estimate the vertical velocity so that
!                              an average vertical velocity can be
!                              calculated to estimate the time required
!                              ascent between model levels.
!
      QEST=0.5*(QTOT+QNEW)
!
!                              Solve for the vertical motion.
!

      G1=WTW+BOTERM-ENTERM-2.*GRAV*DZ*QEST/1.5
      IF(G1.LT.0.0)G1=0.
      WAVG=(SQRT(WTW)+SQRT(G1))/2.
!
!                              Conversion
!
      CONV=RATE*DZ/WAVG
!
!
!                              RATIO3 is the fraction of liquid water
!                              in fresh condensate.
!                              RATIO4 is the fraction of liquid water in
!                              the total amount of condensate involved
!                              in the precipitation process.
!                              Note that only 60% of the fresh condensate
!                              is allowed to participate in the conversion
!                              process.
!
      RATIO3=QNEWLQ/(QNEW+1.E-10)
      QTOT=QTOT+0.6*QNEW
      OLDQ=QTOT
      RATIO4=(0.6*QNEWLQ+QLIQ)/(QTOT+1.E-10)
!
!                             Apply the rate equation.
!

      QTOT=QTOT*EXP(-CONV)
!
!
!                             Determine the amount of precipitation that
!                             falls out of the updraft parcel at this
!                             level.
!
      DQ=OLDQ-QTOT
      QLQOUT=RATIO4*DQ
      QICOUT=(1.-RATIO4)*DQ
!
!
!                             Estimate the mean load of condensate on
!                             the updraft in the layer.
!                             Calculate the vertical velocity.
!
!                             PPTDRG = precipitation drag
!
      PPTDRG=0.5*(OLDQ+QTOT-0.2*QNEW)
!
!                             Apply again equation for vertical motion.
!
      WTW=WTW+BOTERM-ENTERM-2.*GRAV*DZ*PPTDRG/1.5
!
!
!                             Determine the new liquid water and ice
!                             concentrations including lost due to
!                             precipitation and gains from condensation.
!                             Don't forget to include the fraction 0.4qnew
!                             which was not involved in the conversion
!                             process.
!
      QLIQ=RATIO4*QTOT+RATIO3*0.4*QNEW
      QICE=(1.-RATIO4)*QTOT+(1.-RATIO3)*0.4*QNEW
      QNEWLQ=0.
      QNEWIC=0.
!
!
      RETURN
   END SUBROUTINE CONDLOAD


!** S/P CONDLOAD_INFCONV
      SUBROUTINE CONDLOAD_INFCONV(QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ, &
                        QNEWIC,QLQOUT,QICOUT,GRAV)
!
!
      implicit none
!!!#include <arch_specific.hf>
!
!
      REAL BOTERM,CONV,DQ,DZ
      REAL ENTERM,GRAV,G1,OLDQ,PPTDRG
      REAL QEST,QICE,QICOUT,QLIQ,QLQOUT,QNEW
      REAL QNEWIC,QNEWLQ,QTOT,RATE
      REAL RATIO3,RATIO4,WAVG,WTW
!
!AUTHOR
!          Kain and Fritsch   (1988)
!
!REVISION
! 001      Stephane Belair (1994)
! 002      A-M Leduc       (2001)  Adaptation to physics 3.7
!
!OBJECT
!          to perform the precipitation fallout
!
!ARGUMENTS
!
!          - INPUT/OUTPUT -
! QLIQ     liquid water in the cloud layer as input
!          liquid water after the precipitation fallout as output
! QICE     cloud ice in the cloud layer as input
!          cloud ice after the precipitation fallout as output
! WTW      W*W of the updraft as input
!          new vertical velocity with the condensation load
!          effect as output
!
!          - INPUT -
! DZ       row number
! BOTERM   buoyancy term for the updraft
! ENTERM   entrainment term for the updraft
! RATE     rate of precipitation fallout production in the
!          updraft (= 0.01)
! QNEWLQ   fresh liquid water condensate in the updraft
! QNEWIC   fresh cloud ice in the updraft
!
!          - Output -
! QLQOUT   precipitation fallout (liquid water)
! QICOUT   precipitation fallout (cloud ice)
!
!
!Notes
!  9/18/88...This precipitation fallout scheme is based on the scheme
!  used by Ogura and Cho (1973). Liquid water fallout from a parcel
!  is calculated using the equation dq=-rate*q*dt, but to simulate a
!  quasi-continuous process, and to eliminate a dependency on vertical
!  resolution, this is expressed as q=q*exp(-rate*dz).
!*
!
      QTOT=QLIQ+QICE
      QNEW=QNEWLQ+QNEWIC
!
!
!
!                              Estimate the vertical velocity so that
!                              an average vertical velocity can be
!                              calculated to estimate the time required
!                              ascent between model levels.
!
      QEST=0.5*(QTOT+QNEW)
!
!                              Solve for the vertical motion.
!

      G1=WTW+BOTERM-ENTERM-2.*GRAV*DZ*QEST/1.5
      IF(G1.LT.0.0)G1=0.
      WAVG=(SQRT(WTW)+SQRT(G1))/2.
!
!
!                              RATIO3 is the fraction of liquid water
!                              in fresh condensate.
!                              RATIO4 is the fraction of liquid water in
!                              the total amount of condensate involved
!                              in the precipitation process.
!                              Note that only 60% of the fresh condensate
!                              is allowed to participate in the conversion
!                              process.
!
      RATIO3=QNEWLQ/max(QNEW,1.E-10)
      QTOT=QTOT+0.6*QNEW
      OLDQ=QTOT
      RATIO4=(0.6*QNEWLQ+QLIQ)/max(QTOT,1.E-10)
!
!                             Apply the rate equation.
!
      if (abs(wavg) > epsilon(wavg)) then
         CONV = RATE*DZ/WAVG ! Compute conversion
         QTOT = QTOT*EXP(-CONV)
      else
         QTOT = 0.
      endif
!
!
!                             Determine the amount of precipitation that
!                             falls out of the updraft parcel at this
!                             level.
!
      DQ=OLDQ-QTOT
      QLQOUT=RATIO4*DQ
      QICOUT=(1.-RATIO4)*DQ
!
!
!                             Estimate the mean load of condensate on
!                             the updraft in the layer.
!                             Calculate the vertical velocity.
!
!                             PPTDRG = precipitation drag
!
      PPTDRG=0.5*(OLDQ+QTOT-0.2*QNEW)
!
!                             Apply again equation for vertical motion.
!
      WTW=WTW+BOTERM-ENTERM-2.*GRAV*DZ*PPTDRG/1.5
!
!
!                             Determine the new liquid water and ice
!                             concentrations including lost due to
!                             precipitation and gains from condensation.
!                             Don't forget to include the fraction 0.4qnew
!                             which was not involved in the conversion
!                             process.
!
      QLIQ=RATIO4*QTOT+RATIO3*0.4*QNEW
      QICE=(1.-RATIO4)*QTOT+(1.-RATIO3)*0.4*QNEW
      QNEWLQ=0.
      QNEWIC=0.
!
!
      RETURN
   END SUBROUTINE CONDLOAD_INFCONV

 
   subroutine condload_safe(QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ, &
        QNEWIC,QLQOUT,QICOUT,GRAV)
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      use ifport
#endif
      implicit none
!!!#include <arch_specific.hf>
      real BOTERM,CONV,DZ
      real ENTERM,GRAV,G1
      real QEST,QICE,QICOUT,QLIQ,QLQOUT,QNEW
      real QNEWIC,QNEWLQ,QTOT,RATE
      real RATIO3,RATIO4,WAVG,WTW
      integer(2) control,control_nodivzero,status
      !@author Kain and Fritsch   (1988)
      !
      !@revision
      ! 001      Stephane Belair (1994)
      ! 002      A-M Leduc       (2001)  Adaptation to physics 3.7
      !
      !@object perform the precipitation fallout
      !
      !@arguments
      !
      !          - INPUT/OUTPUT -
      ! QLIQ     liquid water in the cloud layer as input
      !          liquid water after the precipitation fallout as output
      ! QICE     cloud ice in the cloud layer as input
      !          cloud ice after the precipitation fallout as output
      ! WTW      W*W of the updraft as input
      !          new vertical velocity with the condensation load
      !          effect as output
      !
      !          - INPUT -
      ! DZ       row number
      ! BOTERM   buoyancy term for the updraft
      ! ENTERM   entrainment term for the updraft
      ! RATE     rate of precipitation fallout production in the
      !          updraft (= 0.01)
      ! QNEWLQ   fresh liquid water condensate in the updraft
      ! QNEWIC   fresh cloud ice in the updraft
      !
      !          - Output -
      ! QLQOUT   precipitation fallout (liquid water)
      ! QICOUT   precipitation fallout (cloud ice)
      !
      !@notes
      !  9/18/88...This precipitation fallout scheme is based on the scheme
      !  used by Ogura and Cho (1973). Liquid water fallout from a parcel
      !  is calculated using the equation dq=-rate*q*dt, but to simulate a
      !  quasi-continuous process, and to eliminate a dependency on vertical
      !  resolution, this is expressed as q=q*exp(-rate*dz).
      !*

      QTOT = QLIQ+QICE
      QNEW = QNEWLQ+QNEWIC

      ! Estimate the vertical velocity so that
      ! an average vertical velocity can be
      ! calculated to estimate the time required
      ! ascent between model levels.

      QEST = 0.5*(QTOT+QNEW)

      ! Solve for the vertical motion.


      G1 = WTW+BOTERM-ENTERM-2.*GRAV*DZ*QEST/1.5
      if (G1 < 0.0) G1=0.
      WAVG = (sqrt(WTW)+sqrt(G1))/2.

      ! Conversion

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      call getcontrolfpqq(control)
      control_nodivzero = ior(control, FPCW$ZERODIVIDE)
      call setcontrolfpqq(control_nodivzero)
      ! Run all possible div-by-zero calculations.
      CONV = RATE*DZ/WAVG
      RATIO3 = QNEWLQ/(QNEW+1.E-10)
      QTOT = QTOT+0.6*QNEW
      RATIO4 = (0.6*QNEWLQ+QLIQ)/(QTOT+1.E-10)
      call getstatusfpqq(status)
      call clearstatusfpqq()
      call setcontrolfpqq(control)
      if (iand(status,FPSW$ZERODIVIDE) > 0) then
         call condload_infconv(QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ, &
              QNEWIC,QLQOUT,QICOUT,GRAV)
      else
         call condload(QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ, &
              QNEWIC,QLQOUT,QICOUT,GRAV)
      endif
#else
      call condload(QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ, &
           QNEWIC,QLQOUT,QICOUT,GRAV)
#endif

      ! End of safe wrapper
      return
   end subroutine condload_safe


end module condload_mod
