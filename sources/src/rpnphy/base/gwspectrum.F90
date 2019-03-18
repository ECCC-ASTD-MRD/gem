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

!/@*
SUBROUTINE gwspectrum4(kproma , kbdim  , klev   &
                         ,s      , sh     , sexpk  &
                         ,shexpk , pressg , th     &
                         ,ptm1   , pum1   , pvm1   &
                         ,ptte   , pvol   , pvom   &
                         ,hflt   , g      , rd     &
                         ,tau    , rmscon , iheatcal &
                         ,kount, trnch               &
                         , std_p_prof, non_oro_pbot)
    USE mo_gwspectrum, ONLY: kstar, naz
    use phy_status, only: phy_error_L
    implicit none
#include <arch_specific.hf>

    ! scalar argument with intent(IN)
    INTEGER ,INTENT(in) :: kproma, kbdim, klev, iheatcal, hflt, kount, trnch 
    REAL,    INTENT(in) :: non_oro_pbot        ! pressure to compute bottom level

    !  Array arguments with intent(IN):
    ! Input 1D
    REAL*8 ,INTENT(IN) :: pressg(kproma)       ! Surface pressure (pascal)
    REAL, INTENT(IN) :: std_p_prof(klev)          ! Standard Pressure Profil (pascal)
    ! Input 2D
    REAL*8 ,INTENT(IN) :: th(kbdim,klev)       ! half level temperature
    ! Input constants
    REAL   ,INTENT(IN) :: g
    REAL   ,INTENT(IN) :: rd
    REAL   ,INTENT(IN) :: tau
    REAL   ,INTENT(IN) :: rmscon

    !  Array arguments with intent(InOut):
    ! - input/output 2d
    REAL*8 ,INTENT(INOUT) :: pum1(kbdim,klev)  ! zonal wind (t-dt)
    REAL*8 ,INTENT(INOUT) :: pvm1(kbdim,klev)  ! meridional wind (t-dt)
    REAL*8 ,INTENT(INOUT) :: ptm1(kbdim,klev)  ! temperature (t-dt)

    !  Array arguments with intent(Out):
    ! - output 2d
    REAL*8 ,INTENT(OUT) :: ptte(kbdim,klev)    ! tendency of temperature
    REAL*8 ,INTENT(OUT) :: pvol(kbdim,klev)    ! tendency of meridional wind
    REAL*8 ,INTENT(OUT) :: pvom(kbdim,klev)    ! tendency of zonal wind
    REAL*8              :: diffco(kbdim,klev)  ! turbulent diffusion coefficient
!
!Authors
!    
!   n. mcfarlane                cccma     may 1995
!   c. mclandress               ists      august 1995
!   m. charron                  mpi       2000-2001
!   e. manzini                  mpi       february 2002 (re-write, based on cccgwd)
!   th.schoenemeyer/h.schmidt   nec/mpi   july 2002 (optimized for vector architecture)
!   h. schmidt                  mpi       march 2003
!   m. charron                  rpn       june 2004
!
!
!Revision
! 001  L. Spacek (Sep 2008) - Density calulation uses th instead ptm1
! 002  A. Zadra & R. McTaggart-Cowan (June 2011) - Provide upper boundary index through interface
! 003  A. Plante (Sept. 2011) - Compute lower boundary index with the standard pressure profil std_p_prof
!
!Object
! 
!   Hines parameterization from ccc/mam (Hines, 1997a,b):       
!   physical tendencies of the prognostic variables u,v 
!   due to vertical transports by a broad band spectrum
!   of gravity waves. 
!
!   Note that diffusion coefficient and  heating rate 
!   only calculated if iheatcal = 1.
!
!   *gwspectrum* is called from *physc*.
! 
!Arguments
! 
!            - Input/Ouput -
! pum1     zonal wind (t-dt)
! pvm1     meridional wind (t-dt)
! ptm1     temperature (t-dt)
!
!          - Output -
! utendgw  zonal tend, gravity wave spectrum (m/s^2)
! pvol     tendency of meridional wind
! pvom     tendency of zonal wind
! diffco   turbulent diffusion coefficient
!
!          - Input -
! pressg   surface pressure (pascal)
! th       half  level temperature
! std_p_prof  STanDard Pressure PRoFil to get emiss_lev
!*@/

    !  Local arrays for ccc/mam hines gwd scheme:

    ! Important local parameter (passed to all subroutines):
    INTEGER, PARAMETER :: nazmth = 8  ! max azimuth array dimension size 

    ! * Vertical positioning arrays and work arrays:                       
    REAL*8, INTENT(in) :: s(kproma,klev), sh(kproma,klev), shexpk(kproma,klev), sexpk(kproma,klev)
    REAL*8 :: dttdsf,dttdzl

    REAL*8 :: utendgw(kproma,klev) ! zonal tend, gravity wave spectrum (m/s^2)
    REAL*8 :: vtendgw(kproma,klev) ! merid tend, gravity wave spectrum (m/s^2)
    REAL*8 :: ttendgw(kproma,klev) ! temperature tend, gravity wave spectrum (K/s)
    REAL*8 ::  flux_u(kproma,klev) ! zonal momentum flux (pascals)
    REAL*8 ::  flux_v(kproma,klev) ! meridional momentum flux (pascals) 

    REAL*8 :: uhs(kproma,klev)      ! zonal wind (m/s), input for hines param
    REAL*8 :: vhs(kproma,klev)      ! merid wind (m/s), input for hines param
    REAL*8 :: bvfreq(kproma,klev)   ! background brunt vassala frequency (rad/s)
    REAL*8 :: density(kproma,klev)  ! background density (kg/m^3)
    REAL*8 :: visc_mol(kproma,klev) ! molecular viscosity (m^2/s) 
    REAL*8 :: alt(kproma,klev)      ! background altitude (m)

    REAL*8 :: rmswind(kproma)        ! rms gravity wave  wind, lowest level (m/s) 
    REAL*8 :: anis(kproma,nazmth)    ! anisotropy factor (sum over azimuths = 1) 
    REAL*8 :: k_alpha(kproma,nazmth) ! horizontal wavenumber of each azimuth (1/m)
    LOGICAL :: lorms(kproma)       ! .true. for rmswind /=0 at launching level 

    REAL*8 :: m_alpha(kproma,klev,nazmth) ! cutoff vertical wavenumber (1/m)
    REAL*8 :: mmin_alpha(kproma,nazmth)   ! minumum value of m_alpha
    REAL*8 :: sigma_t(kproma,klev)        ! total rms gw wind (m/s)
    ! gw variances from orographic sources (for coupling to a orogwd)
    REAL*8 :: sigsqmcw(kproma,klev,nazmth), sigmatm(kproma,klev)

    !
    ! Local scalars:
    INTEGER  :: jk, jl
    INTEGER  :: levbot     ! gravity wave spectrum lowest level
    REAL*8     :: hscal, ratio


       !
       !--  Initialize the ccc/mam hines gwd scheme
       !

       utendgw(:,:) = 0.
       vtendgw(:,:) = 0.
       ttendgw(:,:) = 0.

       diffco(:,:) = 0

       flux_u(:,:) = 0.
       flux_v(:,:) = 0. 

       uhs(:,:) = 0.
       vhs(:,:) = 0.

       ! Wind variances form orographic gravity waves
       ! Note: the code is NOT fully implemeted for this case!

       sigsqmcw(:,:,:) = 0.
       sigmatm(:,:)    = 0.

!     * CALCULATE  B V FREQUENCY EVERYWHERE.

       DO jk=2,klev
         DO jl=1,kproma
           dttdsf=(th(jl,jk)/SHEXPK(jl,jk)-th(jl,jk-1)/SHEXPK(jl,jk-1)) &
                         /(SH(jl,jk)-SH(jl,jk-1))
           dttdsf=MIN(dttdsf, -5./S(jl,jk))
           dttdzl=-dttdsf*S(jl,jk)*g/(rd*ptm1(jl,jk))
           bvfreq(jl,jk)=SQRT(g*dttdzl*SEXPK(jl,jk)/ptm1(jl,jk))
         ENDDO
       ENDDO

       bvfreq(:,1)=bvfreq(:,2)

       DO jk=2,klev
         DO jl=1,kproma
           ratio=5.*LOG(S(jl,jk)/S(jl,jk-1))
           bvfreq(jl,jk) = (bvfreq(jl,jk-1) + ratio*bvfreq(jl,jk))/(1.+ratio)
         END DO
       END DO

       !     * altitude and density at bottom.

       alt(:,klev) = 0.

       DO jl=1,kproma
          hscal = rd * ptm1(jl,klev) / g
          density(jl,klev) = s(jl,klev) * pressg(jl) / (g*hscal)
       END DO

       !     * altitude and density at remaining levels.

       DO jk=klev-1,1,-1
          DO jl=1,kproma
             hscal = rd * th(jl,jk) / g
             alt(jl,jk) = alt(jl,jk+1) + hscal * LOG(s(jl,jk+1)/s(jl,jk))
             density(jl,jk) = s(jl,jk) * pressg(jl) / (rd * ptm1(jl,jk))
          END DO
       END DO

       !
       !     * set molecular viscosity to a very small value.
       !     * if the model top is greater than 100 km then the actual
       !     * viscosity coefficient could be specified here.


       DO jk=1,klev
          DO jl=1,kproma
             visc_mol(jl,jk) = 3.90E-7*ptm1(jl,jk)**.69 / density(jl,jk)
          ENDDO
       ENDDO

       ! use single value for azimuthal-dependent horizontal wavenumber:
       ! kstar = (old latitudinal dependence, introduce here if necessary)

       k_alpha(:,:) = kstar

       !     * defile bottom launch level (emission level of gws)
       levbot=-1
       DO jk=1,klev
          if(std_p_prof(jk)>non_oro_pbot)then
             levbot=jk-1
             exit
          endif
       ENDDO
       IF(levbot.lt.1)then
          write(6,1000)std_p_prof(1),std_p_prof(klev),non_oro_pbot
          call physeterror('gwspectrum', 'Problem with non_oro_pbot values, out of range')
          return
       ENDIF

       !     * initialize switch for column calculation

       lorms(:) = .FALSE.

       !     * background wind minus value at bottom launch level.

       DO jk=1,levbot
         DO jl=1,kproma 
           uhs(jl,jk) = pum1(jl,jk) - pum1(jl,levbot)
           vhs(jl,jk) = pvm1(jl,jk) - pvm1(jl,levbot)
         END DO
       END DO

       !     * specify root mean square wind at bottom launch level.

       DO jl=1,kproma 
          rmswind(jl) = DBLE(rmscon)
          anis(jl,:)   = 1./FLOAT(naz)
       END DO

       DO jl=1,kproma 
          IF (rmswind(jl) .GT. 0.0) THEN
             lorms(jl) = .TRUE.
          ENDIF
       END DO

       !
       !     * calculate gw tendencies (note that diffusion coefficient and
       !     * heating rate only calculated if iheatcal = 1).
       !

       CALL hines_extro4( kproma, klev, nazmth,                          & 
                          utendgw, vtendgw, ttendgw, diffco(1:kproma,:), &
                          flux_u, flux_v,                                & 
                          uhs, vhs, bvfreq, density, visc_mol, alt,      & 
                          rmswind, anis, k_alpha, sigsqmcw,              &
                          m_alpha,  mmin_alpha ,sigma_t, sigmatm,        & 
                          kount, trnch,                                  &
                          levbot, hflt, lorms, iheatcal)
       if (phy_error_L) return

!       DO jk=1, klev
!          DO jl=1,kproma
!             pum1(jl,jk)=pum1(jl,jk)+tau*utendgw(jl,jk)
!             pvm1(jl,jk)=pvm1(jl,jk)+tau*vtendgw(jl,jk)
!             ptm1(jl,jk)=ptm1(jl,jk)+tau*ttendgw(jl,jk)
!          ENDDO
!       ENDDO

       !   update tendencies: 
       !
       DO jk=1, klev
          DO jl=1,kproma
             pvom(jl,jk) = utendgw(jl,jk)
             pvol(jl,jk) = vtendgw(jl,jk)
             ptte(jl,jk) = ttendgw(jl,jk)
          END DO
       END DO
       !

    !     * end of hines calculations.

    !-----------------------------------------------------------------------
1000   FORMAT ( ' *****************************************', &
              / ' *****************************************', &
              / ' *                                       *', &
              / ' ***** ABORT ***** ABORT ***** ABORT *****', &
              / ' *                                       *', &
              / ' * IN S/R GWSPECTRUM                     *', &
              / ' * PROBLEM WITH NON_ORO_PBOT             *', &
              / ' * MUST BE BETWEEN THE FOLLOWING VALUES  *', & 
              / ' * ',E12.5,' AND ',E12.5,' Pa', '        *', & 
              / ' * GOT ',  E12.5,   '                    *', &
              / ' *                                       *', &
              / ' *****************************************', &
              / ' *****************************************')
       
    END SUBROUTINE gwspectrum4
