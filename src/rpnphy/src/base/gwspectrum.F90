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
subroutine gwspectrum6(s, sh, sexpk, shexpk, pressg, th, &
     ptm1, pum1, pvm1, ptte, pvol, pvom, hflt, &
     g, rd, rmscons, iheatcal, std_p_prof, non_oro_pbot, &
     ni, nk)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use mo_gwspectrum, only: kstar, naz
   use phy_status, only: phy_error_L
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   ! scalar argument with intent(IN)
   integer, intent(in) :: ni, nk, iheatcal, hflt
   real,    intent(in) :: non_oro_pbot          ! Pressure to compute bottom level

   !  Array arguments with intent(IN):
   ! Input 1D
   real(REAL64), intent(IN) :: pressg(ni)       ! Surface pressure (Pascals)
   real, intent(IN) :: std_p_prof(nk)           ! Standard Pressure Profil (Pascals)
   real, intent(IN) :: rmscons(ni)              ! RMS of wind speed at departure level (m/s)
   ! Input 2D
   real(REAL64), intent(IN) :: th(ni,nk)        ! Half level temperature
   ! Input constants
   real, intent(IN) :: g
   real, intent(IN) :: rd

   !  Array arguments with intent(InOut):
   ! - input/output 2d
   real(REAL64), intent(INOUT) :: pum1(ni,nk)   ! zonal wind (t-dt)
   real(REAL64), intent(INOUT) :: pvm1(ni,nk)   ! meridional wind (t-dt)
   real(REAL64), intent(INOUT) :: ptm1(ni,nk)   ! temperature (t-dt)

   !  Array arguments with intent(Out):
   ! - output 2d
   real(REAL64), intent(OUT) :: ptte(ni,nk)     ! tendency of temperature
   real(REAL64), intent(OUT) :: pvol(ni,nk)     ! tendency of meridional wind
   real(REAL64), intent(OUT) :: pvom(ni,nk)     ! tendency of zonal wind
   real(REAL64)              :: diffco(ni,nk)   ! turbulent diffusion coefficient

   !@Authors
   !   n. mcfarlane                cccma     may 1995
   !   c. mclandress               ists      august 1995
   !   m. charron                  mpi       2000-2001
   !   e. manzini                  mpi       february 2002 (re-write, based on cccgwd)
   !   th.schoenemeyer/h.schmidt   nec/mpi   july 2002 (optimized for vector architecture)
   !   h. schmidt                  mpi       march 2003
   !   m. charron                  rpn       june 2004

   !@Revision
   ! 001  L. Spacek (Sep 2008) - Density calulation uses th instead ptm1
   ! 002  A. Zadra & R. McTaggart-Cowan (June 2011) - Provide upper boundary index through interface
   ! 003  A. Plante (Sept. 2011) - Compute lower boundary index with the standard pressure profil std_p_prof

   !@Object
   !   Hines parameterization from ccc/mam (Hines, 1997a,b):
   !   physical tendencies of the prognostic variables u,v
   !   due to vertical transports by a broad band spectrum
   !   of gravity waves.
   !   Note that diffusion coefficient and  heating rate
   !   only calculated if iheatcal = 1.
   !   *gwspectrum* is called from *physc*.

   !@Arguments
   !            - Input/Ouput -
   ! pum1     zonal wind (t-dt)
   ! pvm1     meridional wind (t-dt)
   ! ptm1     temperature (t-dt)
   !          - Output -
   ! utendgw  zonal tend, gravity wave spectrum (m/s^2)
   ! pvol     tendency of meridional wind
   ! pvom     tendency of zonal wind
   ! diffco   turbulent diffusion coefficient
   !          - Input -
   ! pressg   surface pressure (pascal)
   ! th       half  level temperature
   ! std_p_prof  STanDard Pressure PRoFil to get emiss_lev
   !*@/

   !  Local arrays for ccc/mam hines gwd scheme:

   ! Important local parameter (passed to all subroutines):
   integer, parameter :: nazmth = 8  ! max azimuth array dimension size

   ! * Vertical positioning arrays and work arrays:
   real(REAL64), intent(in) :: s(ni,nk), sh(ni,nk), shexpk(ni,nk), sexpk(ni,nk)


   real(REAL64) :: dttdsf,dttdzl

   real(REAL64) :: utendgw(ni,nk) ! zonal tend, gravity wave spectrum (m/s^2)
   real(REAL64) :: vtendgw(ni,nk) ! merid tend, gravity wave spectrum (m/s^2)
   real(REAL64) :: ttendgw(ni,nk) ! temperature tend, gravity wave spectrum (K/s)
   real(REAL64) ::  flux_u(ni,nk) ! zonal momentum flux (pascals)
   real(REAL64) ::  flux_v(ni,nk) ! meridional momentum flux (pascals)

   real(REAL64) :: uhs(ni,nk)      ! zonal wind (m/s), input for hines param
   real(REAL64) :: vhs(ni,nk)      ! merid wind (m/s), input for hines param
   real(REAL64) :: bvfreq(ni,nk)   ! background brunt vassala frequency (rad/s)
   real(REAL64) :: density(ni,nk)  ! background density (kg/m^3)
   real(REAL64) :: visc_mol(ni,nk) ! molecular viscosity (m^2/s)
   real(REAL64) :: alt(ni,nk)      ! background altitude (m)

   real(REAL64) :: rmswind(ni)        ! rms gravity wave  wind, lowest level (m/s)
   real(REAL64) :: anis(ni,nazmth)    ! anisotropy factor (sum over azimuths = 1)
   real(REAL64) :: k_alpha(ni,nazmth) ! horizontal wavenumber of each azimuth (1/m)
   logical :: lorms(ni)       ! .true. for rmswind /=0 at launching level

   real(REAL64) :: m_alpha(ni,nk,nazmth) ! cutoff vertical wavenumber (1/m)
   real(REAL64) :: mmin_alpha(ni,nazmth)   ! minumum value of m_alpha
   real(REAL64) :: sigma_t(ni,nk)        ! total rms gw wind (m/s)
   ! gw variances from orographic sources (for coupling to a orogwd)
   real(REAL64) :: sigsqmcw(ni,nk,nazmth), sigmatm(ni,nk)

   ! Local scalars:
   integer  :: jk, jl
   integer  :: levbot     ! gravity wave spectrum lowest level
   real(REAL64) :: hscal, ratio

   !--  Initialize the ccc/mam hines gwd scheme

   utendgw(:,:) = 0.
   vtendgw(:,:) = 0.
   ttendgw(:,:) = 0.

   diffco(:,:) = 0

   flux_u(:,:) = 0.
   flux_v(:,:) = 0.

   uhs(:,:) = 0.
   vhs(:,:) = 0.

   ! Wind variances form orographic gravity waves
   ! Note: the code is NOT fully implemeted for this case

   sigsqmcw(:,:,:) = 0.
   sigmatm(:,:)    = 0.

   !     * CALCULATE  B V FREQUENCY EVERYWHERE.

   do jk=2,nk
      do jl=1,ni
         dttdsf=(th(jl,jk)/SHEXPK(jl,jk)-th(jl,jk-1)/SHEXPK(jl,jk-1)) &
              /(SH(jl,jk)-SH(jl,jk-1))
         dttdsf=min(dttdsf, -5./S(jl,jk))
         dttdzl=-dttdsf*S(jl,jk)*g/(rd*ptm1(jl,jk))
         bvfreq(jl,jk)=sqrt(g*dttdzl*SEXPK(jl,jk)/ptm1(jl,jk))
      enddo
   enddo

   bvfreq(:,1)=bvfreq(:,2)

   do jk=2,nk
      do jl=1,ni
         ratio=5.*log(S(jl,jk)/S(jl,jk-1))
         bvfreq(jl,jk) = (bvfreq(jl,jk-1) + ratio*bvfreq(jl,jk))/(1.+ratio)
      end do
   end do

   !     * altitude and density at bottom.

   alt(:,nk) = 0.

   do jl=1,ni
      hscal = rd * ptm1(jl,nk) / g
      density(jl,nk) = s(jl,nk) * pressg(jl) / (g*hscal)
   end do

   !     * altitude and density at remaining levels.

   do jk=nk-1,1,-1
      do jl=1,ni
         hscal = rd * th(jl,jk) / g
         alt(jl,jk) = alt(jl,jk+1) + hscal * log(s(jl,jk+1)/s(jl,jk))
         density(jl,jk) = s(jl,jk) * pressg(jl) / (rd * ptm1(jl,jk))
      end do
   end do


   !     * set molecular viscosity to a very small value.
   !     * if the model top is greater than 100 km then the actual
   !     * viscosity coefficient could be specified here.


   do jk=1,nk
      do jl=1,ni
         visc_mol(jl,jk) = 3.90E-7*ptm1(jl,jk)**.69 / density(jl,jk)
      enddo
   enddo

   ! use single value for azimuthal-dependent horizontal wavenumber:
   ! kstar = (old latitudinal dependence, introduce here if necessary)

   k_alpha(:,:) = kstar

   !     * defile bottom launch level (emission level of gws)
   levbot=-1
   do jk=1,nk
      if(std_p_prof(jk)>non_oro_pbot)then
         levbot=jk-1
         exit
      endif
   enddo
   if(levbot.lt.1)then
      write(6,1000)std_p_prof(1),std_p_prof(nk),non_oro_pbot
      call physeterror('gwspectrum', 'Problem with non_oro_pbot values, out of range')
      return
   endif

   !     * initialize switch for column calculation

   lorms(:) = .false.

   !     * background wind minus value at bottom launch level.

   do jk=1,levbot
      do jl=1,ni
         uhs(jl,jk) = pum1(jl,jk) - pum1(jl,levbot)
         vhs(jl,jk) = pvm1(jl,jk) - pvm1(jl,levbot)
      end do
   end do

   !     * specify root mean square wind at bottom launch level.

   do jl=1,ni
      rmswind(jl) = dble(rmscons(jl))
      anis(jl,:)   = 1./FLOAT(naz)
   end do

   do jl=1,ni
      if (rmswind(jl) .gt. 0.0) then
         lorms(jl) = .true.
      endif
   end do


   !     * calculate gw tendencies (note that diffusion coefficient and
   !     * heating rate only calculated if iheatcal = 1).


   call hines_extro5(ni, nk, nazmth,             &
        utendgw, vtendgw, ttendgw, diffco(1:ni,:), &
        flux_u, flux_v,                                &
        uhs, vhs, bvfreq, density, visc_mol, alt,      &
        rmswind, anis, k_alpha, sigsqmcw,              &
        m_alpha,  mmin_alpha ,sigma_t, sigmatm,        &
        levbot, hflt, lorms, iheatcal)
   if (phy_error_L) return

   !       DO jk=1, nk
   !          DO jl=1,ni
   !             pum1(jl,jk)=pum1(jl,jk)+tau*utendgw(jl,jk)
   !             pvm1(jl,jk)=pvm1(jl,jk)+tau*vtendgw(jl,jk)
   !             ptm1(jl,jk)=ptm1(jl,jk)+tau*ttendgw(jl,jk)
   !          ENDDO
   !       ENDDO

   !   update tendencies:

   do jk=1, nk
      do jl=1,ni
         pvom(jl,jk) = utendgw(jl,jk)
         pvol(jl,jk) = vtendgw(jl,jk)
         ptte(jl,jk) = ttendgw(jl,jk)
      end do
   end do


   !     * end of hines calculations.

   !-----------------------------------------------------------------------
1000 format ( ' *****************************************', &
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

end subroutine gwspectrum6
