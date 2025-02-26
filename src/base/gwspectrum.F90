
module gwspectrum
   implicit none
   private
   public :: gwspectrum6

contains

!/@*
subroutine gwspectrum6(s, sh, pressg, th, &
     ptm1, pum1, pvm1, vtendgw, utendgw, hflt, &
     rmscons, std_p_prof, non_oro_pbot, &
     ni, nkm1)
   use, intrinsic :: iso_fortran_env, only: REAL32
   use tdpack_const, only: CAPPA, GRAV, RGASD
   use mo_gwspectrum, only: naz, rnaz
   use hines_extro, only: hines_extro5
   use phy_status, only: phy_error_L
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   ! scalar argument with intent(IN)
   integer, intent(in) :: ni, nkm1, hflt
   real,    intent(in) :: non_oro_pbot    ! Pressure to compute bottom level

   !  Array arguments with intent(IN):
   ! Input 1D
   real, intent(IN) :: pressg(ni)         ! Surface pressure (Pascals)
   real, intent(IN) :: std_p_prof(nkm1)   ! Standard Pressure Profil (Pascals)
   real, intent(IN) :: rmscons(ni)        ! RMS of wind speed at departure level (m/s)

   real, intent(IN) :: th(ni,nkm1)        ! Half level temperature

   real, intent(IN) :: pum1(ni,nkm1)      ! zonal wind (t-dt)
   real, intent(IN) :: pvm1(ni,nkm1)      ! meridional wind (t-dt)
   real, intent(IN) :: ptm1(ni,nkm1)      ! temperature (t-dt)

   real, intent(OUT) :: vtendgw(ni,nkm1)  ! tendency of meridional wind
   real, intent(OUT) :: utendgw(ni,nkm1)  ! tendency of zonal wind
   
   ! * Vertical positioning arrays and work arrays:
   real, intent(in) :: s(ni,nkm1)
   real, intent(in) :: sh(ni,nkm1)

   
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
  !          - Input -
   ! pressg   surface pressure (pascal)
   ! th       half  level temperature
   ! std_p_prof  STanDard Pressure PRoFil to get emiss_lev
   !*@/

   !  Local arrays for ccc/mam hines gwd scheme:

   ! Important local parameter (passed to all subroutines):
   real(REAL32) :: dttdsf,dttdzl

   real(REAL32) :: uhs(ni,nkm1)      ! zonal wind (m/s), input for hines param
   real(REAL32) :: vhs(ni,nkm1)      ! merid wind (m/s), input for hines param
   real(REAL32) :: bvfreq(ni,nkm1)   ! background brunt vassala frequency (rad/s)
   real(REAL32) :: density(ni,nkm1)  ! background density (kg/m^3)
   real(REAL32) :: visc_mol(ni,nkm1) ! molecular viscosity (m^2/s)
   real(REAL32) :: alt(ni,nkm1)      ! background altitude (m)

   real(REAL32) :: rmswind(ni)     ! rms gravity wave  wind, lowest level (m/s)

   real(REAL32) :: pressg_8(ni), s_8(ni,nkm1)

   real(REAL32) :: th_8(ni,nkm1), pum1_8(ni,nkm1), pvm1_8(ni,nkm1), ptm1_8(ni,nkm1)
   real(REAL32) :: sh_8(ni,nkm1), shexpk(ni,nkm1), sexpk(ni,nkm1)

   real :: tmpu(ni), tmpv(ni)
   integer :: idx(ni)
   
   ! Local scalars:
   integer  :: jk, jl, j2
   integer  :: nig      ! number of gatthers points (scope of gwspectrum)
   integer  :: levbot   ! gravity wave spectrum lowest level
   real(REAL32) :: hscal, ratio
   !--------------------------------------------------------------------

   !--  Initialize the ccc/mam hines gwd scheme
   
   nig = ni  !# nig: number of gathered points
   IFGATHER: if (any(rmscons <= 0.)) then

      !# Gather
      nig = 0
      idx(:) = 1
      do j2=1,ni
         if (rmscons(j2) > 0.) then
            nig = nig + 1
            idx(nig) = j2
            pressg_8(nig) = pressg(j2)
            rmswind(nig)  = rmscons(j2)
         endif
      enddo
      do jk=1,nkm1
         do jl=1,nig
            j2 = idx(jl)
            th_8(jl,jk)   = th(j2,jk)
            pum1_8(jl,jk) = pum1(j2,jk)
            pvm1_8(jl,jk) = pvm1(j2,jk)
            ptm1_8(jl,jk) = ptm1(j2,jk)
            s_8(jl,jk)    = s(j2,jk)
            sh_8(jl,jk)   = sh(j2,jk)
         enddo
      enddo

   else ! IFGATHER
      
      do jl=1,nig
         pressg_8(jl) = pressg(jl)
         rmswind(jl)  = rmscons(jl)
      enddo
      do jk=1,nkm1
         do jl=1,nig
            th_8(jl,jk)   = th(jl,jk)
            pum1_8(jl,jk) = pum1(jl,jk)
            pvm1_8(jl,jk) = pvm1(jl,jk)
            ptm1_8(jl,jk) = ptm1(jl,jk)
            s_8(jl,jk)    = s(jl,jk)
            sh_8(jl,jk)   = sh(jl,jk)
         enddo
      enddo
      
   endif IFGATHER
      
   !#TODO: use computed values from gwd (bit pattern change?)
   do jk=1,nkm1
      do jl=1,nig
         shexpk(jl,jk) = exp(CAPPA*log(sh_8(jl,jk)))
         sexpk(jl,jk)  = exp(CAPPA*log(s_8(jl,jk)))
!!$         sexpk(jl,jk)  = s_8(jl,jk)**CAPPA
!!$         shexpk(jl,jk) = sh_8(jl,jk)**CAPPA
      enddo
   enddo

   !     * CALCULATE  B V FREQUENCY EVERYWHERE.

   !#TODO: compute in gwd, shared in sgoflx?
   do jk=2,nkm1
      do jl=1,nig
!!$         sh =
!!$         shexpk =
!!$         sexpk = 
         dttdsf=(th_8(jl,jk)/shexpk(jl,jk)-th_8(jl,jk-1)/shexpk(jl,jk-1)) &
              /(sh_8(jl,jk)-sh_8(jl,jk-1))
         dttdsf=min(dttdsf, -5./s_8(jl,jk))
         dttdzl=-dttdsf*s_8(jl,jk)*GRAV/(RGASD*ptm1_8(jl,jk))
         bvfreq(jl,jk)=sqrt(GRAV*dttdzl*sexpk(jl,jk)/ptm1_8(jl,jk))
      enddo
   enddo

   bvfreq(1:nig,1) = bvfreq(1:nig,2)

   do jk=2,nkm1
      do jl=1,nig
         ratio=5.*log(s_8(jl,jk)/s_8(jl,jk-1))
         bvfreq(jl,jk) = (bvfreq(jl,jk-1) + ratio*bvfreq(jl,jk))/(1.+ratio)
      end do
   end do

   !     * altitude and density at bottom.
   alt(1:nig,nkm1) = 0.
   do jl=1,nig
      density(jl,nkm1) = s_8(jl,nkm1) * pressg_8(jl) / (RGASD * ptm1_8(jl,nkm1))
   end do

   !     * altitude and density at remaining levels.
   do jk=nkm1-1,1,-1
      do jl=1,nig
         hscal =  RGASD * th_8(jl,jk) / GRAV   !# TODO: th_8(jl,jk) * rgrav  
         alt(jl,jk) = alt(jl,jk+1) + hscal * log(s_8(jl,jk+1)/s_8(jl,jk))
         density(jl,jk) = s_8(jl,jk) * pressg_8(jl) / (RGASD * ptm1_8(jl,jk))
      end do
   end do


   !     * set molecular viscosity to a very small value.
   !     * if the model top is greater than 100 km then the actual
   !     * viscosity coefficient could be specified here.
   do jk=1,nkm1
      do jl=1,nig
         visc_mol(jl,jk) = 3.90E-7*ptm1_8(jl,jk)**.69 / density(jl,jk)
      enddo
   enddo

   !     * defile bottom launch level (emission level of gws)
   levbot=-1
   do jk=1,nkm1
      if(std_p_prof(jk)>non_oro_pbot)then
         levbot=jk-1
         exit
      endif
   enddo
   if(levbot.lt.1)then
      write(6,1000)std_p_prof(1),std_p_prof(nkm1),non_oro_pbot
      call physeterror('gwspectrum', 'Problem with non_oro_pbot values, out of range')
      return
   endif

   !     * background wind minus value at bottom launch level.
   do jk=levbot+1,nkm1
      do jl=1,nig
         uhs(jl,jk) = 0.
         vhs(jl,jk) = 0.
      enddo
   enddo
   do jk=1,levbot
      do jl=1,nig
         uhs(jl,jk) = pum1_8(jl,jk) - pum1_8(jl,levbot)
         vhs(jl,jk) = pvm1_8(jl,jk) - pvm1_8(jl,levbot)
      end do
   end do

   !     * calculate gw tendencies (note that diffusion coefficient and
   !     * heating rate only calculated if iheatcal = 1).
   call hines_extro5(ni, nig, nkm1, naz,  &
        utendgw, vtendgw,  &
        uhs, vhs, bvfreq, density, visc_mol, alt,  &
        rmswind,  &
        levbot, hflt)
   if (phy_error_L) return

   IFSCATTER: if (ni /= nig) then
      do jk=1,nkm1
         tmpu(:) = 0.
         tmpv(:) = 0.
         do jl=1,nig
            j2 = idx(jl)
            tmpu(j2) = utendgw(jl,jk)
            tmpv(j2) = vtendgw(jl,jk)
         enddo
         do j2=1,ni
            utendgw(j2,jk) = tmpu(j2)
            vtendgw(j2,jk) = tmpv(j2)
         enddo
      enddo
   endif IFSCATTER
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

end module gwspectrum
