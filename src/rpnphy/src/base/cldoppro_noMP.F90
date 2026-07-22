
module cldoppro_noMP
   implicit none
   private
   public :: cldoppro_noMP1

contains


   !/@*
   subroutine cldoppro_noMP1(pvars,taucs, omcs, gcs, taucl, omcl, gcl, &
        liqwcin, icewcin, &
        liqwpin, icewpin, cldfrac, &
        tt, sig, ps, dontneedall, trnch, m, &
        ni, nkm1, nk)
      use debug_mod, only: init2nan
      use tdpack_const
      use phy_options
      use phybusidx
      use phymem, only: phyvar
      use series_mod, only: series_xst
      use ens_perturb, only: ens_spp_get
      use cldop_utils, only: cldop_rew, cldop_rei, cldop_proplw, cldop_propli, &
           cldop_propsw, cldop_propsi
      implicit none
!!!#include <arch_specific.hf>
#include "nbsnbl.cdk"

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, m, nkm1,nk, trnch
      real, intent(out), dimension(ni,nkm1,NBS) :: taucs, omcs, gcs
      real, intent(out), dimension(ni,nkm1,NBL) :: taucl, omcl, gcl
      real, intent(inout), dimension(ni,nkm1) :: liqwcin, icewcin
      real, intent(inout), dimension(ni,nkm1) :: liqwpin, icewpin
      real, intent(inout), dimension(ni,nkm1) :: cldfrac
      real, intent(in), dimension(ni,nkm1) :: tt, sig
      real, intent(in), dimension(ni)    :: ps
      logical, intent(in) :: dontneedall                   ! .true. = need only band1 in SW and band6 in LW

      !@Authors p. vaillancourt, d. talbot, j. li, rpn, cmc, cccma; (may 2006)
      !
      !Object
      !        calculate cloud optical properties for ccc radiative transfer scheme
      !
      !Arguments
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - output -
      ! taucs    cloud solar optical thickness
      ! omcs     cloud solar scattering albedo
      ! gcs      cloud solar asymmetry factor
      ! taucl    cloud longwave optical thickness
      ! omcl     cloud longwave scattering albedo
      ! gcl      cloud longwave asymmetry factor
      ! topthw   total integrated optical thickness of water in the visible
      ! topthi   total integrated optical thickness of ice in the visible
      !          - input -
      ! liqwpin  liquid water path in g/m2
      ! icewpin  solid water path in g/m2
      ! liqwcin  in-cloud liquid water content
      !          in kg water/kg air
      ! icewcin  in-cloud ice water content
      !          in kg water/kg air
      ! cldfrac  layer cloud amount (0. to 1.) (ni,nkm1)
      ! tt       layer temperature (k) (m,nkm1)
      ! sig      sigma levels (0. to 1.) (ni,nkm1; local sigma)
      ! ps       surface pressure (n/m2) (ni)
      ! mg       ground cover (ocean=0.0,land <= 1.)  (ni)
      ! ml       fraction of lakes (0.-1.) (ni)
      ! mrk2     Markov chains for stochastic perturbations   
      ! ni      number of profiles to compute
      ! m        first dimension of temperature (usually ni)
      ! nkm1       number of layers
      !
      !*@/

#include <rmnlib_basics.hf>
#include "phymkptr.hf"

      real, dimension(m,nkm1) :: rew
      real, dimension(m,nkm1) :: rei
      real, dimension(m,nkm1) :: vs1
      real, dimension(ni) :: reifac
      real, dimension(ni) :: rewfac

      include "nocld.cdk"

      integer :: i, j, k, swi, swf, lwi, lwf
      real :: rec_grav
      real :: epsilon, epsilon2, betan, betad
      real :: tausw, omsw, gsw, tausi, omsi, gsi, y1, y2, taulw
      real :: omlw, glw, tauli, omli, gli

      real, pointer, dimension(:), contiguous   :: ztlwp, ztiwp
      real, pointer, dimension(:,:), contiguous :: zlwcrad, ziwcrad, zcldrad, zmrk2
      real, pointer, dimension(:), contiguous :: zmg,zml,ztopthw,ztopthi,zdlat
      real, pointer, dimension(:,:), contiguous   :: zrewx, zreix
     !----------------------------------------------------------------

      MKPTR1D(zdlat, dlat, pvars)
      MKPTR1D(zmg, mg, pvars)
      MKPTR1D(zml, ml, pvars)
      MKPTR1D(ztlwp, tlwp, pvars)
      MKPTR1D(ztiwp, tiwp, pvars)
      MKPTR2D(zlwcrad, lwcrad, pvars)
      MKPTR2D(ziwcrad, iwcrad, pvars)
      MKPTR2D(zcldrad, cldrad, pvars)
      MKPTR1D(ztopthi, topthi, pvars)
      MKPTR1D(ztopthw, topthw, pvars)

      MKPTR2D(zmrk2, mrk2, pvars)
      MKPTR2Dm1(zrewx, rewx, pvars)
      MKPTR2Dm1(zreix, reix, pvars)


      !----------------------------------------------------------------
      call init2nan(rew, rei, vs1)
      call init2nan(reifac, rewfac)

      rec_grav=1./grav
      zrewx    = 0.0
      zreix    = 0.0


      ! impose coherence so that output profile variables are fully consistent with what RT has seen
      ! le premier if change le bit-pattern , je ne comprends pas pourquoi
      do k = 1,nkm1
         do i = 1,ni
            if (liqwpin(i,k) <= wpth .and. icewpin(i,k) <= wpth) &
                 cldfrac(i,k) = 0.0
            if (liqwpin(i,k) <= wpth ) then
               liqwpin(i,k) = 0.0
               liqwcin(i,k) = 0.0
            endif
            if (icewpin(i,k) <= wpth ) then
               icewpin(i,k) = 0.0
               icewcin(i,k) = 0.0
            endif
         end do
      end do

      ! initialize output fields
      do i =1, ni
         ztopthw(i) =  0.0
         ztopthi(i) =  0.0
         ztlwp(i)      = 0.0
         ztiwp(i)      = 0.0
         zlwcrad(i,nk) = 0.0
         ziwcrad(i,nk) = 0.0
         zcldrad(i,nk) = 0.0
      end do

      ! calculate output fields
      do k=1,nkm1
         do i=1,ni
            ztlwp(i)   = ztlwp(i)   + liqwpin(i,k)*cldfrac(i,k)
            ztiwp(i)   = ztiwp(i)   + icewpin(i,k)*cldfrac(i,k)
            zlwcrad(i,k) = liqwcin(i,k)*cldfrac(i,k)
            ziwcrad(i,k) = icewcin(i,k)*cldfrac(i,k)
            zcldrad(i,k) = cldfrac(i,k)
         enddo
      enddo
      !# conversion d'unites : tlwp et tiwp en kg/m2
      do i=1,ni
         ztlwp(i)   = ztlwp(i) * 0.001
         ztiwp(i)   = ztiwp(i) * 0.001
      enddo

      !# extraction pour diagnostics

      call series_xst(ziwcrad, 'iwcr', trnch)
      call series_xst(zlwcrad, 'lwcr', trnch)
      call series_xst(zcldrad, 'cldr', trnch)

      ! Effective radius for liquid water clouds
      rew(:,:) = cldop_rew(tt, liqwcin, sig, ps, zmg, zml, rew_const, rad_cond_rew, ni, nkm1)
      zrewx(:,:) = rew(:,:) * 1e-6

      ! Adjust the effective radius using stochastic perturbations
      rewfac(:) = ens_spp_get('rew_mult', zmrk2, default=1.)
      do k=1, nkm1
         rew(:,k) = rewfac(:) * rew(:,k)
      enddo

      ! Effective radius for ice clouds
      rei(:,:) = cldop_rei(tt, icewcin, sig, ps, zdlat, rei_const, rad_cond_rei, ni, nkm1)
      zreix(:,:) = rei(:,:) * 1e-6
    
      ! Adjust the effective radius using stochastic perturbations
      reifac(:) = ens_spp_get('rei_mult', zmrk2, default=1.)
      do k=1, nkm1
         rei(:,k) = reifac(:) * rei(:,k)
      enddo

! mask output effective radii where there is no cloud
      where (cldfrac(:,:) <= cldfth)
            zreix(:,:)= 0.0
            zrewx(:,:)= 0.0
      endwhere

      !----------------------------------------------------------------------
      !     cloud radiative properties for radiation.
      !     taucs, omcs, gcs (taucl, omcl, gcl): optical depth, single
      !     scattering albedo, asymmetry factor for solar (infrared).
      !     rew: effective radiu (in micrometer) for water cloud
      !     rei: effective radiu (in micrometer) for ice cloud
      !     dg: geometry length for ice cloud
      !     liqwcin  (icewcin): liquid water (ice) content (in gram / m^3)
      !     liqwpin (icewpin): liquid water (ice) path length (in gram / m^2)
      !     cloud: cloud fraction
      !     parameterization for water cloud:
      !     dobbie, etc. 1999, jgr, 104, 2067-2079
      !     lindner, t. h. and j. li., 2000, j. clim., 13, 1797-1805.
      !     parameterization for ice cloud:
      !     fu 1996, j. clim., 9, 2223-2337.
      !     fu et al. 1998 j. clim., 11, 2223-2337.
      ! NOTE:  In two Fu papers, opt prop are tested for DG from 15-130microns or REI from 10-90 microns(using dg=1.5395xrei)
      !        but in practice, rei must be below 70microns or the model crashes
      !----------------------------------------------------------------------

     if(dontneedall) then
      swi=1;swf=1;lwi=6;lwf=6
     else
      swi=1;swf=NBS;lwi=1;lwf=NBL
     endif

      DO_LOOP1: do j = swi, swf
         do k = 1, nkm1
            do i = 1, ni
               if (cldfrac(i,k) <= cldfth) then
                  taucs(i,k,j) =  0.
                  omcs(i,k,j)  =  0.
                  gcs(i,k,j)   =  0.
               else

                  ! Shortwave properties of liquid water clouds
                  call cldop_propsw(tausw, omsw, gsw, liqwpin(i,k), rew(i,k), j, liqwpin(i,k) > wpth)

                  ! Shortwave properties of ice clouds
                  call cldop_propsi(tausi, omsi, gsi, icewpin(i,k), rei(i,k), j, icewpin(i,k) > wpth)
                  
                  ! Combine shortwave properties
                  taucs(i,k,j)  =  tausw + tausi
                  if (taucs(i,k,j) > 0.0) then
                     y1          =  omsw * tausw
                     y2          =  omsi * tausi
                     omcs(i,k,j) = (y1 + y2) / taucs(i,k,j)
                     gcs (i,k,j) = (y1 * gsw + y2 * gsi) / (y1 + y2)
                  else
                     omcs(i,k,j) =  0.
                     gcs (i,k,j) =  0.
                  endif

                  ! calculate the optical depth for water and ice cloud in visible

                  if (j .eq. 1) then
                     ztopthw(i) =  ztopthw(i) + tausw
                     ztopthi(i) =  ztopthi(i) + tausi
                  endif

               endif
            enddo
         enddo
      enddo DO_LOOP1

      DO_LOOP2: do j = lwi, lwf
         do k = 1, nkm1
            do i = 1, ni
               if (cldfrac(i,k) <= cldfth) then
                  taucl(i,k,j) =  0.
                  omcl(i,k,j)  =  0.
                  gcl(i,k,j)   =  0.
               else

                  ! Longwave properties of liquid water clouds
                  call cldop_proplw(taulw, omlw, glw, liqwpin(i,k), rew(i,k), j, liqwpin(i,k) > wpth)

                  ! Longwave properties of ice clouds
                  call cldop_propli(tauli, omli, gli, icewpin(i,k), rei(i,k), j, icewpin(i,k) > wpth)

                  ! Combine longwave properties
                  taucl(i,k,j)   =  taulw + tauli
                  if (taucl(i,k,j) > 0.0) then
                     y1           =  omlw * taulw
                     y2           =  omli * tauli
                     omcl(i,k,j)  = (y1 + y2) / taucl(i,k,j)
                     gcl (i,k,j)  = (glw * y1 + gli * y2) / (y1 + y2)
                  else
                     omcl(i,k,j)  =  0.
                     gcl (i,k,j)  =  0.
                  endif
               endif
            enddo
         enddo
      enddo DO_LOOP2

      !----------------------------------------------------------------
      return
   end subroutine cldoppro_noMP1

end module cldoppro_noMP

