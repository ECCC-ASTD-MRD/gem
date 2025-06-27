
module diagno_clouds
   use debug_mod, only: init2nan
   use tdpack_const
   use phy_options
   use phybusidx
   use phymem, only: phyvar
   use vintage_nt, only: vintage_nt1
   implicit none
   private
   public :: diagno_clouds2

#include "phymkptr.hf"
#include "nocld.cdk"
#include "nbsnbl.cdk"
#include "cldop.cdk"
   
   real, parameter :: THIRD = 0.3333333 !#TODO test with 1./3. (bit pattern change)
   integer, parameter :: TOPC2 = 5000.  !#Top level for cloud tendencies (Pa) (we take 10xtopc in nocld... ad hoc)
   character(len=4), parameter :: OLIST(14) = (/ &
        'ECC ', 'ECCH', 'ECCM', 'ECCL', 'TCC ', 'NT  ', &
        'TCSH', 'TCSM', 'TCSL', 'TCZH', 'TCZM', 'TCZL', &
        'ECTP', 'ECTT'  &
        /)

contains

   !/@*
   subroutine diagno_clouds2(pvars, taucs, taucl,  &
        gz, cloud, &
        tt, sig, ps,  trnch, &
        ni, nkm1, nk)
      implicit none
!!!#include <arch_specific.hf>

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: trnch, ni, nkm1, nk
      real, intent(in), dimension(ni,nkm1,nbs) :: taucs
      real, intent(in), dimension(ni,nkm1,nbl) :: taucl
      real, intent(in), dimension(ni,nkm1) :: gz, cloud, tt, sig
      real, intent(in), dimension(ni)    :: ps

      !Object
      !        calculate true and effective cloud covers, cloud top pressure and temperature, calculate NT
      !
      !Arguments
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - output -
      ! ctp      cloud top pressure
      ! ctt      cloud top temperature
      ! ecc      effective cloud cover (nt)
      !          - input -
      ! taucs    cloud solar optical thickness
      ! taucl    cloud longwave optical thickness
      ! cldfrac  layer cloud amount (0. to 1.) (ni,nkm1)
      ! tt       layer temperature (k) (m,nkm1)
      ! sig      sigma levels (0. to 1.) (ni,nkm1; local sigma)
      ! ps       surface pressure (n/m2) (ni)
      ! trnch    number of the slice
      ! ni      number of profiles to compute
      ! m        first dimension of temperature (usually ni)
      ! nkm1       number of layers
      !
      !*@/

      real, pointer, dimension(:), contiguous :: ztopthw,ztopthi,znt
      real, pointer, dimension(:), contiguous :: ztcc,zecc,zeccl,zeccm,zecch
      real, pointer, dimension(:), contiguous :: ztcsl,ztcsm,ztcsh
      real, pointer, dimension(:), contiguous :: ztczl,ztczm,ztczh
      real, pointer, dimension(:), contiguous :: zctp,zctt
      real, pointer, dimension(:,:), contiguous :: zlwc
      integer :: ktopi(ni), k, i, ktop
      real :: press
      logical :: needoutput
      real, target :: dummy1Dni(ni)
      
      real, dimension(ni,nkm1) :: transmissint, trans_exp, cldfrac, mask
      logical, dimension(ni) :: top
      real, dimension(ni,nkm1) :: aird, rec_cdd, vs1
      real, dimension(ni) :: trmin, tmem, trmin2, tmem2
      real, dimension(ni,nk,nk) :: ff, ff2
      integer, dimension(ni) :: ih, ih2, ih3, ib, ib2, ib3    

      real, parameter :: THIRD = 0.3333333 !#TODO test with 1./3. (bit pattern change)
      integer :: k1, ip, l, km1, lmax
      real :: xnu, xnu2, mask1
      !----------------------------------------------------------------
      call init2nan(transmissint, trans_exp, cldfrac)
      call init2nan(aird, rec_cdd, vs1)
      call init2nan(trmin, tmem, trmin2, tmem2 )
      call init2nan(ff, ff2)

      needoutput = (etccdiagout .or. ISREQSTEPL(OLIST))
      if (.not.needoutput) return

#undef MKPTR1DNI
#undef MKPTR1D
#define MKPTR1DNI(PNAME,IDXV,VARS)       if (IDXV>0) then ; PNAME(1:VARS(IDXV)%meta%ni) => VARS(IDXV)%data(1:VARS(IDXV)%meta%ni) ; else ; PNAME => dummy1Dni ; endif
#define MKPTR1D(PNAME,IDXV,VARS)         MKPTR1DNI(PNAME,IDXV,VARS)
      
      MKPTR1D(zctp, ctp, pvars)
      MKPTR1D(zctt, ctt, pvars)
      MKPTR1D(zecc, ecc, pvars)
      MKPTR1D(zecch, ecch, pvars)
      MKPTR1D(zeccl, eccl, pvars)
      MKPTR1D(zeccm, eccm, pvars)
      MKPTR1D(znt, nt, pvars)
      MKPTR1D(ztcc, tcc, pvars)
      MKPTR1D(ztcsl, tcsl, pvars)
      MKPTR1D(ztcsm, tcsm, pvars)
      MKPTR1D(ztcsh, tcsh, pvars)
      MKPTR1D(ztczl, tczl, pvars)
      MKPTR1D(ztczm, tczm, pvars)
      MKPTR1D(ztczh, tczh, pvars)
      MKPTR1D(ztopthi, topthi, pvars)
      MKPTR1D(ztopthw, topthw, pvars)

      nullify(zlwc)
      if (stcond(1:3) /= 'MP_') then
         MKPTR2D(zlwc, lwc, pvars)
      endif

      ktopi = 1
      do k=1,nkm1-1
         do i=1,ni
            press = sig(i,k)*ps(i)
            if (press < TOPC) ktopi(i) = k
         enddo
      enddo
      ktop = minval(ktopi)  !#Note: this may cause bit pattern change with change of ptopo or runlgt
      
      !     diagnostics: cloud top pressure (ctp) and temperature (ctt)
      !     using the cloud optical depth at window region (band 6) to
      !     calculate the emissivity

      !     calcul des indices IH et IB pour nuages 2-D
      !     IH = niveau le plus pres de sigma=0.4
      !     IB = niveau le plus pres de sigma=0.7
      !     IH2 = niveau le plus pres de sigma=0.45 pour compatibilite avec meme variable pour era5
      !     IB2 = niveau le plus pres de sigma=0.8  pour compatibilite avec meme variable pour era5
      !     IH3 = niveau le plus pres de height=6.5km  for compatibility with the same variable for Calipso-GOCCP
      !     IB3 = niveau le plus pres de height=3.2km  for compatibility with the same variable for Calipso-GOCCP 

      do k = 1, nkm1
         do i = 1, ni
            trans_exp(i,k) = exp(- 1.64872 * taucl(i,k,6))
            if (sig(i,k) <= rad_siglim(4)) ih(i) = k
            if (sig(i,k) <= rad_siglim(3)) ib(i) = k
         enddo
      enddo

      if (etccdiag) then
         do k = 1, nkm1
            do i = 1, ni
               if (sig(i,k) <= rad_siglim(2)) ih2(i) = k
               if (sig(i,k) <= rad_siglim(1)) ib2(i) = k
               if (gz(i,k) >= rad_zlim(2)) ih3(i) = k
               if (gz(i,k) >= rad_zlim(1)) ib3(i) = k
            enddo
         enddo
      endif

      do i = 1, ni
         zctp (i)   = 110000.
         zctt (i)   = 310.
         top(i) = .true.
         transmissint(i,1) = 1. - cloud(i,1) * (1.-trans_exp(i,1) )
         if ( (1. - transmissint(i,1)) > 0.99 .and. top(i) ) then
            zctp(i) = sig(i,1)*ps(i)
            zctt(i) = tt(i,1)
            top(i) = .false.
         end if
      end do

      do k = 2, nkm1
         do i = 1, ni
            ! transmissint(i,k)=transmissint(i,k-1) * (1. - cloud(i,k) *
            !    1              exp (- 1.64872 * taucl(i,k,6)))
            transmissint(i,k)=transmissint(i,k-1) * (1. - cloud(i,k) * &
                 (1.-trans_exp(i,k) ) )
            if ( (1. - transmissint(i,k)) > 0.99 .and. top(i) ) then
               zctp(i) = sig(i,k)*ps(i)
               zctt(i) = tt(i,k)
               top(i) = .false.
            end if
         end do
      end do

      !...  compute total, high cloud, middle cloud and low cloud effective cloud cover (nt) as in radir7
      !     using the cloud optical depth at window region (band 6) to
      !     calculate the emissivity

      do k=2,nk
         k1 = max(k-2,1)
         do i=1,ni
            mask(i,k-1) = MASK_GE(cloud(i,k1), 0.01)
         enddo
      enddo
      if (etccdiag) then
         lmax = min(max(maxval(IH), maxval(IB), maxval(IH2), maxval(IB2), maxval(IH3), maxval(IB3)), nk-1)
      else
         lmax = min(max(maxval(IH), maxval(IB)), nk-1)
      endif
      
      do l=1,lmax
         ip = max(ktop, l+1)
         do i=1,ni
            tmem(i)=1.
            trmin(i)=1.
            tmem2(i)=1.
            trmin2(i)=1.
            ff(i,ip-1,l) = 1.
            ff2(i,ip-1,l) = 1.
         enddo
         do k=ip,nk
            km1 = k-1
            do i=1,ni
               xnu   = 1. - cloud(i,km1)*(1.-trans_exp(i,km1))
               xnu2  = 1. - cloud(i,km1)
               mask1 = 1. - mask(i,km1)
               tmem(i)   = mask(i,km1)*tmem(i)  + mask1*ff(i,km1,l)
               tmem2(i)  = mask(i,km1)*tmem2(i) + mask1*ff2(i,km1,l)
               trmin(i)  = mask(i,km1)*min(trmin(i),xnu)   + mask1*xnu
               trmin2(i) = mask(i,km1)*min(trmin2(i),xnu2) + mask1*xnu2
               ff(i,k,l) = tmem(i)  * trmin(i)
               ff2(i,k,l)= tmem2(i) * trmin2(i)
            enddo
         enddo
      enddo
      do i=1,ni
         zecc(i)  = 1. - ff(i,nk,1)
         zecch(i) = 1. - ff(i,IH(i),1)
         zeccm(i) = 1. - ff(i,IB(i),IH(i))
         zeccl(i) = 1. - ff(i,nk,IB(i))
      enddo
      
      if (ISREQSTEPL((/"TCC","NT "/)) .or. ISREQOUTL((/"TCCM", "NF  "/))) &
           ztcc = 1. - ff2(:,nk,1)

      if (ISREQSTEP("NT") .or. ISREQOUT("NF")) then
         if (stcond(1:3)=='MP_') then
            do i=1,ni
               znt(i) = ztcc(i)*(1.-exp(-0.1*(ztopthw(i) + ztopthi(i))))
            enddo
         else
            !#Note: vintage_nt modify cldfrac, need a copy
            cldfrac = cloud
            call vintage_nt1( &
                 tt, ps, sig, zlwc, &
                 cldfrac, znt, trnch, ni, nk, nkm1)
         endif
      endif


      if (etccdiag) then
         do i=1,ni
            ztcsh(i) = 1. - ff2(i,IH2(i),1)
            ztcsm(i) = 1. - ff2(i,IB2(i),IH2(i))
            ztcsl(i) = 1. - ff2(i,nk,IB2(i))
            ! to calculate 2-d cloud fraction with calipso-GOCCP height criteria
            ztczh(i) = 1. - ff2(i,IH3(i),1)   
            ztczm(i) = 1. - ff2(i,IB3(i),IH3(i))
            ztczl(i) = 1. - ff2(i,nk,IB3(i))
         enddo
      endif

      !----------------------------------------------------------------
      return
   end subroutine diagno_clouds2
end module diagno_clouds

