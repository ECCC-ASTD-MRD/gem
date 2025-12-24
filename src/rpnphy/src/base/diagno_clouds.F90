
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
   ! liste complète des variables calculées par diagno_clouds 
   character(len=4), parameter :: ILIST(17) = (/ &
        'ECC ', 'ECCH', 'ECCM', 'ECCL', 'TCC ', 'NT  ', &
        'TCSH', 'TCSM', 'TCSL', 'TCZH', 'TCZM', 'TCZL', &
        'ECTP', 'ECTT','TCCF', 'HCBR', 'HCBC'           &
        /)
   ! liste complète des variables moyennes ou accum associées à la liste ci-haut
   character(len=4), parameter :: ALIST(8) = (/ &
        'TCCM', 'NF  ',                                 &
        'TSHM', 'TSMM', 'TSLM', 'TZHM', 'TZMM', 'TZLM' &
        /)

   !liste des couvertures par couches (exclus les total)
   character(len=4), parameter :: TILIST(7) = (/ &
        'TCSH', 'TCSM', 'TCSL', 'TCZH', 'TCZM', 'TCZL', 'TCCF' &
        /)
   !liste des couvertures MOYHR
   character(len=4), parameter :: TMLIST(6) = (/ &
        'TSHM', 'TSMM', 'TSLM', 'TZHM', 'TZMM', 'TZLM' &
        /)

contains

   !/@*
   subroutine diagno_clouds2(pvars, taucsv, taucl6,  &
        gz, cloud, &
        tt, sig, ps,  trnch, &
        ni, nkm1, nk)
      implicit none
!!!#include <arch_specific.hf>

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: trnch, ni, nkm1, nk
      real, intent(in), dimension(ni,nkm1) :: taucsv
      real, intent(in), dimension(ni,nkm1) :: taucl6
      real, intent(in), dimension(ni,nkm1) :: gz, cloud, tt, sig
      real, intent(in), dimension(ni)      :: ps

      !Object
      !        calculate true and effective cloud covers, cloud top pressure and temperature, calculate NT
      !...Effective cloud covers: cloud fraction is weigted by transmittance calculated from the cloud optical depth in window region (band 6) 
      !
      !Arguments
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - output -
      ! tcc      true cloud cover (total, low, middle, high) 
      ! tccf     true cloud cover with filtering
      ! tcsl,tcsm,tcsh      true cloud covers (low, middle, high with sigma limits) 
      ! tczl,tczm,tczh      true cloud covers (low, middle, high with height limits) 
      ! ecc/l,m,h      effective cloud covers (total, low, middle, high) 
      ! hcbr     cloud base height 
      ! hcbc     cloud base height with filtering compatible with ceilometer observations
      ! ctp      cloud top pressure
      ! ctt      cloud top temperature
      !          - input -
      ! taucsv    cloud solar optical thickness in visible
      ! taucl    cloud longwave optical thickness
      ! cloud    layer cloud amount (0. to 1.) (ni,nkm1)
      ! tt       layer temperature (k) (m,nkm1)
      ! sig      sigma levels (0. to 1.) (ni,nkm1; local sigma)
      ! ps       surface pressure (n/m2) (ni)
      ! trnch    number of the slice
      ! ni      number of profiles to compute
      ! nkm1       number of layers
      !
      !*@/

      real, pointer, dimension(:), contiguous :: ztopthw, ztopthi, znt
      real, pointer, dimension(:), contiguous :: ztcc, zecc, zeccl, zeccm, zecch
      real, pointer, dimension(:), contiguous :: ztcsl, ztcsm, ztcsh
      real, pointer, dimension(:), contiguous :: ztczl, ztczm, ztczh
      real, pointer, dimension(:), contiguous :: zctp, zctt, zhcbc, zhcbr
      real, pointer, dimension(:,:), contiguous :: zlwc, ztccf

      integer :: ktopi(ni), k, i, ktop
      real :: press
      logical :: needoutput
      real, target :: dummy1Dni(ni)
      
      real, dimension(ni,nkm1) :: transmissint, trans_exp,  mask
      logical, dimension(ni) :: top
      real, dimension(ni,nkm1) :: aird, rec_cdd, vs1
      real, dimension(ni,nk) :: tausum
      real, dimension(ni) :: trmin, tmem, trmin2, tmem2
      real, dimension(ni,nk,nk) :: ff, ff2
      integer, dimension(ni) :: ih, ih2, ih3, ib, ib2, ib3    
      integer, dimension(ni) :: filterkc
      integer, dimension(ni,6) :: filterkh
      integer, dimension(ni,3) :: filterzk
      logical :: keeplooking, foundbase

      integer :: k1, ip, l, km1, lmax, m
      real :: xnu, xnu2, mask1, tauth
      !----------------------------------------------------------------
      call init2nan(transmissint, trans_exp)
      call init2nan(aird, rec_cdd, vs1)
      call init2nan(trmin, tmem, trmin2, tmem2 )
      call init2nan(ff, ff2)

      needoutput = (ISREQOUTL(ALIST) .or. ISREQSTEPL(ILIST))
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
      MKPTR2D(ztccf, tccf, pvars)
      MKPTR1D(ztcsl, tcsl, pvars)
      MKPTR1D(ztcsm, tcsm, pvars)
      MKPTR1D(ztcsh, tcsh, pvars)
      MKPTR1D(ztczl, tczl, pvars)
      MKPTR1D(ztczm, tczm, pvars)
      MKPTR1D(ztczh, tczh, pvars)
      MKPTR1D(ztopthi, topthi, pvars)
      MKPTR1D(ztopthw, topthw, pvars)
      MKPTR1D(zhcbc, hcbc, pvars)
      MKPTR1D(zhcbr, hcbr, pvars)

      nullify(zlwc)
      if (stcond(1:3) /= 'MP_') then
         MKPTR2D(zlwc, lwc, pvars)
      endif

      !     calcul des indices IH et IB pour nuages 2-D
      !     IH = niveau le plus pres de sigma=0.4 (coincé avec ça parce que eccl,eccm,ecch sort en ops)
      !     IB = niveau le plus pres de sigma=0.7 (coincé avec ça parce que eccl,eccm,ecch sort en ops)
      !     IH2 = niveau le plus pres de sigma=0.45 pour compatibilite avec meme variable pour era5 (sert en mheep)
      !     IB2 = niveau le plus pres de sigma=0.8  pour compatibilite avec meme variable pour era5 (sert en mheep)
      !     IH3 = niveau le plus pres de height=6.5km  for compatibility with the same variable for Calipso-GOCCP (sert en mheep)
      !     IB3 = niveau le plus pres de height=3.2km  for compatibility with the same variable for Calipso-GOCCP (sert en mheep)

      ktopi = 1
      do k=1,nkm1-1
         do i=1,ni
            press = sig(i,k)*ps(i)
            if (press < TOPC) ktopi(i) = k
         enddo
      enddo
      ktop = minval(ktopi)  !#Note: this may cause bit pattern change with change of ptopo or runlgt

      ! could use ktop as starting index for k everywhere but must manage k-1 

      ! need following for ctp,ctt, and effective cloud covers
      do k = 1, nkm1
         do i = 1, ni
            trans_exp(i,k) = exp(-1.64872 * taucl6(i,k))
         enddo
      enddo

      filterzk = 1               
      IF_TCFF: if (ISREQSTEPL((/"TCCF","HCBC"/))) then
         do k = 2, nkm1
            do i = 1, ni
               if (gz(i,k) >= 7000.) filterzk(i,1) = k+1 ! find first level below rad_zth 
               if (gz(i,k) >= 5800.) filterzk(i,2) = k+1 ! find first level below rad_zth 
               if (gz(i,k) >= 3660.) filterzk(i,3) = k+1 ! find first level below rad_zth 
            enddo
         enddo

         ! filter lower cloud layers with optical thickness < threshold (approximates detectability threshold by human eye or ceilometer)
         ! filter clouds above a height threshold(approximates height range limit of CL31 ceilometers)
         ! filterkc: filtering with optical thrershold of 0.01, k level of clear level just above highest cloud layer that must be filtered 
         ! filterzk: k level corresponding to chosen height threshold , filter above this value
         ! filterkh: filtering only as a fct of optical thickness, k level of clear level just above highest cloud layer that must be filtered 

         ! CEILOMETER: with hard-coded rad_tauth=.01 
         filterkc = nkm1+1
         zhcbc    = 20000. ! if no cloud base
         tausum   = 0.0
         do i = 1, ni
            keeplooking = .true.
            foundbase   = .false.
            do k = nkm1,1,-1
               if (taucsv(i,k) > tiny(1.0)) then                            ! level k is cloudy
                  tausum(i,k) = tausum(i,k+1)+ taucsv(i,k)                  ! sum optical thickness 
                  if (.not.foundbase) then 
                     zhcbc(i)  = gz(i,k)                                    ! intercept first base height before knowing if layer will be filtered
                     foundbase = .true.                         
                  endif
               else                                                         ! level k is clear
                  tausum(i,k) = 0.                                            ! clear layer found so reinit tausum   
                  if (tausum(i,k+1) >= .01) then                            ! cloud layer below is visible
                     keeplooking = .false.                                  ! therefore, stop looking
                  elseif (keeplooking .and. tausum(i,k+1) > tiny(1.0)) then ! cloud layer is NOT visible and had not found a visible one yet
                     filterkc(i) = k                                        ! cloud layer must be filtered, therefore increase filterk
                     foundbase = .false.                                    ! start looking for base again
                  endif
               endif
            enddo
            if (zhcbc(i) > 3660.) zhcbc(i) = 20000.
         enddo

         ! Find filterkh for  5 thresholds of (0.01-0.05) every 0.01 
         tauth    = 0.0
         filterkh = nkm1+1
         do l = 2,6  ! loop on 5 thresholds
            tausum = 0.0  
            tauth = tauth+.01
            do i = 1, ni  !#TODO: invert loop order
               keeplooking = .true.
               foundbase   = .false. 
               do k = nkm1,1,-1
                  if (taucsv(i,k) > tiny(1.0)) then                   ! inside a cloud layer
                     tausum(i,k) = tausum(i,k+1) + taucsv(i,k)        ! sum optical thickness
                     if (.not.foundbase) then 
                        foundbase = .true.                         
                     endif
                  else                                       
                     tausum(i,k) = 0.                            ! clear layer so reinit tausum   
                     if (tausum(i,k+1) >= tauth) then            ! cloud layer below is visible
                        keeplooking = .false.                    ! therefore, stop looking
                     elseif (keeplooking .and. tausum(i,k+1) > tiny(1.0)) then   ! cloud layer is NOT visible and had not found a visible one yet
                        filterkh(i,l) = k                                        ! cloud layer must be filtered, therefore increase filterk
                        foundbase     = .false.                                  ! start looking for base again
                     endif
                  endif
               enddo
            enddo
         enddo
      endif IF_TCFF

      do k=2,nk
         k1 = max(k-2,1)
         do i=1,ni
            mask(i,k-1) = MASK_GE(cloud(i,k1), 0.01)
         enddo
      enddo

      ! combine all code for effective cloud covers
      ! if only ecc is asked, lmax=ktop
      lmax = ktop
      if (ISREQSTEPL((/"ECCL",'ECCM','ECCH'/))) then
         do k = 1, nkm1
            do i = 1, ni
               if (sig(i,k) <= rad_siglim(4)) ih(i) = k
               if (sig(i,k) <= rad_siglim(3)) ib(i) = k
            enddo
         enddo
         lmax = min(max(maxval(IH), maxval(IB)), nk-1)
!!$      elseif (ISREQSTEP("ECC")) then 
!!$         lmax = ktop
      endif
      
      IF_FF: if (ISREQSTEPL((/"ECC ","ECCL",'ECCM','ECCH'/))) then
         do l=ktop,lmax
            ip = l+1
            do i=1,ni
               tmem(i)  = 1.
               trmin(i) = 1.
               ff(i,ip-1,l) = 1.
            enddo
            do k=ip,nk
               km1 = k-1
               do i=1,ni
                  xnu   = 1. - cloud(i,km1)*(1.-trans_exp(i,km1))
                  mask1 = 1. - mask(i,km1)
                  tmem(i)   = mask(i,km1)*tmem(i)  + mask1*ff(i,km1,l)
                  trmin(i)  = mask(i,km1)*min(trmin(i),xnu)   + mask1*xnu
                  ff(i,k,l) = tmem(i)  * trmin(i)
               enddo
            enddo
         enddo
      endif IF_FF
      
      if (ISREQSTEP("ECC")) &
           zecc = 1. - ff(:,nk,ktop)
      if (ISREQSTEPL((/"ECCL",'ECCM','ECCH'/))) then
         do i = 1, ni
            zecch(i) = 1. - ff(i,IH(i),ktop)
            zeccm(i) = 1. - ff(i,IB(i),IH(i))
            zeccl(i) = 1. - ff(i,nk,IB(i))
         enddo
      endif

      ! combine all code for true cloud covers
      lmax = ktop
      if (ISREQSTEPL(TILIST) .or. ISREQOUTL(TMLIST)) then
         do k = 1, nkm1
            do i = 1, ni
               if (sig(i,k) <= rad_siglim(2)) ih2(i) = k
               if (sig(i,k) <= rad_siglim(1)) ib2(i) = k
               if (gz(i,k)  >= rad_zlim(2))   ih3(i) = k
               if (gz(i,k)  >= rad_zlim(1))   ib3(i) = k
            enddo
         enddo
         lmax = min(max(maxval(filterzk),maxval(IH2), maxval(IB2), maxval(IH3), maxval(IB3)), nk-1)
!!$      elseif (ISREQSTEPL((/"TCC","NT "/)) .or. ISREQOUTL((/"TCCM","NF  "/))) then 
!!$         lmax = ktop
      endif
      
      IF_FF2: if (ISREQSTEPL(TILIST) .or. ISREQOUTL(TMLIST) .or. ISREQSTEPL((/"TCC","NT "/)) .or. ISREQOUTL((/"TCCM","NF  "/))) then
         do l=ktop,lmax
            ip = l+1
            do i=1,ni
               tmem2(i) = 1.
               trmin2(i) = 1.
               ff2(i,ip-1,l) = 1.
            enddo
            do k=ip,nk
               km1 = k-1
               do i=1,ni
                  xnu2  = 1. - cloud(i,km1)
                  mask1 = 1. - mask(i,km1)
                  tmem2(i)  = mask(i,km1)*tmem2(i) + mask1*ff2(i,km1,l)
                  trmin2(i) = mask(i,km1)*min(trmin2(i),xnu2) + mask1*xnu2
                  ff2(i,k,l)= tmem2(i) * trmin2(i)
               enddo
            enddo
         enddo
      endif IF_FF2
      
      if (ISREQSTEPL(TILIST) .or. ISREQOUTL(TMLIST)) then
         do i = 1, ni
            ztcsh(i) = 1. - ff2(i,IH2(i),ktop)
            ztcsm(i) = 1. - ff2(i,IB2(i),IH2(i))
            ztcsl(i) = 1. - ff2(i,nk,IB2(i))
            ! to calculate 2-d cloud fraction with calipso-GOCCP height criteria
            ztczh(i) = 1. - ff2(i,IH3(i),ktop)   
            ztczm(i) = 1. - ff2(i,IB3(i),IH3(i))
            ztczl(i) = 1. - ff2(i,nk,IB3(i))
         enddo
      endif
      
      if (ISREQSTEPL((/"TCC","NT "/)) .or. ISREQOUTL((/"TCCM","NF  "/))) &
           ztcc = 1. - ff2(:,nk,ktop)

      IF_TCCF: if (ISREQSTEP("TCCF")) then
         do l=1,3
            m = l+6
            do i=1,ni
               if (filterkc(i)>filterzk(i,l)) then
                  ztccf(i,m) = 1. - ff2(i,filterkc(i),filterzk(i,l))
               elseif (filterkc(i) == filterzk(i,l)) then
                  ztccf(i,m) = cloud(i,filterkc(i))
               else
                  ztccf(i,m) = 0.0 
               endif
            enddo
         enddo
         do i=1,ni !#TODO: invert loop order
            ztccf(i,1) = 1. - ff2(i,nk,ktop)
            do l=2,6
               ztccf(i,l) = 1. - ff2(i,filterkh(i,l),ktop)
            enddo
         enddo
      endif IF_TCCF

      IF_NT: if (ISREQSTEP("NT") .or. ISREQOUT("NF")) then
         if (stcond(1:3) == 'MP_') then
            do i=1,ni
               znt(i) = ztcc(i)*(1.-exp(-0.1*(ztopthw(i) + ztopthi(i))))
            enddo
         else
            call vintage_nt1( &
                 tt, ps, sig, zlwc, &
                 cloud, znt, trnch, ni, nk, nkm1)
         endif

         if (ISREQSTEPL((/"ECTP","ECTT"/))) then
            !     diagnostics: cloud top pressure (ctp) and temperature (ctt)
            !     using the cloud optical depth at window region (band 6) to
            !     calculate the emissivity
            do i = 1, ni
               zctp(i) = 110000.
               zctt(i) = 310.
               top(i)  = .true.
               transmissint(i,1) = 1. - cloud(i,1) * (1.-trans_exp(i,1) )
               if ((1. - transmissint(i,1)) > 0.99 .and. top(i)) then
                  zctp(i) = sig(i,1)*ps(i)
                  zctt(i) = tt(i,1)
                  top(i)  = .false.
               endif
            enddo
            do k = 2, nkm1
               do i = 1, ni
                  transmissint(i,k)=transmissint(i,k-1) * (1. - cloud(i,k) * &
                       (1.-trans_exp(i,k) ) )
                  if ((1. - transmissint(i,k)) > 0.99 .and. top(i)) then
                     zctp(i) = sig(i,k)*ps(i)
                     zctt(i) = tt(i,k)
                     top(i)  = .false.
                  endif
               enddo
            enddo
         endif

         if (ISREQSTEP("HCBR")) then
            ! cloud base defined as first level with non zero taucsv (should correspond to first level with CLDR)
            zhcbr= 20000. ! if no cloud base
            do i = 1, ni  !#TODO: invert loop order
               keeplooking = .true.
               do k = nkm1,1,-1
                  if (taucsv(i,k) > tiny(1.0) .and. keeplooking) then                  !inside a cloud layer
                     zhcbr(i) = gz(i,k)
                     keeplooking = .false.
                  endif
               enddo
            enddo
         endif
      endif IF_NT

      !----------------------------------------------------------------
      return
   end subroutine diagno_clouds2
end module diagno_clouds

