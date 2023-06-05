!-------------------------------------- LICENCE BEGIN ------------------------
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

module energy_budget
   use, intrinsic :: iso_fortran_env, only: REAL64
   use tdpack, only: CHLC, CHLF, CPD, CPI, CPV, CPW, GRAV, TCDK, TRPL
   use phy_status, only: PHY_OK, PHY_ERROR
   implicit none
   private
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   ! Public API
   public :: eb_init                                            !Initialize energy budget package
   public :: eb_en                                              !Vertically integrated thermal energy
   public :: eb_den                                             !Vertically integrated energy tendency
   public :: eb_pw                                              !Vertically integrated total water
   public :: eb_dpw                                             !Vertically integrated tendency of total water
   public :: eb_residual_en                                     !Normalized energy residual due to a physical process
   public :: eb_residual_pw                                     !Normalized total water residual due to a physical process
   public :: eb_conserve_en                                     !Impose energy conservation by adjusting temperature tendencies
   public :: eb_conserve_pw                                     !Impose total water conservation by adjusting humidity tendencies

   ! Local generic interfaces
   interface compute_lv
      module procedure compute_lv2d
      module procedure compute_lv3d
   end interface compute_lv
   interface ice_partition
      module procedure ice_partition2d
      module procedure ice_partition3d
   end interface ice_partition

   ! Local parameters
   integer, parameter :: LONG_CHAR=16                           !Length of standard character strings
   character(len=LONG_CHAR), parameter :: DEFAULT_INTTYPE='pchip'!Default type of integration strategy

   ! Local variables
   integer :: ksfc=1                                            !Near-surface level index
   real :: norm_en=1.                                           !Normalizing factor for energy unit conversion (to W/m2)
   real :: norm_pw=1./CHLC                                      !Normalizing factor for total water unit conversion (to W/m2)
   logical :: initialized=.false.                               !Initialization status of package
   logical :: lv_temp_dependent=.true.                          !Use temperature-dependent latent heat of vapourization and sublimation

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_init(F_lv_temp_dependent,F_norm_en,F_norm_pw,F_nearsfc_k) result(F_stat)
      ! Initialize the energy budget package
      implicit none

      ! Arguments
      logical, intent(in), optional :: F_lv_temp_dependent      !Use temperature-dependent latent heats [.true.]
      real, intent(in), optional :: F_norm_en                   !Normalizing factor for energy unit conversion [1.]
      real, intent(in), optional :: F_norm_pw                   !Normalizing factor for total water unit conversion [1/CHLC]
      integer, intent(in), optional :: F_nearsfc_k              !Index of the lowest model level [1]
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Set return values
      F_stat = PHY_ERROR

      ! Check initialization status
      if (initialized) then
         call msg(MSG_ERROR,'(energy_budget::eb_init) multiple initializations of EB package are not allowed')
         return
      endif
      initialized = .true.

      ! Set global variables for EB module
      if (present(F_lv_temp_dependent)) lv_temp_dependent = F_lv_temp_dependent
      if (present(F_norm_en)) norm_en = F_norm_en
      if (present(F_norm_pw)) norm_pw = F_norm_pw
      if (present(F_nearsfc_k)) ksfc = F_nearsfc_k

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function eb_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_en(F_eni,F_tt,F_hu,F_qc,F_sig,F_p0,F_nk,F_qi,F_inttype,F_en) result(F_stat)
      ! Compute the vertically integrated thermal energy
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_hu                  !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_qc                  !Liquid condensate (kg/kg)
      real, dimension(:,:), intent(in) :: F_sig                 !Sigma levels ()
      real, dimension(:), intent(in) :: F_p0                    !Surface pressure (Pa)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:,:), intent(in), optional :: F_qi        !Solid condensate (kg/kg) [diagnosed]
      character(len=*), intent(in), optional :: F_inttype       !Integral type ['pchip']
      real(REAL64), dimension(:), intent(out) :: F_eni          !Integrated energy (m3/s2)
      real, dimension(:,:), intent(out), optional :: F_en       !Profile of energy (m3/s2)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: istat
      real, dimension(size(F_tt,dim=1),F_nk) :: cpm,lv,ls,myqc,myqi,en,qw
      character(len=LONG_CHAR) :: myInttype

      ! Set return values
      F_stat = PHY_ERROR
      F_eni = 0.

      ! Handle optional arguments
      myInttype = DEFAULT_INTTYPE
      if (present(F_inttype)) myInttype = F_inttype

      ! Estimate ice fraction if no solid ice is provided
      if (present(F_qi)) then
         myqc = F_qc(:,1:F_nk)
         myqi = F_qi(:,1:F_nk)
      else
         istat = ice_partition(F_tt(:,1:F_nk),F_qc(:,1:F_nk),myqc,myqi)
      endif

      ! Compute total water
      qw = F_hu(:,1:F_nk) + myqc + myqi

      ! Compute thermodynamic parameters
      istat = thermo_param(F_tt(:,1:F_nk),qw,cpm,lv,ls)
      en = cpm*F_tt(:,1:F_nk) - lv*myqc - ls*myqi

      ! Integrate energy values in height (dp/g)
      if (integrate_profile(F_eni,en,F_sig(:,1:F_nk),F_p0,myInttype) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_en) error in integral calculation')
         return
      endif

      ! Set optional return values
      if (present(F_en)) F_en = en

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function eb_en

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_den(F_deni,F_dtt,F_dhu,F_tt,F_hu,F_qc,F_sig,F_p0,F_nk,F_dqc,F_dqi,F_qi,F_inttype,F_abs,F_den) result(F_stat)
      ! Compute the change of vertically integrated thermal energy
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_dtt                 !Dry air temperature tendency (K/s)
      real, dimension(:,:), intent(in) :: F_dhu                 !Specific humidity tendency (kg/kg/s)
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_hu                  !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_qc                  !Liquid condensate (kg/kg)
      real, dimension(:,:), intent(in) :: F_sig                 !Sigma levels ()
      real, dimension(:), intent(in) :: F_p0                    !Surface pressure (Pa)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:,:), intent(in), optional :: F_dqc       !Liquid condensate tendency (kg/kg/s) [0.]
      real, dimension(:,:), intent(in), optional :: F_dqi       !Solid condensate tendency (kg/kg/s) [diagnosed]
      real, dimension(:,:), intent(in), optional :: F_qi        !Solid condensate (kg/kg) [diagnosed]
      character(len=*), intent(in), optional :: F_inttype       !Integral type ['pchip']
      logical, intent(in), optional :: F_abs                    !Computed integrated absolute value of energy [.false.]
      real(REAL64), dimension(:), intent(out) :: F_deni   !Integrated energy tendency (m3/s3)
      real, dimension(:,:), intent(out), optional :: F_den      !Profile of energy tendency (m3/s3)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: istat
      real, dimension(size(F_tt,dim=1),F_nk) :: mydqc,mydqi,myqc,myqi,qw,dqw,cpm,lv,ls,dcpm,dlv,dls,den
      character(len=LONG_CHAR) :: myInttype
      logical :: myAbs

      ! Set return values
      F_stat = PHY_ERROR
      F_deni = 0.

      ! Handle optional arguments
      mydqc = 0.
      if (present(F_dqc)) mydqc = F_dqc(:,1:F_nk)
      myInttype = DEFAULT_INTTYPE
      if (present(F_inttype)) myInttype = F_inttype
      myAbs = .false.
      if (present(F_abs)) myAbs = F_abs

      ! Estimate tendencies for solid/liquid components if no ice tendency is provided
      if (present(F_dqi)) then
         mydqi = F_dqi(:,1:F_nk)
      else
         istat = ice_partition(F_tt(:,1:F_nk),F_dqc(:,1:F_nk),mydqc,mydqi)
      endif
      if (present(F_qi)) then
         myqc = F_qc(:,1:F_nk)
         myqi = F_qi(:,1:F_nk)
      else
         istat = ice_partition(F_tt(:,1:F_nk),F_qc(:,1:F_nk),myqc,myqi)
      endif

      ! Compute total water quantities
      qw = F_hu(:,1:F_nk) + myqc + myqi
      dqw = F_dhu(:,1:F_nk) + mydqc + mydqi

      ! Compute thermodynamic parameters
      istat = thermo_param(F_tt(:,1:F_nk),qw,cpm,lv,ls)
      istat = dthermo_param(dcpm,dlv,dls)

      ! Compute profiles of energy tendency
      den = cpm*F_dtt(:,1:F_nk) + F_tt(:,1:F_nk)*dcpm*dqw - (myqc*dlv*F_dtt(:,1:F_nk) + lv*mydqc + myqi*dls*F_dtt(:,1:F_nk) + ls*mydqi)

      ! Convert to absolute value on request
      if (myAbs) den = abs(den)

      ! Integrate energy tendency values in height (dp/g)
      if (integrate_profile(F_deni,den,F_sig(:,1:F_nk),F_p0,myInttype) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_den) error in integral calculation')
         return
      endif

      ! Set optional return values
      if (present(F_den)) F_den = den

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function eb_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_conserve_en(F_dttc,F_dtt,F_dhu,F_tt,F_hu,F_qc,F_sig,F_p0,F_nk,F_dqc,F_dqi,F_qi, &
        F_rain,F_snow,F_shf,F_wvf,F_rad,F_mask,F_inttype) result(F_stat)
      ! Apply corrective offset to tendencies for energy conservation
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_dtt                 !Dry air temperature tendency (K/s)
      real, dimension(:,:), intent(in) :: F_dhu                 !Specific humidity tendency (kg/kg/s)
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_hu                  !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_qc                  !Liquid condensate (kg/kg)
      real, dimension(:,:), intent(in) :: F_sig                 !Sigma levels ()
      real, dimension(:), intent(in) :: F_p0                    !Surface pressure (Pa)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:,:), intent(in), optional :: F_dqc       !Liquid condensate tendency (kg/kg/s)
      real, dimension(:,:), intent(in), optional :: F_dqi       !Solid condensate tendency (kg/kg/s)
      real, dimension(:,:), intent(in), optional :: F_qi        !Solid condensate (kg/kg)
      real, dimension(:), intent(in), optional :: F_rain        !Surface liquid precipitation flux (kg/m2/s) [0.]
      real, dimension(:), intent(in), optional :: F_snow        !Surface solid precipitation flux (kg/m2/s) [0.]
      real, dimension(:), intent(in), optional :: F_shf         !Surface turbulent sensible heat flux (W/m2) [0.]
      real, dimension(:), intent(in), optional :: F_wvf         !Surface water vapour flux (kg/m2/s) [0.]
      real, dimension(:), intent(in), optional :: F_rad         !Net radiation flux TOA-surface (W/m2) [0.]
      logical, intent(in), optional :: F_mask                   !Apply offset only where tendencies are non-zero [true]
      character(len=*), intent(in), optional :: F_inttype       !Integral type ['pchip']
      real, dimension(:,:), intent(out) :: F_dttc               !Conservative temperature humidity tendency (K/s)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: k,istat
      real, dimension(size(F_tt,dim=1)) :: myrain,mysnow,myshf,mywvf,myrad,src,mult
      real, dimension(size(F_tt,dim=1),F_nk) :: myqc,myqi,mydqc,mydqi,qw,cpm,ls,lv,dcpm,dls,dlv,dqw,den,denc
      real(REAL64), dimension(size(F_tt,dim=1)) :: deni,dabseni
      character(len=LONG_CHAR) :: myInttype
      logical :: mymask

      ! Set return values
      F_stat = PHY_ERROR
      F_dttc(:,1:F_nk) = F_dtt(:,1:F_nk)

      ! Handle optional arguments
      mydqc = 0.
      if (present(F_dqc)) mydqc = F_dqc(:,1:F_nk)
      myrain = 0.
      if (present(F_rain)) myrain = F_rain
      mysnow = 0.
      if (present(F_snow)) mysnow = F_snow
      myshf = 0.
      if (present(F_shf)) myshf = F_shf
      mywvf = 0.
      if (present(F_wvf)) mywvf = F_wvf
      myrad = 0.
      if (present(F_rad)) myrad = F_rad
      mymask = .true.
      if (present(F_mask)) mymask = F_mask
      myInttype = DEFAULT_INTTYPE
      if (present(F_inttype)) myInttype = F_inttype

      ! Estimate tendencies for solid/liquid components if no ice tendency is provided
      if (present(F_dqi)) then
         mydqi = F_dqi(:,1:F_nk)
      else
         istat = ice_partition(F_tt(:,1:F_nk),mydqc,mydqc,mydqi)
      endif
      if (present(F_qi)) then
         myqc = F_qc(:,1:F_nk)
         myqi = F_qi(:,1:F_nk)
      else
         istat = ice_partition(F_tt(:,1:F_nk),F_qc(:,1:F_nk),myqc,myqi)
      endif

      ! Compute total water
      qw = F_hu(:,1:F_nk) + myqc + myqi
      dqw = F_dhu(:,1:F_nk) + mydqc + mydqi

      ! Compute thermodynamic properties
      istat = thermo_param(F_tt(:,1:F_nk),qw,cpm,lv,ls)
      istat = dthermo_param(dcpm,dlv,dls)

      ! Compute source term from provided boundary values
      istat = compute_source_en(src,F_tt(:,1:F_nk),F_hu(:,1:F_nk),myqc,myqi,myrain,mysnow,myshf,mywvf,myrad)

      ! Compute tendency-based energy tendency
      if (eb_den(deni,F_dtt(:,1:F_nk),F_dhu(:,1:F_nk),F_tt(:,1:F_nk),F_hu(:,1:F_nk),myqc,F_sig(:,1:F_nk), &
           F_p0,F_nk,F_dqc=mydqc(:,1:F_nk),F_dqi=mydqi(:,1:F_nk),F_qi=myqi,F_inttype=myInttype,F_den=den) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_conserve_en) error in tendency integral')
         return
      endif

      ! Compute absolute value of energy tendency
      if (eb_den(dabseni,F_dtt(:,1:F_nk),F_dhu(:,1:F_nk),F_tt(:,1:F_nk),F_hu(:,1:F_nk),myqc,F_sig(:,1:F_nk), &
           F_p0,F_nk,F_dqc=mydqc(:,1:F_nk),F_dqi=mydqi(:,1:F_nk),F_qi=myqi,F_inttype=myInttype,F_abs=.true.) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_conserve_en) error in tendency integral')
         return
      endif

      ! Compute correction factor for conservation
      ! Compute correction factor for conservation
      where (dabseni > epsilon(dabseni)*norm_en)
         mult = (deni - src) / dabseni
      elsewhere
         mult = 0.
      endwhere

      ! Apply multiplicative correction to dry air temperature
      do k=1,F_nk
         denc(:,k) = den(:,k) * (1.-sign(1.,den(:,k))*mult(:))
         F_dttc(:,k) = F_dtt(:,k) + (denc(:,k)-den(:,k)) / (cpm(:,k) - dlv(:,k)*myqc(:,k) - dls(:,k)*myqi(:,k))
      enddo

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function eb_conserve_en

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_residual_en(F_resid,F_fld0,F_fld1,F_tt,F_hu,F_qc,F_dt,F_nk,F_qi,F_rain,F_snow,F_shf,F_wvf,F_rad) result(F_stat)
      ! Compute the budget residual for thermal energy
      implicit none

      ! Arguments
      real(REAL64), dimension(:), intent(in) :: F_fld0          !Energy before process
      real(REAL64), dimension(:), intent(in) :: F_fld1          !Energy after process
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_hu                  !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_qc                  !Liquid condensate (kg/kg)
      real, intent(in) :: F_dt                                  !Time step (s)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:,:), intent(in), optional :: F_qi        !Solid condensate (kg/kg) [automatic]
      real, dimension(:), intent(in), optional :: F_rain        !Surface liquid precipitation flux (kg/m2/s) [0.]
      real, dimension(:), intent(in), optional :: F_snow        !Surface solid precipitation flux (kg/m2/s) [0.]
      real, dimension(:), intent(in), optional :: F_shf         !Surface turbulent sensible heat flux (W/m2) [0.]
      real, dimension(:), intent(in), optional :: F_wvf         !Surface turbulent water vapour flux (kg/m2/s) [0.]
      real, dimension(:), intent(in), optional :: F_rad         !Net radiation flux TOA-surface (W/m2) [0.]
      real(REAL64), dimension(:), intent(out) :: F_resid        !Energy residual(Wm-2 if norm_en=1.) (non-conservation)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: istat
      real, dimension(size(F_fld0,dim=1)) :: myrain,mysnow,myshf,mywvf,myrad,src
      real, dimension(size(F_tt,dim=1),F_nk) :: myqc,myqi

      ! Set return values
      F_stat = PHY_ERROR
      F_resid = 0.

      ! Handle optional arguments
      myrain = 0.
      if (present(F_rain)) myrain = F_rain
      mysnow = 0.
      if (present(F_snow)) mysnow = F_snow
      myshf = 0.
      if (present(F_shf)) myshf = F_shf
      mywvf = 0.
      if (present(F_wvf)) mywvf = F_wvf
      myrad = 0.
      if (present(F_rad)) myrad = F_rad

      ! Estimate solid/liquid components if no ice tendency is provided
      if (present(F_qi)) then
         myqc = F_qc(:,1:F_nk)
         myqi = F_qi(:,1:F_nk)
      else
         istat = ice_partition(F_tt(:,1:F_nk),F_qc(:,1:F_nk),myqc,myqi)
      endif

      ! Compute source term from provided boundary values
      istat = compute_source_en(src,F_tt(:,1:F_nk),F_hu(:,1:F_nk),myqc,myqi,myrain,mysnow,myshf,mywvf,myrad)

      ! Use generic residual function with specified normalization factor
      F_stat = compute_residual(F_resid,F_fld0,F_fld1,F_dt,norm_en,src)

      ! End of subprogram
      return
   end function eb_residual_en

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function compute_source_en(F_src,F_tt,F_qv,F_qc,F_qi,F_rain,F_snow,F_shf,F_wvf,F_rad) result(F_stat)
      ! Compute the boundary source of energy
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_qv                  !Water vapour specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_qc                  !Liquid condensate specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_qi                  !Solid condensate specific humidity (kg/kg)
      real, dimension(:), intent(in) :: F_rain                  !Surface liquid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in) :: F_snow                  !Surface solid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in) :: F_shf                   !Surface turbulent sensible heat flux (W/m2)
      real, dimension(:), intent(in) :: F_wvf                   !Surface turbulent water vapour flux (kg/m2/s)
      real, dimension(:), intent(in) :: F_rad                   !Net radiation flux TOA-surface (W/m2)
      real, dimension(:), intent(out) :: F_src                  !Energy boundary source (W/m2)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: istat
      real, dimension(size(F_tt,dim=1)) :: lv,ls,qw

      ! Set return values
      F_stat = PHY_ERROR
      F_src = 0.

      ! Compute latent heats
      istat = compute_lv(F_tt(:,ksfc),F_lv=lv,F_ls=ls)

      ! Compute total water at the near-surface level
      qw = F_qv(:,ksfc) + F_qc(:,ksfc) + F_qi(:,ksfc)

      ! Precipitation energy source
      F_src = F_src + (lv+(CPD-CPV)*F_tt(:,ksfc))*F_rain + (ls+(CPD-CPV)*F_tt(:,ksfc))*F_snow

      ! Turbulent flux energy source
      F_src = F_src + (1.+(CPV/CPD-1.)*qw)*F_shf + ((CPV-CPD)*F_tt(:,ksfc))*F_wvf

      ! Radiative energy source
      F_src = F_src + F_rad

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function compute_source_en

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_pw(F_pwi,F_hu,F_qc,F_sig,F_p0,F_nk,F_qi,F_inttype,F_pw) result(F_stat)
      ! Compute the vertically integrated total water
      use integrals, only: int_profile,INT_OK
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_hu                  !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_qc                  !Liquid condensate (kg/kg)
      real, dimension(:,:), intent(in) :: F_sig                 !Sigma levels ()
      real, dimension(:), intent(in) :: F_p0                    !Surface pressure (Pa)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:,:), intent(in), optional :: F_qi        !Solid condensate (kg/kg)
      character(len=*), intent(in), optional :: F_inttype       !Integral type ['pchip']
      real(REAL64), dimension(:), intent(out) :: F_pwi          !Integrated total water (Kg/m2)
      real, dimension(:,:), intent(out), optional :: F_pw       !Profile of total water (kg/kg)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      real, dimension(size(F_hu,dim=1),F_nk) :: qw
      character(len=LONG_CHAR) :: myInttype

      ! Set return values
      F_stat = PHY_ERROR
      F_pwi = 0.

      ! Handle optional arguments
      myInttype = DEFAULT_INTTYPE
      if (present(F_inttype)) myInttype = F_inttype

      ! Compute total water
      qw = F_hu(:,1:F_nk) + F_qc(:,1:F_nk)
      if (present(F_qi)) qw = qw + F_qi(:,1:F_nk)

      ! Integrate total water values in height (dp/g)
      if (integrate_profile(F_pwi,qw,F_sig(:,1:F_nk),F_p0,myInttype) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_pw) error in integral calculation')
         return
      endif

      ! Set optional return values
      if (present(F_pw)) F_pw = qw

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function eb_pw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_dpw(F_dpwi,F_dhu,F_hu,F_sig,F_p0,F_nk,F_dqc,F_dqi,F_inttype,F_abs,F_dpw) result(F_stat)
      ! Compute the vertically integrated total water tendency
      use integrals, only: int_profile,INT_OK
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_dhu                 !Specific humidity tendency (kg/kg/s)
      real, dimension(:,:), intent(in) :: F_hu                  !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_sig                 !Sigma levels ()
      real, dimension(:), intent(in) :: F_p0                    !Surface pressure (Pa)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:,:), intent(in), optional :: F_dqc       !Liquid condensate tendency (kg/kg/s) [0.]
      real, dimension(:,:), intent(in), optional :: F_dqi       !Solid condensate tendency (kg/kg/s) [0.]
      character(len=*), intent(in), optional :: F_inttype       !Integral type ['pchip']
      logical, intent(in), optional :: F_abs                    !Computed integrated absolute value of total water [.false.]
      real(REAL64), dimension(:), intent(out) :: F_dpwi         !Integrated total water (kg/m2)
      real, dimension(:,:), intent(out), optional :: F_dpw      !Profile of total water (kg/kg)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      real, dimension(size(F_hu,dim=1),F_nk) :: dqw,mydqc,mydqi
      character(len=LONG_CHAR) :: myInttype
      logical :: myAbs

      ! Set return values
      F_stat = PHY_ERROR
      F_dpwi = 0.

      ! Handle optional arguments
      mydqc = 0.
      if (present(F_dqc)) mydqc = F_dqc(:,1:F_nk)
      mydqi = 0.
      if (present(F_dqi)) mydqi = F_dqi(:,1:F_nk)
      myInttype = DEFAULT_INTTYPE
      if (present(F_inttype)) myInttype = F_inttype
      myAbs = .false.
      if (present(F_abs)) myAbs = F_abs

      ! Compute total water tendency
      dqw = F_dhu(:,1:F_nk) + mydqc + mydqi

      ! Convert to absolute value on request
      if (myAbs) dqw = abs(dqw)

      ! Integrate total water values in height (dp/g)
      if (integrate_profile(F_dpwi,dqw,F_sig(:,1:F_nk),F_p0,myInttype) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_dpw) error in integral calculation')
         return
      endif

      ! Set optional return values
      if (present(F_dpw)) F_dpw = dqw

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function eb_dpw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_conserve_pw(F_dhuc,F_dhu,F_tt,F_hu,F_sig,F_p0,F_nk,F_dqc,F_dqi,F_rain,F_snow,F_wvf,F_inttype) result(F_stat)
      ! Apply corrective offset to tendencies for total water conservation
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_dhu                 !Specific humidity tendency (kg/kg/s)
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_hu                  !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: F_sig                 !Sigma levels ()
      real, dimension(:), intent(in) :: F_p0                    !Surface pressure (Pa)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:,:), intent(in), optional :: F_dqc       !Liquid condensate tendency (kg/kg/s) [0.]
      real, dimension(:,:), intent(in), optional :: F_dqi       !Solid condensate tendency (kg/kg/s) [0.]
      real, dimension(:), intent(in), optional :: F_rain        !Surface liquid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in), optional :: F_snow        !Surface solid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in), optional :: F_wvf         !Surface turbulent water vapour flux (kg/m2/s)
      character(len=*), intent(in), optional :: F_inttype       !Integral type ['pchip']
      real, dimension(:,:), intent(out) :: F_dhuc               !Conservative tendency (kg/kg/s)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: k,istat
      real, dimension(size(F_tt,dim=1)) :: myrain,mysnow,mywvf,src,mult
      real, dimension(size(F_tt,dim=1),F_nk) :: mydqc,mydqi,dqw,dqwc
      real(REAL64), dimension(size(F_tt,dim=1)) :: dpw,dabspw
      character(len=LONG_CHAR) :: myInttype

      ! Set return values
      F_stat = PHY_ERROR
      F_dhuc(:,1:F_nk) = F_dhu(:,1:F_nk)

      ! Handle optional arguments
      mydqc = 0.
      if (present(F_dqc)) mydqc = F_dqc(:,1:F_nk)
      mydqi = 0.
      if (present(F_dqi)) mydqi = F_dqi(:,1:F_nk)
      myrain = 0.
      if (present(F_rain)) myrain = F_rain
      mysnow = 0.
      if (present(F_snow)) mysnow = F_snow
      mywvf = 0.
      if (present(F_wvf)) mywvf = F_wvf
      myInttype = DEFAULT_INTTYPE
      if (present(F_inttype)) myInttype = F_inttype

      ! Compute source term from provided boundary values
      istat = compute_source_pw(src,F_tt(:,1:F_nk),myrain,mysnow,mywvf)

      ! Compute tendency-based total water tendency
      if (eb_dpw(dpw,F_dhu(:,1:F_nk),F_hu(:,1:F_nk),F_sig(:,1:F_nk),F_p0,F_nk, &
           F_dqc=mydqc,F_dqi=mydqi,F_inttype=myInttype,F_dpw=dqw) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_conserve_pw) error in tendency integral')
         return
      endif

      ! Compute absolute value of total water tendency
      if (eb_dpw(dabspw,F_dhu(:,1:F_nk),F_hu(:,1:F_nk),F_sig(:,1:F_nk),F_p0,F_nk, &
           F_dqc=mydqc,F_dqi=mydqi,F_inttype=myInttype,F_abs=.true.) /= PHY_OK) then
         call msg(MSG_ERROR,'(energy_budget::eb_conserve_pw) error in absolute tendency integral')
         return
      endif

      ! Compute correction factor for conservation
      where (dabspw > epsilon(dabspw)*norm_pw)
         mult = (dpw - src) / dabspw
      elsewhere
         mult = 0.
      endwhere

      ! Apply correction to specific humidity tendency
      do k=1,F_nk
         dqwc(:,k) = dqw(:,k) * (1.-sign(1.,dqw(:,k))*mult(:))
         F_dhuc(:,k) = F_dhu(:,k) + (dqwc(:,k)-dqw(:,k))
      enddo

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function eb_conserve_pw

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function eb_residual_pw(F_resid,F_fld0,F_fld1,F_tt,F_dt,F_nk,F_rain,F_snow,F_wvf) result(F_stat)
      ! Compute the budget residual for total water
      implicit none

      ! Arguments
      real(REAL64), dimension(:), intent(in) :: F_fld0          !Total water before process
      real(REAL64), dimension(:), intent(in) :: F_fld1          !Total water after process
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, intent(in) :: F_dt                                  !Time step (s)
      integer, intent(in) :: F_nk                               !Number of prognostic levels
      real, dimension(:), intent(in), optional :: F_rain        !Surface liquid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in), optional :: F_snow        !Surface solid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in), optional :: F_wvf         !Surface turbulent water vapour flux (kg/m2/s)
      real(REAL64), dimension(:), intent(out) :: F_resid        !Residual of field (kg/m2/s if norm_pw=1. or W/m2 if norm_pw=1./CHLC)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: istat
      real, dimension(size(F_fld0,dim=1)) :: myrain,mysnow,mywvf,src

      ! Set return values
      F_stat = PHY_ERROR
      F_resid = 0.

      ! Handle optional arguments
      myrain = 0.
      if (present(F_rain)) myrain = F_rain
      mysnow = 0.
      if (present(F_snow)) mysnow = F_snow
      mywvf = 0.
      if (present(F_wvf)) mywvf = F_wvf

      ! Compute source term from provided boundary values
      istat = compute_source_pw(src,F_tt(:,1:F_nk),myrain,mysnow,mywvf)

      ! Use generic residual function with specified normalization factor
      F_stat = compute_residual(F_resid,F_fld0,F_fld1,F_dt,norm_pw,src)

      ! End of subprogram
      return
   end function eb_residual_pw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function compute_source_pw(F_src,F_tt,F_rain,F_snow,F_wvf) result(F_stat)
      ! Compute the boundary source of total water
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:), intent(in) :: F_rain                  !Surface liquid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in) :: F_snow                  !Surface solid precipitation flux (kg/m2/s)
      real, dimension(:), intent(in) :: F_wvf                   !Surface turbulent water vapour flux (kg/m2/s)
      real, dimension(:), intent(out) :: F_src                  !Boundary source of total water (kg/m2/s)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: istat
      real, dimension(size(F_tt,dim=1)) :: lv

      ! Set return values
      F_stat = PHY_ERROR
      F_src = 0.

      ! Compute latent heats
      istat = compute_lv(F_tt(:,ksfc),F_lv=lv)

      ! Precipitation water source
      F_src = F_src - F_rain - F_snow

      ! Turbulent flux water source
      F_src = F_src + F_wvf

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function compute_source_pw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function compute_lv3d(F_tt,F_lv,F_ls) result(F_stat)
      ! Compute the latent heat of vapourization and sublimation
      implicit none

      ! Arguments
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(out), optional :: F_lv       !Latent heat of vapourization (J/kg)
      real, dimension(:,:), intent(out), optional :: F_ls       !Latent heat of sublimation (J/kg)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local varibles
      real, dimension(size(F_tt,dim=1),size(F_tt,dim=2)) :: mylv,myls
      integer :: k

      ! Set return values
      F_stat = PHY_ERROR

      ! Call 2D form of function for calculation
      F_stat = PHY_OK
      do k=1,size(F_tt,dim=2)
         F_stat = min(F_stat,compute_lv(F_tt(:,k),F_lv=mylv(:,k),F_ls=myls(:,k)))
      enddo

      ! Provide requested outputs
      if (present(F_lv)) F_lv = mylv
      if (present(F_ls)) F_ls = myls

      ! End of subprogram
      return
   end function compute_lv3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function compute_lv2d(F_tt,F_lv,F_ls,F_lv_temp_dependent) result(F_stat)
      ! Compute the latent heats of vapourization and sublimation
      implicit none

      ! Arguments
      real, dimension(:), intent(in) :: F_tt                    !Dry air temperature (K)
      real, dimension(:), intent(out), optional :: F_lv         !Latent heat of vapourization (J/kg)
      logical, intent(in), optional :: F_lv_temp_dependent      !Use temperature-dependent latent heats [as set in eb_init()]
      real, dimension(:), intent(out), optional :: F_ls         !Latent heat of sublimation (J/kg)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      logical :: myLv_temp_dependent

      ! Set return values
      F_stat = PHY_ERROR

      ! Handle command line arguments
      myLv_temp_dependent = lv_temp_dependent
      if (present(F_lv_temp_dependent)) myLv_temp_dependent = F_lv_temp_dependent

      ! Compute temperature-dependent Lv that is internally consistent for the package
      if (myLv_temp_dependent) then
         if (present(F_lv)) F_lv = CHLC + (CPV - CPW) * (F_tt - TRPL)
         if (present(F_ls)) F_ls = (CHLC + CHLF) + (CPV - CPI) * (F_tt - TRPL)
      else
         if (present(F_lv)) F_lv = CHLC
         if (present(F_ls)) F_ls = CHLC+CHLF
      endif

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function compute_lv2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function mixing_from_specific(F_mix,F_spec) result(F_stat)
      ! Compute mixing ratio from specific value
      implicit none
      real, dimension(:,:), intent(in) :: F_spec                !Specfic value
      real, dimension(:,:), intent(out) :: F_mix                !Mixing ratio value
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)
      F_stat = PHY_OK
      F_mix = F_spec / (1.-F_spec)
      return
   end function mixing_from_specific

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function dmixing_from_dspecific(F_dmix,F_dspec,F_spec) result(F_stat)
      ! Compute mixing ratio from specific value
      implicit none
      real, dimension(:,:), intent(in) :: F_dspec               !Tendency of specfic value
      real, dimension(:,:), intent(in) :: F_spec                !Specfic value
      real, dimension(:,:), intent(out) :: F_dmix               !Tendency of mixing ratio value
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)
      real, dimension(size(F_spec,dim=1),size(F_spec,dim=2)) :: mix
      F_stat = mixing_from_specific(mix,F_spec)
      F_dmix = (1.+mix)**2 * F_dspec
      return
   end function dmixing_from_dspecific

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function thermo_param(F_tt,F_qw,F_cpm,F_lv,F_ls) result(F_stat)
      ! Compute thermodynamic parameters as in Bechtold (2001)
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_qw                  !Total water specific humidity (kg/kg)
      real, dimension(:,:), intent(out) :: F_cpm                !Specific heat of moist air
      real, dimension(:,:), intent(out) :: F_lv                 !Latent heat of vaporization
      real, dimension(:,:), intent(out) :: F_ls                 !Latent heat of sublimation
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)
      integer :: istat
      F_stat = PHY_OK
      F_cpm = CPD + (CPV-CPD) * F_qw
      istat = compute_lv(F_tt,F_lv=F_lv,F_ls=F_ls)
      return
   end function thermo_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function dthermo_param(F_dcpm,F_dlv,F_dls) result(F_stat)
      ! Compute thermodynamic parameters as in Bechtold (2001)
      real, dimension(:,:), intent(out) :: F_dcpm               !Specific heat of moist air tendency
      real, dimension(:,:), intent(out) :: F_dlv                !Latent heat of vaporization tendency
      real, dimension(:,:), intent(out) :: F_dls                !Latent heat of sublimation tendency
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)
      F_stat = PHY_OK
      F_dcpm = CPV-CPD
      F_dlv = CPV - CPW
      F_dls = CPV - CPI
      return
   end function dthermo_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ice_partition2d(F_tts,F_qcs,F_qliqs,F_qices) result(F_stat)
      ! Compute fractional ice content (Rockel et al. Beitr. Atmos. Phy. 1991)
      implicit none
      real, dimension(:), intent(in) :: F_tts                   !Dry air temperature (K)
      real, dimension(:), intent(in) :: F_qcs                   !Total condensate (kg/kg)
      real, dimension(:), intent(out) :: F_qliqs                !Liquid condensate (kg/kg)
      real, dimension(:), intent(out) :: F_qices                !Solid condensate (kg/kg)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      real, dimension(size(F_tts)) :: wfrac

      ! Set return value
      F_stat = PHY_OK

      ! Partition condensate field based on environmental temperature
      where (F_tts >= TCDK)
         wfrac = 1.
      elsewhere
         wfrac = 0.0059 + 0.9941*exp(-0.003102*(F_tts-TCDK)**2)
      endwhere
      where (wfrac < 0.01) wfrac = 0.
      F_qices = (1.-wfrac) * F_qcs
      F_qliqs = wfrac * F_qcs

      ! End of subprogram
      return
   end function ice_partition2d

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ice_partition3d(F_tt,F_qc,F_qliq,F_qice) result(F_stat)
      ! Compute fractional ice content (Rockel et al. Beitr. Atmos. Phy. 1991)
      implicit none
      real, dimension(:,:), intent(in) :: F_tt                  !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: F_qc                  !Total condensate (kg/kg)
      real, dimension(:,:), intent(out) :: F_qliq               !Liquid condensate (kg/kg)
      real, dimension(:,:), intent(out) :: F_qice               !Solid condensate (kg/kg)
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: k,istat

      ! Set return value
      F_stat = PHY_OK

      ! Loop to compute condensate fields for all levels
      do k=1,size(F_tt,dim=2)
         istat = ice_partition2d(F_tt(:,k),F_qc(:,k),F_qliq(:,k),F_qice(:,k))
      enddo

      ! End of subprogram
      return
   end function ice_partition3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function integrate_profile(F_fldi,F_fld,F_sig,F_p0,F_inttype) result(F_stat)
      use integrals, only: int_profile,INT_OK
      implicit none
      real, dimension(:,:), intent(in) :: F_fld                 !Field to integrate
      real, dimension(:,:), intent(in) :: F_sig                 !Sigma levels ()
      real, dimension(:), intent(in) :: F_p0                    !Surface pressure (Pa)
      character(len=*), intent(in), optional :: F_inttype       !Integral type
      real(REAL64), dimension(:), intent(out) :: F_fldi   !Integrated value
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Local variables
      integer :: k,nk
      real, dimension(size(F_fld,dim=1)) :: pbot,ptop
      real, dimension(size(F_fld,dim=1),size(F_fld,dim=2)) :: pinv,fldinv

      ! Set return status
      F_stat = PHY_ERROR

      ! Invert pressure slab and integrate
      nk = size(F_fld,dim=2)
      do k=1,nk
         pinv(:,nk-(k-1)) = F_sig(:,k)*F_p0
         fldinv(:,nk-(k-1)) = F_fld(:,k)
      enddo
      ptop = pinv(:,nk)
      pbot = F_p0
      if (int_profile(F_fldi,fldinv/GRAV,pinv,ptop,pbot,type=F_inttype) /= INT_OK) then
         call msg(MSG_ERROR,'(energy_budget::integrate_profile) error in integral calculation')
         return
      endif

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function integrate_profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function compute_residual(F_resid,F_fld0,F_fld1,F_dt,F_norm,F_src) result(F_stat)
      ! Compute the budget residual
      implicit none

      ! Arguments
      real(REAL64), dimension(:), intent(in) :: F_fld0    !Field before process
      real(REAL64), dimension(:), intent(in) :: F_fld1    !Field after process
      real, intent(in) :: F_dt                                  !Time step (s)
      real, intent(in) :: F_norm                                !Normalizing factor for conserved field
      real, dimension(:), intent(in) :: F_src                   !Sources for field (boundary conditions /s)
      real(REAL64), dimension(:), intent(out) :: F_resid  !Residual of field
      integer :: F_stat                                         !Return status (PHY_OK or PHY_ERROR)

      ! Set return values
      F_stat = PHY_ERROR
      F_resid = 0.

      ! Compute residual of conserved quantity
      F_resid = F_fld1 - (F_fld0 + dble(F_dt*F_src))

      ! Normalize residual by expected "noise" level and convert to tendency /s
      F_resid = F_resid / F_dt / F_norm

      ! End of subprogram
      F_stat = PHY_OK
      return
   end function compute_residual

end module energy_budget
