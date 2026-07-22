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
module sfc_nml_mod
  implicit none
  public
contains
  
function sfc_nml2(F_namelist) result(F_istat)
   use clib_itf_mod, only: clib_toupper, clib_isreadok
   use wb_itf_mod
   use str_mod, only: str_concat, str_toreal
   use sfc_options
   use sfcbus_mod
#ifdef HAVE_NEMO
   use cpl_itf, only: cpl_nml
#endif
   use sfclayer, only: sl_put, SL_OK
   implicit none
!!!#include <arch_specific.hf>
   !@objective Set defaults values and read surface namelist
   !@Arguments
   ! F_namelist    File name containing the namelists to read
   character(len=*), intent(in) :: F_namelist
   !@return
   integer :: F_istat
   !*@/

#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   include "tebcst.cdk"

   integer, parameter :: SFC_NML_ERR = RMN_ERR
   integer, parameter :: SFC_NML_OK  = RMN_OK + 1
   integer, parameter :: CPL_NML_OK  = 1

   integer :: err, unout
   !-------------------------------------------------------------------
   F_istat = SFC_NML_ERR

   err = sfc_options_init()
   if (.not.RMN_IS_OK(err)) return
   err = sfc_nml_read(F_namelist)
   if (.not.RMN_IS_OK(err)) return
   err = sfc_nml_check()
   if (.not.RMN_IS_OK(err)) return
   call sfc_nml_print()
   err = sfc_nml_post_init()
   if (.not.RMN_IS_OK(err)) return

   unout = msg_getUnit(MSG_INFO)
#ifdef HAVE_NEMO
   err   = cpl_nml(F_namelist, unout)
#else
   err   = 0
#endif
   if (.not.RMN_IS_OK(err)) then
      call msg(MSG_ERROR,'(sfc_nml) Probleme in cpl_nml')
      return
   endif
   cplocn = (err ==  CPL_NML_OK)
#ifdef HAVE_NEMO
   if (cplocn) err = cpl_nml('print', unout)
#endif

   F_istat = SFC_NML_OK
   !-------------------------------------------------------------------
   return


contains
 

   function sfc_nml_read(m_namelist) result(m_istat)
      implicit none
      character(len=*), intent(in) :: m_namelist
      integer :: m_istat
      character(len=1024) :: msg_S, namelist_S, name_S
      integer :: istat, unf
      !----------------------------------------------------------------
      m_istat = RMN_ERR
      name_S = '(sfc_nml)'
      namelist_S = '&surface_cfgs'

      unf = 0
      istat = clib_isreadok(m_namelist)
      if (RMN_IS_OK(istat)) istat = fnom(unf, m_namelist, 'SEQ+OLD', 0)
      if (.not.RMN_IS_OK(istat)) then
         write(msg_S,'(a,a,a,a)') trim(name_S), &
              ' Using default config, ', &
              'Namelist file Not found/readable: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
         return
      endif

      read(unf, nml=surface_cfgs, iostat=istat)
      if (istat == 0) then     !# Read ok
         m_istat = RMN_OK + 1
      else if (istat < 0) then !# EOF, nml not found
         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
              ' Using default config, Namelist ',&
              trim(namelist_S), ' not available in file: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
      else !# Error
         write(msg_S,'(a,a,a,a,a)') trim(name_S), &
              ' Namelist', trim(namelist_S), &
              ' invalid in file: ', trim(m_namelist)
         call msg(MSG_ERROR, msg_S)
      endif
      istat = fclos(unf)
      !----------------------------------------------------------------
      return
   end function sfc_nml_read


   subroutine sfc_nml_print()
      implicit none
      integer, external :: msg_getUnit
      integer :: unout
      !----------------------------------------------------------------
      unout = msg_getUnit(MSG_INFO)
      if (unout > 0) then
         write(unout, nml=surface_cfgs)
      endif
      !----------------------------------------------------------------
      return
   end subroutine sfc_nml_print


   function sfc_nml_check() result(m_istat)
      use sfc_options, only : nl_svs_default, dp_svs_default, &
                       nb_teb_snowpar, teb_snroof_default, teb_snroad_default
      use svs_configs, only : nl_svs,dl_svs,ntypel,ntypeh,vl_type,vh_type
      use mode_crodebug
      implicit none
      integer :: m_istat

      integer :: istat,k,nk
      character(len=1024) :: msg_S
      !----------------------------------------------------------------
      m_istat = RMN_ERR

      istat = clib_toupper(diusst)
      istat = clib_toupper(ice_emiss)
      istat = clib_toupper(schmsol)
      istat = clib_toupper(schmurb)
      istat = clib_toupper(sl_diag_type)
      istat = clib_toupper(sl_lmin_type)
      istat = clib_toupper(sl_func_stab)
      istat = clib_toupper(sl_func_unstab)
      istat = clib_toupper(snow_emiss)
      istat = clib_toupper(isba_soil_emiss)
      istat = clib_toupper(isba_snowfrac_bare)
      istat = clib_toupper(soiltext)
      istat = clib_toupper(soil_cond)
      istat = clib_toupper(soil_ksat_ice)
      istat = clib_toupper(svs_hrsurf_method)
      istat = clib_toupper(isba_hrsurf_method)
      istat = clib_toupper(svs_snow_rain)
      istat = clib_toupper(vf_type)
      istat = clib_toupper(water_emiss)
      istat = clib_toupper(z0mtype)
      istat = clib_toupper(z0ttype)
      istat = clib_toupper(z0tevol)
      istat = clib_toupper(schmlake)
      istat = clib_toupper(schmriver)

      if (.not.any(schmsol == SCHMSOL_OPT)) then
         call str_concat(msg_S, SCHMSOL_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) schmsol = '//trim(schmsol)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(schmurb == SCHMURB_OPT)) then
         call str_concat(msg_S, SCHMURB_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) schmurb = '//trim(schmurb)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(schmlake == SCHMLAKE_OPT)) then
         call str_concat(msg_S, SCHMLAKE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) schmlake = '//trim(schmlake)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(schmriver == SCHMRIVER_OPT)) then
         call str_concat(msg_S, SCHMRIVER_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) schmriver = '//trim(schmriver)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(sfc_poids == SFC_POIDS_OPT)) then
         call str_concat(msg_S, SFC_POIDS_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) sfc_poids = '//trim(sfc_poids)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(diusst == DIUSST_OPT)) then
         call str_concat(msg_S, DIUSST_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) diusst = '//trim(diusst)//' : Should be one of: '//trim(msg_S))
         return
      endif
      
      if (.not.any(isba_snowfrac_bare == ISBA_SNOWFRAC_BARE_OPT)) then
         call str_concat(msg_S, ISBA_SNOWFRAC_BARE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) diusst = '//trim(isba_snowfrac_bare)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(z0mtype == Z0MTYPE_OPT)) then
         call str_concat(msg_S, Z0MTYPE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) z0mtype = '//trim(z0mtype)//' : Should be one of: '//trim(msg_S))
         return
      endif

     if (.not.any(z0ttype == Z0TTYPE_OPT)) then
         call str_concat(msg_S, Z0TTYPE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) z0ttype = '//trim(z0ttype)//' : Should be one of: '//trim(msg_S))
         return
      endif

     if (.not.any(z0tevol == Z0TEVOL_OPT)) then
         call str_concat(msg_S, Z0TEVOL_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) z0tevol = '//trim(z0tevol)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if ( (min(z0tlat(1), z0tlat(2)) <  0.0)    .or. &
           (max(z0tlat(1), z0tlat(2)) > 90.0)    .or. &
           (z0tlat(1) > z0tlat(2))         )    then
         write(msg_S, '(2f0.3)') z0tlat(1), z0tlat(2)
         call msg(MSG_ERROR, '(sfc_nml_check) zz0tlat = '//trim(msg_S)//' : should be in 0. to 90. and z0tlat(1) < z0tlat(2)')
         return
      endif

      if (.not.any(sl_diag_type == SL_DIAG_TYPE_OPT)) then
         call str_concat(msg_S, SL_DIAG_TYPE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) sl_diag_type = '//trim(sl_diag_type)//' : Should be one of: '//trim(msg_S))
         return
      endif

      if (.not.any(sl_lmin_type == SL_LMIN_TYPE_OPT)) then
         call str_concat(msg_S, SL_LMIN_TYPE_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) sl_lmin_type = '//trim(sl_lmin_type)//' : Should be one of: '//trim(msg_S))
         return
      endif
      
      if (.not.any(sl_func_stab == SL_FUNC_STAB_OPT)) then
         call str_concat(msg_S, SL_FUNC_STAB_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) sl_func_stab = '//trim(sl_func_stab)//' : Should be one of: '//trim(msg_S))
         return
      endif
      if (sl_func_stab /= 'DELAGE97' .and. .not.sl_z0ref) then
         call msg(MSG_ERROR,'(sfc_nml_check) sl_z0ref must be .TRUE. for any class of stability functions (sl_func_stab) other than DELAGE97')
         return
      endif

      if (.not.any(sl_func_unstab == SL_FUNC_UNSTAB_OPT)) then
         call str_concat(msg_S, SL_FUNC_UNSTAB_OPT,', ')
         call msg(MSG_ERROR,'(sfc_nml_check) sl_func_unstab = '//trim(sl_func_unstab)//' : Should be one of: '//trim(msg_S))
         return
      endif
      if (sl_func_unstab /= 'DELAGE92' .and. .not.sl_z0ref) then
         call msg(MSG_ERROR,'(sfc_nml_check) sl_z0ref must be .TRUE. for any class of stability functions (sl_func_unstab) other than DELAGE92')
         return
      endif

      if (sl_re2 > 0.) then
         if (tdiaglim) then
            call msg(MSG_ERROR,'(sfc_nml_check) tdiaglim must be .FALSE. if sl_re2 is used (>0).')
            return
         endif
         if (sl_diag_type /= 'INTERP') then
            call msg(MSG_ERROR,'(sfc_nml_check) sl_diag_type must be INTERP if sl_re2 is used (>0).')
            return
         endif
         if (sl_ximax <= 0.) then
            call msg(MSG_ERROR,'(sfc_nml_check) sl_ximax must be >0 if sl_re2 is used (>0).')
            return
         endif
      endif

      if (sl_ximax > 0. .and. (sl_Lmin_glacier > 0. .or. sl_Lmin_seaice > 0. .or. sl_Lmin_soil > 0. &
           .or. sl_Lmin_water > 0. .or. sl_Lmin_town > 0.)) then
         call msg(MSG_ERROR,'(sfc_nml_check) no sl_Lmin* limits should be set if sl_ximax is used (>0).')
         return
      endif
         
      if (.not.RMN_IS_OK(str_toreal(ice_emiss_const,ice_emiss))) then
         call msg(MSG_ERROR,'(sfc_nml_check) ice_emiss = '//trim(ice_emiss)//' : Should be a floating point number between 0 and 1')
         return
      endif

      if (.not.RMN_IS_OK(str_toreal(snow_emiss_const,snow_emiss))) then
         call msg(MSG_ERROR,'(sfc_nml_check) snow_emiss = '//trim(snow_emiss)//' : Should be a floating point number between 0 and 1')
         return
      endif

      if (.not.RMN_IS_OK(str_toreal(isba_soil_emiss_const,isba_soil_emiss))) then
         if (.not.any(isba_soil_emiss == ISBA_SOIL_EMISS_OPT)) then
            call str_concat(msg_S, ISBA_SOIL_EMISS_OPT,', ')
            call msg(MSG_ERROR,'(sfc_nml_check) isba_soil_emiss = '//trim(isba_soil_emiss)//' : Should be a number or one of: '//trim(msg_S))
            return
         endif
      endif

      if (.not.RMN_IS_OK(str_toreal(water_emiss_const,water_emiss))) then
         call msg(MSG_ERROR,'(sfc_nml_check) water_emiss = '//trim(water_emiss)//' : Should be a floating point number between 0 and 1')
         return
      endif

      ! ------ CHECK SVS OPTIONS --------------
      IF_SVS: if (schmsol == 'SVS' .or. schmsol == 'SVS2') then

         ! check number of SOIL LEVELS
         nk = 0
         do k=1,size(dp_svs)
            if (dp_svs(k) < 0.) exit
            nk = nk+1
         enddo

         ! KDP option for SVS now OBSOLETE, USERS SHOULD USE KHYD instead
         !--- MAKE MODEL ABORT IF USER SPECIFIES A VALUE WHEN USING SVS
         if (kdp /= -1) then
            write(msg_S, '(a,1x,i0,1x,a)') '(sfc_nml_check) KDP=',kdp,'should not be used with SVS. '//&
                 'This sfc config variable is OBSOLETE and should be removed from surface_cfgs. '//&
                 'To change last/deepest soil layer considered during the accumulation '//&
                 'of lateral flow and drainage, use surface_cfg variable KHYD instead'
            call msg(MSG_ERROR, msg_S)
            return
         endif

         ! USE DEFAULT NUMBER & DEPTH of LAYERS if NOT SPECIFIED
         ! FOR KHYD: the "last soil layer considered during the accumulation of lateral flow and drainage"
         ! (1) if using default layers, and khyd not specified then:
         ! set khyd = last SVS layer
         ! (2) if using default layer and khyd specified, use khyd
         ! (3) if users specifies soil layers, then forced user to specify khyd also.
         !#TODO: should keep only the nml check here, var setting should be move to sfc_nml_post_init
         if (nk == 0) then
            call msg(MSG_INFO,'(sfc_nml_check) SVS soil layer depth dp_svs'//&
                 ' NOT SPECIFIED by user, will use DEFAULT VALUE')
            ! default number of soil layers
            nk = nl_svs_default
            ! By default "last soil layer for accumulation of lateral flow and drainage" LEVEL = nk,
            if (khyd <= 1) then
               call msg(MSG_INFO,'(sfc_nml_check) Using default SVS configuration with '//&
                    '"last hydro drainage/lateral flow layer" level KHYD = last soil layer')
               khyd = nl_svs_default
            else
               call msg(MSG_INFO, '(sfc_nml_check) Using default SVS configuration for soil '//&
                    'layer but NOT for "last hydro drainage/lateral flow layer" KHYD')
            endif
            ! default values for depth of SVS values [in meters]
            do k=1,nk
               dp_svs(k) = dp_svs_default(k)
            enddo
         else if (nk < 3) then
            write(msg_S,*) '(sfc_nml_check) SVS requires a minimum of 3 soil layers. '//&
                 'Only ',nk, ' depths were specified in dp_svs=' , dp_svs(1:nk)
            call msg(MSG_ERROR, msg_S)
            return
         else if (nk >= 3 .and. khyd <= 1) then
            msg_S = '(sfc_nml_check)  Last hydro drainage/lateral flow layer level: '//&
                 'KHYD needs to be specified by USER because not using DEFAULT SVS set-up'
            call msg(MSG_ERROR, msg_S)
            return
         endif

         ! CHECK THAT DEPTH IS INCREASING and SAVE NUMBER OF LAYERS
         do k=1,nk-1
            if (dp_svs(k) >= dp_svs(k+1)) then
               write(msg_S, '(a,1x,f0.3,1x,a,1x,i0,1x,a,1x,f0.3,1x,a,1x,i0,1x,a)') '(sfc_nml_check) SVS layer depths badly specified: '//&
                    'dp_svs(k)=', dp_svs(k), ' (k=', k, ') should be shallower than '//&
                    'dl_svs(k+1)=', dp_svs(k+1), ' (k+1=', k+1 ,')'
               call msg(MSG_ERROR, msg_S)
               return
            endif
         enddo
         ! Save depth in dl_svs ... to be used by svs
         !#TODO: this is not a check per se, should be move to sfc_nml_post_init
         nl_svs = nk
         if (.not.allocated(dl_svs)) allocate(dl_svs(nl_svs))
         do k=1,nl_svs
            dl_svs(k) = dp_svs(k)
         enddo

         if (khyd > nl_svs) then
            write(msg_S, '(a,i3,a,i3)') '(sfc_nml_check) Last hydro drainage/lateral '//&
                 'flow layer in SVS khyd=',khyd,' cannot exceed number of SVS soil layers =', nl_svs
            call msg(MSG_ERROR, msg_S)
            return
         endif

         !Check if macropores are activated only if soil freezing is activated
         !If soil freezing is activated, make sure that the soil thermal conductivity is either the model from PL1998 or TIAN2016
         !If soil freezing is activated, make sure that the soil-snow heat flux is set to either DST_HD, DST_FD, DST_MAXD (default) or ST_D_DD
         !If soil freezing is activated, make sure that the option for setting the skin conductivity/resistance for bare-ground is set to either RTH_GRND or LAM_BOU2021
         !If soil freezing is activated, make sure that the dbtm option is set to either MID (default), DEEP ort NOFL
         if (lsoil_freezing_svs1) then
            if (.not.any(soil_cond == SOIL_COND_OPT)) then
                call str_concat(msg_S, SOIL_COND_OPT,', ')
                call msg(MSG_ERROR,'(sfc_nml_check) soil_cond = '//trim(soil_cond)//' : Should be one of: '//trim(msg_S))
                return
            endif
            if (.not.any(soilgrndhf_svs1 == SOILGRNDHF_SVS1_OPT)) then
                call str_concat(msg_S, SOILGRNDHF_SVS1_OPT,', ')
                call msg(MSG_ERROR,'(sfc_nml_check) soilgrndhf_svs1 = '//trim(soilgrndhf_svs1)//' : Should be one of: '//trim(msg_S))
                return
            endif
            if (.not.any(soilsnowhf_svs1 == SOILSNOWHF_SVS1_OPT)) then
                call str_concat(msg_S, SOILSNOWHF_SVS1_OPT,', ')
                call msg(MSG_ERROR,'(sfc_nml_check) soilsnowhf_svs1 = '//trim(soilsnowhf_svs1)//' : Should be one of: '//trim(msg_S))
                return
            endif
            if (.not.any(soildbtm_svs1 == SOILDBTM_SVS1_OPT)) then
                call str_concat(msg_S, SOILDBTM_SVS1_OPT,', ')
                call msg(MSG_ERROR,'(sfc_nml_check) soildbtm_svs1 = '//trim(soildbtm_svs1)//' : Should be one of: '//trim(msg_S))
                return
            endif
            if (lmacropores_svs) then
                if (mp_alpha < 0 .or. mp_alpha > 1) then
                    call msg(MSG_ERROR, '(sfc_nml_check) mp_alpha must be within the 0 and 1 interval')
                    return
                endif
            endif
            if (lphase_change_eff_svs1) then
                if (chi_min < 0 .or. chi_min > 1) then
                    call msg(MSG_ERROR, '(sfc_nml_check) chi_min must be within the 0 and 1 interval')
                    return
                endif
            endif
        endif
         
         ! Check consistency of kplough and ktdrain values with number of svs soil layers
         ! when the option "svs_tdrains_plough" is activated. Otherwise, disregard 
         ! kplough and ktdrains values and use -1 for these.
         if (.not.svs_tdrains_plough) then
           !--- FORCE kplough and ktdrains to have non-working values
           kplough = -1
           ktdrains = -1
         else 
           !--- MAKE MODEL ABORT IF svs_tdrains_plough ACTIVATED AND kplough OR ktdrains
           !--- ARE NOT SPECIFIED
           if (kplough < 1 .or. ktdrains < 1) then
            write(msg_S, '(a,i3,a,i3,a)') '(sfc_nml_check) KPLOUGH=',kplough,' and '//&
                 'KTDRAINS=',ktdrains,' should be defined when the option '//&
                 '"SVS_TDRAINS_PLOUGH" is set to TRUE;'//&
                 ' double-check content of section "&surface_cfgs" in "sps.cfg"'
            call msg(MSG_ERROR, msg_S)
            return
           else
            !--- CHECK CONSISTENCY OF KPLOUGH AND KTDRAINS VALUES: 
            !--- THEY CANNOT BE GREATER THAN THE NUMBER OF SVS SOIL LAYERS
            if (kplough > nl_svs .or. ktdrains > nl_svs) then
              write(msg_S, '(a,i3,a,i3,a,i3,a)') '(sfc_nml_check) KPLOUGH=',kplough,' and '//&
                 'KTDRAINS=',ktdrains,' can not be greater than the number of '//&
                 'SVS SOIL LAYERS=',nl_svs,'; double-check'//&
                 ' content of section "&surface_cfgs" in "sps.cfg"'
              call msg(MSG_ERROR, msg_S)
              return
            endif
           endif
         endif

         !# check that svs_local_z0m=True when using tofd i.e, z0veg_only=.true.
         if (z0veg_only .and. .not.svs_local_z0m) then
            call msg(MSG_ERROR, '(sfc_nml_check) svs_local_z0m must be .TRUE. when using tofd=/NIL in &physics_cfgs')
            return
         endif

         if (.not.any(svs_snow_rain == SVS_SNOW_RAIN_OPT)) then
            call str_concat(msg_S, SVS_SNOW_RAIN_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) svs_snow_rain = '//trim(svs_snow_rain)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif
         
         if (.not.any(svs_snowalb == SVS_SNOWALB_OPT)) then
            call str_concat(msg_S, SVS_SNOWALB_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) svs_snowalb = '//trim(svs_snowalb)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif  

         if (.not.any(svs_snowfrac_ground == SVS_SNOWFRAC_GROUND_OPT)) then
            call str_concat(msg_S, SVS_SNOWFRAC_GROUND_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) svs_snowfrac_ground = '//trim(svs_snowfrac_ground)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif  

         
         if (.not.any(soiltext == SOILTEXT_OPT)) then
            call str_concat(msg_S, SOILTEXT_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) soiltext = '//trim(soiltext)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif
         
         if (.not.any(soil_ksat_ice == SOIL_KSAT_ICE_OPT)) then
            call str_concat(msg_S, SOIL_KSAT_ICE_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) soil_ksat_ice = '//trim(soil_ksat_ice)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif
         
         if (.not.any(svs_hrsurf_method == SVS_HRSURF_METHOD_OPT)) then
            call str_concat(msg_S, SVS_HRSURF_METHOD_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) svs_hrsurf_method = '//trim(svs_hrsurf_method)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif

          if (.not.any(isba_hrsurf_method == ISBA_HRSURF_METHOD_OPT)) then
            call str_concat(msg_S, ISBA_HRSURF_METHOD_OPT, ', ')
            call msg(MSG_ERROR, '(sfc_nml_check) isba_hrsurf_method = '//trim(svs_hrsurf_method)//&
                 ' : Should be one of: '//trim(msg_S))
            return
         endif        
         
         if (.not.any(vf_type == VFTYPE_OPT)) then
            call str_concat(msg_S, VFTYPE_OPT,', ')
            call msg(MSG_ERROR,'(sfc_nml_check) vf_type = '//trim(vf_type)//' : Should be one of: '//trim(msg_S))
            return
         endif

         
         if (vf_type == "CCILCECO" .and. .not.(use_photo)) then 
            ! Abort, need to use photosynthesis with this option'
            call msg(MSG_ERROR,'(sfc_nml_check) vf_type = '//trim(vf_type)//' can ONLY be used with USE_PHOTO=.true. ')
            return
         endif


         if(schmsol == 'SVS2') then

            if (.not.any(hsnowscheme == HSNOWSCHEME_OPT)) then
               call str_concat(msg_S, HSNOWSCHEME_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowscheme = '//trim(hsnowscheme)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(hsnowdrift_es == HSNOWDRIFT_ES_OPT)) then
               call str_concat(msg_S, HSNOWDRIFT_ES_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowdrift_es = '//trim(hsnowdrift_es)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(hsnowdrift_cro == HSNOWDRIFT_CRO_OPT)) then
               call str_concat(msg_S, HSNOWDRIFT_CRO_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowdrift_cro = '//trim(hsnowdrift_cro)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(hsnowmetamo == HSNOWMETAMO_OPT)) then
               call str_concat(msg_S, HSNOWMETAMO_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowmetamo = '//trim(hsnowmetamo)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(hsnowcond == HSNOWCOND_OPT)) then
               call str_concat(msg_S, HSNOWCOND_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowcond = '//trim(hsnowcond)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(hsnowrad == HSNOWRAD_OPT)) then
               call str_concat(msg_S, HSNOWRAD_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowrad = '//trim(hsnowrad)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(hsnowcomp == HSNOWCOMP_OPT)) then
               call str_concat(msg_S, HSNOWCOMP_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowcomp = '//trim(hsnowcomp)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(hsnowfall == HSNOWFALL_OPT)) then
               call str_concat(msg_S, HSNOWFALL_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowfall = '//trim(hsnowfall)//' : Should be one of: '//trim(msg_S))
               return
            endif

             if (.not.any(hsnowres == HSNOWRES_OPT)) then
               call str_concat(msg_S, HSNOWRES_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowres = '//trim(hsnowres)//' : Should be one of: '//trim(msg_S))
               return
            endif           

            if (.not.any(hsnowhold == HSNOWHOLD_OPT)) then
               call str_concat(msg_S, HSNOWHOLD_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) hsnowhold = '//trim(hsnowhold)//' : Should be one of: '//trim(msg_S))
               return
            endif


            if (vf_type == "CCILCECO" ) then 
               ! Abort, CCILCECO not revised for SVS2 for now'
               call msg(MSG_ERROR,'(sfc_nml_check) vf_type = '//trim(vf_type)//' cannot be used with SVS2 for now ')
               return
            endif
          
            if (svs_dynamic_z0h) then 
               ! Abort, svs_dynamic_z0h not revised for SVS2 for now'
               call msg(MSG_ERROR,'(sfc_nml_check) svs_dynamic_z0h = .true. cannot be used with SVS2 for now ')
               return
            endif

            if (.not.any(z0snow_svs2 == Z0SNOW_SVS2_OPT)) then
               call str_concat(msg_S, Z0SNOW_SVS2_OPT,', ')
               call msg(MSG_ERROR,'(sfc_nml_check) z0snow_svs2 = '//trim(z0snow_svs2)//' : Should be one of: '//trim(msg_S))
               return
            endif

            if (.not.any(lbcheat_svs2 == LBCHEAT_SVS2_OPT)) then
               call str_concat(msg_S, LBCHEAT_SVS2_OPT, ', ')
               call msg(MSG_ERROR, '(sfc_nml_check) lbcheat_svs2 = '//trim(lbcheat_svs2)//&
                    ' : Should be one of: '//trim(msg_S))
               return
            endif
         

            if (.not.any(cano_ref_forcing == CANO_REF_FORCING_OPT)) then
               call str_concat(msg_S, CANO_REF_FORCING_OPT, ', ')
               call msg(MSG_ERROR, '(sfc_nml_check) cano_ref_forcing = '//trim(cano_ref_forcing)//&
                    ' : Should be one of: '//trim(msg_S))
               return
            endif
            
            if (read_oc .and. .not. (soiltext=='SOILGRIDSV2' .or. soiltext=='GSDE') ) then
               call str_concat(msg_S, SOILTEXT_OPT, ', ')
               call msg(MSG_ERROR, '(sfc_nml_check) read_ocshould be used with SOILGRIDSV2')
               return
            endif

         endif


         ! save high and low vegetation split

          if (vf_type == "CCILCECO") then    
             ntypel = 11
             ntypeh = 10
          else
             ntypel = 13
             ntypeh = 8
          endif
          if ( .not. allocated(vl_type) ) allocate( vl_type(ntypel) )
          if ( .not. allocated(vh_type) ) allocate( vh_type(ntypeh) )
          if (vf_type == "CCILCECO") then    
            vl_type = (/ 10, 11, 12, 13, 14, 15, 16, 17, 20, 22 , 23  /)
            vh_type = (/ 4, 5, 6, 7, 8, 9, 18, 19, 25, 26 /)
         else
            vl_type = (/ 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22 , 23  /)
            vh_type = (/ 4, 5, 6, 7, 8, 9, 25, 26 /)
         endif


         if (svs_read_aldat) then
            if (any(svs_aldat < 0.) .or. any(svs_aldat > 10.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_aldat out of range [0., 10.]')
               return
            endif
         endif
         if (svs_read_d2dat) then
            if (any(svs_d2dat < 0.) .or. any(svs_d2dat > 10.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_d2dat out of range [0., 10.]')
               return
            endif
         endif
         if (svs_read_d50dat) then
            if (any(svs_d50dat < 0.) .or. any(svs_d50dat > 10.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_d50dat out of range [0., 10.]')
               return
            endif
         endif
         if (svs_read_d95dat) then
            if (any(svs_d95dat < 0.) .or. any(svs_d95dat > 10.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_d95dat out of range [0., 10.]')
               return
            endif
         endif
         if (svs_read_vegdat) then
            if (any(svs_vegdat < -99.) .or. any(svs_vegdat > 1.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_vegdat out of range [-99.., 1.]')
               return
            endif
         endif
         if (svs_read_z0mdat) then
            if (any(svs_z0mdat < 0.) .or. any(svs_z0mdat > 10.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_z0mdat out of range [0., 10.]')
               return
            endif
         endif
         if (svs_read_vf2ctemdat) then
            do k=1,nclass
               if ( (k .lt. 4) .or. (k .eq. 21) .or. (k .gt. 23)) then
                  if ((svs_vf2ctemdat(k) .ge. 1) .and. (svs_vf2ctemdat(k) .le. 9)) then       
                     write(msg_S, '(a)') '(sfc_nml_check) Values of svs_vf2ctemdat for classes '//&
                                         'that can not be remapped have to be different from [1, 9]. '//&
                                         'All classes can be remapped except 1, 2, 3, 21, 24, 25 and 26'
                     call msg(MSG_ERROR, msg_S)
                     return
                  endif
               else
                  if ((svs_vf2ctemdat(k) .lt. 1) .or. (svs_vf2ctemdat(k) .gt. 9)) then     
                     write(msg_S, '(a)') '(sfc_nml_check) Value of svs_vf2ctemdat for classes '//&
                                         'that can be remapped is out of range [1, 9]. '//&
                                         'All classes can be remapped except 1, 2, 3, 21, 24, 25 and 26'
                     call msg(MSG_ERROR, msg_S)
                     return
                  endif
               endif
            enddo
         endif
         if (svs_read_vegcrops) then
            if (any(svs_vegcrops < 0.) .or. any(svs_vegcrops > 1.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_vegcrops out of range [0., 1.]')
               return
            endif
         endif
         if (svs_read_vegdat14) then
            if (any(svs_vegdat14 < 0.) .or. any(svs_vegdat14 > 1.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_vegdat14 out of range [0., 1.]')
               return
            endif
         endif
         if (svs_read_vegdat16) then
            if (any(svs_vegdat16 < 0.) .or. any(svs_vegdat16 > 1.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_vegdat16 out of range [0., 1.]')
               return
            endif
         endif
         if (svs_read_vegdat17) then
            if (any(svs_vegdat17 < 0.) .or. any(svs_vegdat17 > 1.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_vegdat17 out of range [0., 1.]')
               return
            endif
         endif
         if (svs_read_lai11) then
            if (any(svs_lai11 < 0.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_lai11 out of range < 0.')
               return
            endif
         endif
         if (svs_read_lai14) then
            if (any(svs_lai14 < 0.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_lai14 out of range < 0.')
               return
            endif
         endif
         if (svs_read_lai15) then
            if (any(svs_lai15 < 0.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_lai15 out of range < 0.')
               return
            endif
         endif
         if (svs_read_lai16) then
            if (any(svs_lai16 < 0.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_lai16 out of range < 0.')
               return
            endif
         endif
         if (svs_read_lai17) then
            if (any(svs_lai17 < 0.)) then
               call msg(MSG_ERROR,'(sfc_nml_check) svs_lai17 out of range < 0.')
               return
            endif
         endif

      endif IF_SVS


      ! ------ CHECK TEB OPTIONS --------------
      IF_TEB: if (schmurb == 'TEB' ) then

         if (.not.any(teb_snow == TEB_SNOW_OPT)) then
               call str_concat(msg_S, TEB_SNOW_OPT, ', ')
               call msg(MSG_ERROR, '(sfc_nml_check) teb_snow = '//trim(teb_snow)//&
                    ' : Should be one of: '//trim(msg_S))
               return
            endif

         if (teb_snow == 'CUSTOM') then
            do k=1, nb_teb_snowpar
             if (teb_snroof(k) == -1) then
                 call msg(MSG_INFO,'(sfc_nml_check) TEB snow parameters over ROOF teb_snroof'//&
                 ' NOT SPECIFIED by user, will use DEFAULT VALUE')
                  teb_snroof = TEB_SNROOF_DEFAULT
                  EXIT
             endif
            enddo
            do k=1, nb_teb_snowpar
              if (teb_snroad(k) == -1) then
                 call msg(MSG_INFO,'(sfc_nml_check) TEB snow parameters over ROAD teb_snroad'//&
                 ' NOT SPECIFIED by user, will use DEFAULT VALUE')
                  teb_snroad = TEB_SNROAD_DEFAULT
                  EXIT
               endif
            enddo
         endif

         if (.not.any(teb_hydro == TEB_HYDRO_OPT)) then
               call str_concat(msg_S, TEB_HYDRO_OPT, ', ')
               call msg(MSG_ERROR, '(sfc_nml_check) teb_hydro = '//trim(teb_hydro)//&
                    ' : Should be one of: '//trim(msg_S))
               return
            endif
        if (teb_hydro == 'MA00  ') then
            do k=1, nb_teb_hydropar
                  teb_hydropar = TEB_HYDROPAR_DEFAULT
            enddo
        endif

        if (teb_hydro == 'CUSTOM') then
            do k=1, nb_teb_hydropar
             if (teb_hydropar(k) == -1) then
                 call msg(MSG_INFO,'(sfc_nml_check) TEB hydrological parameters over ROOF and ROAD '//&
                 ' NOT SPECIFIED by user, will use DEFAULT VALUE')
                  teb_hydropar = TEB_HYDROPAR_DEFAULT
                  EXIT
             endif
            enddo
         endif

      endif IF_TEB


      ! check that use crodebug only with SV2
      ! read in environement variables here if ok to avoid looping on TRNCH
      if( svs2_crodebug ) then
         if ( schmsol == 'SVS2' ) then
            ! Activate debug option for Crocus (SVS2)
            call INIT_CRODEBUG(HSNOWSCHEME)
         else
            write(msg_S,'(a, l1, a)') '(sfc_nml_check) svs2_crodebug = ', svs2_crodebug, ' can only be used with schmsol = SVS2'
            call msg(MSG_ERROR,msg_S)
            return
         endif
      endif
      

      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function sfc_nml_check


   function sfc_nml_post_init() result(m_istat)
      use tdpack_const, only: PI
      use sfclayer_compz0, only: compz0_set
      implicit none
      integer :: m_istat

      integer :: istat, options, iverb, idxmax
      !----------------------------------------------------------------
      m_istat = RMN_ERR
 
      !# nsurf nb of surface types (schemes)
      !# phybusinit / gesdict defines nrow = nagg = nsurf+1
      idxmax = max(indx_soil, indx_glacier, indx_water, indx_ice, indx_agrege)
      if (schmurb /= 'NIL') idxmax = max(idxmax, indx_urb)
      if (schmlake  /= 'NIL') idxmax = max(idxmax, indx_lake)
      if (schmriver /= 'NIL') idxmax = max(idxmax, indx_river)
      nsurf = idxmax - 1

      !# Z0tlat to radians
      z0tlat(:) = z0tlat(:) * PI/180.

      !# Put some sfc options value on th WB for the main physic to access
      options = WB_REWRITE_NONE+WB_IS_LOCAL
      iverb   = wb_verbosity(WB_MSG_INFO)
      istat = WB_OK
      istat = min(wb_put('sfc/icelac', icelac, options), istat)
      istat = min(wb_put('sfc/nclass', nclass, options), istat)
      istat = min(wb_put('sfc/nicelvl', nl, options), istat)
      istat = min(wb_put('sfc/nsurf', nsurf, options), istat)
      istat = min(wb_put('sfc/schmsol', schmsol, options), istat)
      istat = min(wb_put('sfc/schmurb', schmurb, options), istat)
      istat = min(wb_put('sfc/schmlake', schmlake, options), istat)
      istat = min(wb_put('sfc/schmriver', schmriver, options), istat)
      istat = min(wb_put('sfc/snoalb_anl', snoalb_anl, options), istat)
      iverb = wb_verbosity(iverb)
      if (istat /= WB_OK) then
         call msg(MSG_ERROR, '(sfc_nml) Problem in wb_put')
         return
      endif

      !# Surface Layer module init
      istat = SL_OK
      if (istat == SL_OK) istat = sl_put('afd', sl_afd)
      if (istat == SL_OK) istat = sl_put('beta', beta)
      if (istat == SL_OK) istat = sl_put('ximax',sl_ximax)
      if (istat == SL_OK) istat = sl_put('re2', sl_re2)
      if (istat == SL_OK) istat = sl_put('rineutral',sl_rineutral)
      if (istat == SL_OK) istat = sl_put('tdiaglim',tdiaglim)
      if (istat == SL_OK) istat = sl_put('tdlrate',tdlrate)
      if (istat == SL_OK) istat = sl_put('sl_stabfunc_stab',sl_func_stab)
      if (istat == SL_OK) istat = sl_put('sl_stabfunc_unstab',sl_func_unstab)
      if (istat == SL_OK) istat = sl_put('z0ref',sl_z0ref)
      if (istat == SL_OK) istat = sl_put('diag_type',sl_diag_type)
      if (istat == SL_OK) istat = sl_put('lmin_type',sl_lmin_type)
      if (istat /= SL_OK) then
         call msg(MSG_ERROR,'(sfc_nml) Problem initializing the surface layer module')
         return
      endif

      istat = RMN_OK
      if (istat == RMN_OK) istat = compz0_set('Z0MIN', z0min)
      if (istat == RMN_OK) istat = compz0_set('Z0GLA', Z0GLA)
      if (istat == RMN_OK) istat = compz0_set('Z0HCON', z0hcon)
      if (istat == RMN_OK) istat = compz0_set('Z0LAT', z0tlat(1), z0tlat(2))
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(sfc_nml) Problem initializing the compz0 module')
         return
      endif

      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function sfc_nml_post_init

 end function sfc_nml2
end module sfc_nml_mod
