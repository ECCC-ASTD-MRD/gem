module svs_configs
  !!!!! ------------- IMPORTANT OPEN MP INFORMATION --- PLEASE READ   ---------------  !!!!!
  ! The physics of GEM is in open-mp by default. SVS_CONFIG is mostly called in the PHYSICS and should thus be coded
  ! with open-mp "safe" code.
  ! In open-mp, "save" and  conditional allocates can lead to problems (race condition) 
  ! The save/allocate below are OK because they are ACTUALLY used in the DYNAMICS in areas with no open-mp. Specifically :
  ! 1) the config reading section of GEM: sfc_nml --> phy_nml --> itf_phy_nml --> set_world_view --> gemdm
  ! 2) the physics/surface bus initialization of GEM : sfc_businit --> phybusinit --> phydebu --> phy_init --> itf_phy_init --> gemdm
  ! >>>> If you plan to add save or allocates... please make sure they are done in NO OPENMP regions <<<<< 
  !----------------------------------------------------------------------------------  !!!!!

  ! MODULE to SAVE SVS INFORMATION SPECIFIED in NML, to initialize SVS variables that are NML dependent 
  ! ALSO contains 2 function used in SVS to aggregate the SVS surface tiles (bare ground, vegetation and 2 snowpacks)


  implicit none

  ! INDEXES FOR SVS SURFACES 
  integer, parameter :: svs_tilesp1  = 5 ! SVS tile + 1 for aggregated value
  !
  integer, parameter :: indx_svs_ag  = 1 ! AGGREGATED VALUE 
  integer, parameter :: indx_svs_bg  = 2 ! BARE GROUND
  integer, parameter :: indx_svs_vg  = 3 ! VEGETATION (NOT COVERED BY SNOW)  
  integer, parameter :: indx_svs_sn  = 4 ! SNOW ON BARE GROUND & LOW VEG
  integer, parameter :: indx_svs_sv  = 5 ! SNOW ON HIGH VEG 
 


  ! Number of soil ("GROUND") layers for SVS land scheme
  integer, save :: nl_svs

  ! Depth of deep soil layers in METERS   
  real, allocatable, save :: dl_svs(:)
  !=  (/ 0.05, 0.1, 0.2, 0.4, 1.0, 2.0, 3.0 /) 


  ! Thicknesses of layers in METERS
  real,  allocatable, save :: delz(:) !nl_svs

  !---------------------------------------------------------------
  ! SOIL TEXTURE OPTIONS: 
  !  1:   "GSDE"      --  8 layers of soil texture from CHINESE DATASET !!!
  !  2:   "SLC"       --  5 layers of soil texture:  SOIL LANDSCAPE of CANADA
  !  3:   "SOILGRIDS" --  7  layers of soil texture:  ISRIC ? World Soil Information
  ! ENTRY  bus number of levels for clay & sand variables
  integer,  save :: nl_ste
  ! PERMANENT PHYSICS bus number of levels for clay & sand variables
  integer, save :: nl_stp

  ! GSDE NUMBER OF LAYERS 
  integer,  parameter :: nl_gsde = 8 
  ! GSDE SOIL DEPTH in METERS... 
  real, parameter , dimension(nl_gsde):: dl_gsde =  (/ 0.045, 0.091, 0.166, 0.289, 0.493, 0.829, 1.38, 2.229 /) 
  

   ! SLC NUMBER OF LAYERS 
  integer,  parameter :: nl_slc = 6
  ! SLC SOIL DEPTH in METERS... ( 0-5cm, 15-30cm, 30-60cm, 60-100cm, 100-200cm ) 
  real, parameter , dimension(nl_slc):: dl_slc =  (/ 0.05, 0.15, 0.3, 0.6, 1.0, 2.0 /) 
  
  ! SOILGRIDS NUMBER OF LAYERS 
  integer,  parameter :: nl_soilgrids = 7 
  ! SOILGRIDS SOIL DEPTH in METERS... 
  real, parameter , dimension(nl_soilgrids):: dl_soilgrids =  (/ 0.025, 0.1, 0.225, 0.45, 0.8, 1.5, 2.0 /) 


  !  WEIGHTS TO MAP soil parameters calculated on SOIL TEXTURE LAYERS  layers unto MODEL SOIL LAYERS 
  real, allocatable, save :: weights(:,:)       !(max_nl_svs,nl_soil_texture)


  !------------------------------------------------------------------
  ! NUMERICAL METHOD HYDRO
  ! Numerical method used to solve Richards' equations in hydro_svs
  ! hydro_svs_method=0: Euler, forward, 1st order
  ! hydro_svs_method=1: Runge-Kutta 4th order (RK4)
  integer, parameter :: hydro_svs_method = 1

  !---------------------------------------------------------
  ! OPTION FOR CALCULATION of AVERAGE LAND SURFACE TEMPERATURE AND HUMIDITY
  ! D. Deacu 06/2018
  !!!!!  use_eff_surf_tq  --- moved to surface namelist 
  ! use_eff_surf_tq=.FALSE. :  Area-average only calculation for sfc T and Hum.
  ! use_eff_surf_tq=.TRUE.  :  NEW OPTION that uses effective surface temperature
  !                            and specific humidity instead of composite (area-averaged only) 
  !                            counterparts in surface flux calculations
  !           Reference: Chehbouni, A., E.G. Njoku, J. Lhomme, and Y.H. Kerr, 1995: Approaches for Averaging 
  !                      Surface Parameters and Fluxes over Heterogeneous Terrain. J. Climate, 8, 1386?1393
  !---------------------------------------------------------

  !-----------------------------------------------------------------
  ! VEG HIGH/LOW SPLIT 
  integer, save :: ntypel 
  integer, save :: ntypeh
  integer, allocatable, save :: vl_type(:)
  integer, allocatable, save :: vh_type(:)

  ! Urban Class: 21, split between high and low veg
  ! fraction of urban in low veg (+ fraction of urban in high veg = 1.0)
  real, parameter :: furb_vl = 0.5
  

  !---------------------------------------------------------------  
  ! other constants...
  ! M_TO_MM : CONVERT METERS TO MM
  real, parameter :: M_TO_MM = 1000.00
  ! FRACTION OF URBAN CLASS CONSIDERED IMPERVIOUS FOR RUNOFF CALC. IN HYDRO...
  REAL, PARAMETER :: IMP_URB = 0.33
  ! TRACE AMOUNT OF SNOW MASS [kg/m3] (measured=0.1 ... but in model allow smaller threshold)
  REAL, PARAMETER :: CRITSNOWMASS=0.01 ! model integration snow mass criteria
  REAL, PARAMETER :: LAI0 = 1.0  ! DEFAULT LAI VALUE
  ! espilon to be used to avoid division by zero/testing on zero (defined here 
  ! so that consistent between subroutines
  REAL, PARAMETER :: EPSILON_SVS = 1.00E-7
  ! espilon to be used for testing temperature (in Kelvin)
  REAL, PARAMETER :: EPSILON_SVS_TK = 1.00E-4

  ! Constants related to roughness
  !---------------------------------------------------------
!  ROUGNESS
!

!  GENERAL COMMENT:
!   want z0m *with orography* at resolution coarser than 3km... 
!  and local onl_svsy at higher resolution. Makes sure that the local z0 
!  calculated based on veg. height database, and NOT
!  on look-up table (values are quite weak...)
!
!     Conversion factor to convert from momemtum roughness to thermal roughness
  REAL, PARAMETER :: Z0M_TO_Z0H = 0.2
!
!     Thermal roughness for snow
  REAL, PARAMETER :: Z0HSNOW = 0.010
!
!     Momentum roughness for bare ground 
  REAL, PARAMETER :: Z0MBG = 0.05
!---------------------------------------------------------
!
 

  
  private :: weights_soil_texture


contains


  subroutine init_soil_text_levels()

    use sfc_options
    implicit none

    ! OBJECT:  Initialize number of levels of GEOPHY soil texture based on sfc nml input
    !          Calculate weights to map GEOPHY soil texture unto SVS layers
    !          Also calculate thickness of SVS soil layers here for OPEN MP reasons
    nl_ste = 3
    nl_stp = 3

    ! number of levels of entry and bus clay, sand variables
    if ( soiltext == "GSDE" ) then

       ! ENTRY  bus number of levels for clay & sand variables
       nl_ste = nl_gsde
       ! PERMANENT PHYSICS bus number of levels for clay & sand variables
       nl_stp = nl_gsde

       call weights_soil_texture()

    else  if ( soiltext == "SLC" ) then

       ! ENTRY  bus number of levels for clay & sand variables
       nl_ste = nl_slc
       ! PERMANENT PHYSICS bus number of levels for clay & sand variables
       nl_stp = nl_slc
       call weights_soil_texture()

    else if ( soiltext == "SOILGRIDS" ) then

       ! ENTRY  bus number of levels for clay & sand variables
       nl_ste = nl_soilgrids
       ! PERMANENT PHYSICS bus number of levels for clay & sand variables
       nl_stp = nl_soilgrids
       call weights_soil_texture()

    endif
    ! Calculate SVS soil layer thickness from layer Depths specified by sfc nml
    call layer_thickness()
    
    return
  end subroutine INIT_SOIL_TEXT_LEVELS


  subroutine weights_soil_texture()

#include <rmn/msg.h>
    use sfc_options
    implicit none
    integer :: unout
    integer, external :: msg_getUnit
    integer k, kk
    real d_svs(nl_svs+1), d_soil_texture(nl_stp+1)

  
    if ( .not. allocated(weights) ) allocate( weights(nl_svs,nl_stp) )

    ! defined array of depth with value starting at 0m, to make looping easier
    ! for same reason , set maximum depth of SOIL_TEXTURE equal to that of MODEL layers... i.e., assume
    ! texture constant past maximum depth of SOIL_TEXTURE texture
    d_svs(1) = 0
    do k=2,nl_svs+1
       d_svs(k)=dl_svs(k-1)
    enddo

    ! first depth set to 0
    d_soil_texture(1) = 0
    if ( soiltext == "GSDE" ) then
    
       do k=2,nl_stp
          d_soil_texture(k)=dl_gsde(k-1)
       enddo
       ! last depth of soil texture database set to max depth of SVS 
       ! i.e, deepest soil texture measured extends to the bottom of last SVS layer
       d_soil_texture(nl_stp+1)=max( dl_svs(nl_svs) , dl_gsde(nl_stp) )
    else  if ( soiltext == "SLC" ) then
         
       do k=2,nl_stp
          d_soil_texture(k)=dl_slc(k-1)
       enddo
       ! last depth of soil texture database set to max depth of SVS 
       ! i.e, deepest soil texture measured extends to the bottom of last SVS layer
       d_soil_texture(nl_stp+1)=max( dl_svs(nl_svs) , dl_slc(nl_stp) )
    else if ( soiltext == "SOILGRIDS" ) then
     
       do k=2,nl_stp
          d_soil_texture(k)=dl_soilgrids(k-1)
       enddo
       ! last depth of soil texture database set to max depth of SVS 
       ! i.e, deepest soil texture measured extends to the bottom of last SVS layer
       d_soil_texture(nl_stp+1)=max( dl_svs(nl_svs) , dl_soilgrids(nl_stp) )    


    endif
    
    
    ! calculate weights

    do k = 1,nl_svs ! model layers
       do kk = 1, nl_stp ! database layers

          weights(k,kk) = max( min(d_svs(k+1),d_soil_texture(kk+1)) - max(d_svs(k),d_soil_texture(kk) ) , 0.0 ) &
               / (d_svs(k+1) - d_svs(k)) 

       end do
    end do


    ! PRINT INFO in LISTING
    unout         = msg_getUnit(MSG_INFO)
    if (unout > 0) then
       if ( soiltext == "GSDE" ) then
          write(unout, *) ' ****** GSDE SOIL TEXTURE ******* '
       else if (soiltext == "SLC" ) then
          write(unout, *) ' ****** SLC SOIL TEXTURE ******* '
       else if (soiltext == "SOILGRIDS" ) then
          write(unout, *) ' ****** SOILGRIDS SOIL TEXTURE ******* '
       endif
          
       write(unout, *) ' ****** SOIL MAPPING WEIGHTS [METERS] ******* '
       do k = 1,nl_svs ! model layers

          write(unout, *) '--- for SVS layer k=', k, ' depth=', dl_svs(k), ' ---'

          if ( soiltext == "GSDE" ) then
             do kk = 1, nl_stp ! database layers
                write(unout, *) 'for GSDE layer kk=', kk,' depth=', dl_gsde(kk),' weight=', weights(k,kk)
             enddo
          else if (soiltext == "SLC" ) then
             do kk = 1, nl_stp ! database layers
                write(unout, *) 'for SLC layer kk=', kk,' depth=', dl_slc(kk),' weight=', weights(k,kk)
             enddo
          else if (soiltext == "SOILGRIDS" ) then
             do kk = 1, nl_stp ! database layers
                write(unout, *) 'for SOILGRIDS layer kk=', kk,' depth=', dl_soilgrids(kk),' weight=', weights(k,kk)
             enddo
          endif
       enddo
    endif

    return
  end subroutine weights_soil_texture

  subroutine layer_thickness()
    implicit none

    integer k
    if ( .not. allocated(delz) ) allocate( delz(nl_svs) )
    DELZ(1)=DL_SVS(1)            
    DO K=2,NL_SVS        
       DELZ(K)=DL_SVS(K)-DL_SVS(K-1)
    ENDDO

  end subroutine layer_thickness

     real function AG(FBG,FVG,FSNO,FSNV,BG,VG,SNO,SNV)
    !     FUNCTION TO AGGREGATE SURFACE VARIABLES OVER SVS SURFACE TYPES
    !     (Bare Ground, Vegetation, low-veg-and-bare-ground snow, under-high-veg snow)

    !          - Input -
    ! FBG      Fraction of BARE GROUND
    ! FVG      Fraction of VEGETATION
    ! FSNO     Fraction of SNOW ON LOW VEG + BARE GROUND
    ! FSNV     Fraction of SNOW UNDER HIGH VEG
    ! BG       Bare Ground Variable 
    ! VG       Vegetation  Variable
    ! SNO      Snow over low veg & bare ground Variable 
    ! SNV      Snow-under-vegetation variable
    !
    !         - Output -
    ! AG       Aggregated Variable
    REAL, intent(in) :: FBG,FVG,FSNO,FSNV,BG,VG,SNO,SNV
     
    AG =   FBG             * BG   &
         + FVG             * VG   &
         + FSNO            * SNO  &
         + FSNV            * SNV  
    
    return
  end function AG

  
  subroutine weights_svs(FVH,FVL,FSN,FSNVH,NI, WT)
    implicit none
    !     FUNCTION TO CALCULATE WEIGHTS FOR AGGREGATION OF SVS VARIABLE FOR the 4 SVS SURFACES
    !     (Bare Ground, Vegetation, low-veg-and-bare-ground snow, under-high-veg snow)
    !         - Input -
    ! FVH      Fraction of grid cell covered by HIGH vegetation
    ! FVL      Fraction of grid cell covered by LOW  vegetation
    ! FSN      Fraction of LOW  vegetation or bare ground covered by snow
    ! FSNVH   Fraction of HIGH vegetation covered by snow (from ground or atm)*** NOTE THE USER MUST PROVIDE THE APPROPRIATE FRACTION !
    ! NI       Horizontal dimension
    !
    !         - Output -
    ! WT       WEIGHTS FOR 4 SVS SURFACES + AGGREGATED VALUE
    !          LEVEL 1 : BARE GROUND
    !          LEVEL 2 : VEGETATION
    !          LEVEL 3 : SNOW ON BARE GROUND AND LOW VEGETATION
    !          LEVEL 4 : SNOW "ON" HIGH VEGETATION (
    !          LEVEL 5 : SUM of LEVEL 1-4 
    !  
    INTEGER I, NI
    REAL, dimension(ni) :: FVH,FVL,FSNVH,FSN
    REAL, dimension(ni,svs_tilesp1) :: WT
    REAL, parameter ::  min_wt=0.001 ! [1.=100%] minimum value of weight, otherwise set to zero.

    DO I=1,NI
      
       WT(i,indx_svs_sn)  =  (1. - FVH(i)          ) *       FSN(i)
       if( WT(i,indx_svs_sn) .lt. min_wt )  WT(i,indx_svs_sn) = 0.0

       WT(i,indx_svs_sv)  =        FVH(i)            *       FSNVH(i)
       if( WT(i,indx_svs_sv) .lt. min_wt )  WT(i,indx_svs_sv) = 0.0

       WT(i,indx_svs_vg)  =        FVH(i)            * (1. - FSNVH(i)) &
                                 + FVL(i)            * (1. - FSN(i)  )  
       if( WT(i,indx_svs_vg) .lt. min_wt )  WT(i,indx_svs_vg) = 0.0

       !
       !WT(i,indx_svs_bg)  =  (1. - FVH(i) - FVL(i) ) * (1. - FSN(i) )
       WT(i,indx_svs_bg) = 1.0 - WT(i,indx_svs_sn) -  WT(i,indx_svs_sv) - WT(i,indx_svs_vg)

 
       ! CLEANUP -- NOT DONE FOR NOW ...

       ! Aggregate
       WT(i,indx_svs_ag) =  WT(i,indx_svs_bg) + WT(i,indx_svs_vg) + &
                            WT(i,indx_svs_sn) + WT(i,indx_svs_sv)
    ENDDO
  
    return
  end subroutine weights_svs


end module svs_configs
  
