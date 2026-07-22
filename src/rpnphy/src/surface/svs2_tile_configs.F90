module svs2_tile_configs

  ! Store functions used in SVS2 to aggregate the SVS2 surface tiles (bare ground, vegetation and 2 snowpacks)

  implicit none

  ! INDEXES FOR SVS SURFACES 
  integer, parameter :: svs2_tilesp1  = 7 ! SVS tile + 1 for aggregated value
  !
  integer, parameter :: indx_svs2_ag  = 1 ! AGGREGATED VALUE 
  integer, parameter :: indx_svs2_bg  = 2 ! EXPOSED BARE GROUND (NOT COVERED BY SNOW)
  integer, parameter :: indx_svs2_vl  = 3 ! LOW VEGETATION (NOT COVERED BY SNOW)
  integer, parameter :: indx_svs2_vh  = 4 ! HIGH VEGETATION (NOT COVERED BY SNOW)
  integer, parameter :: indx_svs2_sn  = 5 ! SNOW ON BARE GROUND & LOW VEG
  integer, parameter :: indx_svs2_sv  = 6 ! SNOW BELOW HIGH VEG  
  integer, parameter :: indx_svs2_gv  = 7 ! GROUND BELOW HIGH VEG  (NOT COVERED BY SNOW)

contains

     real function AG_SVS2(FBG,FVL, FVH, FSNO,FSNV, FGVH, BG,VL,VH,SNO,SNV, GVH)
    !     FUNCTION TO AGGREGATE SURFACE VARIABLES OVER SVS2 SURFACE TYPES
    !     (Bare Ground, low vegetation, high vegetation, 
    !     low-veg-and-bare-ground snow, under-high-veg snow, under-high-veg ground 

    !          - Input -
    ! FBG      Fraction of BARE GROUND
    ! FVL      Fraction of LOW VEGETATION
    ! FVH      Fraction of HIGH VEGETATION
    ! FSNO     Fraction of SNOW ON LOW VEG + BARE GROUND
    ! FSNV     Fraction of SNOW UNDER HIGH VEG
    ! FGVH     Fraction of GROUND UNDER  HIGH VEG
    ! BG       Bare Ground Variable 
    ! VL       Low vegetation  Variable
    ! VH       High vegetation  Variable
    ! SNO      Snow over low veg & bare ground Variable 
    ! SNV      Snow-under-high vegetation variable
    ! GVH      Ground-under-high vegetation variable
    
    !         - Output -
    ! AG       Aggregated Variable
    REAL, intent(in) :: FBG,FVL,FVH,FSNO,FSNV,FGVH
    REAL, intent(in) :: BG,VL,VH,SNO,SNV, GVH
     
    AG_SVS2 =   FBG             * BG   &
         + FVL             * VL   &
         + FVH             * VH   &
         + FSNO            * SNO  &
         + FSNV            * SNV  &  
         + FGVH            * GVH
   
    return
  end function AG_SVS2

  
  subroutine weights_svs2(FVH,FVL,FSN,FSNVH,FSURFVH, HTYPE, NI, WT)
    implicit none
    !     FUNCTION TO CALCULATE WEIGHTS FOR AGGREGATION OF SVS2 VARIABLE FOR the 6 SVS2 SURFACES
    !     (Bare Ground, low vegetation, high vegetation, 
    !     low-veg-and-bare-ground snow, under-high-veg snow, under-high-veg ground     
    !         - Input -
    ! FVH     Fraction of grid cell covered by HIGH vegetation
    ! FVL     Fraction of grid cell covered by LOW  vegetation
    ! FSN     Fraction of LOW vegetation or bare ground covered by snow
    ! FSNVH   Fraction of ground covered by snow under high vegetation 
    ! FSURFVH Fraction of surface seen from the atmosphere through sparse high veg. 
    ! HTYPE   Fraction seen from the ground ot the atmosphere
    ! NI       Horizontal dimension
    ! 
    !
    !         - Output -
    ! WT       WEIGHTS FOR 6 SVS2 SURFACES + AGGREGATED VALUE
    !          LEVEL 1 : SUM of LEVEL 2-7
    !          LEVEL 2 : EXPOSED BARE GROUND (NOT COVERED BY SNOW)
    !          LEVEL 3 : LOW VEGETATION (NOT COVERED BY SNOW)
    !          LEVEL 4 : HIGH VEGETATION (NOT COVERED BY SNOW)
    !          LEVEL 5 : SNOW ON BARE GROUND AND LOW VEGETATION
    !          LEVEL 6 : SNOW BELOW HIGH VEGETATION 
    !          LEVEL 7 : GROUND BELOW HIGH VEG  (NOT COVERED BY SNOW)
    !  
    INTEGER I, NI
    REAL, dimension(ni) :: FVH,FVL,FSNVH,FSN, FSURFVH
    REAL, dimension(ni,svs2_tilesp1) :: WT
    CHARACTER(LEN=*)          :: HTYPE 
    REAL, parameter ::  min_wt=0.001 ! [1.=100%] minimum value of weight, otherwise set to zero.

    DO I=1,NI

      !  LOW VEGETATION (NOT COVERED BY SNOW)
      WT(i,indx_svs2_vl)  =      FVL(i)            * (1. - FSN(i)  )  
      if( WT(i,indx_svs2_vl) .lt. min_wt )  WT(i,indx_svs2_vl) = 0.0  

      !   SNOW ON BARE GROUND AND LOW VEGETATION
      WT(i,indx_svs2_sn)  =  (1. - FVH(i)          ) *       FSN(i)
      if( WT(i,indx_svs2_sn) .lt. min_wt )  WT(i,indx_svs2_sn) = 0.0
       
      ! Case of high vegetation 
      IF(HTYPE =='ATM') THEN ! Fraction seen from the atmosphere 
          !    HIGH VEGETATION (NOT COVERED BY SNOW)
          WT(i,indx_svs2_vh)  =      FVH(I) * (1. - FSURFVH(I)) 
          !   SNOW BELOW HIGH VEGETATION 
          WT(i,indx_svs2_sv)  =      FVH(I) * FSURFVH(I)   * FSNVH(i) 
          !    GROUND BELOW HIGH VEGETATION 
          WT(i,indx_svs2_gv)  =     FVH(i)  * FSURFVH(I) *  (1- FSNVH(i))     
      ELSE IF(HTYPE =='GROUND') THEN   ! Fraction seen from the ground
          !    HIGH VEGETATION (NOT COVERED BY SNOW)
          ! Here we assume that WTG(vh)=VGH to use this weight in the fraction of turbulent fluxes towards the atmosphere
          ! All turbulent fluxes are initially scaled with VEG_DENS to account for canopy sparseness
          WT(i,indx_svs2_vh)  =    FVH(i) 
          !   SNOW BELOW HIGH VEGETATION 
          WT(i,indx_svs2_sv)  =        FVH(i) *  FSNVH(i)
          !    GROUND BELOW HIGH VEGETATION 
          WT(i,indx_svs2_gv)  =        FVH(i)  *  (1- FSNVH(i)) 
      ENDIF

      if( WT(i,indx_svs2_vl) .lt. min_wt )  WT(i,indx_svs2_vl) = 0.0
      if( WT(i,indx_svs2_sv) .lt. min_wt )  WT(i,indx_svs2_sv) = 0.0             
      if( WT(i,indx_svs2_gv) .lt. min_wt )  WT(i,indx_svs2_gv) = 0.0  
      if( (WT(i,indx_svs2_gv)+WT(i,indx_svs2_sv)) .lt. min_wt )  WT(i,indx_svs2_vh) = 0.0  

      IF(HTYPE =='ATM') THEN ! Fraction seen from the atmosphere 
          !     EXPOSED BARE GROUND (NOT COVERED BY SNOW)
          WT(i,indx_svs2_bg) = 1.0 - WT(i,indx_svs2_sn) -  WT(i,indx_svs2_sv) - WT(i,indx_svs2_vl) &
                               -  WT(i,indx_svs2_vh) -  WT(i,indx_svs2_gv)
          ! Aggregate
          WT(i,indx_svs2_ag) =  WT(i,indx_svs2_bg) + WT(i,indx_svs2_vl) + WT(i,indx_svs2_vh) + &
                             WT(i,indx_svs2_sn) + WT(i,indx_svs2_sv) + WT(i,indx_svs2_gv)

      ELSE IF(HTYPE =='GROUND') THEN   ! Fraction seen from the ground
          !     EXPOSED BARE GROUND (NOT COVERED BY SNOW)
          WT(i,indx_svs2_bg) = 1.0 - WT(i,indx_svs2_sn)  - WT(i,indx_svs2_vl) &
                               -  WT(i,indx_svs2_vh) 
          ! Aggregate
          WT(i,indx_svs2_ag) =  WT(i,indx_svs2_bg) + WT(i,indx_svs2_vl) + WT(i,indx_svs2_vh) + &
                                 WT(i,indx_svs2_sn) 
      ENDIF


    ENDDO
  
    return
  end subroutine weights_svs2


end module svs2_tile_configs  
