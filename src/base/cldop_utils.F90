module cldop_utils
  private

  ! API functions
  public :: cldop_propli                                !Longwave radiative properties for ice clouds
  public :: cldop_proplw                                !Longwave radiative properties for water clouds
  public :: cldop_propsi                                !Shortwave radiative properties for ice clouds
  public :: cldop_propsw                                !Shortwave radiative properties for water clouds
  public :: cldop_rei                                   !Effective ice radius calculations
  public :: cldop_rew                                   !Effective water radius calculations

  ! Internal constants
#include "nbsnbl.cdk"
  real :: aws(4,NBS), bws(4,NBS), cws(4,NBS)
  real :: awl(5,NBL), bwl(4,NBL), cwl(4,NBL)
  real :: ais(2,NBS), bis(4,NBS), cis(4,NBS)
  real :: ail(3,NBL), bil(4,NBL), cil(4,NBL)

  ! Internal variables
  integer :: i, j

  !     --------------------------------------------------------
  !     new water properties for sw (dobbie, li, and chylek 1999, jgr)
  !     --------------------------------------------------------

  data ((aws(i,j), i = 1, 4), j = 1, NBS)           / &
       4.554e-04,   1.500e+00,   7.190e-01,  -9.419e-01, &
       3.859e-04,   1.508e+00,   9.512e-01,  -1.053e+00, &
       -3.946e-05,   1.538e+00,   1.035e+00,   2.638e-01, &
       2.936e-04,   1.541e+00,   1.698e-00,   1.521e+00  /
  data ((bws(i,j), i = 1, 4), j = 1, NBS)           / &
       6.481e-08,   1.553e-07,  -7.755e-10,   7.616e-12, &
       1.072e-05,   1.345e-05,  -1.799e-08,  -3.146e-11, &
       4.078e-04,   2.169e-03,  -2.177e-05,   1.506e-07, &
       2.013e-01,   1.109e-02,  -2.897e-04,   3.055e-06  /
  data ((cws(i,j), i = 1, 4), j = 1, NBS)           / &
       8.069e-01,   6.188e-03,  -2.065e-04,   2.352e-06, &
       7.685e-01,   9.337e-03,  -3.101e-04,   3.527e-06, &
       7.471e-01,   9.440e-03,  -2.616e-04,   2.614e-06, &
       7.956e-01,   8.138e-03,  -1.861e-04,   1.611e-06  /

  !     --------------------------------------------------------
  !     new water properties for lw (lindner and li, 2000 jcl)
  !     --------------------------------------------------------

  data ((awl(i,j), i = 1, 5), j = 1, NBL)                      / &
       -.21671e-01, .79578e-03, .14899e+01, .62606e+01,-.12705e+02, &
       -.14126e+00, .28208e-02, .35125e+01,-.34541e+01,-.22679e+01, &
       -.18829e+00, .34065e-02, .46731e+01,-.11664e+02, .87105e+01, &
       -.16383e+00, .26574e-02, .48670e+01,-.16442e+02, .16128e+02, &
       -.20294e-01, .85110e-04, .28650e+01,-.11202e+02, .12047e+02, &
       .28752e-01,-.37315e-03, .14591e+01,-.48788e+01, .49725e+01, &
       -.40386e-01, .80822e-03, .25318e+01,-.64641e+01, .55609e+01, &
       -.48716e-01, .81275e-03, .30390e+01,-.97845e+01, .95101e+01, &
       .64794e-01,-.98530e-03, .12797e+01,-.55272e+01, .62599e+01 /
  data ((bwl(i,j), i = 1, 4), j = 1, NBL)          / &
       .36899e-02,-.54184e-03, .14561e-01,-.18451e-03, &
       .62141e-02, .61190e-01, .21127e-01,-.29731e-03, &
       .87326e-01, .29908e+00, .22928e-01,-.35569e-03, &
       -.37551e-01, .70237e+00, .26945e-01,-.37999e-03, &
       .51671e-01, .10199e+01, .18296e-01,-.21209e-03, &
       .52184e+00, .72352e+00,-.48090e-02, .10414e-03, &
       .57688e+00, .63008e+00,-.56325e-02, .87852e-04, &
       .50346e+00, .79407e+00,-.13179e-02, .25467e-04, &
       .67792e+00, .68259e+00,-.12136e-01, .20941e-03 /
  data ((cwl(i,j), i = 1, 4), j = 1, NBL)          / &
       .73147e+00, .11761e+00, .86402e-02,-.10761e-03, &
       .81284e+00,-.60287e-01, .45367e-02,-.33372e-04, &
       .92468e+00,-.39653e+00, .30494e-03, .20980e-04, &
       .10006e+01,-.71422e+00,-.46784e-02, .10114e-03, &
       .10635e+01,-.10097e+01,-.58726e-02, .97485e-04, &
       .10762e+01,-.12482e+01,-.40343e-02, .54330e-04, &
       .97445e+00,-.13875e+01, .79204e-03,-.27995e-04, &
       .79053e+00,-.13566e+01, .10452e-01,-.18111e-03, &
       .35512e+00,-.80671e+00, .30384e-01,-.47204e-03 /

  !----------------------------------------------------------------------
  !    ice fu 1997 jcl, fu et al. 1998 jcl
  !----------------------------------------------------------------------

  data ((ais(i,j), i = 1, 2), j = 1, NBS) / &
       -0.24276e-04, 2.51884e+00, &
       -0.48500e-04, 2.52275e+00, &
       -0.98503e-05, 2.52048e+00, &
       0.24435e-03, 2.49116e+00                /

  data ((bis(i,j), i = 1, 4), j = 1, NBS)            / &
       0.13031e-06, 0.94102e-07,-0.75971e-10, 0.33977e-12, &
       -0.77603e-06, 0.73420e-05, 0.11514e-09,-0.90818e-12, &
       0.10007e-02, 0.10992e-02,-0.45043e-05, 0.12637e-07, &
       0.21201e+00, 0.25713e-02,-0.19228e-04, 0.62183e-07 /

  data ((cis(i,j), i = 1, 4), j = 1, NBS)            / &
       0.74821e+00, 0.92318e-03,-0.72862e-06,-0.95642e-08, &
       0.75227e+00, 0.10653e-02,-0.24930e-05,-0.29114e-08, &
       0.75553e+00, 0.17297e-02,-0.87585e-05, 0.19201e-07, &
       0.84323e+00, 0.20925e-02,-0.18302e-04, 0.60381e-07 /

  data ((ail(i,j), i = 1, 3), j = 1, NBL)    / &
       -.8839455e-03,  .2662598e+01,  .2196338e+01, &
       -.2066995e-02,  .2787904e+01,  .1397838e+01, &
       -.3085730e-02,  .2906257e+01, -.1911363e+01, &
       -.6968920e-02,  .3284275e+01, -.6973825e+01, &
       -.8372696e-02,  .3455018e+01, -.1516692e+02, &
       -.1691632e-02,  .2765756e+01, -.8331033e+01, &
       -.7098616e-02,  .3343404e+01, -.8144649e+01, &
       -.1041746e-01,  .3824226e+01, -.2255945e+02, &
       .5689700e-02,  .2285636e+01, -.1430752e+02 /

  data ((bil(i,j), i = 1, 4), j = 1, NBL)                   / &
       .5723611e+00,  .1627863e-01, -.1684272e-03,  .6061332e-06, &
       .4402328e+00,  .1736939e-01, -.1656608e-03,  .5709622e-06, &
       .8802908e+00,  .1249744e-01, -.1550609e-03,  .6105065e-06, &
       .6351165e+00,  .1781519e-01, -.1979682e-03,  .6000892e-06, &
       .5409536e-00,  .1949649e-01, -.2050908e-03,  .7364680e-06, &
       .1195515e+01,  .3350616e-02, -.5266996e-04,  .2233377e-06, &
       .1186334e+01,  .6213290e-02, -.1044277e-03,  .2233377e-06, &
       .2279562e+00,  .2017007e-01, -.1756872e-03,  .5703918e-06, &
       .7718967e+00,  .2120626e-01, -.2587649e-03,  .9878070e-06 /

  data ((cil(i,j), i = 1, 4), j = 1, NBL)                   / &
       .7975757e+00,  .3843973e-02, -.3540463e-04,  .1179791e-06, &
       .7947997e+00,  .3190423e-02, -.2386042e-04,  .6691811e-07, &
       .8737279e+00,  .2465886e-02, -.2468764e-04,  .8686448e-07, &
       .8577221e+00,  .2321034e-02, -.1897764e-04,  .8641223e-07, &
       .8906280e-00,  .1903269e-02, -.1733552e-04,  .5855071e-07, &
       .8663385e-00,  .2797934e-02, -.3187011e-04,  .1217209e-06, &
       .7644037e+00,  .4427001e-02, -.4494615e-04,  .1217209e-06, &
       .7200100e+00,  .3220301e-02, -.2195542e-04,  .6604318e-07, &
       .5355918e+00,  .1127081e-01, -.1234705e-03,  .4567953e-06 /

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cldop_propli(tauli, omli, gli, icewpin, rei, j, is_cloud, rlim)
    implicit none

    ! Input arguments
    real, intent(in) :: icewpin                         !In-cloud ice water path
    real, intent(in) :: rei                             !Effective radius of ice crystals (um)
    integer, intent(in) :: j                            !Radiative band
    logical, intent(in) :: is_cloud                     !Threshold for calculations
    real, dimension(2), intent(in), optional :: rlim    !Radius limits (/min, max/) (um)

    ! Output arguments
    real, intent(out) :: tauli                          !Longwave transmittance
    real, intent(out) :: omli                           !
    real, intent(out) :: gli                            !

    ! Internal variables
    real :: dg, dg2, dg3

    ! No-cloud case
    if (.not. is_cloud) then
       tauli = 0.
       omli = 0.
       gli = 0.
       return
    endif

    ! Compute radiatve properties, noting that because in fu etc. the param
    ! is for absorptance, a factor icewpin / tauli is needed for single
    ! scattering albedo
    dg   =  1.5396 * rei
    if (present(rlim)) dg = min(max(dg, rlim(1)), rlim(2))
    dg2  =  dg  * dg
    dg3  =  dg2 * dg
    tauli =  icewpin * (ail(1,j) + ail(2,j) / dg + &
         ail(3,j) / dg2)
    omli  =  1.0 - (bil(1,j) / dg + bil(2,j) + &
         bil(3,j) * dg + bil(4,j) * dg2) * &
         icewpin / tauli
    gli   =  cil(1,j) + cil(2,j) * dg + cil(3,j) * dg2 + &
         cil(4,j) * dg3

    ! End of subprogram
    return
  end subroutine cldop_propli

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cldop_proplw(taulw, omlw, glw, liqwpin, rew, j, is_cloud)
    use phy_status, only: physeterror
    implicit none

    ! Input arguments
    real, intent(in) :: liqwpin                         !In-cloud liquid water path
    real, intent(in) :: rew                             !Effective radius of water droplets (um)
    integer, intent(in) :: j                            !Radiative band
    logical, intent(in) :: is_cloud                     !Threshold for calculations

    ! Output arguments
    real, intent(out) :: taulw                          !Transmittance
    real, intent(out) :: omlw                           !
    real, intent(out) :: glw                            !

    ! Internal variables
    real :: rew2, rew3

    ! No-cloud case
    if (.not. is_cloud) then
       taulw = 0.
       omlw = 0.
       glw = 0.
       return
    endif

    ! Compute radiatve properties
    rew2 =  rew * rew
    rew3 =  rew2 * rew
    taulw =  liqwpin * (awl(1,j) + awl(2,j) * rew+ &
         awl(3,j) / rew + awl(4,j) / rew2 + &
         awl(5,j) / rew3)
    if (taulw < 0.) call physeterror('cldop_proplw', 'Found negative taulw value')
    omlw  =  1.0 - (bwl(1,j) + bwl(2,j) / rew + &
         bwl(3,j) * rew + bwl(4,j) * rew2)
    glw   =  cwl(1,j) + cwl(2,j) / rew + &
         cwl(3,j) * rew + cwl(4,j) * rew2

    ! End of subprogram
    return
  end subroutine cldop_proplw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cldop_propsi(tausi, omsi, gsi, icewpin, rei, j, is_cloud, rlim)
    implicit none

    ! Input arguments
    real, intent(in) :: icewpin                         !In-cloud ice water path
    real, intent(in) :: rei                             !Effective radius of ice crystals (um)
    integer, intent(in) :: j                            !Radiative band
    logical, intent(in) :: is_cloud                     !Threshold for calculations
    real, dimension(2), intent(in), optional :: rlim    !Radius limits (/min, max/) (um)

    ! Output arguments
    real, intent(out) :: tausi                          !Longwave transmittance
    real, intent(out) :: omsi                           !
    real, intent(out) :: gsi                            !

    ! Internal variables
    real :: dg, dg2, dg3

    ! No-cloud case
    if (.not. is_cloud) then
       tausi = 0.
       omsi = 0.
       gsi = 0.
       return
    endif

    ! Compute radiatve properties, noting that because in fu etc. the param
    ! is for absorptance, a factor icewpin / tausi is needed for single
    ! scattering albedo
    dg   =  1.5396 * rei
    if (present(rlim)) dg = min(max(dg, rlim(1)), rlim(2))
    dg2  =  dg  * dg
    dg3  =  dg2 * dg
    tausi =  icewpin * ( ais(1,j) + ais(2,j) / dg )
    omsi  =  1.0 - (bis(1,j) + bis(2,j) * dg + &
         bis(3,j) * dg2 + bis(4,j) * dg3)
    gsi   =  cis(1,j) + cis(2,j) * dg + cis(3,j) * dg2 + &
         cis(4,j) * dg3

    ! End of subprogram
    return
  end subroutine cldop_propsi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cldop_propsw(tausw, omsw, gsw, liqwpin, rew, j, is_cloud)
    implicit none

    ! Input arguments
    real, intent(in) :: liqwpin                         !In-cloud liquid water path
    real, intent(in) :: rew                             !Effective radius of water droplets (um)
    integer, intent(in) :: j                            !Radiative band
    logical, intent(in) :: is_cloud                     !Threshold for calculations

    ! Output arguments
    real, intent(out) :: tausw                          !Transmittance
    real, intent(out) :: omsw                           !
    real, intent(out) :: gsw                            !

    ! Internal variables
    real :: rew2, rew3

    ! No-cloud case
    if (.not. is_cloud) then
       tausw = 0.
       omsw = 0.
       gsw = 0.
       return
    endif

    ! Compute radiatve properties
    rew2 =  rew * rew
    rew3 =  rew2 * rew
    tausw =  liqwpin * &
         (aws(1,j) + aws(2,j) / rew + &
         aws(3,j) / rew2 + aws(4,j) / rew3)
    omsw  =  1.0 - (bws(1,j) + bws(2,j) * rew + &
         bws(3,j) * rew2 + bws(4,j) * rew3)
    gsw   =  cws(1,j) + cws(2,j) * rew + &
         cws(3,j) * rew2 + cws(4,j) * rew3

    ! End of subprogram
    return
  end subroutine cldop_propsw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function cldop_rew(tt, liqwcin, sigt, ps, mg, ml, rew_const, formula, ni, nkm1, remp) result(rew)
    use tdpack_const, only: RGASD
    use phy_status, only: physeterror
    implicit none

    ! Input arguments
    integer, intent(in) :: ni                           !Horizontal dimension
    integer, intent(in) :: nkm1                         !Vertical dimension
    real, dimension(ni,nkm1), intent(in) :: tt          !Dry air temperature (K)  
    real, dimension(ni,nkm1), intent(in) :: liqwcin     !In-cloud liquid water content (kg/kg)
    real, dimension(ni,nkm1), intent(in) :: sigt        !Sigma values of thermo levels
    real, dimension(ni), intent(in) :: ps               !Surface pressure (Pa)
    real, dimension(ni), intent(in) :: mg               !Land-sea mask
    real, dimension(ni), intent(in) :: ml               !
    real, intent(in) :: rew_const                       !User-specified constant for rew (um)
    character(len=*), intent(in) :: formula             !Formula for calcluation
    real, dimension(ni,nkm1), intent(in), optional :: remp  !Effective radius from microphysics (m)

    ! Output arguments
    real, dimension(ni,nkm1) :: rew                     !Effective radius of water droplets (um)

    ! Internal parameters
    real, parameter :: THIRD = 1./3.

    ! Internal variables
    integer :: i, k
    real :: epsilon, epsilon2, betad, betan
    real, dimension(ni,nkm1) :: rec_cdd, aird, vs1

    ! Precompute formula inputs
    if (any(formula == (/'BARKER    ', 'NEWRAD    ', 'ROTSTAYN03'/))) then
       do k = 1, nkm1
          do i = 1, ni
             ! Set cloud droplet concentration per cm^3 to 100 over oceans and 500 over land
             if (mg(i) <= 0.5 .and. ml(i) <= 0.5) then
                rec_cdd(i,k) = 0.01
             else
                rec_cdd(i,k) = 0.002
             endif
             aird(i,k) = sigt(i,k) * ps(i) / ( tt(i,k) * RGASD )  !aird is air density in kg/m3
          end do
       end do
    endif

    select case (formula)
    case ('BARKER')
       ! Radius as in newrad: from H. Barker based on aircraft data (range 4-17um from Slingo)
       rew(:,:) = min(max(4., 754.6 * (liqwcin*aird*rec_cdd)**THIRD), 17.0)
    case ('NEWRAD')
       ! Radius as in newrad: corresponds to so called new optical properties
       vs1(:,:) = (1.0 + liqwcin(:,:) * 1.e4) &
            * liqwcin(:,:) * aird(:,:) * rec_cdd(:,:)
       rew(:,:) =  min(max(2.5, 3000. * vs1**THIRD), 50.0)
    case ('ROTSTAYN03')
       ! Radius according to Rotstayn and Liu (2003)
       do k = 1, nkm1
          do i = 1, ni
             epsilon =  1.0 - 0.7 * exp(- 0.001 / rec_cdd(i,k))
             epsilon2 =  epsilon * epsilon
             betad =  1.0 + epsilon2
             betan =  betad + epsilon2
             rew(i,k) = 620.3504944*((betan*betan*liqwcin(i,k)*aird(i,k)) &
                  / (betad / rec_cdd(i,k)) )**THIRD
             rew(i,k) =  min (max (2.5, rew(i,k)), 17.0)
          end do
       end do
    case ('MP')
       if (present(remp)) then
          rew(:,:) = remp(:,:) * 1.e6
       else
          rew = rew_const  !# make sure rew has value not crashing model
          call physeterror('cldop_rew', 'Microphysics radius estimate missing')
       endif
    case DEFAULT
       ! Radius is a user-specified constant (in microns)
       rew = rew_const
    end select

    ! End of subprogram
    return
  end function cldop_rew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function cldop_rei(tt, icewcin, sigt, ps, lat, rei_const, formula, ni, nkm1, remp) result(rei)
    use tdpack_const, only: RGASD, TCDK
    use phy_status, only: physeterror
    implicit none

    ! Input arguments
    integer, intent(in) :: ni                           !Horizontal dimension
    integer, intent(in) :: nkm1                         !Vertical dimension
    real, dimension(ni,nkm1), intent(in) :: tt          !Dry air temperature (K)  
    real, dimension(ni,nkm1), intent(in) :: icewcin     !In-cloud liquid water content (kg/kg)
    real, dimension(ni,nkm1), intent(in) :: sigt        !Sigma values of thermo levels
    real, dimension(ni), intent(in) :: ps               !Surface pressure (Pa)
    real, dimension(ni), intent(in) :: lat              !Latitude (rad)
    real, intent(in) :: rei_const                       !User-specified constant rei (um)
    character(len=*), intent(in) :: formula             !Formula for calcluation
    real, dimension(ni,nkm1), intent(in), optional :: remp  !Effective radius from microphysics (m)

    ! Output arguments
    real, dimension(ni,nkm1) :: rei                         !Effective radius of ice crystals (um)

    ! Internal variables
    integer :: i, k
    real, dimension(ni, nkm1) :: aird

    ! Precompute formula inputs
    if (any(formula == (/'CCCMA', 'ECMWF'/))) then
       do k = 1, nkm1
          do i = 1, ni
             aird(i,k) = sigt(i,k) * ps(i) / ( tt(i,k) * RGASD )  !aird is air density in kg/m3
          end do
       end do
    endif
    
    ! Effective radius of crystals in ice clouds
    select case (formula)
    case ('CCCMA')
       ! Units of icewcin must be in g/m3 for this parameterization of rei (in microns)
       rei(:,:) = (1000. * icewcin * aird)**0.216
       where (icewcin(:,:) >= 1.e-9)
          rei(:,:) = 83.8 * rei(:,:)
       elsewhere
          rei(:,:) = 20.
       endwhere
       rei(:,:) =  max(min(rei(:,:), 50.0), 20.0)
    case ('SIGMA')
       ! Radius varies from 60um (near-surface) to 15um (upper-troposphere)
       rei(:,:) = max(sigt(:,:)-0.25, 0.0)*60. + 15.
    case ('ECMWF')
       ! see IFS documentation for Cy47R3 -eqns 2.74 and 2.75 - beware of parenthesis error for first term
       do k = 1, nkm1
          do i = 1, ni
             rei(i,k) = 1000. * icewcin(i,k) * aird(i,k) ! convert to gm-3
             rei(i,k) = (1.2351 + 0.0105*(tt(i,k) - TCDK)) * (45.8966*rei(i,k)**0.2214 + 0.7957*rei(i,k)**0.2535*(tt(i,k) - 83.15))
             rei(i,k) =  max(min(rei(i,k), 155.0), (20.+40.*cos(lat(i)))) ! impose a lat dependent min
             rei(i,k) = 0.64952*rei(i,k)
             rei(i,k) =  min(rei(i,k), 70.0) ! necessary to avoid crashes
          enddo
       enddo
    case ('MP')
       ! Use what microphy provides, some filters were applied before and are applied below !!!!coded only for mpcat=1
       if (present(remp)) then
          rei(:,:) = remp(:,:) * 1.e6
       else
          rei(:,:) = rei_const  !# make sure rew has value not crashing model
          call physeterror('cldop_rei', 'Microphysics radius estimate missing')
       endif
    case DEFAULT
       ! Radius is a user-specified constant (in microns)
       rei(:,:) = rei_const
    end select

    ! End of subprogram
    return
  end function cldop_rei

end module cldop_utils
