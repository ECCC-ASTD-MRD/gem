!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------
!
!!if_on
subroutine chm_load_metvar(busdyn, busper, busvol, metvar2d, metvar3d)
  use chm_metvar_mod
  use chm_ptopo_grid_mod, only: chm_ni, chm_nk
!!if_off
  use chm_utils_mod,      only: ik, chm_lun_out, global_debug
  use chm_consphychm_mod, only: grav, delta, rgasd, pi
  use chm_nml_mod,        only: chm_pblh_min_l, chm_indirect_l, chm_cffeps_online_l
  use chm_phyvar_mod
! From rpnphy
  use sfcbus_mod,         only: indx_agrege
  use phy_options,        only: conv_mid

  implicit none
!!if_on
 real(kind=4),    dimension(:), pointer, contiguous :: busdyn
 real(kind=4),    dimension(:), pointer, contiguous :: busper
 real(kind=4),    dimension(:), pointer, contiguous :: busvol
 real(kind=4),    intent(out) :: metvar2d(chm_ni, SIZE_MV2D)
 real(kind=4),    intent(out) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
!!if_off
!
!  Local variables
!
  integer(kind=4)         :: i, k, this_ik, this_ikp, this_ikm
!  integer :: iomp_alloc_lock
  logical(kind=4)         :: local_dbg
  real(kind=4)            :: p_kfc_m, p_kfcl_m, p_kfcs_m, p_kfc, p_kfcprod, &
                             p_kfcevap, p_cs_p, p_cs, p_csprod, p_csevap,   &
                             p_kfml_m, p_kfms_m, wdir
  real(kind=4), parameter :: smf = 1.0e-15, rad2deg = 180.0 / pi
!
!  BEGIN CODE
!
  local_dbg = ((.false. .or. global_debug) .and. ( chm_lun_out > 0 ))

  if (local_dbg) then
     write(chm_lun_out, *) "Entre dans chm_load_metvar "
  endif

  if (local_dbg) then
     write(chm_lun_out, *) "chm_ni, chm_nk        : " , chm_ni, chm_nk
     write(chm_lun_out, *) "============ 2D METVAR ================"
     write(chm_lun_out, *) "---------------------------------------"
     write(chm_lun_out, *) "SIZE_MV2D                  : " , SIZE_MV2D
     write(chm_lun_out, *) "MV2D_PPLUS      , p0_plus  : " , MV2D_PPLUS   , p0_plus
     write(chm_lun_out, *) "MV2D_DXDY       , dxdy     : " , MV2D_DXDY    , dxdy
     write(chm_lun_out, *) "MV2D_TSURF      , tsurf    : " , MV2D_TSURF   , tsurf
     write(chm_lun_out, *) "MV2D_WSDIAG     , u/vdiag  : " , MV2D_WSDIAG  , udiag, " + ", vdiag
     write(chm_lun_out, *) "MV2D_TDIAG      , tdiag    : " , MV2D_TDIAG   , tdiag
     write(chm_lun_out, *) "MV2D_QDIAG      , qdiag    : " , MV2D_QDIAG   , qdiag
     write(chm_lun_out, *) "MV2D_GLSEA      , glsea    : " , MV2D_GLSEA   , glsea
     write(chm_lun_out, *) "MV2D_SNODP      , snodp    : " , MV2D_SNODP   , snodp
     write(chm_lun_out, *) "MV2D_H          , h        : " , MV2D_H       , h
     write(chm_lun_out, *) "MV2D_DLAT       , dlat     : " , MV2D_DLAT    , dlat
     write(chm_lun_out, *) "MV2D_DLON       , dlon     : " , MV2D_DLON    , dlon
     write(chm_lun_out, *) "MV2D_FLUSOLIS   , flusolis : " , MV2D_FLUSOLIS, flusolis
     write(chm_lun_out, *) "MV2D_MT         , me_moins : " , MV2D_MT      , me_moins
     write(chm_lun_out, *) "MV2D_ILMO       , ilmo     : " , MV2D_ILMO    , ilmo
     write(chm_lun_out, *) "MV2D_WSOIL      , wsoil    : " , MV2D_WSOIL   , wsoil
     write(chm_lun_out, *) "MV2D_UE         , ue       : " , MV2D_UE      , ue
     write(chm_lun_out, *) "MV2D_CANG       , cang     : " , MV2D_CANG    , cang
     write(chm_lun_out, *) "MV2D_RAINRATE   , rainrate : " , MV2D_RAINRATE, rainrate
     write(chm_lun_out, *) "MV2D_SNOF       , psn      : " , MV2D_SNOF    , psn
     write(chm_lun_out, *) "MV2D_MG         , mg       : " , MV2D_MG      , mg
     write(chm_lun_out, *) "MV2D_AL5        , alvis    : " , MV2D_AL5     , alvis
     write(chm_lun_out, *) "---------------------------------------"
     write(chm_lun_out, *) "============ 3D METVAR ================"
     write(chm_lun_out, *) "---------------------------------------"
     write(chm_lun_out, *) "SIZE_MV3D                  : " , SIZE_MV3D
     write(chm_lun_out, *) "MV3D_TPLUS      , tplus    : " , MV3D_TPLUS   , tplus
     write(chm_lun_out, *) "MV3D_WS         , u/vplus  : " , MV3D_WS      , uplus, " + ",vplus
     write(chm_lun_out, *) "MV3D_HUPLUS     , huplus   : " , MV3D_HUPLUS  , huplus
     write(chm_lun_out, *) "MV3D_QCPLUS     , qcplus   : " , MV3D_QCPLUS  , qcplus
     write(chm_lun_out, *) "MV3D_SIGM       , sigm     : " , MV3D_SIGM    , sigm
     write(chm_lun_out, *) "MV3D_SIGT       , sigt     : " , MV3D_SIGT    , sigt
     write(chm_lun_out, *) "MV3D_WPLUS      , wplus    : " , MV3D_WPLUS   , wplus
     write(chm_lun_out, *) "MV3D_FTOT       , ftot     : " , MV3D_FTOT    , ftot
     write(chm_lun_out, *) "MV3D_ZMOM       , gzmom    : " , MV3D_ZMOM    , gzmom
     write(chm_lun_out, *) "MV3D_ZPLUS      , gztherm  : " , MV3D_ZPLUS   , gztherm
     write(chm_lun_out, *) "MV3D_KT         , kt       : " , MV3D_KT      , kt
     write(chm_lun_out, *) "MV3D_RNFLX      , rnflx    : " , MV3D_RNFLX   , rnflx , " + ",kfcrf
     write(chm_lun_out, *) "MV3D_SNOFLX     , snoflx   : " , MV3D_SNOFLX  , snoflx  , " + ",kfcsf
     write(chm_lun_out, *) "MV3D_PEVP       , pevp     : " , MV3D_PEVP    , fevp , " + eqn. "
     write(chm_lun_out, *) "MV3D_PPRO       , ppro     : " , MV3D_PPRO    , f12  , " + ", qrkfc
     write(chm_lun_out, *) "MV3D_LWC        , lwc      : " , MV3D_LWC     , lwc
     write(chm_lun_out, *) "MV3D_CLDRAD     , cldrad   : " , MV3D_CLDRAD  , cldrad
     write(chm_lun_out, *) "MV3D_NCPLUS     , ncplus   : " , MV3D_NCPLUS  , ncplus
     write(chm_lun_out, *) "MV3D_O3L        , o3lplus  : " , MV3D_O3L     , o3lplus
     write(chm_lun_out, *) "---------------------------------------"
  endif
!
! Load the 2d arrays
!
!     Dynamic bus
!
  do i = 1, chm_ni
     metvar2d(i,MV2D_MT)       = busdyn(me_moins + i - 1) / grav
     metvar2d(i,MV2D_PPLUS)    = busdyn(p0_plus  + i - 1)
!
!     Permanent bus
!
     metvar2d(i,MV2D_DLAT)     = busper(dlat     + i - 1)
     metvar2d(i,MV2D_DLON)     = busper(dlon     + i - 1)
     metvar2d(i,MV2D_DXDY)     = busper(dxdy     + i - 1)
     metvar2d(i,MV2D_FLUSOLIS) = busper(flusolis + i - 1)
     metvar2d(i,MV2D_GLSEA)    = busper(glsea    + i - 1)
     metvar2d(i,MV2D_QDIAG)    = busper(qdiag    + i - 1)
     metvar2d(i,MV2D_TDIAG)    = busper(tdiag    + i - 1)

     metvar2d(i,MV2D_WSDIAG)   = sqrt(busper(udiag + i - 1)**2 + &
                                      busper(vdiag + i - 1)**2)

     metvar2d(i,MV2D_WSOIL)    = busper(wsoil    + i - 1)
     metvar2d(i,MV2D_MG)       = busper(mg       + i - 1)
! impose minimum to the PBL height
     if (chm_pblh_min_l) then
        metvar2d(i,MV2D_H)     = max(100.,busper(h + i - 1))
     else
        metvar2d(i,MV2D_H)     = busper(h        + i - 1)
     end if
! aggregated Inverse of Monin-Obukhov, snow depth, and surface temperature
     this_ik = ik(i, indx_agrege, chm_ni)
     metvar2d(i,MV2D_ILMO)     = busper(ilmo  + this_ik)
     metvar2d(i,MV2D_SNODP)    = busper(snodp + this_ik)
     metvar2d(i,MV2D_TSURF)    = busper(tsurf + this_ik)
     metvar2d(i,MV2D_AL5)      = busper(alvis + this_ik)
!
!     Volatile bus
!
     metvar2d(i,MV2D_CANG)     = amax1(-1.0, amin1(1.0, busvol(cang  + i - 1)))
     metvar2d(i,MV2D_RAINRATE) = busvol(rainrate + i - 1)
     metvar2d(i,MV2D_SNOF)     = busvol(psn      + i - 1)
     metvar2d(i,MV2D_UE)       = busvol(ue       + i - 1)
  end do
  if (chm_cffeps_online_l) then
     do i = 1, chm_ni
        if (metvar2d(i,MV2D_WSDIAG) == 0.0) then
           wdir = 0.0
        else
           if (busper(udiag + i - 1) == 0.0) then
              if (busper(vdiag + i - 1) >= 0.0) then
                 wdir = rad2deg * metvar2d(i,MV2D_DLON) - 90.0
              else
                 wdir = rad2deg * metvar2d(i,MV2D_DLON) + 90.0
              end if
           else
              wdir = rad2deg * (metvar2d(i,MV2D_DLON) - &
                     atan(busper(vdiag + i - 1), busper(udiag + i - 1)))
           end if
        end if
        metvar2d(i,MV2D_WDDIAG) = amod(amod(wdir, 360.0) + 360.0, 360.0)
     end do
  else
     metvar2d(:, MV2D_WDDIAG) = 0.0
  end if
!
! Load the 3D arrays
!
  do k = 1, chm_nk
     do i = 1, chm_ni
        this_ik = ik(i, k, chm_ni)
!
!           Dynamic bus
!
        metvar3d(i, k, MV3D_HUPLUS)   = max(0.,busdyn(huplus    + this_ik))
        metvar3d(i, k, MV3D_QCPLUS)   = max(0.,busdyn(qcplus    + this_ik))
        metvar3d(i, k, MV3D_SIGM)     = busdyn(sigm      + this_ik)
        metvar3d(i, k, MV3D_SIGT)     = busdyn(sigt      + this_ik)
        metvar3d(i, k, MV3D_TPLUS)    = busdyn(tplus     + this_ik)

        metvar3d(i, k, MV3D_WS)       = sqrt(busdyn(uplus + this_ik)**2 + &
                                             busdyn(vplus + this_ik)**2)

        metvar3d(i, k, MV3D_WPLUS)    = busdyn(wplus     + this_ik)

        if (chm_indirect_l) then
           metvar3d(i, k, MV3D_NCPLUS)= busdyn(ncplus     + this_ik)
        else
           metvar3d(i, k, MV3D_NCPLUS)= 0.0
        end if
        if (o3lplus > 0) then
           metvar3d(i, k, MV3D_O3L)   = busdyn(o3lplus   + this_ik)
        else
           metvar3d(i, k, MV3D_O3L)   = 0.0
        end if
!
!           Permanent bus
!
        metvar3d(i, k, MV3D_FTOT)     = busper(ftot      + this_ik)
        metvar3d(i, k, MV3D_LWC)      = busper(lwc       + this_ik)
!
!           Volatile bus
!
        metvar3d(i, k, MV3D_ZMOM)     = busvol(gzmom     + this_ik)
        metvar3d(i, k, MV3D_ZPLUS)    = busvol(gztherm   + this_ik)
        metvar3d(i, k, MV3D_KT)       = busvol(kt        + this_ik)
        metvar3d(i, k, MV3D_CLDRAD)   = busvol(cldrad    + this_ik)
!
!           Mixed sources


        if (qrkfc > 0 ) then                ! KFC used

!       Condensation fluxes, production and evaporation (CONSUN or MY) to be combined with those of KFC
          this_ikp = ik(i, k+1, chm_ni)
!
          p_cs_p = busper(rnflx + this_ikp) + busper(snoflx + this_ikp)
          p_cs   = busper(rnflx + this_ik ) + busper(snoflx + this_ik )
          p_csevap = p_cs
          if (busvol(fevp + this_ik) < 1.) p_csevap =  busvol(fevp + this_ik) * p_cs_p /(1.- busvol(fevp + this_ik))
          p_csprod = p_cs_p - p_cs + p_csevap
!
!         KFC fluxes, production and evaporation
          this_ikm = ik(i, k-1, chm_ni)
!
          p_kfc_m =  0.0
          p_kfcl_m =  0.0
          p_kfcs_m =  0.0
          p_kfml_m = 0.0
          p_kfms_m = 0.0
          if (k > 1) then
             if (conv_mid /= 'NIL') then
                p_kfml_m = busvol(kfmrf + this_ikm)
                p_kfms_m = busvol(kfmsf + this_ikm)
             end if
             p_kfcl_m = busper(kfcrf + this_ikm) + p_kfml_m
             p_kfcs_m = busper(kfcsf + this_ikm) + p_kfms_m
             p_kfc_m  = p_kfcl_m + p_kfcs_m
          end if
          p_kfc   =  busper(kfcrf + this_ik ) + busper(kfcsf + this_ik )
          if (conv_mid /= 'NIL') then
             p_kfc = p_kfc + busvol(kfmrf + this_ik) + busvol(kfmsf + this_ik)
          end if
          p_kfcevap =  max(0.,p_kfc_m - p_kfc)
          p_kfcprod =  p_kfc - p_kfc_m + p_kfcevap
          metvar3d(i, k, MV3D_RNFLX)    = busper(rnflx     + this_ik) + p_kfcl_m
          metvar3d(i, k, MV3D_SNOFLX)   = busper(snoflx    + this_ik) + p_kfcs_m
          metvar3d(i, k, MV3D_PEVP)     = min(1.,max(0., (p_csevap + p_kfcevap) / &
                                          (p_cs + p_csprod + p_kfc_m + p_kfcprod + smf)))
          metvar3d(i, k, MV3D_PPRO)     = max(0., (busvol(f12 + this_ik) + busper(qrkfc + this_ik)))
        else
          metvar3d(i, k, MV3D_RNFLX)    = busper(rnflx     + this_ik)
          metvar3d(i, k, MV3D_SNOFLX)   = busper(snoflx    + this_ik)
          metvar3d(i, k, MV3D_PEVP)     = busvol(fevp      + this_ik)
          metvar3d(i, k, MV3D_PPRO)     = max(0., busvol(f12 + this_ik))
        end if
!
!           Evaluated fields
!
! Air density in kg/m3
        metvar3d(i, k, MV3D_RHO)      = metvar3d(i, k, MV3D_SIGT) * metvar2d(i, MV2D_PPLUS) / &
                                        (rgasd * metvar3d(i, k, MV3D_TPLUS) *                 &
                                        (1.0 + delta * metvar3d(i, k, MV3D_HUPLUS)))

     end do
  end do

  if (local_dbg) then
     write(chm_lun_out, *) "============ END METVAR PRINT ========="
     write(chm_lun_out, *) " And now some values from metvar2d: "
     write(chm_lun_out, *) " metvar2d(1,1),metvar2d(chm_ni,SIZE_MV2D)         :  ", metvar2d(1,1),metvar2d(chm_ni,SIZE_MV2D)
     write(chm_lun_out, *) " And now some values from metvar3d: "
     write(chm_lun_out, *) " metvar3d(1,1,1),metvar3d(chm_ni,chm_nk,SIZE_MV3D):  ", metvar3d(1,1,1),metvar3d(chm_ni,chm_nk,SIZE_MV3D)
     write(chm_lun_out, *) "EXIT chm_load_metvar"
     write(chm_lun_out, *) " "
  end if

  return
end subroutine chm_load_metvar

