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
subroutine chm_load_metvar(pvars, metvar2d, metvar3d)
  use phymem,             only: phyvar
  use chm_metvar_mod
  use chm_ptopo_grid_mod, only: chm_ni, chm_nk
!!if_off
  use chm_utils_mod,      only: ik, chm_lun_out, chm_msg_debug, chm_error_l
  use chm_consphychm_mod, only: grav, delta, rgasd, rad2deg, tcdk
  use chm_nml_mod,        only: chm_pblh_min_l, chm_indirect_l, &
                                chm_cffeps_online_l
  use chm_species_info_mod, only: print_phymeta
  use chm_phyvar_mod
! From rpnphy
  use phymem,             only: phymeta, phymem_find, PHY_NAMELEN
  use sfcbus_mod,         only: indx_agrege
  use phy_options,        only: conv_mid
  implicit none
!!if_on
  type(phyvar),    pointer, contiguous :: pvars(:)
  real(kind=4),    intent(out) :: metvar2d(chm_ni, SIZE_MV2D)
  real(kind=4),    intent(out) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
!!if_off
!
! Local variables
!
  integer(kind=4)         :: i, k, this_ik, this_ikp, this_ikm, istat, idxv1(1)
  real(kind=4)            :: p_kfc_m, p_kfcl_m, p_kfcs_m, p_kfc, p_kfcprod, &
                             p_kfcevap, p_cs_p, p_cs, p_csprod, p_csevap,   &
                             p_kfml_m, p_kfms_m, wdir
  real(kind=4), parameter :: smf = 1.0e-15

! for debug
  logical(kind=4)         :: local_dbg
  integer, parameter      :: nvnamedbg=5
  character(len=PHY_NAMELEN), parameter:: vnamedbg(nvnamedbg) = &
                              (/'dlat','wsoil','snodp','ftot','vegf'/)
! Declaration of external subroutines
  external :: msg_toall
!
! BEGIN CODE
!
  local_dbg = (.false. .and. (chm_lun_out > 0))
  call msg_toall(chm_msg_debug, 'chm_load_metvar [BEGIN]')

! Load the 2d variable with vmeta%fmul=1 (or vmeta%size=vmeta%ni)
  metvar2d(:,MV2D_DLAT)     = pvars(dlat)%data
  metvar2d(:,MV2D_DLON)     = pvars(dlon)%data
  metvar2d(:,MV2D_DXDY)     = pvars(dxdy)%data
  metvar2d(:,MV2D_MT)       = pvars(me_moins)%data(:) / grav
  metvar2d(:,MV2D_PPLUS)    = pvars(p0_plus)%data
  metvar2d(:,MV2D_DLAT)     = pvars(dlat)%data
  metvar2d(:,MV2D_DLON)     = pvars(dlon)%data
  metvar2d(:,MV2D_DXDY)     = pvars(dxdy)%data
  metvar2d(:,MV2D_FLUSOLIS) = pvars(flusolis)%data
  metvar2d(:,MV2D_GLSEA)    = pvars(glsea)%data
  metvar2d(:,MV2D_QDIAG)    = pvars(qdiag)%data
  metvar2d(:,MV2D_TDIAG)    = pvars(tdiag)%data
  metvar2d(:,MV2D_WSDIAG)   = sqrt(pvars(udiag)%data(:)**2 + &
                                   pvars(vdiag)%data(:)**2)
  metvar2d(:,MV2D_MG)       = pvars(mg)%data(:)
  metvar2d(:,MV2D_CANG)     = amax1(-1.0, amin1(1.0, pvars(cang)%data(:)))
  metvar2d(:,MV2D_RAINRATE) = pvars(rainrate)%data
  metvar2d(:,MV2D_SNOF)     = pvars(psn)%data
  metvar2d(:,MV2D_UE)       = pvars(ue)%data
  ! Sea-surf temp in K, converted to C
  metvar2d(:,MV2D_TWATER)   = pvars(twater)%data - tcdk

! impose minimum to the PBL height
  if (chm_pblh_min_l) then
     metvar2d(:,MV2D_H)     = max(100.,pvars(h)%data(:))
  else
     metvar2d(:,MV2D_H)     = pvars(h)%data(:)
  end if

! 2d variable with vmeta%fmul=2 -- use only first layer
!! CHECK THIS (from phphyvar.hf)
!IF_ISBA: PHYVAR3D1(wsoil,        'VN=wsoil        ;ON=I1  ;VD=soil volumetric water contents                 ;VS=A*2          ;VB=p1        ;MIN=0')
!IF_SVS:  PHYVAR3D1(wsoil,        'VN=wsoil        ;ON=WSOL;VD=soil volm water content per layer              ;VS=A*'//ngl//'  ;VB=p1')
!         PHYVAR2D1(wsoilm,       'VN=wsoilm       ;ON=WSLM;VD=mean soil volm watr cont for the whole column                   ;VB=p0')
  metvar2d(:,MV2D_WSOIL)    = pvars(wsoil)%data(1:chm_ni)

! variable with vmeta%fmul=5 (use 5th aggregate)
  do i = 1, chm_ni
     this_ik = ik(i, indx_agrege, chm_ni) + 1
     metvar2d(i,MV2D_ILMO)     = pvars(ilmo)%data(this_ik)
     metvar2d(i,MV2D_SNODP)    = pvars(snodp)%data(this_ik)
     metvar2d(i,MV2D_TSURF)    = pvars(tsurf)%data(this_ik)
     metvar2d(i,MV2D_AL5)      = pvars(alvis)%data(this_ik)
  end do

! screen wind direction in deg
  if (chm_cffeps_online_l) then
     do i = 1, chm_ni
        if (metvar2d(i,MV2D_WSDIAG) == 0.0) then
           wdir = 0.0
        else
           if (pvars(udiag)%data(i) == 0.0) then
              if (pvars(vdiag)%data(i) >= 0.0) then
                 wdir = rad2deg * metvar2d(i,MV2D_DLON) - 90.0
              else
                 wdir = rad2deg * metvar2d(i,MV2D_DLON) + 90.0
              end if
           else
              wdir = rad2deg * (metvar2d(i,MV2D_DLON) - &
                     atan(pvars(vdiag)%data(i), pvars(udiag)%data(i)))
           end if
        end if
        metvar2d(i,MV2D_WDDIAG) = amod(amod(wdir, 360.0) + 360.0, 360.0)
     end do
  else
     metvar2d(:, MV2D_WDDIAG) = 0.0
  end if
!
! 3D variable with vmeta%fmul=1 (or vmeta%size=vmeta%ni)
  metvar3d(:,:,MV3D_HUPLUS) = reshape(max(0.,pvars(huplus)%data), (/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_QCPLUS) = reshape(max(0.,pvars(qcplus)%data), (/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_SIGM)   = reshape(pvars(sigm)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_SIGT)   = reshape(pvars(sigt)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_TPLUS)  = reshape(pvars(tplus)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_WS)     = sqrt(reshape(pvars(uplus)%data,(/chm_ni,chm_nk/))**2 + &
                                   reshape(pvars(vplus)%data,(/chm_ni,chm_nk/))**2 )
  metvar3d(:,:,MV3D_WPLUS)  = reshape(pvars(wplus)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_FTOT)   = reshape(pvars(ftot)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_ZMOM)   = reshape(pvars(gzmom)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_ZPLUS)  = reshape(pvars(gztherm)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_KT)     = reshape(pvars(kt)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_LWCRAD) = reshape(pvars(lwcrad)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_CLDRAD) = reshape(pvars(cldrad)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_RNFLX)  = reshape(pvars(rnflx)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_SNOFLX) = reshape(pvars(snoflx)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_PEVP)   = reshape(pvars(fevp)%data,(/chm_ni,chm_nk/))
  metvar3d(:,:,MV3D_PPRO)   = max(0.,reshape(pvars(f12)%data,(/chm_ni,chm_nk/)))

! Account for rain water from the MP schemes in the total cloud water content
  if (qrplus > 0) then
     metvar3d(:,:,MV3D_QCPLUS) = metvar3d(i, k, MV3D_QCPLUS) + &
                          max(0.,reshape(pvars(qrplus)%data,(/chm_ni,chm_nk/)))
  endif
  if (chm_indirect_l) then
     metvar3d(:,:,MV3D_NCPLUS)= reshape(pvars(ncplus)%data,(/chm_ni,chm_nk/))
  else
     metvar3d(:,:,MV3D_NCPLUS)= 0.0
  endif
  if (o3lplus > 0) then
     metvar3d(:,:,MV3D_O3L) = reshape(pvars(o3lplus)%data,(/chm_ni,chm_nk/))
  else
     metvar3d(:,:,MV3D_O3L) = 0.0
  endif

! if KFC used condensation fluxes, production and evaporation (CONSUN or MY) to be combined
  if (qrkfc > 0) then
     do k = 1, chm_nk
        do i = 1, chm_ni
          this_ik = ik(i, k, chm_ni) + 1  ! to match pvars idx start with 1
          this_ikp = ik(i, k+1, chm_ni) + 1 !  match idx with k+1
          this_ikm = ik(i, k-1, chm_ni) + 1 !  match idx with k-1

          p_cs_p = pvars(rnflx)%data(this_ikp) + pvars(snoflx)%data(this_ikp)
          p_cs   = pvars(rnflx)%data(this_ik) + pvars(snoflx)%data(this_ik)
          p_csevap = p_cs
          if (pvars(fevp)%data(this_ik) < 1.) &
              p_csevap = pvars(fevp)%data(this_ik) * &
                         p_cs_p /(1.- pvars(fevp)%data(this_ik))
          p_csprod = p_cs_p - p_cs + p_csevap
!
!         KFC fluxes, production and evaporation
          p_kfc_m =  0.0
          p_kfcl_m = 0.0
          p_kfcs_m = 0.0
          p_kfml_m = 0.0
          p_kfms_m = 0.0
          if (k > 1) then
             if (conv_mid /= 'NIL') then
                p_kfml_m = pvars(kfmrf)%data(this_ikm)
                p_kfms_m = pvars(kfmsf)%data(this_ikm)
             end if
             p_kfcl_m = pvars(kfcrf)%data(this_ikm) + p_kfml_m
             p_kfcs_m = pvars(kfcsf)%data(this_ikm) + p_kfms_m
             p_kfc_m  = p_kfcl_m + p_kfcs_m
          end if
          p_kfc = pvars(kfcrf)%data(this_ik) + pvars(kfcsf)%data(this_ik)
          if (conv_mid /= 'NIL') then
             p_kfc = p_kfc + pvars(kfmrf)%data(this_ik) + pvars(kfmsf)%data(this_ik)
          end if
          p_kfcevap = max(0.,p_kfc_m - p_kfc)
          p_kfcprod = p_kfc - p_kfc_m + p_kfcevap
          metvar3d(i, k, MV3D_RNFLX)  = pvars(rnflx)%data(this_ik) + p_kfcl_m
          metvar3d(i, k, MV3D_SNOFLX) = pvars(snoflx)%data(this_ik) + p_kfcs_m
          metvar3d(i, k, MV3D_PEVP)   = min(1.,max(0., (p_csevap + p_kfcevap) / &
                               (p_cs + p_csprod + p_kfc_m + p_kfcprod + smf)))
          metvar3d(i, k, MV3D_PPRO)   = max(0., (pvars(f12)%data(this_ik) + &
                                        pvars(qrkfc)%data(this_ik)))
        end do
     end do
  end if

! Air density in kg/m3
  do i = 1, chm_ni
  do k = 1, chm_nk
     metvar3d(i, k, MV3D_RHO) = metvar3d(i, k, MV3D_SIGT) * metvar2d(i, MV2D_PPLUS) / &
                                (rgasd * metvar3d(i, k, MV3D_TPLUS) *                 &
                                (1.0 + delta * metvar3d(i, k, MV3D_HUPLUS)))
  enddo
  enddo

!  Debug print
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
     write(chm_lun_out, *) "MV3D_QCPLUS     , qcplus   : " , MV3D_QCPLUS  , qcplus, " + ",qrplus
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
     write(chm_lun_out, *) "MV3D_LWCRAD     , lwcrad   : " , MV3D_LWCRAD  , lwcrad
     write(chm_lun_out, *) "MV3D_CLDRAD     , cldrad   : " , MV3D_CLDRAD  , cldrad
     write(chm_lun_out, *) "MV3D_NCPLUS     , ncplus   : " , MV3D_NCPLUS  , ncplus
     write(chm_lun_out, *) "MV3D_O3L        , o3lplus  : " , MV3D_O3L     , o3lplus
     write(chm_lun_out, *) "---------------------------------------"
     write(chm_lun_out,*) 'Listing pvars element for:',vnamedbg(1:nvnamedbg)
     do i = 1, nvnamedbg
        write(chm_lun_out,*) 'phymem_find: i, vname:',i, vnamedbg(i)
        istat = phymem_find(idxv1, vnamedbg(i), F_npath='VOI', &
                            F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
        if (istat < 0 .or. idxv1(1) < 0) then
           write(chm_lun_out,*) 'phymem_find error istat', istat
           chm_error_l = .true.
           return
        end if
        istat = print_phymeta(idxv1(1), chm_lun_out)
        write(chm_lun_out,*) '****** '
     enddo

     write(chm_lun_out, *) "============"
     write(chm_lun_out, *) " And now some values from metvar2d: "
     write(chm_lun_out, *) " metvar2d(1,1),metvar2d(chm_ni,SIZE_MV2D)         :  ", metvar2d(1,1),metvar2d(chm_ni,SIZE_MV2D)
     write(chm_lun_out, *) " And now some values from metvar3d: "
     write(chm_lun_out, *) " metvar3d(1,1,1),metvar3d(chm_ni,chm_nk,SIZE_MV3D):  ", metvar3d(1,1,1),metvar3d(chm_ni,chm_nk,SIZE_MV3D)
     write(chm_lun_out, *) " ============"
     write(chm_lun_out, *) " METVAR2D(i), MAXVAL,MINVAL:"
     do i = 1, SIZE_MV2D
        write(chm_lun_out, *)"METVAR2D:",i,maxval(metvar2d(:,i)),minval(metvar2d(:,i))
     enddo
     write(chm_lun_out, *) " ============"
     write(chm_lun_out, *) " METVAR3D(i), MAXVAL,MINVAL:"
     do i = 1, SIZE_MV3D
        write(chm_lun_out, *)"METVAR3D:",i,maxval(metvar3d(:,:,i)),minval(metvar3d(:,:,i))
     enddo
     write(chm_lun_out, *) "============ END METVAR PRINT ========="
  end if

  call msg_toall(chm_msg_debug, 'chm_load_metvar [END]')

  return
end subroutine chm_load_metvar

