!---------------------------------- LICENCE BEGIN ------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END --------------------------------

module series_init_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use mu_jdate_mod, only: jdate_to_cmc
   use ptopo_utils, only: ptopo_init_var, ptopo_isblocmaster_L
   use phygridmap, only: ijphy, phydim_nj, drv_glb_ni, drv_glb_nj
   use series_options
   use series_write_mod, only: series_write
   use rmn_gmm
   implicit none
   private

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   public :: series_init

contains

   !/@*
   function series_init(phy_lcl_gid, drv_glb_gid, &
        phydim_nk , phy_lcl_ni, phy_lcl_nj, &
        F_moyhr, F_delt, F_jdate, F_satuco_L, F_master_pe) result(F_istat)
      implicit none
      !@arguments
      !#TODO: never used: phydim_ni, phydim_nj, p_runlgt
      integer, intent(in) :: phy_lcl_gid, drv_glb_gid
      integer, intent(in) :: phydim_nk
      integer, intent(in) :: phy_lcl_ni, phy_lcl_nj
      integer, intent(in) :: F_moyhr
      integer, intent(in) :: F_master_pe
      logical, intent(in) :: F_satuco_L
      real,    intent(in) :: F_delt
      integer(INT64), intent(in) :: F_jdate
      !@returns
      integer :: F_istat
      !@author: ?
      !@revisions
      ! 2017-01, S.Chamberland - Major revision
      !@description
      !      This routine initializes the physics variables
      !      related to time series extraction: variable names to
      !      extract (profil:3D, surface:2D), grid point indicies where
      !      to extract from, number or vertical levels...
      !
      !      It also performs memory allocation for buffers based
      !      on: the number of 2D and 3D variables, the number of
      !      vertical level and the number of grid point where to
      !      extract from.
      !*@/
      integer :: istat, nn, cmcdateo, bidon
      type(gmm_metadata) :: mymeta
      !---------------------------------------------------------------
      !#TODO: protect from multi thread

      series_master_pe = F_master_pe
      series_moyhr     = F_moyhr
      series_delt      = F_delt
      series_jdate     = F_jdate
      series_satuco_L  = F_satuco_L

      F_istat = RMN_OK
      if (series_initok_L .or. .not.series_on_L) return
      F_istat = RMN_ERR

      cmcdateo = jdate_to_cmc(series_jdate)
      istat = newdate(cmcdateo, series_cmcdateo14, bidon, RMN_DATE_STAMP2OLD)
      series_interval_sec = int(series_delt * P_serg_srwri)

      !# Find stn i,j coor (from provided lat,lon) and trim to local stn list
      istat = priv_stnij(phy_lcl_gid, drv_glb_gid, phy_lcl_ni, phy_lcl_nj)
      if (.not.RMN_IS_OK(istat)) return

      !# Add Mandatory P0 for profiles
      if (series_nprof > 0) then
         if (.not.any(p_serg_srsrf_s(1:series_nsurf) == 'P0')) then
            if (series_nsurf >= NVARMAX) then
               call msg(MSG_WARNING, PKGNAME_S//'Too many surface var')
               return
            endif
            call msg(MSG_INFO, PKGNAME_S//'Adding P0 to surface var for profiles')
            series_nsurf = series_nsurf + 1
            p_serg_srsrf_s(series_nsurf) = 'P0'
         endif
      endif

      !# Allocate buffers
      call ptopo_init_var()
 
      !#TODO: review index order
      mymeta = GMM_NULL_METADATA
      nn = max(1,series_ngeo)  ; mymeta%l(1) = gmm_layout(1,nn,0,0,nn)
      istat = gmm_create('series_gdone', series_gdone, mymeta, GMM_FLAG_RSTR)
      nn = max(1,series_nstnb) ; mymeta%l(2) = gmm_layout(1,nn,0,0,nn)
      istat = gmm_create('series_gdata', series_gdata, mymeta, GMM_FLAG_RSTR)
      if (ptopo_isblocmaster_L) then
         istat = gmm_create('series_gdata_out', series_gdata_out,mymeta)
      else
         series_gdata_out => series_gdata
      endif

      mymeta = GMM_NULL_METADATA
      nn = max(1,series_nsurf)      ; mymeta%l(1) = gmm_layout(1,nn,0,0,nn)
      istat = gmm_create('series_sdone', series_sdone, mymeta, GMM_FLAG_RSTR)
      if (ptopo_isblocmaster_L) then
         allocate(series_sdone_out(nn))
      else
         series_sdone_out => series_sdone
      endif
      allocate(series_sdone2(nn,max(1,series_nstnl)))
      nn = max(1,series_nstnb)      ; mymeta%l(2) = gmm_layout(1,nn,0,0,nn)
      nn = max(1,series_out_nsteps) ; mymeta%l(3) = gmm_layout(1,nn,0,0,nn)
      istat = gmm_create('series_sdata', series_sdata, mymeta, GMM_FLAG_RSTR)
      if (ptopo_isblocmaster_L) then
         istat = gmm_create('series_sdata_out', series_sdata_out, mymeta)
      else
         series_sdata_out => series_sdata
      endif
      mymeta = GMM_NULL_METADATA
      nn = max(1,series_nprof)      ; mymeta%l(1) = gmm_layout(1,nn,0,0,nn)
      istat = gmm_create('series_pdone', series_pdone, mymeta, GMM_FLAG_RSTR)
      if (ptopo_isblocmaster_L) then
         allocate(series_pdone_out(nn))
      else
         series_pdone_out => series_pdone
      endif
      allocate(series_pdone2(nn,max(1,series_nstnl)))
      nn = max(1,series_nstnb)      ; mymeta%l(2) = gmm_layout(1,nn,0,0,nn)
      nn = max(1,phydim_nk)         ; mymeta%l(3) = gmm_layout(1,nn,0,0,nn)
      nn = max(1,series_out_nsteps) ; mymeta%l(4) = gmm_layout(1,nn,0,0,nn)
      istat = gmm_create('series_pdata', series_pdata, mymeta, GMM_FLAG_RSTR)
      if (ptopo_isblocmaster_L) then
         istat = gmm_create('series_pdata_out', series_pdata_out, mymeta)
      else
         series_pdata_out => series_pdata
      endif
      series_gdone = 0
      series_sdone = 0
      series_pdone = 0
      series_sdone2 = 0
      series_pdone2 = 0
      series_gdata = 0.
      series_sdata = 0.
      series_pdata = 0.

      F_istat = RMN_OK
      series_initok_L = .true.
      !---------------------------------------------------------------
      return
   end function series_init


   !#---- Private functions ------------------------------------------
 
   !/@*
   function priv_stnij(phy_lcl_gid, drv_glb_gid, phy_lcl_ni, phy_lcl_nj) &
        result(F_istat)
      implicit none
!!!#include <arch_specific.hf>
      !@object: compute i,j from lat,lon
      !@params
      integer, intent(in) :: phy_lcl_gid, drv_glb_gid, phy_lcl_ni, phy_lcl_nj
      !@returns
      integer :: F_istat
      !@author: V. Lee - May 2000
      !@revisions
      ! v4_50 - Desgagne M.    - Major revision
      ! 2017-01, S.Chamberland - Major revision
      !*@/
      character(len=128) :: str128
      integer :: i, k, k2, il, jl, ni, nj, istat, pos(1)
      real, dimension(:), allocatable :: xl, yl, xg, yg, lat, lon
      type(SER_STN_T), allocatable :: stn_tmp(:)
      integer :: stn_islocal(NSTATMAX), stn_isbloc(NSTATMAX)
      !---------------------------------------------------------------
      F_istat = RMN_OK
      if (series_nstng < 1) return
      F_istat = RMN_ERR

      call msg(MSG_INFO, PKGNAME_S//'Processing time-series grid points')

      allocate(stn_tmp(series_nstng), &
           xl(series_nstng), yl(series_nstng), &
           xg(series_nstng), yg(series_nstng), &
           lat(series_nstng), lon(series_nstng))

      !# Get grid points i,j
      lat(1:series_nstng) = series_stng(1:series_nstng)%lat
      lon(1:series_nstng) = series_stng(1:series_nstng)%lon

      istat = gdxyfll(phy_lcl_gid, xl, yl, lat, lon, series_nstng)
      istat = min(gdxyfll(drv_glb_gid, xg, yg, lat, lon, series_nstng), istat)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_WARNING, PKGNAME_S//'Problem getting x,y coor from lat,lon in gdxyfll')
         return
      endif

      !# Keep a list of stations on Grid
      k2 = 0
      do k= 1,series_nstng
         series_stng(k)%ig = nint(xg(k))
         series_stng(k)%jg = nint(yg(k))
         IF_GRID: if ( &
              series_stng(k)%ig > 0 .and. &
              series_stng(k)%jg > 0 .and. &
              series_stng(k)%ig < drv_glb_ni .and. &
              series_stng(k)%jg < drv_glb_nj) then
            k2 = k2 + 1
            il = nint(xl(k))
            jl = nint(yl(k))
            if (il > 0 .and. il <= phy_lcl_ni .and. &
                 jl > 0 .and. jl <= phy_lcl_nj) then
               !# Note: i,j are not folded values, ser_init converts them
               series_stng(k)%il = il
               series_stng(k)%jl = jl
            else
               series_stng(k)%il = STN_MISSING
               series_stng(k)%jl = STN_MISSING
            endif
            series_stng(k2) = series_stng(k)
         endif IF_GRID
      enddo
      series_nstng = k2

      !# Put the stations in increasing order of index in a list
      ni = maxval(series_stng(1:series_nstng)%ig)
      nj = maxval(series_stng(1:series_nstng)%jg)
      do k = 1, series_nstng
         series_stng(k)%gidx = series_stng(k)%ig+(series_stng(k)%jg-1)*ni
      enddo

      i = 1
      do k = 1, series_nstng
         pos = minloc(series_stng(1:series_nstng)%gidx)
         stn_tmp(i) = series_stng(pos(1))
         series_stng(pos(1))%gidx = ni*nj + 1
         i = i + 1
      enddo
      series_stng(1:series_nstng) = stn_tmp(1:series_nstng)

      deallocate(stn_tmp, xl, yl, xg, yg, lat, lon, stat=istat)

      !# Print list of globally available stations
      call msg(MSG_INFO, &
           '_________________________________________________________________')
      call msg(MSG_INFO,  PKGNAME_S// &
           'Reordered grid points with ACTUAL lat-lon values and short names')
      call msg(MSG_INFO, &
           '_________________________________________________________________')
      call msg(MSG_INFO, &
           '    N    |        NAME        |   I    |   J    |  LAT   |  LON  |')
      call msg(MSG_INFO, &
           '_________________________________________________________________')

      do k = 1, series_nstng
         write(str128, 912) k, series_stng(k)%name, series_stng(k)%ig, &
              series_stng(k)%jg, series_stng(k)%lat, series_stng(k)%lon
912      format(1x,I5,'    ',a18,'   ',I5,'    ',I5,'    ',f8.3,' ',f8.3,' ')
         call msg(MSG_INFO, str128)
      enddo
      call msg(MSG_INFO, &
           '_________________________________________________________________')

      !# Keep a list of local PE stations
      !# Get folded i,j for local PE stations
      stn_isbloc(:)  = 0
      stn_islocal(:) = 0
      do i = 1, series_nstng
         if ((series_stng(i)%il /= STN_MISSING) .and. &
              (series_stng(i)%jl /= STN_MISSING)) then
             stn_islocal(i) = 1
         endif
      end do

      call RPN_COMM_allreduce(stn_islocal, stn_isbloc, series_nstng, &
           RPN_COMM_INTEGER, RPN_COMM_MAX, RPN_COMM_BLOC_COMM, istat)

      series_nstnb = 0
      series_nstnl = 0
      do i = 1, series_nstng
         if (stn_isbloc(i) > 0) then
            series_nstnb = series_nstnb + 1
            if (stn_islocal(i) > 0) then
               series_nstnl = series_nstnl + 1
               series_stnl_idxb(series_nstnl) = series_nstnb
               series_stng(i)%iphy = ijphy(1, series_stng(i)%il, series_stng(i)%jl)
               series_stng(i)%jphy = ijphy(2, series_stng(i)%il, series_stng(i)%jl)
            endif
            series_stnb(series_nstnb) = series_stng(i)
         endif
      end do
 
      F_istat = RMN_OK
      !---------------------------------------------------------------
      return
   end function priv_stnij


end module series_init_mod
