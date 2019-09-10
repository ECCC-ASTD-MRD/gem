!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!--------------------------------------------------------------------------

!/@*
module drv_grid_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use wb_itf_mod
   use config_mod
   use ptopo_utils
   use hgrid_wb
   use ezgrid_mod
   implicit none
   private
   !@objective Manage Grid configuration
   !@author
   !  Michel Desgagne, Feb 2008
   !  Ron McTaggart-Cowan, Feb 2008
   !  Stephane Chamberland, Feb 2008
   !@revisions
   !  2012-02, Stephane Chamberland: RPNPhy offline
   !@public_functions
   public :: drv_grid_config, drv_grid_init
   !@public_vars
   !
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <mu_gmm.hf>
#include <msg.h>

   include "drv_dyn_itf.inc"

   integer,save   :: l_ni,l_nj,l_nk
   integer,save   :: l_minx,l_maxx,l_miny,l_maxy
   integer,save   :: l_i0,l_j0


contains


   !/@*
   function drv_grid_config(F_cfg_basename_S) result(F_istat)
      implicit none
      !@objective Read grid config from file to WB
      !@arguments
      character(len=*),intent(in) :: F_cfg_basename_S  !-
      !@returns
      integer :: F_istat
      !@author
      !  Stephane Chamberland, Feb 2008
   !*@/
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] grid_config')
      F_istat = config_read(F_cfg_basename_S,'grid_cfgs')
      call msg(MSG_DEBUG,'[END] grid_config')
      !---------------------------------------------------------------------
      return
   end function drv_grid_config


   !/@*
   function drv_grid_init() result(F_istat)
      implicit none
      !@objective Initialize grid configuration
      !@returns
      integer :: F_istat
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
      !@revisions
      !  2012-02, Stephane Chamberland: RPNPhy offline
   !*@/
      integer,external :: gdll
      integer  :: G_ni,G_nj,G_grid_id,istat,gid,i0,j0,in,jn,hx,hy,grid_id2
      integer  :: G_halox,G_haloy,iperiodx,iperiody
      logical  :: G_periodx,G_periody
      character(len=256) :: tmp_S
      real,pointer,dimension(:,:) :: dxdy !#,lat,lon
      type(gmm_metadata) :: mymeta
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] grid_init')
      F_istat =  dyn_grid_init(G_ni,G_nj,G_halox,G_haloy,G_periodx,G_periody,G_grid_id)
      F_istat =  min(priv_init(G_ni,G_nj,G_halox,G_haloy),F_istat)

      if (RMN_IS_OK(F_istat)) then
         write(tmp_S,'(4(A,2I6),a,2l,a)') &
              '(drv_grid) Initialization OK: [G_ni,G_nj=',G_ni,G_nj, &
              "] [l_ni,l_nj=",l_ni,l_nj,"] [l_i0,l_j0=",l_i0,l_j0, &
              "] [hx,hy=",G_halox,G_haloy,"] [periodx,y=", &
              G_periodx,G_periody,"]"
!!$      call msg(MSG_INFO,tmp_S)
         call msg_toall(MSG_INFO,tmp_S)
      else
         call msg(MSG_ERROR,'(drv_grid) Problem in Initialisation')
      endif

      F_istat = min(dyn_grid_post_init(l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj),F_istat)

      !TODO: Define l_haloxy potentially different from G_haloxy

      iperiodx=0 ; if (G_periodx) iperiodx=1
      iperiody=0 ; if (G_periody) iperiody=1
      hx=G_halox ; hy=G_haloy
      istat = hgrid_wb_put('full+h',G_grid_id,1+hx,1+hy,G_ni,G_nj,hx,hy, &
           F_periodX=iperiodx,F_periodY=iperiody)
      istat = hgrid_wb_put('local#+h',G_grid_id,l_i0+hx,l_i0+hy,l_ni,l_nj, &
           hx,hy,F_dieze=HGRID_DIEZE,F_periodX=iperiodx,F_periodY=iperiody)

      i0=l_i0+hx ; j0=l_j0+hy
      in=l_i0+l_ni-1+hx ; jn=l_j0+l_nj-1+hy
      gid = ezgrid_sub(G_grid_id,i0-hx,j0-hy,in+hx,jn+hy)
      istat = hgrid_wb_put('localz+h',gid,1+hx,1+hy,l_ni,l_nj,hx,hy)

      grid_id2 = G_grid_id
      if (G_halox > 0 .or. G_haloy > 0) then
         hx=G_halox ; hy=G_haloy
         i0=1+hx ; j0=1+hy
         in=i0+G_ni-1 ; jn=j0+G_nj-1
         grid_id2 = ezgrid_sub(G_grid_id,i0,j0,in,jn)

         in=l_i0+l_ni-1 ; jn=l_j0+l_nj-1
         gid = ezgrid_sub(grid_id2,l_i0,l_j0,in,jn)
      endif
      hx=0 ; hy=0 ; i0=1 ; j0=1
      istat = hgrid_wb_put('full',grid_id2,i0,j0,G_ni,G_nj,hx,hy, &
           F_periodX=iperiodx,F_periodY=iperiody)
      istat = hgrid_wb_put('local#',grid_id2,l_i0+hx,l_j0+hx,l_ni,l_nj,hx,hy, &
           F_dieze=HGRID_DIEZE,F_periodX=iperiodx,F_periodY=iperiody)
      istat = hgrid_wb_put('localz',gid,i0,j0,l_ni,l_nj,hx,hy)

      !TODO-later: save G_periodx,G_periody,l_minx,l_maxx,l_miny,l_maxy,... somewhere (wb?)

      F_istat = min(hgrid_wb_gmmmeta('localz',mymeta),F_istat)
      nullify(dxdy)
      istat = gmm_create('DXDY',dxdy,mymeta)
      if (.not.associated(dxdy)) F_istat = RMN_ERR
      if (RMN_IS_OK(F_istat)) &
           F_istat = min(dyn_dxdy(dxdy,l_i0,l_j0,l_ni,l_nj),F_istat)

!!$      nullify(lat,lon)
!!$      F_istat = min(gmm_create('LAT',lat,mymeta),F_istat)
!!$      F_istat = min(gmm_create('LON',lon,mymeta),F_istat)
!!$      if (RMN_IS_OK(gid)) then
!!$         F_istat = min(gdll(gid,lat,lon),F_istat)
!!$      else
!!$         F_istat = min(dyn_lat_lon(lat,lon),F_istat)
!!$      endif
      call msg(MSG_DEBUG,'[BEGIN] grid_init')
      !---------------------------------------------------------------------
      return
   end function drv_grid_init


   !/@*
   function priv_init(G_ni,G_nj,G_halox,G_haloy) result(F_istat)
      implicit none
      !@objective Initialize Grid configuration
      !@arguments
      integer  :: G_ni,G_nj
      integer  :: G_halox,G_haloy
      !@returns
      integer :: F_istat
      !@author
      !  Michel Desgagne, Feb 2008
      !  Ron McTaggart-Cowan, Feb 2008
      !  Stephane Chamberland, Feb 2008
   !*@/
      integer :: G_lnimax,G_lnjmax
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'[BEGIN] grid_init')
      F_istat = RMN_OK
      !TODO-later: check that ptopo_init [or at least rpn_comm_init] is already done, if not error

      !- Domain Decomposition
      F_istat = min(rpn_comm_topo(G_ni,l_minx,l_maxx,l_ni,G_lnimax, &
           G_halox,l_i0,RPN_COMM_TOPO_X,RPN_COMM_TOPO_FILL),F_istat)
      F_istat = min(rpn_comm_topo(G_nj,l_miny,l_maxy,l_nj,G_lnjmax, &
           G_haloy,l_j0,RPN_COMM_TOPO_Y,RPN_COMM_TOPO_FILL),F_istat)

      !TODO-later: check_partition

      call msg(MSG_DEBUG,'[END] grid_init')
      !---------------------------------------------------------------------
      return
   end function priv_init

end module drv_grid_mod
