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

!/@
module hgrid_wb
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use wb_itf_mod
   use ezgrid_mod
   use rmn_gmm
   implicit none
   private
   !@objective Whiteboard (data store) for ezscint id + halo, ij0, stag,...
   !@author Stephane Chamberland, 2012-01
   !@description
   ! Public functions
   public :: hgrid_wb_put, hgrid_wb_get, hgrid_wb_gmmmeta, &
        HGRID_DIEZE,HGRID_PERIODIC
   ! Public constants
   !
!@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   integer,parameter :: HGRID_DIEZE = 1
   integer,parameter :: HGRID_PERIODIC = 1

   character(len=*),parameter :: PREFIX_S = 'HGWB/'
   integer,parameter :: IDX_ID = 1
   integer,parameter :: IDX_I0 = 2
   integer,parameter :: IDX_J0 = 3
   integer,parameter :: IDX_LNI = 4
   integer,parameter :: IDX_LNJ = 5
   integer,parameter :: IDX_HX = 6
   integer,parameter :: IDX_HY = 7
   integer,parameter :: IDX_DIEZE = 8
   integer,parameter :: IDX_PERIODX = 9
   integer,parameter :: IDX_PERIODY = 10
   integer,parameter :: NMAXMETA = 10

!!$   !GMM has/need:
!!$      integer,dimension(4) :: low,high,halo,halomax,n
!!$      integer   :: flags
!!$        GMM_FLAG_STAG_X,GMM_FLAG_STAG_Y,GMM_FLAG_STAG_Z
!!$   !RPN_COMM has:
!!$      gni,gnj,lni,lnj,minx,maxx,miny,maxy,halox,haloy,i0,j0

   !TODO-later: include more meta - minx,maxx,miny,maxy

contains


   !/@*
   function hgrid_wb_put(F_name_S,F_ez_id,F_i0,F_j0,F_lni,F_lnj,F_hx,F_hy, &
        F_dieze,F_periodx,F_periody,F_quiet_L,F_rewrite_L) result(F_id)
      implicit none
      !@objective Store a new grid
      !@arguments
      character(len=*),intent(in) :: F_name_S !- Key (internal var name)
      integer,intent(in) :: F_ez_id
      integer,intent(in),optional :: F_i0,F_j0,F_lni,F_lnj,F_hx,F_hy,F_dieze,F_periodx,F_periody
      logical,intent(in),optional :: F_quiet_L, F_rewrite_L
      !@return
      integer :: F_id !- hgrid id or RMN_ERR
      !@author  S. Chamberland, 2012-01
      !*@/
      character(len=256) :: msg_S
      character(len=WB_MAXNAMELENGTH) :: name_S
      integer :: istat,istat2,grid_meta(NMAXMETA),nij(2),ig14(4),ij0(2),myhgrid,n,v0
      character(len=2) :: grtyp_S,grref_S,dieze_S
      logical :: msgok_L, rewrite_L
      !---------------------------------------------------------------------
      F_id = RMN_ERR
      msgok_L = .true.
      if (present(F_quiet_L)) msgok_L = .not.F_quiet_L
      rewrite_L = .false.
      if (present(F_rewrite_L)) rewrite_L = F_rewrite_L
      if (len_trim(F_name_S) == 0) then
         if (msgok_L) call msg(MSG_ERROR,'(hgrid_wb_put) need to provide a internal name')
         return
      endif
      if (len_trim(PREFIX_S)+len_trim(F_name_S) > WB_MAXNAMELENGTH) then
         if (msgok_L) call msg(MSG_WARNING,'(hgrid_wb_put) name too long, will be trimed: '//trim(F_name_S))
      endif
      if (.not.rewrite_L) then
         v0 = wb_verbosity(WB_MSG_FATAL)
         istat = wb_get(trim(PREFIX_S)//trim(F_name_S),myhgrid)
         istat2 = wb_verbosity(v0)
         if (RMN_IS_OK(istat)) then
            if (msgok_L) call msg(MSG_ERROR,'(hgrid_wb_put) hgrid already exists: '//trim(F_name_S))
            return
         endif
      endif
      istat = ezgrid_params(F_ez_id,nij,grtyp_S,grref_S,ig14,ij0)
      if (.not.RMN_IS_OK(istat)) then
         if (msgok_L) call msg(MSG_ERROR,'(hgrid_wb_put) Need to provide a valid ezscint grid id')
         return
      endif

      grid_meta(IDX_ID) = F_ez_id
      grid_meta(IDX_I0:NMAXMETA) = (/ij0(1),ij0(2),nij(1),nij(2),0,0,0,0,0/)
      if (present(F_dieze)) grid_meta(IDX_DIEZE) = F_dieze
      if (present(F_periodx)) grid_meta(IDX_PERIODX) = F_periodx
      if (present(F_periody)) grid_meta(IDX_PERIODY) = F_periody
      if (grid_meta(IDX_DIEZE) == HGRID_DIEZE) then
         if (present(F_i0)) grid_meta(IDX_I0) = min(max(ij0(1),F_i0),nij(1))
         if (present(F_j0)) grid_meta(IDX_J0) = min(max(ij0(2),F_j0),nij(2))
         if (present(F_lni)) grid_meta(IDX_LNI) = min(F_lni,nij(1)-grid_meta(IDX_I0)+1)
         if (present(F_lnj)) grid_meta(IDX_LNJ) = min(F_lnj,nij(2)-grid_meta(IDX_J0)+1)
      else
         if (present(F_i0)) grid_meta(IDX_I0) = F_i0
         if (present(F_j0)) grid_meta(IDX_J0) = F_j0
         if (present(F_lni)) grid_meta(IDX_LNI) = min(F_lni,nij(1))
         if (present(F_lnj)) grid_meta(IDX_LNJ) = min(F_lnj,nij(2))
      endif
      if (present(F_hx)) grid_meta(IDX_HX) = max(0,F_hx)
      if (present(F_hy)) grid_meta(IDX_HY) = max(0,F_hy)

      istat = wb_put(trim(PREFIX_S)//trim(F_name_S),grid_meta,WB_REWRITE_MANY)
      if (.not.RMN_IS_OK(istat)) then
         if (msgok_L) call msg(MSG_ERROR,'(hgrid_wb_put) problem storing data: '//trim(F_name_S))
         return
      endif
!!$      write(id_S,'(a,I6.6,a)') 'id/',F_ez_id,'/'//trim(F_name_S)
!!$      istat = wb_put(trim(PREFIX_S)//trim(id_S),trim(F_name_S))
      F_id = F_ez_id
      name_S = F_name_S
      n = max(8,len_trim(name_S))
      dieze_S = ' '
      if (grid_meta(IDX_DIEZE) == HGRID_DIEZE) dieze_S = ' #'
      !TODO: print periodx/y info
      write(msg_S,'(a,4(a,2I6,a),a,i6)') &
           '(hgrid_wb) Put: '//name_S(1:n)//dieze_S, &
           ' [G_ni,G_nij=',nij(1),nij(2),']', &
           ' [l_ni,l_nj=',grid_meta(IDX_LNI),grid_meta(IDX_LNJ),']', &
           ' [l_i0,l_j0=',grid_meta(IDX_I0),grid_meta(IDX_J0),']', &
           ' [hx,hy=',grid_meta(IDX_HX),grid_meta(IDX_HY),']', &
           ' id=',F_id
!!$      call msg_toall(MSG_INFO,msg_S)
      if (msgok_L) call msg(MSG_INFO,msg_S)
      !---------------------------------------------------------------------
      return
   end function hgrid_wb_put


   !/@*
   function hgrid_wb_get(F_name_S,F_ez_id,F_i0,F_j0,F_lni,F_lnj,F_hx,F_hy,F_dieze,F_periodx,F_periody,F_quiet_L) result(F_istat)
      implicit none
      !@objective  Retrieve grid info
      !@arguments
      character(len=*),intent(in) :: F_name_S !- Key (internal var name)
      integer,intent(out) :: F_ez_id
      integer,intent(out),optional :: F_i0,F_j0,F_lni,F_lnj,F_hx,F_hy,F_dieze,F_periodx,F_periody
      logical,intent(in),optional :: F_quiet_L
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-01
      !*@/
      integer :: grid_meta(NMAXMETA),n
      character(len=WB_MAXNAMELENGTH) :: name_S
      character(len=256) :: msg_S
      character(len=2)   :: dieze_S
      logical :: msgok_L
      !---------------------------------------------------------------------
      msgok_L = .true.
      if (present(F_quiet_L)) msgok_L = .not.F_quiet_L
      F_istat = wb_get(trim(PREFIX_S)//trim(F_name_S),grid_meta,n)
      if (.not.RMN_IS_OK(F_istat)) then
         if (msgok_L) call msg(MSG_INFO,'(hgrid_wb_get) hgrid not found: '//trim(F_name_S))
         return
      endif
      F_istat = grid_meta(IDX_ID)
      F_ez_id = grid_meta(IDX_ID)
      if (present(F_i0)) F_i0 = grid_meta(IDX_I0)
      if (present(F_j0)) F_j0 = grid_meta(IDX_J0)
      if (present(F_lni)) F_lni = grid_meta(IDX_LNI)
      if (present(F_lnj)) F_lnj = grid_meta(IDX_LNJ)
      if (present(F_hx)) F_hx = grid_meta(IDX_HX)
      if (present(F_hy)) F_hy = grid_meta(IDX_HY)
      if (present(F_dieze)) F_dieze = grid_meta(IDX_DIEZE)
      if (present(F_periodx)) F_periodx = grid_meta(IDX_PERIODX)
      if (present(F_periody)) F_periody = grid_meta(IDX_PERIODY)
      name_S = trim(F_name_S)//'      '
      n = max(8,len_trim(name_S))
      dieze_S = ' '
      if (grid_meta(IDX_DIEZE) == HGRID_DIEZE) dieze_S = ' #'
      !TODO: print periodx/y info
      write(msg_S,'(a,3(A,2I6,a),a,i6)') &
           '(hgrid_wb) Get: '//name_S(1:n)//dieze_S, &
!!$           ' [G_ni,G_nij=',nij(1),nij(2),']', &
           ' [l_ni,l_nj=',grid_meta(IDX_LNI),grid_meta(IDX_LNJ),']', &
           ' [l_i0,l_j0=', grid_meta(IDX_I0), grid_meta(IDX_J0),']', &
           ' [hx,hy=', grid_meta(IDX_HX), grid_meta(IDX_HY),']', &
           ' id=',grid_meta(IDX_ID)
!!$      call msg_toall(MSG_INFOPLUS,msg_S)
      if (msgok_L) call msg(MSG_INFOPLUS,msg_S)
      !---------------------------------------------------------------------
      return
   end function hgrid_wb_get


   !/@*
   function hgrid_wb_gmmmeta(F_name_S,F_meta,F_kn,F_k0) result(F_istat)
      implicit none
      !@objective Create GMM metadata from hgrid_wb entry
      !@arguments
      character(len=*),intent(in) :: F_name_S !- Key (internal var name)
      type(gmm_metadata),intent(out) :: F_meta
      integer,intent(in),optional :: F_kn,F_k0
      !@return
      integer :: F_istat !- exit status
      !@author  S. Chamberland, 2012-03
      !*@/
      integer :: grid_meta(NMAXMETA),n,kn,k0,nk,lowx,highx,hx,nx,lowy,highy,hy,ny
      character(len=WB_MAXNAMELENGTH) :: name_S
      character(len=256) :: msg_S
      !---------------------------------------------------------------------
      F_istat = wb_get(trim(PREFIX_S)//trim(F_name_S),grid_meta,n)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(hgrid_wb_gmmmeta) hgrid not found: '//trim(F_name_S))
         return
      endif
      F_istat = grid_meta(IDX_ID) !grid_id

      if (.not.RMN_IS_OK(F_istat)) return

      k0 = 0 ; kn = 0 ; nk = 0
      if (present(F_kn)) then
         kn = F_kn
         k0 = min(1,kn)
      endif
      if (present(F_k0)) k0 = F_k0
      kn = max(k0,kn)
      if (k0/=0 .or. kn/=0) nk = kn-k0+1 

      nx = grid_meta(IDX_LNI) ; ny = grid_meta(IDX_LNJ)
      hx = grid_meta(IDX_HX)  ; hy = grid_meta(IDX_HY)
      lowx = 1-hx ; highx = nx+hx
      lowy = 1-hy ; highy = ny+hy

      F_meta = GMM_NULL_METADATA
      F_meta%l(1) = gmm_layout(lowx,highx,hx,hx,nx)
      F_meta%l(2) = gmm_layout(lowy,highy,hy,hy,ny)
      F_meta%l(3) = gmm_layout(k0,kn,0,0,nk)
      name_S = F_name_S
      n = max(8,len_trim(name_S))
      write(msg_S,'(a,3(a,i6,a,i6,a,i6,a,i6,a),a,i6)') &
           '(hgrid_wb) GMM: '//name_S(1:n), &
           ' [i0:in,ni,hx=',lowx,':',highx,',',nx,',',hx,']', &
           ' [j0:jn,nj,hy=',lowy,':',highy,',',ny,',',hy,']', &
           ' [k0:kn,nk,hz=',k0,':',kn,',',nk,',',0,']', &
           ' id=',grid_meta(IDX_ID)
!!$      call msg_toall(MSG_INFOPLUS,msg_S)
      call msg(MSG_INFOPLUS,msg_S)
!!$      meta%a%initmode = initmode
!!$      meta%a%flags = flags
      !---------------------------------------------------------------------
      return
   end function hgrid_wb_gmmmeta

end module hgrid_wb
