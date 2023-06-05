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

!/@
module ptopo_utils
   use iso_c_binding
   use rpn_comm_itf_mod
   use wb_itf_mod
   implicit none
   private
   !@objective
   !@author  Stephane Chamberland, 2011-03
   !@description
   !   WORLD     COMM is split in pe (not matrix), DOMAIN (not matrix)
   !   DOMAIN    COMM is split in pe (not matrix), MULTIGRID (not matrix)
   !   MULTIGRID COMM is split in pe (not matrix), GRID (not matrix)
   !   GRID      COMM is split in pe (2d-matrix),  BLOC (2d-matrix)
   !   BLOC      COMM is split in pe (2d-matrix)
   !   I/O       COMM is split in pe (not matrix), Sub set of GRID PE used for I/O, optimally placed

   ! Public functions
   public :: ptopo_init_var, ptopo_bloc_set, ptopo_io_set, ptopo_openMP_set, &
        ptopo_ismaster_L, ptopo_collect_dims, ptopo_collect_dims_ij0, &
        ptopo_collect, ptopo_copyfirst2last, ptopo_iprocxy, &
        ptopo_comm_pe_info, ptopo_get_io_params

   ! Public constants
!!$   integer, parameter, public :: PTOPO_WORLD = 1
!!$   integer, parameter, public :: PTOPO_DOMAIN = 2
!!$   integer, parameter, public :: PTOPO_MULTIGRID = 3
!!$   integer, parameter, public :: PTOPO_GRID = 4
   integer, parameter, public :: PTOPO_BLOC = 5
   integer, parameter, public :: PTOPO_IO = 6
   integer, parameter, public :: PTOPO_IODIST = 6

   ! Public var
   integer,public,save :: &
        ptopo_iotype = PTOPO_BLOC, &

        ptopo_npeOpenMP_resv = -1, &
        ptopo_npeOpenMP = -1, &
        ptopo_smt_fact = 1, &
        ptopo_world_ndoms = 1, &      !- WORLD, nb of DOMAINs
        ptopo_world_idom = 0, &       !- WORLD, DOMAIN idx
        ptopo_dom_ngrids = 1, &       !- DOMAIN, nb of GRIDs
        ptopo_dom_igrid = 0, &        !- DOMAIN, GRID idx

        ptopo_grid_npe = 1, &         !- GRID, nb of PE in GRID
        ptopo_grid_npex = 1, &
        ptopo_grid_npey = 1, &        !- GRID, nb of PE along axes(ncol,nrow)
        ptopo_grid_ipe = 0, &         !- GRID, PE idx (rel. to GRID master)
        ptopo_grid_ipex = 0, &
        ptopo_grid_ipey = 0, &        !- GRID, PE idx along axes (col,row)

        ptopo_grid_nbloc = 1, &       !- GRID, nb of BLOC in GRID
        ptopo_grid_nblocx = 1, &
        ptopo_grid_nblocy = 1, &         !- GRID, nb of BLOC along axes
        ptopo_grid_ipe_blocmaster = 0, & !- GRID, idx of PE's blocmaster
        ptopo_grid_ibloc = 0, &          !- GRID, BLOC idx
        ptopo_grid_iblocx = 0, &
        ptopo_grid_iblocy = 0, &      !- GRID,  BLOC idx along axes

        ptopo_bloc_npe = 1, &         !- BLOC,  nb of PE in BLOC
        ptopo_bloc_npex = 1, &
        ptopo_bloc_npey = 1,&         !- BLOC,  nb of PE along axes(ncol,nrow)
        ptopo_bloc_ipe = 0, &         !- BLOC,  PE idx (rel. to bloc master)
        ptopo_bloc_ipex = 0, &
        ptopo_bloc_ipey = 0, &        !- BLOC,  PE idx along axes (col,row)

        ptopo_io_npe = 1, &           !- I/O, Nb of PE for I/O on the GRID
        ptopo_io_setno = -1, &        !- I/O, io_set number
        ptopo_io_ipe = -1, &          !- I/O, PE idx in io_set
        ptopo_io_comm_id = -1, &      !-
        ptopo_grid_ipe_io_master = 0  !- GRID, idx of PE's I/O master

   logical,public,save :: &
        ptopo_isblocmaster_L = .true., &  !- .T. if in communicator "BLOCMASTER"
        ptopo_isiope_L = .true., &        !- .T. if ptopo_io_ipe >= 0
        ptopo_isiomaster_L = .true.       !- .T. if ptopo_grid_ipe == ptopo_grid_ipe_io_master
   character(len=16),public,save :: &
        ptopo_domname_S = ' '             !- Domain name of the current PE


!!$   character(len=1024),public,save :: &
!!$        ptopo_grid_basedir_S, &      !- basdir
!!$        ptopo_pe_basedir_S,   &      !- basdir
   !@Description
   !@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface ptopo_get_io_params
      module procedure ptopo_get_io_params_0
      module procedure ptopo_get_io_params_1
   end interface

   interface ptopo_collect
      module procedure ptopo_collect_r4_3d
      module procedure ptopo_collect_r4_2d
   end interface

   interface ptopo_copyfirst2last
      module procedure ptopo_copyfirst2last_r4_2d
      module procedure ptopo_copyfirst2last_r4_3d
   end interface

   integer, parameter :: TAG = 210
   integer, parameter :: IDX_I0  = 1
   integer, parameter :: IDX_J0  = 2
   integer, parameter :: IDX_NIL = 3
   integer, parameter :: IDX_NJL = 4
   integer, parameter :: N_PARAM = 4

contains


   !/@
   subroutine ptopo_init_var(F_ndomains,F_idomain,F_ngrids,F_igrid)
      implicit none
      !@objective
      !@argument
      integer,intent(in),optional :: F_ndomains,F_idomain,F_ngrids,F_igrid
      !@/
      logical,save :: is_init_L = .false.
      integer :: npex,npey,istat,npew,ipew,nped,iped
      !---------------------------------------------------------------------
      if (.not.is_init_L) call msg(MSG_DEBUG,'(ptopo) init_var [BEGIN]')

      !#Note: Allow re-init bloc topo since it may change
      call priv_update_bloc(npex,npey)
      ptopo_grid_npex = npex
      ptopo_grid_npey = npey

      if (is_init_L) return
      is_init_L = .true.

      istat = ptopo_openMP_set()

      call rpn_comm_size(RPN_COMM_GRID,ptopo_grid_npe,istat)
      istat = rpn_comm_mype(ptopo_grid_ipe,ptopo_grid_ipex,ptopo_grid_ipey)

!!$      RPN_COMM_WORLD,RPN_COMM_ALLGRIDS, RPN_COMM_MULTIGRID
!!$      call RPN_COMM_size('MULTIGRID',Glb_numproc,ierr)
!!$      call RPN_COMM_rank('MULTIGRID',Glb_myproc ,ierr)

      call rpn_comm_size(RPN_COMM_WORLD,npew,istat)
      call rpn_comm_rank(RPN_COMM_WORLD,ipew,istat)
      call rpn_comm_size(RPN_COMM_MULTIGRID,nped,istat)
      call rpn_comm_rank(RPN_COMM_MULTIGRID,iped,istat)

      !TODO-later: n/idom, n/igrid needs to be tested
      ptopo_world_ndoms = npew/max(1,nped)
      ptopo_world_idom = ipew/max(1,nped)
      ptopo_dom_ngrids = nped/max(1,ptopo_grid_npe)
      ptopo_dom_igrid = iped/max(1,ptopo_grid_npe)

      if (present(F_ndomains)) ptopo_world_ndoms = min(max(1,F_ndomains),npew)
      if (present(F_idomain)) ptopo_world_idom = max(0,F_idomain)
      if (present(F_ngrids)) ptopo_dom_ngrids = min(max(1,F_ngrids),npew)
      if (present(F_igrid)) ptopo_dom_igrid = min(max(0,F_igrid),ptopo_dom_ngrids-1)

      istat = wb_put('ptopo/ndomains',ptopo_world_ndoms, WB_REWRITE_AT_RESTART)
      istat = wb_put('ptopo/idomain',ptopo_world_idom, WB_REWRITE_AT_RESTART)
      istat = wb_put('ptopo/ngrids',ptopo_dom_ngrids, WB_REWRITE_AT_RESTART)
      istat = wb_put('ptopo/igrid',ptopo_dom_igrid, WB_REWRITE_AT_RESTART)

      istat = ptopo_io_set(ptopo_io_npe)

      call msg(MSG_DEBUG,'(ptopo) init_var [END]')
      !---------------------------------------------------------------------
      return
   end subroutine ptopo_init_var


   !/@
   function ptopo_bloc_set(F_nblocx,F_nblocy,F_bcast_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(inout) :: F_nblocx,F_nblocy
      logical,intent(in),optional :: F_bcast_L
      !@return
      integer :: F_istat
      !@/
      integer :: istat,npex,npey,nblocxy(2)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) bloc_set [BEGIN]')
      F_istat = RMN_OK
      nblocxy(1:2) = (/max(1,min(F_nblocx,ptopo_grid_npex)), &
           max(1,min(F_nblocy,ptopo_grid_npey))/)
      if (present(F_bcast_L)) then
         if (F_bcast_L) then
            call rpn_comm_bcast(nblocxy,size(nblocxy),RPN_COMM_INTEGER, &
                 RPN_COMM_MASTER,RPN_COMM_GRID,istat)
         endif
      endif
      F_nblocx = nblocxy(1)
      F_nblocy = nblocxy(2)
      if (nblocxy(1) == ptopo_grid_nblocx .and. &
           nblocxy(2) == ptopo_grid_nblocy) return
      F_istat = rpn_comm_bloc(nblocxy(1),nblocxy(2))
      call priv_update_bloc(npex,npey)
      call msg(MSG_DEBUG,'(ptopo) bloc_set [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_bloc_set


   !/@
   function ptopo_io_set(F_io_npe) result(F_istat)
      implicit none
      !@objective
      !@arguments
      integer,intent(in) :: F_io_npe
      !@return
      integer :: F_istat
      !@/
      logical, parameter :: DIAG_L = .false.
      integer, parameter :: FILL_METHODE = 0
      integer, parameter :: ORIGIN_PE = 0
      integer :: istat
      integer :: ptopo_io_pe_xcoord(F_io_npe), ptopo_io_pe_ycoord(F_io_npe)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) io_set [BEGIN]')
      F_istat = RMN_ERR
      call ptopo_init_var()

      if (ptopo_io_setno >= 0 .and. F_io_npe /= ptopo_io_npe) then
!!$         print *,ptopo_grid_ipe,'rpn_comm_free_io_set',ptopo_io_setno, F_io_npe, ptopo_io_npe
         istat = rpn_comm_free_io_set(ptopo_io_setno)
         ptopo_io_setno = -1
      endif

      ptopo_io_npe = max(1, min(F_io_npe, &
                                min(ptopo_grid_npex, ptopo_grid_npey)**2))

      if (ptopo_io_setno < 0) then
!!$         print *,ptopo_grid_ipe,'rpn_comm_io_pe_valid_set', F_io_npe, ptopo_io_npe
         istat = rpn_comm_io_pe_valid_set( &
              ptopo_io_pe_xcoord, ptopo_io_pe_ycoord, &
              ptopo_io_npe, ptopo_grid_npex, ptopo_grid_npey, &
              DIAG_L, FILL_METHODE)
         if (.not.RMN_IS_OK(istat)) then
            ptopo_io_npe = 1
            ptopo_io_setno = -1
            ptopo_io_ipe = -1
            ptopo_grid_ipe_io_master = 0
            ptopo_isiope_L = .true.
            ptopo_isiomaster_L = .true.
            call msg(MSG_WARNING,'(ptopo) Unable to set requested I/O PE distribution')
            return
         endif
         ptopo_io_setno = rpn_comm_create_io_set(ptopo_io_npe, FILL_METHODE)
!!$         print *,ptopo_grid_ipe,'rpn_comm_create_io_set', ptopo_io_setno, ptopo_io_npe
      endif

      ptopo_io_ipe   = rpn_comm_is_io_pe(ptopo_io_setno)
      ptopo_grid_ipe_io_master = rpn_comm_io_pe_gridid(ptopo_io_setno, ORIGIN_PE)
      ptopo_io_comm_id   = rpn_comm_io_pe_comm(ptopo_io_setno)
      ptopo_isiope_L     = (ptopo_io_ipe >= 0)
      ptopo_isiomaster_L = (ptopo_grid_ipe_io_master == ptopo_grid_ipe)

!!$      print *,ptopo_grid_ipe,'ptopo_io_set', ptopo_io_setno, ptopo_io_npe, ptopo_io_ipe, ptopo_grid_ipe_io_master
!!$      F_grid_id    = rpn_comm_create_2dgrid (  G_ni, G_nj, &
!!$                              l_minx, l_maxx, l_miny, l_maxy )

      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(ptopo) io_set [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_io_set


   !/@
   function ptopo_openMP_set(F_npeOpenMP,F_smt_fact) result(F_npeOpenMPout)
      implicit none
      !@objective
      !@arguments
      integer,intent(in),optional :: F_npeOpenMP,F_smt_fact
      !@return
      integer :: F_npeOpenMPout
      !@/
      integer,external :: gnthread
      integer,parameter :: MAX_SMT = 2 !TODO-later: this is machine dependent
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) openMP_set [BEGIN]')
      if (ptopo_npeOpenMP_resv < 1) ptopo_npeOpenMP_resv = max(1,gnthread())
      if (ptopo_npeOpenMP < 1) ptopo_npeOpenMP = ptopo_npeOpenMP_resv
      if (present(F_npeOpenMP)) ptopo_npeOpenMP = F_npeOpenMP
      if (present(F_smt_fact)) ptopo_smt_fact = F_smt_fact
      ptopo_npeOpenMP = min(max(1,ptopo_npeOpenMP),ptopo_npeOpenMP_resv)
      ptopo_smt_fact = min(max(1,ptopo_smt_fact),MAX_SMT)
      F_npeOpenMPout = ptopo_npeOpenMP * ptopo_smt_fact
      !NOTE: pe_rebind is only usefull on AIX p5 (AIX p7?)
!!$      ptopo_npeOpenMP_resv = max(1,gnthread())
!!$      ptopo_npeOpenMP = 1 !- number of processors requested for OpenMp
!!$      ptopo_smtphy    = 1 !- number of threads around PHYSICS
!!$      ptopo_smtglb    = 1 !- number of threads for the global run
      call msg(MSG_DEBUG,'(ptopo) openMP_set [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_openMP_set


   !/@*
   function ptopo_ismaster_L(F_comm_S) result(F_ismaster_L)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_comm_S
      !@return
      logical :: F_ismaster_L
      !@author S. Chamberland, 2010-03
      !*@/
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) ismaster [BEGIN]')
      call ptopo_init_var()
      F_ismaster_L = .false.
      select case(F_comm_S)
      case(RPN_COMM_GRID)
         F_ismaster_L = (ptopo_grid_ipe == RPN_COMM_MASTER)
      case(RPN_COMM_BLOC_COMM)
         F_ismaster_L = ptopo_isblocmaster_L
      end select
      call msg(MSG_DEBUG,'(ptopo) ismaster [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_ismaster_L

   !/@*
   subroutine ptopo_comm_pe_info(F_comm_S,F_ismaster_L,F_npe,F_npx,F_npy, &
        F_ipe,F_ipx,F_ipy)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_comm_S
      logical,intent(out) :: F_ismaster_L
      integer,intent(out) :: F_npe,F_npx,F_npy,F_ipe,F_ipx,F_ipy
      !@author S. Chamberland, 2014-12
      !*@/
      !---------------------------------------------------------------------
      F_ismaster_L = .false.
      F_npe = -1; F_npx = -1; F_npy = -1
      F_ipe = -1; F_ipx = -1; F_ipy = -1
      call ptopo_init_var()
      select case(F_comm_S)
      case(RPN_COMM_GRID)
         F_ismaster_L = (ptopo_grid_ipe == RPN_COMM_MASTER)
         F_npe = ptopo_grid_npe
         F_npx = ptopo_grid_npex
         F_npy = ptopo_grid_npey
         F_ipe = ptopo_grid_ipe
         F_ipx = ptopo_grid_ipex
         F_ipy = ptopo_grid_ipey
      case(RPN_COMM_BLOC_COMM)
         F_ismaster_L = ptopo_isblocmaster_L
         F_npe = ptopo_grid_nbloc
         F_npx = ptopo_grid_nblocx
         F_npy = ptopo_grid_nblocy
         F_ipe = ptopo_grid_ibloc
         F_ipx = ptopo_grid_iblocx
         F_ipy = ptopo_grid_iblocy
      end select
      call msg(MSG_DEBUG,'(ptopo) ismaster [END]')
      !---------------------------------------------------------------------
      return
   end subroutine ptopo_comm_pe_info

  !/@*
    function ptopo_get_io_params_1(F_isiomaster_L, F_isiope_L, &
        F_comm_ipe_io_master, F_communicator_S) result(F_istat)
      implicit none
      !@objective Select ptopo/rpn_comm params depending if blocking or not
      !@arguments
      logical, intent(out) :: F_isiomaster_L, F_isiope_L
      integer, intent(out) :: F_comm_ipe_io_master
      character(len=*) :: F_communicator_S
      !@return
      integer :: F_istat
      !*@/
      !---------------------------------------------------------------------
      F_istat = ptopo_get_io_params_0(ptopo_iotype, F_isiomaster_L, F_isiope_L,&
           F_comm_ipe_io_master, F_communicator_S)
      !---------------------------------------------------------------------
      return
   end function ptopo_get_io_params_1


   !/@*
    function ptopo_get_io_params_0(F_iotype, F_isiomaster_L, F_isiope_L, &
        F_comm_ipe_io_master, F_communicator_S) result(F_istat)
      implicit none
      !@objective Select ptopo/rpn_comm params depending if blocking or not
      !@arguments
      integer, intent(in)  :: F_iotype
      logical, intent(out) :: F_isiomaster_L, F_isiope_L
      integer, intent(out) :: F_comm_ipe_io_master
      character(len=*) :: F_communicator_S
      !@return
      integer :: F_istat
      !*@/
      !---------------------------------------------------------------------
      call ptopo_init_var()
      F_istat = RMN_OK
      F_isiomaster_L = .false.
      F_comm_ipe_io_master = RMN_ERR
      F_communicator_S = ''
      if (F_iotype == PTOPO_BLOC) then
         F_isiomaster_L = ptopo_isblocmaster_L
         F_isiope_L = ptopo_isblocmaster_L
         F_comm_ipe_io_master = RPN_COMM_MASTER !#TODO: check this
         !#TODO: = ptopo_grid_ipe_blocmaster
         F_communicator_S = RPN_COMM_BLOC_COMM
      elseif (F_iotype == PTOPO_IODIST) then
         F_isiomaster_L = ptopo_isiomaster_L
         F_isiope_L = ptopo_isiope_L
         F_comm_ipe_io_master = ptopo_grid_ipe_io_master
         F_communicator_S = RPN_COMM_GRID
!!$      elseif (F_iotype == PTOPO_GRID) then  !#TODO:
!!$         F_isiomaster_L = (ptopo_grid_ipe == RPN_COMM_MASTER)
!!$         F_isiope_L = RPN_COMM_MASTER
!!$         F_comm_ipe_io_master = RPN_COMM_MASTER
!!$         F_communicator_S = RPN_COMM_GRID
!!$      elseif (F_iotype == ) then !#TODO: other WORLD, DOMAIN, MULTIGRID
      else
          F_istat = RMN_ERR
      endif
      !---------------------------------------------------------------------
      return
   end function ptopo_get_io_params_0


   !/@*
   function ptopo_collect_dims(F_comm_S,F_nil,F_njl,F_nic,F_njc,F_i0c,F_j0c) &
        result(F_istat)
      implicit none
      !@objective Get dims of bloc
      !@arguments
      character(len=*),intent(in) :: F_comm_S
      integer,intent(in) :: F_nil,F_njl    !dims of the local data (w/o halo)
      integer,intent(out) :: F_nic,F_njc   !dims of the bloc data
      integer,intent(out) :: F_i0c,F_j0c   !pos of F_fld_l(1,1) in bloc
      !@return
      integer :: F_istat
      !@author S. Chamberland, 2010-03
      !*@/
      integer :: gi0c,gj0c
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) collect_dims [BEGIN]')
      call ptopo_init_var()
      F_istat = RMN_ERR
      select case(F_comm_S)
      case(RPN_COMM_GRID)
         F_istat = priv_collect_dims_grid(ptopo_grid_npex,ptopo_grid_npey, &
              ptopo_grid_ipex,ptopo_grid_ipey,F_nil,F_njl,1,1, &
              F_nic,F_njc,F_i0c,F_j0c,gi0c,gj0c)
      case(RPN_COMM_BLOC_COMM)
         F_istat = priv_collect_dims(F_comm_S,ptopo_bloc_npex,ptopo_bloc_npey, &
              ptopo_bloc_ipex,ptopo_bloc_ipey,F_nil,F_njl,1,1, &
              F_nic,F_njc,F_i0c,F_j0c,gi0c,gj0c)
      end select
!!$      print *,'(ptopo1) '//trim(F_comm_S),ptopo_grid_ipex,ptopo_grid_ipey,F_nic,F_njc ; call flush(6)
      call msg(MSG_DEBUG,'(ptopo) collect_dims [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_collect_dims


   !/@*
   function ptopo_collect_dims_ij0(F_comm_S,F_nil,F_njl,F_gi0,F_gj0, &
        F_nic,F_njc,F_i0c,F_j0c,F_gi0c,F_gj0c) result(F_istat)
      implicit none
      !@objective Get dims of bloc
      !@arguments
      character(len=*),intent(in) :: F_comm_S
      integer,intent(in)  :: F_nil,F_njl    !dims of the local data (w/o halo)
      integer,intent(in)  :: F_gi0,F_gj0    !pos of local data(1,1) in full grid
      integer,intent(out) :: F_nic,F_njc    !dims of the bloc data
      integer,intent(out) :: F_i0c,F_j0c    !pos of F_fld_l(1,1) in bloc
      integer,intent(out) :: F_gi0c,F_gj0c  !pos of bloc data(1,1) in full grid
      !@return
      integer :: F_istat
      !@author S. Chamberland, 2010-03
      !*@/
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) collect_dims [BEGIN]')
      call ptopo_init_var()
      F_istat = RMN_ERR
      select case(F_comm_S)
      case(RPN_COMM_GRID)
         F_istat = priv_collect_dims_grid(ptopo_grid_npex,ptopo_grid_npey, &
              ptopo_grid_ipex,ptopo_grid_ipey,F_nil,F_njl,F_gi0,F_gj0, &
              F_nic,F_njc,F_i0c,F_j0c,F_gi0c,F_gj0c)
      case(RPN_COMM_BLOC_COMM)
         F_istat = priv_collect_dims(F_comm_S,ptopo_bloc_npex,ptopo_bloc_npey, &
              ptopo_bloc_ipex,ptopo_bloc_ipey,F_nil,F_njl,F_gi0,F_gj0, &
              F_nic,F_njc,F_i0c,F_j0c,F_gi0c,F_gj0c)
      end select
      call msg(MSG_DEBUG,'(ptopo) collect_dims [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_collect_dims_ij0


   !/@*
   function ptopo_collect_r4_3d(F_fld_b,F_fld_l,F_comm_S,F_i0l,F_j0l, &
        F_lni,F_lnj) result(F_istat)
      implicit none
      !@objective Merge MPI data-tiles belonging to the same bloc
      !@arguments
      real,pointer :: F_fld_b(:,:,:)  !O, F_comm_S data
      real,pointer :: F_fld_l(:,:,:)  !I, local data
      integer,intent(in) :: F_i0l,F_j0l !pos of F_fld_l(1,1) in F_comm_S
      integer,intent(in),optional :: F_lni,F_lnj !Size of fld_l w/o halo
      character(len=*),intent(in) :: F_comm_S
      !@return
      integer :: F_istat
      !@author
      !@description
      !*@/
      logical :: ismaster_L
      integer :: lni,lnj
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) collect_r4_3d [BEGIN]')
      call ptopo_init_var()
      F_istat = RMN_ERR
      if (.not.associated(F_fld_l)) then
         print *,'ERROR: (ptopo_collect) .not.associated(F_fld_l)';call flush(6)
      endif
      lni = ubound(F_fld_l,1)
      lnj = ubound(F_fld_l,2)
      if (present(F_lni))  lni = min(lni,F_lni)
      if (present(F_lnj))  lnj = min(lnj,F_lnj)

      select case(F_comm_S)
      case(RPN_COMM_GRID)
         ismaster_L = (ptopo_grid_ipe == RPN_COMM_MASTER)
         if (ismaster_L .and. .not.associated(F_fld_b)) &
              print *,'ERROR: (ptopo_collect) .not.associated(F_fld_l)'

         F_istat = priv_collect_r4_3d(F_fld_b,F_fld_l,F_i0l,F_j0l,lni,lnj, &
              F_comm_S,ismaster_L,ptopo_grid_npe)
      case(RPN_COMM_BLOC_COMM)
         if (ptopo_isblocmaster_L .and. .not.associated(F_fld_b)) then
            print *,'ERROR: (ptopo_collect) .not.associated(F_fld_b)' ; call flush(6)
         endif
         F_istat = priv_collect_r4_3d(F_fld_b,F_fld_l,F_i0l,F_j0l,lni,lnj, &
              F_comm_S,ptopo_isblocmaster_L,ptopo_bloc_npe)
      end select
      call msg(MSG_DEBUG,'(ptopo) collect_r4_3d [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_collect_r4_3d


   !/@*
   function ptopo_collect_r4_2d(F_fld2d_b,F_fld2d_l,F_comm_S,F_i0l,F_j0l, &
        F_lni,F_lnj) result(F_istat)
      implicit none
      !@objective Merge MPI data-tiles belonging to the same bloc
      !@arguments
      real,pointer :: F_fld2d_b(:,:)  !O, F_comm_S data
      real,pointer :: F_fld2d_l(:,:)  !I, local data
      integer,intent(in) :: F_i0l,F_j0l !pos of F_fld_l(1,1) in F_comm_S
      integer,intent(in),optional :: F_lni,F_lnj !Size of fld_l w/o halo
      character(len=*),intent(in) :: F_comm_S
      !@return
      integer :: F_istat
      !@author
      !@description
      !*@/
      logical :: ismaster_L
      integer :: lni,lnj
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) collect_r4_2d [BEGIN]')
      call ptopo_init_var()
      F_istat = RMN_ERR
      if (.not.associated(F_fld2d_l)) then
         print *,'(ptopo_collect) ERROR: .not.associated(F_fld2d_l)' ; call flush(6)
      endif
      lni = ubound(F_fld2d_l,1)
      lnj = ubound(F_fld2d_l,2)
      if (present(F_lni)) lni = min(lni,F_lni)
      if (present(F_lnj)) lnj = min(lnj,F_lnj)

      select case(F_comm_S)
      case(RPN_COMM_GRID)
         ismaster_L = (ptopo_grid_ipe == RPN_COMM_MASTER)
         if (ismaster_L .and. .not.associated(F_fld2d_b)) then
            print *,'(ptopo_collect) ERROR: .not.associated(F_fld2d_b)' ; call flush(6)
         endif
         F_istat = priv_collect_r4_2d(F_fld2d_b,F_fld2d_l,F_i0l,F_j0l,lni,lnj, &
              F_comm_S,ismaster_L,ptopo_grid_npe)
      case(RPN_COMM_BLOC_COMM)
         if (ptopo_isblocmaster_L .and. .not.associated(F_fld2d_b)) then
            print *,'(ptopo_collect) ERROR: .not.associated(F_fld2d_b)' ; call flush(6)
         endif
         F_istat = priv_collect_r4_2d(F_fld2d_b,F_fld2d_l,F_i0l,F_j0l,lni,lnj, &
              F_comm_S,ptopo_isblocmaster_L,ptopo_bloc_npe)
      end select
      call msg(MSG_DEBUG,'(ptopo) collect_r4_2d [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_collect_r4_2d


   !/@*
   function ptopo_copyfirst2last_r4_2d(F_fld,F_copyx,F_copyy,F_bloc_L) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_fld(:,:)  !I/O, local data on F_comm_S
      integer,intent(in) :: F_copyx,F_copyy !# number of points to copy along x/y
      logical,intent(in) :: F_bloc_L !# true if data is on bloc instead of grid
      !@return
      integer :: F_istat
      !@author S. Chamberland, 2014-12
      !@description
      !*@/
      integer :: copyx,copyy,npe,npx,npy,ipe,ipx,ipy,lijk(2),uijk(2),iproc, &
           istat,ierr,fr_in,fr_i0,to_in,to_i0,bnpe,bnpx,bnpy,bipe,bipx,bipy
      logical :: ismaster_L
      real,pointer :: datarow1(:,:)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) copyfirst2last_r4_2d [BEGIN]')
      F_istat = RMN_OK
      copyx = max(0,F_copyx) ; copyy = max(0,F_copyy)
      if (copyx+copyy == 0) return
      call ptopo_init_var()

      !# Only 1st and last PE row/col should go forward
      call ptopo_comm_pe_info(RPN_COMM_GRID,ismaster_L,npe,npx,npy,ipe,ipx,ipy)
      if (.not.( &
           (copyx > 0 .and. (ipx==RPN_COMM_MASTER .or. ipx==npx-1)) .or. &
           (copyy > 0 .and. (ipy==RPN_COMM_MASTER .or. ipy==npy-1)) )) return

      if (copyy > 0) then
         print *,'(ptopo_copyfirst2last) ERROR: F_copyy not yet supported' ; call flush(6)
         F_istat = RMN_ERR
         return
      endif

      if (F_bloc_L) then
         call ptopo_comm_pe_info(RPN_COMM_BLOC_COMM,ismaster_L,bnpe,bnpx,bnpy, &
              bipe,bipx,bipy)
      else
         call ptopo_comm_pe_info(RPN_COMM_GRID,ismaster_L,bnpe,bnpx,bnpy,bipe, &
              bipx,bipy)
         ismaster_L = .true.
      endif

      if (.not.ismaster_L) return
      F_istat = RMN_ERR

      !# Copy fist copyxy rows to the last ones
      if (.not.associated(F_fld)) then
         print *,'(ptopo_copyfirst2last) ERROR: .not.associated(F_fld)' ; call flush(6)
         return
      endif
      lijk = lbound(F_fld)
      uijk = ubound(F_fld)
      fr_in=lijk(1) ; fr_i0=lijk(1)+(copyx-1)
      to_in=uijk(1) ; to_i0=to_in-(copyx-1)


      IF_COPYX: if (copyx > 0) then

         if (bnpx == 1) then
            F_fld(to_i0:to_in,:) = F_fld(fr_i0:fr_in,:)
         else if (bipy==RPN_COMM_MASTER) then

            nullify(datarow1)
            if (ipx==RPN_COMM_MASTER .or. ipx==npx-1) &
                 allocate(datarow1(1:copyx,lijk(2):uijk(2)))
            if (ipx==RPN_COMM_MASTER) then
               datarow1(1:copyx,:) = F_fld(fr_i0:fr_in,:)
               iproc = ptopo_iprocxy((npx-1),ipy,npx)
               call rpn_comm_send(datarow1,size(datarow1),RPN_COMM_REAL,iproc, &
                    TAG,RPN_COMM_GRID,ierr)
            elseif (ipx==npx-1) then
               iproc = ptopo_iprocxy(0,ipy,npx)
               call rpn_comm_recv(datarow1,size(datarow1),RPN_COMM_REAL,iproc, &
                    TAG,RPN_COMM_GRID,istat,ierr)
               F_fld(to_i0:to_in,:) = datarow1(1:copyx,:)
            endif
            if (associated(datarow1)) deallocate(datarow1,stat=istat)

         endif

      endif IF_COPYX
      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(ptopo) copyfirst2last_r4_2d [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_copyfirst2last_r4_2d


   !/@*
   function ptopo_copyfirst2last_r4_3d(F_fld,F_copyx,F_copyy,F_bloc_L) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_fld(:,:,:)  !I/O, local data on F_comm_S
      integer,intent(in) :: F_copyx,F_copyy !# number of points to copy along x/y
      logical,intent(in) :: F_bloc_L !# true if data is on bloc instead of grid
      !@return
      integer :: F_istat
      !@author S. Chamberland, 2014-12
      !@description
      !*@/
      integer :: copyx,copyy,npe,npx,npy,ipe,ipx,ipy,lijk(3),uijk(3),iproc, &
           istat,ierr,fr_in,fr_i0,to_in,to_i0,bnpe,bnpx,bnpy,bipe,bipx,bipy
      logical :: ismaster_L
      real,pointer :: datarow1(:,:,:)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) copyfirst2last_r4_3d [BEGIN]')
      F_istat = RMN_OK
      copyx = max(0,F_copyx) ; copyy = max(0,F_copyy)
      if (copyx+copyy == 0) return
      call ptopo_init_var()

      !# Only 1st and last PE row/col should go forward
      call ptopo_comm_pe_info(RPN_COMM_GRID,ismaster_L,npe,npx,npy,ipe,ipx,ipy)
      if (.not.( &
           (copyx > 0 .and. (ipx==RPN_COMM_MASTER .or. ipx==npx-1)) .or. &
           (copyy > 0 .and. (ipy==RPN_COMM_MASTER .or. ipy==npy-1)) )) return

      if (copyy > 0) then
         print *,'(ptopo_copyfirst2last) ERROR: F_copyy not yet supported' ; call flush(6)
         F_istat = RMN_ERR
         return
      endif

      if (F_bloc_L) then
         call ptopo_comm_pe_info(RPN_COMM_BLOC_COMM,ismaster_L,bnpe,bnpx,bnpy, &
              bipe,bipx,bipy)
      else
         call ptopo_comm_pe_info(RPN_COMM_GRID,ismaster_L,bnpe,bnpx,bnpy,bipe, &
              bipx,bipy)
         ismaster_L = .true.
      endif

      if (.not.ismaster_L) return
      F_istat = RMN_ERR

      !# Copy fist copyxy rows to the last ones
      if (.not.associated(F_fld)) then
         print *,'(ptopo_copyfirst2last) ERROR: .not.associated(F_fld)' ; call flush(6)
         return
      endif
      lijk = lbound(F_fld)
      uijk = ubound(F_fld)
      fr_in=lijk(1) ; fr_i0=lijk(1)+(copyx-1)
      to_in=uijk(1) ; to_i0=to_in-(copyx-1)


      IF_COPYX: if (copyx > 0) then

         if (bnpx == 1) then
            F_fld(to_i0:to_in,:,:) = F_fld(fr_i0:fr_in,:,:)
         else if (bipy==RPN_COMM_MASTER) then

            nullify(datarow1)
            if (ipx==RPN_COMM_MASTER .or. ipx==npx-1) &
                 allocate(datarow1(1:copyx,lijk(2):uijk(2),lijk(3):uijk(3)))
            if (ipx==RPN_COMM_MASTER) then
               datarow1(1:copyx,:,:) = F_fld(fr_i0:fr_in,:,:)
               iproc = ptopo_iprocxy((npx-1),ipy,npx)
               call rpn_comm_send(datarow1,size(datarow1),RPN_COMM_REAL, &
                    iproc,TAG,RPN_COMM_GRID,ierr)
            elseif (ipx==npx-1) then
               iproc = ptopo_iprocxy(0,ipy,npx)
               call rpn_comm_recv(datarow1,size(datarow1),RPN_COMM_REAL, &
                    iproc,TAG,RPN_COMM_GRID,istat,ierr)
               F_fld(to_i0:to_in,:,:) = datarow1(1:copyx,:,:)
            endif
            if (associated(datarow1)) deallocate(datarow1,stat=istat)

         endif

      endif IF_COPYX
      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(ptopo) copyfirst2last_r4_3d [END]')
      !---------------------------------------------------------------------
      return
   end function ptopo_copyfirst2last_r4_3d


   !/@
   function ptopo_iprocxy(F_ipx,F_ipy,F_npx) result(F_iproc)
      implicit none
      integer,intent(in) :: F_ipx,F_ipy,F_npx
      integer :: F_iproc
      !@/
      !------------------------------------------------------------------
      F_iproc = F_ipx + F_ipy*F_npx
      !------------------------------------------------------------------
      return
   end function ptopo_iprocxy


   !==== Private Functions =================================================

   !/@
   subroutine priv_update_bloc(F_npex,F_npey)
      implicit none
      !@arguments
      integer,intent(out) :: F_npex,F_npey
      !@/
      integer :: me,medomm,mex,mey,blocsizex,blocsizey,ismaster,mymaster, &
           mybloc,myblocx,myblocy,blocme,blocmex,blocmey
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_update_bloc [BEGIN]')
      call rpn_comm_carac(F_npex,F_npey,me,medomm,mex,mey,blocsizex, &
           blocsizey,ismaster,mymaster,mybloc,myblocx,myblocy,blocme, &
           ptopo_domname_S)

      if (all((/F_npex,F_npey/) > 0)) then
         call rpn_comm_bloctopo(blocme,blocmex,blocmey,blocsizex,blocsizey)
      else
         F_npex = 1
         F_npey = 1
         blocme = 0
         blocmex = 0
         blocmey = 0
         blocsizex = 1
         blocsizey = 1
         ismaster = 1
      endif

      ptopo_grid_nblocx = blocsizex
      ptopo_grid_nblocy = blocsizey
      ptopo_grid_nbloc = blocsizex * blocsizey
      ptopo_grid_ipe_blocmaster = mymaster
      ptopo_grid_ibloc = mybloc
      ptopo_grid_iblocx = myblocx
      ptopo_grid_iblocy = myblocy
      ptopo_bloc_npex = F_npex / ptopo_grid_nblocx
      ptopo_bloc_npey = F_npey / ptopo_grid_nblocy
      ptopo_bloc_npe = ptopo_bloc_npex * ptopo_bloc_npey
      ptopo_bloc_ipe = blocme
      ptopo_bloc_ipex = blocmex
      ptopo_bloc_ipey = blocmey
      ptopo_isblocmaster_L = (ismaster == 1)
      call msg(MSG_DEBUG,'(ptopo) priv_update_bloc [END]')
      !------------------------------------------------------------------
      return
   end subroutine priv_update_bloc


   !/@
   function priv_recv_collect_r4_3d(F_fld_b,F_nil2,F_njl2,F_nk2,F_i0l,F_j0l, &
        F_comm_S,F_iproc) result(F_istat)
      implicit none
      !@arguments
      real,pointer :: F_fld_b(:,:,:)            !O, bloc data
      integer,intent(in) :: F_nil2,F_njl2,F_nk2 !dims of the local data
      integer,intent(in) :: F_i0l,F_j0l         !F_fld_l(1,1) pos in F_fld_b
      character(len=*),intent(in) :: F_comm_S   !communicator name
      integer,intent(in) :: F_iproc             !sending ipe in comm_S
      !@return
      integer :: F_istat
      !@/
      integer :: istat,ierr, datalen
      real,target :: buf(F_nil2,F_njl2,F_nk2)
      real,pointer :: p_buf(:,:,:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_recv_collect_r4_3d [BEGIN]')
      F_istat = RMN_ERR
      datalen = F_nil2 * F_njl2 * F_nk2
      call rpn_comm_recv(buf,datalen,RPN_COMM_REAL,F_iproc,TAG,F_comm_S, &
           istat,ierr)
      if (.not.RMN_IS_OK(ierr)) return
      p_buf => buf
      F_istat = priv_copy_r4_3d(F_fld_b,p_buf,F_i0l,F_j0l,F_nil2,F_njl2)
      call msg(MSG_DEBUG,'(ptopo) priv_recv_collect_r4_3d [END]')
      !------------------------------------------------------------------
      return
   end function priv_recv_collect_r4_3d


   !/@
   function priv_recv_collect_r4_2d(F_fld_b,F_nil2,F_njl2,F_i0l,F_j0l, &
        F_comm_S,F_iproc) result(F_istat)
      implicit none
      !@arguments
      real,pointer :: F_fld_b(:,:)              !O, bloc data
      integer,intent(in) :: F_nil2,F_njl2       !dims of the local data
      integer,intent(in) :: F_i0l,F_j0l         !F_fld_l(1,1) pos in F_fld_b
      character(len=*),intent(in) :: F_comm_S   !communicator name
      integer,intent(in) :: F_iproc             !sending ipe in comm_S
      !@return
      integer :: F_istat
      !@/
      integer :: istat,ierr, datalen
      real,target :: buf(F_nil2,F_njl2)
      real,pointer :: p_buf(:,:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_recv_collect_r4_2d [BEGIN]')
      F_istat = RMN_ERR
      datalen = F_nil2 * F_njl2
      call rpn_comm_recv(buf,datalen,RPN_COMM_REAL,F_iproc,TAG,F_comm_S, &
           istat,ierr)
      if (.not.RMN_IS_OK(ierr)) return
      p_buf => buf
      F_istat = priv_copy_r4_2d(F_fld_b,p_buf,F_i0l,F_j0l,F_nil2,F_njl2)
      call msg(MSG_DEBUG,'(ptopo) priv_recv_collect_r4_2d [END]')
      !------------------------------------------------------------------
      return
   end function priv_recv_collect_r4_2d


   !/@
   function priv_copy_r4_3d(F_fld_b,F_fld_l,F_i0l,F_j0l,F_lni,F_lnj) &
        result(F_istat)
      implicit none
      !@arguments
      real,pointer :: F_fld_b(:,:,:)       !O, bloc data
      real,pointer :: F_fld_l(:,:,:)       !I, local data
      integer,intent(in) :: F_i0l,F_j0l    !F_fld_l(1,1) pos in F_fld_b
      integer,intent(in) :: F_lni,F_lnj    !Size of fld_l w/o halo
      !@return
      integer :: F_istat
      !@/
      integer :: i, j, k, i2,j2,l_ijk(3),u_ijk(3),lni,lnj
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_copy_r4_3d [BEGIN]')
      l_ijk = lbound(F_fld_b)
      u_ijk = ubound(F_fld_b)
      lni = min(F_lni,ubound(F_fld_l,1))
      lnj = min(F_lnj,ubound(F_fld_l,2))
      do k = 1, ubound(F_fld_l,3)
         do j = 1, lnj
            j2 =  max(l_ijk(2),min(j + F_j0l-1,u_ijk(2)))
            do i = 1, lni
               i2 = max(l_ijk(1),min(i + F_i0l-1,u_ijk(1)))
               F_fld_b(i2,j2,k) = F_fld_l(i,j,k)
            enddo
         enddo
      enddo
      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(ptopo) priv_copy_r4_3d [END]')
      !------------------------------------------------------------------
      return
   end function priv_copy_r4_3d


   !/@
   function priv_copy_r4_2d(F_fld_b,F_fld_l,F_i0l,F_j0l,F_lni,F_lnj) &
        result(F_istat)
      implicit none
      !@arguments
      real,pointer :: F_fld_b(:,:)       !O, bloc data
      real,pointer :: F_fld_l(:,:)       !I, local data
      integer,intent(in) :: F_i0l,F_j0l  !F_fld_l(1,1) pos in F_fld_b
      integer,intent(in) :: F_lni,F_lnj  !Size of fld_l w/o halo
      !@return
      integer :: F_istat
      !@/
      integer :: i, j, i2,j2,l_ijk(2),u_ijk(2),lni,lnj
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_copy_r4_2d [BEGIN]')
      l_ijk = lbound(F_fld_b)
      u_ijk = ubound(F_fld_b)
      lni = min(F_lni,ubound(F_fld_l,1))
      lnj = min(F_lnj,ubound(F_fld_l,2))
      do j = 1, lnj
         j2 =  max(l_ijk(2),min(j + F_j0l-1,u_ijk(2)))
         do i = 1, lni
            i2 = max(l_ijk(1),min(i + F_i0l-1,u_ijk(1)))
            F_fld_b(i2,j2) = F_fld_l(i,j)
         enddo
      enddo
      F_istat = RMN_OK
      call msg(MSG_DEBUG,'(ptopo) priv_copy_r4_2d [END]')
      !------------------------------------------------------------------
      return
   end function priv_copy_r4_2d

   
   !/@*
   function priv_collect_dims(F_comm_S,F_npex,F_npey,F_ipex,F_ipey, &
        F_nil,F_njl,F_gi0,F_gj0,F_nic,F_njc,F_i0c,F_j0c,F_gi0c,F_gj0c) &
        result(F_istat)
      implicit none
      !@objective Get dims of bloc
      !@arguments
      character(len=*),intent(in) :: F_comm_S
      integer,intent(in) :: F_npex,F_npey,F_ipex,F_ipey !nb of pe and rank in comm_S
      integer,intent(in)  :: F_nil,F_njl   !dims of the local data (w/o halo)
      integer,intent(in)  :: F_gi0,F_gj0   !pos of local data(1,1) in full grid
      integer,intent(out) :: F_nic,F_njc   !dims of the comm data
      integer,intent(out) :: F_i0c,F_j0c   !pos of F_fld_l(1,1) in comm
      integer,intent(out) :: F_gi0c,F_gj0c !pos of bloc data(1,1) in full grid
      !@return
      integer :: F_istat
      !@author S. Chamberland, 2012-01
      !*@/
      integer, parameter :: IDX_GI0 = N_PARAM+1
      integer, parameter :: IDX_GJ0 = N_PARAM+2
      integer :: i,j,ni,nj
      integer :: params(N_PARAM+2),params2d(N_PARAM+2,F_npex,F_npey)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_collect_dims [BEGIN]')
      F_istat = RMN_OK
      F_nic = F_nil
      F_njc = F_njl
      F_i0c = 1
      F_j0c = 1
      F_gi0c = F_gi0
      F_gj0c = F_gj0
      if (F_npex == 1 .and. F_npey == 1) return

      if (.not.priv_alongX()) then
         F_istat = RMN_ERR
         call msg(MSG_ERROR,'(ptopo) priv_collect_dims for comm='//trim(F_comm_S)//' Not yet supported when topology is not alongX')
         return
      endif

      params(IDX_I0)  = F_ipex
      params(IDX_J0)  = F_ipey
      params(IDX_NIL) = F_nil
      params(IDX_NJL) = F_njl
      params(IDX_GI0) = F_gi0
      params(IDX_GJ0) = F_gj0

      params2d = 0
      call rpn_comm_allgather( &
           params, N_PARAM+2, RPN_COMM_INTEGER, &
           params2d, N_PARAM+2, RPN_COMM_INTEGER, &
           F_comm_S, F_istat)
     if (RMN_IS_OK(F_istat)) then
         i = F_ipex+1
         j = F_ipey+1
         ni = F_npex
         nj = F_npey
         F_nic = sum(params2d(IDX_NIL,1:ni ,j))
         F_njc = sum(params2d(IDX_NJL,i    ,1:nj))
         F_i0c = sum(params2d(IDX_NIL,1:i-1,j))     + 1
         F_j0c = sum(params2d(IDX_NJL,i    ,1:j-1)) + 1
         F_gi0c = minval(params2d(IDX_GI0,:,:))
         F_gj0c = minval(params2d(IDX_GJ0,:,:))
      endif
      call msg(MSG_DEBUG,'(ptopo) priv_collect_dims [END]')
      !---------------------------------------------------------------------
      return
   end function priv_collect_dims


   !/@*
   function priv_collect_dims_grid(F_npex,F_npey,F_ipex,F_ipey, &
        F_nil,F_njl,F_gi0,F_gj0,F_nic,F_njc,F_i0c,F_j0c,F_gi0c,F_gj0c) &
        result(F_istat)
      implicit none
      !@objective Get dims of bloc
      !@arguments
      integer,intent(in) :: F_npex,F_npey,F_ipex,F_ipey !nb of pe and rank in comm_S
      integer,intent(in)  :: F_nil,F_njl   !dims of the local data (w/o halo)
      integer,intent(in)  :: F_gi0,F_gj0   !pos of local data(1,1) in full grid
      integer,intent(out) :: F_nic,F_njc   !dims of the comm data
      integer,intent(out) :: F_i0c,F_j0c   !pos of F_fld_l(1,1) in comm
      integer,intent(out) :: F_gi0c,F_gj0c !pos of bloc data(1,1) in full grid
      !@return
      integer :: F_istat
      !@author S. Chamberland, 2012-01
      !*@/
      integer :: params_row(F_npex)
      integer :: params_col(F_npey)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_collect_dims_grid [BEGIN]')
      F_istat = RMN_OK
      F_nic = F_nil
      F_njc = F_njl
      F_i0c = 1
      F_j0c = 1
      F_gi0c = F_gi0
      F_gj0c = F_gj0
      if (F_npex == 1 .and. F_npey == 1) return

      call rpn_comm_allreduce(F_gi0, F_gi0c, 1, RPN_COMM_INTEGER, &
           &    RPN_COMM_MIN, RPN_COMM_EW, F_istat)
      if (RMN_IS_OK(F_istat)) &
           call rpn_comm_allreduce(F_gj0, F_gj0c, 1, RPN_COMM_INTEGER, &
           &    RPN_COMM_MIN, RPN_COMM_NS, F_istat)

      if (RMN_IS_OK(F_istat)) &
           call rpn_comm_allgather(F_nil, 1, RPN_COMM_INTEGER, &
           &    params_row, 1, RPN_COMM_INTEGER, RPN_COMM_EW, F_istat)
      if (RMN_IS_OK(F_istat)) &
           call rpn_comm_allgather(F_njl, 1, RPN_COMM_INTEGER, &
           &    params_col, 1, RPN_COMM_INTEGER, RPN_COMM_NS, F_istat)
            
      F_i0c = sum(params_row(1:F_ipex)) + 1
      F_j0c = sum(params_col(1:F_ipey)) + 1
      F_nic = sum(params_row)
      F_njc = sum(params_col)
 
      call msg(MSG_DEBUG,'(ptopo) priv_collect_dims_grid [END]')
      !---------------------------------------------------------------------
      return
   end function priv_collect_dims_grid

   
   !/@*
   function priv_collect_r4_3d(F_fld_c,F_fld_l,F_i0l,F_j0l,F_lni,F_lnj, &
        F_comm_S,F_ismaster_L,F_npe) result(F_istat)
      implicit none
      !@objective Merge MPI data-tiles belonging to the same comm_S
      !@arguments
      real,pointer :: F_fld_c(:,:,:)    !O, comm_S data
      real,pointer :: F_fld_l(:,:,:)    !I, local data
      integer,intent(in) :: F_i0l,F_j0l !pos of F_fld_l(1,1) in comm_S
      integer,intent(in) :: F_lni,F_lnj !Size of fld_l w/o halo
      character(len=*),intent(in) :: F_comm_S
      logical,intent(in) :: F_ismaster_L
      integer,intent(in) :: F_npe
      !@return
      integer :: F_istat
      !@revision
      !  2010-03, S. Chamberland - split out of gem 4.1.2
      !@description
      !*@/
      integer :: iproc, istat,ierr, datalen, nijk_l(3)
      integer :: params(N_PARAM)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_collect_r4_3d [BEGIN]')
      call ptopo_init_var()
      nijk_l = ubound(F_fld_l)
      nijk_l(1) = min(F_lni,nijk_l(1))
      nijk_l(2) = min(F_lnj,nijk_l(2))

      IF_MASTER: if (F_ismaster_L) then

         istat = priv_copy_r4_3d(F_fld_c,F_fld_l,F_i0l,F_j0l,F_lni,F_lnj)

         DO_IPROC: do iproc = 1, F_npe - 1

            call rpn_comm_recv(params,N_PARAM,RPN_COMM_INTEGER,iproc,TAG, &
                 F_comm_S,istat,ierr)
            if (params(IDX_NIL) * params(IDX_NJL) > 0) then
               istat = priv_recv_collect_r4_3d(F_fld_c,params(IDX_NIL), &
                    params(IDX_NJL),nijk_l(3),params(IDX_I0),params(IDX_J0), &
                    F_comm_S,iproc)
            endif

         enddo DO_IPROC

      else !IF_MASTER

         params(IDX_I0)  = F_i0l
         params(IDX_J0)  = F_j0l
         params(IDX_NIL) = nijk_l(1)
         params(IDX_NJL) = nijk_l(2)

         datalen = nijk_l(1)*nijk_l(2)*nijk_l(3)
         call rpn_comm_send(params,N_PARAM,RPN_COMM_INTEGER,RPN_COMM_MASTER, &
              tag,F_comm_S,ierr)
         if (datalen > 0) call rpn_comm_send(F_fld_l(1:nijk_l(1), &
              1:nijk_l(2),:),datalen,RPN_COMM_REAL,RPN_COMM_MASTER,TAG, &
              F_comm_S,ierr)

      endif IF_MASTER
      F_istat = RMN_OK !TODO-later: check rpn_comm_send/recv status
      call msg(MSG_DEBUG,'(ptopo) priv_collect_r4_3d [END]')
      !---------------------------------------------------------------------
      return
   end function priv_collect_r4_3d


   !/@*
   function priv_collect_r4_2d(F_fld_c,F_fld_l,F_i0l,F_j0l,F_lni,F_lnj, &
        F_comm_S,F_ismaster_L,F_npe) result(F_istat)
      implicit none
      !@objective Merge MPI data-tiles belonging to the same comm_S
      !@arguments
      real,pointer :: F_fld_c(:,:)    !O, comm_S data
      real,pointer :: F_fld_l(:,:)    !I, local data
      integer,intent(in) :: F_i0l,F_j0l !pos of F_fld_l(1,1) in comm_S
      integer,intent(in) :: F_lni,F_lnj !Size of fld_l w/o halo
      character(len=*),intent(in) :: F_comm_S
      logical,intent(in) :: F_ismaster_L
      integer,intent(in) :: F_npe
      !@return
      integer :: F_istat
      !@revision
      !  2010-03, S. Chamberland - split out of gem 4.1.2
      !@description
      !*@/
      integer :: iproc, istat,ierr, datalen, nijk_l(2)
      integer :: params(N_PARAM)
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG,'(ptopo) priv_collect_r4_2d [BEGIN]')
      call ptopo_init_var()
      nijk_l = ubound(F_fld_l)
      nijk_l(1) = min(F_lni,nijk_l(1))
      nijk_l(2) = min(F_lnj,nijk_l(2))

      IF_MASTER: if (F_ismaster_L) then

         istat = priv_copy_r4_2d(F_fld_c,F_fld_l,F_i0l,F_j0l,F_lni,F_lnj)

         DO_IPROC: do iproc = 1, F_npe - 1

            call rpn_comm_recv(params,N_PARAM,RPN_COMM_INTEGER,iproc,TAG, &
                 F_comm_S,istat,ierr)
            if (params(IDX_NIL) * params(IDX_NJL) > 0) then
               istat = priv_recv_collect_r4_2d(F_fld_c,params(IDX_NIL), &
                    params(IDX_NJL),params(IDX_I0),params(IDX_J0), &
                    F_comm_S,iproc)
            endif

         enddo DO_IPROC

      else !IF_MASTER

         params(IDX_I0)  = F_i0l
         params(IDX_J0)  = F_j0l
         params(IDX_NIL) = nijk_l(1)
         params(IDX_NJL) = nijk_l(2)

         datalen = nijk_l(1)*nijk_l(2)
         call rpn_comm_send(params,N_PARAM,RPN_COMM_INTEGER,RPN_COMM_MASTER, &
              tag,F_comm_S,ierr)
         if (datalen > 0) call rpn_comm_send(F_fld_l(1:nijk_l(1),1:nijk_l(2)), &
              datalen,RPN_COMM_REAL,RPN_COMM_MASTER,TAG,F_comm_S,ierr)

      endif IF_MASTER
      F_istat = RMN_OK  !TODO-later: check rpn_comm_send/recv status
      call msg(MSG_DEBUG,'(ptopo) priv_collect_r4_2d [END]')
      !---------------------------------------------------------------------
      return
   end function priv_collect_r4_2d

   
   !/@*
   function priv_alongX() result(F_alongx_L)
      implicit none
      !@objective Check that the topology is compatible with alongX
      !           assumption done in this module
      !@arguments
      !@return
      logical :: F_alongx_L
      !@revision
      !@description
      !*@/
      logical, save :: isinit_L = .false.
      logical, save :: alongx_L = .false.
      integer :: params_row(ptopo_grid_npex,ptopo_grid_npey), istat
      !---------------------------------------------------------------------
      F_alongx_L = alongx_L
      if (isinit_L) return
      isinit_L =  .true.
      
      if (ptopo_grid_npex == 1 .or. ptopo_grid_npey == 1) then
         alongx_L = .true.
      else
         call rpn_comm_allgather(ptopo_grid_ipex, 1, RPN_COMM_INTEGER, &
              &    params_row, 1, RPN_COMM_INTEGER, RPN_COMM_GRID, istat)
         if (RMN_IS_OK(istat)) then
            alongx_L = (params_row(1,1) < params_row(2,1) .and. &
                 params_row(2,1) <= params_row(ptopo_grid_npex,1) .and. &
                 all(params_row(1,1) == params_row(1,:)))
         endif
      endif
      
      F_alongx_L = alongx_L
      return
   end function priv_alongX

end module ptopo_utils
