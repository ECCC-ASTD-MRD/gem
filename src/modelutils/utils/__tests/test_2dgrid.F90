

subroutine test_2dgrid()
   use iso_c_binding
   implicit none
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"

   external :: test_get_doms, test_get_npxy
!!$   integer, external :: rpn_comm_init_multigrid, rpn_comm_topo, rpn_comm_create_2dgrid

   integer,parameter :: GNI = 153
   integer,parameter :: GNJ = 180
   integer, parameter :: HALO = 0
   logical, parameter :: ALONGX = .true.
   logical, parameter :: FILL = .false. !.true.

   integer :: istat, mini, maxi, lni, lnimax, li0, minj, maxj, lnj, lnjmax, lj0, rpncomm_gridid, mydomain, igrid, myproc, numproc, npex, npey, ngrids
   ! ---------------------------------------------------------------------
   call rpn_comm_mydomain(test_get_doms, mydomain)

   ngrids = 1
   npex = 0
   npey = 0
   igrid = rpn_comm_init_multigrid(test_get_npxy, myproc, numproc, &
        npex, npey, ngrids)

   istat = rpn_comm_topo(GNI, mini, maxi, lni, lnimax, HALO, li0, &
        ALONGX, FILL)
   istat = rpn_comm_topo(GNJ, minj, maxj, lnj, lnjmax, HALO, lj0, &
        .not.ALONGX, FILL)

   maxi = lni
   rpncomm_gridid = rpn_comm_create_2dgrid(GNI, GNJ, mini, maxi, minj, maxj)

   print *, myproc, ":", rpncomm_gridid, ":", mini, maxi, lni, ":", minj, maxj, lnj
   ! ---------------------------------------------------------------------
   return
end subroutine test_2dgrid


subroutine test_get_doms(F_ndomains, F_offset, F_istat)
   implicit none
   integer, intent(out) :: F_ndomains,F_offset,F_istat
   ! ---------------------------------------------------------------------
   F_ndomains = 1
   F_offset = 0
   F_istat = 0
   ! ---------------------------------------------------------------------
   return
end subroutine test_get_doms


subroutine test_get_npxy(F_npex, F_npey)
   implicit none
   integer,intent(out) :: F_npex,F_npey
   !---------------------------------------------------------------------
   F_npex = 2
   F_npey = 1
   !---------------------------------------------------------------------
   return
end subroutine test_get_npxy
