!---------------------------------- LICENCE BEGIN -------------------------------
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
!---------------------------------- LICENCE END ---------------------------------

      subroutine gem_init_appmpi (F_COMM)
      use app_mpmd
      implicit none
      
      integer, intent(OUT) :: F_COMM
      
      include 'mpif.h'

      integer :: GEM_comp_id , INs_comp_id ,&
                 OUTs_comp_id, IRIS_comp_id,&
                 gemgrid_id, checkdmpart_id
      integer :: i,ierr
      integer, dimension(:), allocatable :: gem_and_ios_ids
      type :: participants
         character(len=12) :: name
         integer :: id
      end type participants
      type (participants) :: players(6)
!     
!--------------------------------------------------------------------
!
      players(1)%name= 'gem'
      players(2)%name= 'IN-server'
      players(3)%name= 'OUT-server'
      players(4)%name= 'IRIS'
      players(5)%name= 'gemgrid'
      players(6)%name= 'checkdmpart'
      
      ierr= App_MPMD_Init()
      call App_Start()

      do i=1,6
         players(i)%id= App_MPMD_GetComponentId(players(i)%name)
      end do

      GEM_comp_id   = players(1)%id
      INs_comp_id   = players(2)%id
      OUTs_comp_id  = players(3)%id
      IRIS_comp_id  = players(4)%id
      gemgrid_id    = players(5)%id
      checkdmpart_id= players(6)%id

      if ( (App_MPMD_GetSelfComponentId() == GEM_comp_id) .and. (App_MPMD_GetSelfComponentRank() == 0) ) then
         call App_MPMD_PrintSummary()
      end if
      
      ! Create intracommunicators with all the processes of the listed components

      F_COMM = MPI_COMM_NULL
      if ((gemgrid_id >=0) .or. (checkdmpart_id >=0)) then
         F_COMM= app_mpmd_getselfcomm()
      else

         gem_and_ios_ids = [GEM_comp_id]
         if(INs_comp_id >= 0) then
             gem_and_ios_ids = [gem_and_ios_ids, INs_comp_id]
         endif
         if(OUTs_comp_id >= 0) then
             gem_and_ios_ids = [gem_and_ios_ids, OUTs_comp_id]
         endif

         if(size(gem_and_ios_ids) == 1) then
             F_COMM = App_MPMD_GetSelfComm()
         else
             F_COMM = App_MPMD_GetSharedComm(gem_and_ios_ids, include_pes0only=.false.)
         endif

      endif
!
!--------------------------------------------------------------------
!
      return
      end subroutine gem_init_appmpi
