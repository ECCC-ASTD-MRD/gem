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

      logical function sol_numa_space (F_checkparti_L)
      use ISO_C_Binding
      use dyn_fisl_options
      use glb_ld
      use ldnh
      use lun
      use trp
      use numa
      use sol
      use ptopo
      implicit none

      logical, intent(in) :: F_checkparti_L

      logical alongY_L
      integer :: i, mpx,irest, err
      integer(KIND=MPI_ADDRESS_KIND) :: msize

      integer, dimension(4,0:Numa_cores_per_socket-1) :: gindx,socki

      ! Pointer to the shared-memory region as returned by numa_space
      ! Note its integer type
      integer, dimension(:), pointer :: tapon 

      ! Pointer to the shared memory region as a REAL64, three-dimesnional
      ! array, used for assigning the final arrays with proper lower bounds
      real(kind=REAL64), dimension(:,:,:), contiguous, pointer :: temp_ptr

      ! C-pointer as a translation variable, used like casts to/from (void *)
      type(C_PTR) :: cptr

      !! Working variables to compute array bounds
      integer jbound, sol_sock_nj, diviser_gni
!
!----------------------------------------------------------------------
!
      sol_numa_space= .false.
      
      if (.not.Numa_uniform_L) then
         if (Lun_out>0) then
            write (Lun_out, '(/3x,"WRONG CHOICE of NPEY= ",i4," for this configuration that uses: ",i4," cores per socket")') Ptopo_npey,Numa_cores_per_socket
            write (Lun_out, '( 3x,"BEST CHOICE WHEN NPEY is a divider or a multiple of the number of cores per socket")')
         endif
         return
      endif
      
      if (Ptopo_npey > 1) then
         alongY_L = (Ptopo_gindx(1,1) == Ptopo_gindx(1,2))
      else
         alongY_L = .true.
      endif

      sol_numa_space = alongY_L .and. Ptopo_alongY_L

      if (sol_numa_space .and. .not.F_checkparti_L) then

         gindx = 0
         gindx(1,Numa_sockrank) = ldnh_j0
         gindx(2,Numa_sockrank) = ldnh_nj
         gindx(3,Numa_sockrank) = trp_12sn0
         gindx(4,Numa_sockrank) = trp_12sn
         
         socki = 0
         call MPI_Allreduce(gindx,socki,4*Numa_cores_per_socket,&
                         MPI_INTEGER, MPI_BOR,Numa_sockcomm,err)
                               
         Sol_miny= 2*G_nj ; Sol_maxy= -2*G_nj
         do i=0,Numa_cores_per_socket-1
            if (socki(2,i)>0) then
               Sol_miny=min(Sol_miny,socki(1,i))
               Sol_maxy=max(Sol_maxy,socki(1,i)+socki(2,i)-1)
            endif
         end do
         Sol_mink= 2*G_nk ; Sol_maxk= -2*G_nk
         do i=0,Numa_cores_per_socket-1
            if (socki(4,i)>0) then
               Sol_mink=min(Sol_mink,socki(3,i))
               Sol_maxk=max(Sol_maxk,socki(3,i)+socki(4,i)-1)
            endif
         end do
         Sol_sock_nk= max(Sol_maxk-Sol_mink+1,0)
         if (Sol_sock_nk==0) then
            Sol_mink=-10
            Sol_maxk=-10
         endif

         ! If previously associated, clear the pointers to the shared-
         ! memory region
         nullify (Sol_fft,Sol_a,Sol_b,Sol_c)
         
         ! Matrix size: G_ni * (numa_ny + 2) * socket_nk elements
         ! As-used in sol_fft_numa, each shared-memory array has bounds
         ! of (Sol_mink:Sol_maxk) and (1:G_ni) in k and i respectively.
         ! Sol_fft (F_fft) has j-bounds of (Sol_miny-1:Sol_maxy+1),
         ! Sol_c (F_c) has j-bounds of Sol_miny-1:Sol_maxy
         ! Sol_{a,b} (F_{a,b}) have j-bounds of Sol_miny:Sol_maxy

         ! Overall, this gives a shared-memory size of:
         sol_sock_nj = (Sol_maxy - Sol_miny + 1)
         msize = max((4*sol_sock_nj + 3)*G_ni*Sol_sock_nk,1)
        
         ! We allocate the shared memory region via the numa_space procedure.
         ! This procedure takes a pointer to an integer array (which will be
         ! reassociated) and the required number of integer-sized words.
         ! REAL64 values are twice the size (8 bytes) of integer words (4
         ! bytes), so we must use a scale factor:


         call numa_space ( tapon, msize*2, err)
         
         if (err<0) then
            sol_numa_space= .false.
            if (Lun_out>0) write (Lun_out, '(/" WARNING -- UNABLE to use the ONE-TRANSPOSE NUMA SOLVER")')
            return
         endif
        
         ! To perform the proper pointer assignment, we must first have
         ! a view of this region as a REAL64 array:
         cptr = C_LOC(tapon) ! Get the memory address (void *)
         call c_f_pointer(cptr, temp_ptr, & ! Associate that address with temp_ptr
                        [sol_sock_nk,G_ni,4*sol_sock_nj+3]) ! Overall extents

         ! Sol_fft is now associated to the powest portion of this array.
         ! Since the association now exists excusively in the Fortran world, we
         ! can use proper lower bounds

         ! Note that Sol_fft has padding in j both above and below the solved region
         jbound = 1 ! Lowest valid j-index of temp_ptr
         Sol_fft(Sol_mink:,1:,Sol_miny-1:) => temp_ptr(:,:,jbound:(jbound + sol_sock_nj + 2 - 1))
         Sol_fft = 0.0d0
         jbound = jbound + sol_sock_nj + 2 ! Update j-bound

         ! Assign Sol_c next, which has one row of padding in j
         Sol_c(Sol_mink:,1:,Sol_miny-1:) => temp_ptr(:,:,jbound:(jbound + sol_sock_nj + 1 - 1))
         Sol_c = 0.0d0
         jbound = jbound + sol_sock_nj + 1

         ! Sol_a and Sol_b have the same shapes and no padding
         Sol_a(Sol_mink:,1:,Sol_miny:) => temp_ptr(:,:,jbound:(jbound + sol_sock_nj - 1))
         Sol_a = 0.0d0
         jbound = jbound + sol_sock_nj
         Sol_b(Sol_mink:,1:,Sol_miny:) => temp_ptr(:,:,jbound:(jbound + sol_sock_nj - 1))
         Sol_b = 0.0d0
         jbound = jbound + sol_sock_nj

         diviser_gni= Numa_active_cores_per_socket
         mpx= mod( Numa_sockrank, diviser_gni )
         Sol_ldni = G_ni / diviser_gni
         irest  = G_ni  - Sol_ldni * diviser_gni
         Sol_istart = mpx * Sol_ldni + 1
         if ( mpx < irest ) then
            Sol_ldni   = Sol_ldni + 1
            Sol_istart = Sol_istart + mpx
         else
            Sol_istart = Sol_istart + irest
         endif
         Sol_iend= Sol_istart + Sol_ldni - 1

      endif
!     
!     ---------------------------------------------------------------
!
      return
      end function sol_numa_space
      
