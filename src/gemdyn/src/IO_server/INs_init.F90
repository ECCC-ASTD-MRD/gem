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

!**s/r INs_init

      subroutine INs_init ()
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use IOs
      use INs
      implicit none

#include <clib_interface_mu.hf>
      character(len=64) dumc_S
      logical :: bool
      integer :: n1,n2,dim,tag,item(4),ierr
!     
!--------------------------------------------------------------------
!
      if (INs_1o1_L) then
         tag=121
         call MPI_recv ( item, size(item),&
               MPI_INTEGER,INs_gem1o1,tag,INs_GEM_COMM,&
               MPI_STATUSES_IGNORE,ierr)
      endif
      call MPI_bcast (item,size(item),MPI_INTEGER,0,&
                      MY_WORLD_COMM,ierr)
      INs_maxreqs = item(1)
      INs_CMCdate0= item(2)
      INs_Iau_listfst_L= TRANSFER(item(3), bool)
      INs_Spn_listfst_L= TRANSFER(item(4), bool)
      
      call datf2p(INs_runstrt_S,INs_CMCdate0)

      Path_input_S= IOs_path_input_S
                      
      if ( clib_getenv ('INS_MAXNKA', dumc_S) > 0 ) then
         read(dumc_S,*) INs_maxNKA
      else
         INs_maxNKA = 100
      endif
      if (INs_1o1_L) then
         tag= tag+1
         call MPI_send (INs_maxNKA,1,MPI_INTEGER,&
               INs_gem1o1,tag,INs_GEM_COMM,ierr )
      endif
      
! prepare INs_host communicator for later shared memory allocation
      call MPI_Comm_split_type(MY_WORLD_COMM, MPI_COMM_TYPE_SHARED,0,&
                               MPI_INFO_NULL, INs_host, ierr)
      call MPI_COMM_rank (INs_host,INs_hostmyproc ,ierr)

      INs_nid = G_ni+2*G_halox
      INs_njd = G_nj+2*G_haloy
      INs_hord= INs_nid*INs_njd

! Memory allocation here is set to a constant maximum value
! which must exactly match the client allocation. It cannot be
! dynamic because of the non-blocking nature of the isend/irecv
      allocate (SRL(INs_maxreqs))
      allocate (Nes_cBUF(INs_maxreqs*2),&
                Iau_cBUF(INs_maxreqs*2),&
                Spn_cBUF(INs_maxreqs*2) )

      dim= 2*INs_maxreqs + INs_maxNKA*INs_maxreqs
      allocate (Nes_iBUF(dim),Iau_iBUF(dim),Spn_iBUF(dim))

      allocate (Nes_VGD_tbl(10*INs_maxNKA),&
                Iau_VGD_tbl(10*INs_maxNKA),&
                Spn_VGD_tbl(10*INs_maxNKA) )
                
      allocate (INs_Nes_isend(2*client_pelocal+3),&
                INs_Iau_isend(2*client_pelocal+3),&
                INs_Spn_isend(2*client_pelocal+3))

      INs_Nes_isend= MPI_REQUEST_NULL
      INs_Iau_isend= MPI_REQUEST_NULL
      INs_Spn_isend= MPI_REQUEST_NULL
      
      n1= model_gindx(2,1)-model_gindx(1,1)+1+2*G_haloy
      n2= model_gindx(4,1)-model_gindx(3,1)+1+2*G_haloy      
      dim= n1*n2*(INs_maxNKA*6 + 1)
      allocate (GZbuf(dim,client_pestart:client_peend))
      dim= n1*n2*INs_maxNKA*INs_maxreqs
      allocate (Nesbuf(dim,client_pestart:client_peend))
      allocate (Iaubuf(dim,client_pestart:client_peend))
      allocate (Spnbuf(dim,client_pestart:client_peend))
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_init

