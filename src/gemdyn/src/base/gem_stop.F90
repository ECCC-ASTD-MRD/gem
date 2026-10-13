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
!----------------------------------LICENCE END ---------------------------------

      subroutine gem_stop ()
      use svri_mod
      use svro_mod
      use inp_mod
      implicit none

      include 'mpif.h'
      integer tag,send(5),err
!
!     ---------------------------------------------------------------
!
      if (Inp_window/=0) call MPI_Win_free (Inp_window,err)
      
      send= -1
      if (INs_comm_L) then
         tag= 801
         call MPI_send (send,1,MPI_INTEGER,&
                        INs_1on1,tag,INs_COMM,err )
      endif
      if (OUTs_comm_L) then
         tag= 601
         call MPI_send ( send,size(send),MPI_INTEGER,&
                         OUTs_1on1,tag,OUTs_COMM,err )
      endif
      
      call MPI_finalize (err)
      stop
!
!     ---------------------------------------------------------------
!
      end subroutine gem_stop
