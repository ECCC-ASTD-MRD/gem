!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */

      SUBROUTINE RPN_COMM_xch_halox(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,g2,minx2,maxx2,miny2,maxy2,gni,npol_row)
      use rpn_comm
      implicit none
!
!	exchange a halo with neighbours
!
      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      integer minx2,maxx2,miny2,maxy2
      integer gni,npol_row
      logical periodx,periody
      integer g(minx:maxx,miny:maxy,nk)
      real g_adj(minx:maxx,miny:maxy,nk)
      integer g2(minx2:maxx2,miny2:maxy2,nk)
      real g2_adj(minx2:maxx2,miny2:maxy2,nk)
      pointer (g_adj_,g_adj)
      pointer (g2_adj_,g2_adj)
      logical adjoint
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
      integer i,j,k,borneminx,bornemaxx,borneminy,bornemaxy
      integer nimax,ierr
!
      adjoint=.false.

 1      continue
      g_adj_=loc(g)
      g2_adj_=loc(g2)
      if (adjoint) then
        call RPN_COMM_adj_halo(g2,minx2,maxx2,miny2,maxy2,  &
                   ni,nj,nk,halox,haloy,periodx,periody,    &
                   gni,npol_row)
        nimax = (gni + pe_nx - 1)/pe_nx
        if (min(pe_mey+1,pe_ny-pe_mey).le.npol_row) then
          do k=1,nk
          do j=1,nj
!VDIR NODEP
          do i=1,ni
            g_adj(i,j,k)=g2_adj(i+pe_mex*nimax,j,k)+g_adj(i,j,k)
          enddo
          enddo
          enddo      	
        else
          do k=1,nk
          do j=1,nj
          do i=1,ni
            g(i,j,k)=g2(i,j,k)
          enddo
          enddo
          enddo      	
        endif
        else
        do k=1,nk
        do j=1,nj
        do i=1,ni
          g2(i,j,k)=g(i,j,k)
        enddo
        enddo
        enddo      

        call RPN_COMM_xch_halo(g2,minx2,maxx2,miny2,maxy2,  &
                   ni,nj,nk,halox,haloy,periodx,periody,    &
                   gni,npol_row)
        endif

!
      return

      entry RPN_COMM_adj_halox(g,minx,maxx,miny,maxy,      &
                   ni,nj,nk,halox,haloy,periodx,periody,   &
      	           g2,minx2,maxx2,miny2,maxy2,gni,npol_row)
      adjoint=.true.
      goto 1
      end
