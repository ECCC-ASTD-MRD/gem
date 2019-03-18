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

      SUBROUTINE RPN_COMM_xch_halon(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,gni,npol_row,n)
      use rpn_comm
      implicit none
!
!	Calls RPN_COMM_xch_halo for n-length words
!       n=1 -> real*4, n=2 -> real*8 and so on...
!        include 'rpn_comm.h'
!        include 'mpif.h'

      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      integer gni,npol_row,n
      logical periodx,periody
      integer g(n*(minx-1)+1:n*maxx,miny:maxy,nk)


      call RPN_COMM_xch_halo(g,n*(minx-1)+1,n*maxx,miny,maxy,  &
                   n*ni,nj,nk,n*halox,haloy,periodx,periody,   &
                   n*gni,npol_row)
      
        return
        end

        SUBROUTINE RPN_COMM_xch_haloxn(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,g2,minx2,maxx2,miny2,maxy2,gni,npol_row,n)
      use rpn_comm
      implicit none
!
!	Calls RPN_COMM_xch_halox for n-length words
!       n=1 -> real*4, n=2 -> real*8 and so on...
!
!        include 'rpn_comm.h'
!        include 'mpif.h'

      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      integer minx2,maxx2,miny2,maxy2
      integer gni,npol_row,n
      logical periodx,periody
      integer g(n*(minx-1)+1:n*maxx,miny:maxy,nk)
      integer g2(n*(minx2-1)+1:n*maxx2,miny2:maxy2)



        call RPN_COMM_xch_halox(g,n*(minx-1)+1,n*maxx,miny,maxy,   &
                   n*ni,nj,nk,n*halox,haloy,periodx,periody,       &
                   g2,n*(minx2-1)+1,n*maxx2,miny2,maxy2,n*gni,npol_row)


        return
        end
