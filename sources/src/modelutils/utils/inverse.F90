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

!**s/r inverse - compute inverses of F_len  matrices of the same order
!
       subroutine inverse(F_B_8,F_A_8,F_dim,F_len)
!
      implicit none
#include <arch_specific.hf>
!
       integer F_dim, F_len
       real*8 F_A_8(F_dim,F_dim,F_len), F_B_8(F_dim,F_dim,F_len)
!
!author
!     Abdessamad Qaddouri
!
!revision
! v2_10 - Qaddouri A.       - initial version
!
!
!arguments
!  Name        I/O                 Description
!--------------------------------------------------------------
!  F_A_8           I                 matrices to be inverted
!  F_B_8           O                 inverses
!----------------------------------------------------------------
!
!*
       real*8  det_8
       integer kk
!
!     __________________________________________________________________
!     
       if (F_dim .eq. 1) then
        do kk=1,F_len
        F_B_8(1,1,kk) =1/F_A_8(1,1,kk)
        enddo
        endif
       
       if (F_dim .eq. 2) then
       do kk=1,F_len
       det_8= F_A_8(1,1,kk)*F_A_8(2,2,kk)- F_A_8(1,2,kk)*F_A_8(2,1,kk)
         F_B_8(1,1,kk) =   F_A_8(2,2,kk)/det_8
         F_B_8(1,2,kk) = - F_A_8(1,2,kk)/det_8
         F_B_8(2,1,kk) = - F_A_8(2,1,kk)/det_8
         F_B_8(2,2,kk) =   F_A_8(1,1,kk)/det_8
       enddo
       endif
 
       if (F_dim .eq. 3) then
        do kk=1,F_len
        det_8= F_A_8(1,1,kk)*F_A_8(2,2,kk)*F_A_8(3,3,kk)-  &
           F_A_8(1,1,kk)*F_A_8(2,3,kk)*F_A_8(3,2,kk)-F_A_8(2,1,kk) &
           *F_A_8(1,2,kk)*F_A_8(3,3,kk) + F_A_8(2,1,kk)*F_A_8(1,3,kk)*F_A_8(3,2,kk)+F_A_8(3,1,kk) &
           *F_A_8(1,2,kk)*F_A_8(2,3,kk)-F_A_8(3,1,kk)*F_A_8(1,3,kk)*F_A_8(2,2,kk)

         F_B_8(1,1,kk)= (F_A_8(2,2,kk)*F_A_8(3,3,kk)-F_A_8(2,3,kk)*F_A_8(3,2,kk))/det_8
         F_B_8(1,2,kk)=-( F_A_8(1,2,kk)*F_A_8(3,3,kk)-F_A_8(1,3,kk)*F_A_8(3,2,kk))/det_8
         F_B_8(1,3,kk)= (F_A_8(1,2,kk)*F_A_8(2,3,kk)-F_A_8(1,3,kk)*F_A_8(2,2,kk))/det_8
         F_B_8(2,1,kk)= -(F_A_8(2,1,kk)*F_A_8(3,3,kk)-F_A_8(2,3,kk)*F_A_8(3,1,kk))/det_8
         F_B_8(2,2,kk)=(F_A_8(1,1,kk)*F_A_8(3,3,kk)-F_A_8(1,3,kk)*F_A_8(3,1,kk))/det_8
         F_B_8(2,3,kk)=-(F_A_8(1,1,kk)*F_A_8(2,3,kk)-F_A_8(1,3,kk)*F_A_8(2,1,kk))/det_8
         F_B_8(3,1,kk)=(F_A_8(2,1,kk)*F_A_8(3,2,kk)-F_A_8(2,2,kk)*F_A_8(3,1,kk))/det_8
         F_B_8(3,2,kk)=-(F_A_8(1,1,kk)*F_A_8(3,2,kk)-F_A_8(1,2,kk)*F_A_8(3,1,kk))/det_8
         F_B_8(3,3,kk)=(F_A_8(1,1,kk)*F_A_8(2,2,kk)-F_A_8(1,2,kk)*F_A_8(2,1,kk))/det_8
        enddo
        endif

        if (F_dim .eq. 4) then
       do kk=1,F_len
       det_8= F_A_8(1,1,kk)*F_A_8(2,2,kk)*F_A_8(3,3,kk)*F_A_8(4,4,kk) &
               - F_A_8(1,1,kk)*F_A_8(2,2,kk)*F_A_8(3,4,kk)*F_A_8(4,3,kk) &
               -F_A_8(1,1,kk)*F_A_8(3,2,kk)*F_A_8(2,3,kk)*F_A_8(4,4,kk) + &
                F_A_8(1,1,kk)*F_A_8(3,2,kk)*F_A_8(2,4,kk)*F_A_8(4,3,kk) &
               +F_A_8(1,1,kk)*F_A_8(4,2,kk)*F_A_8(2,3,kk)*F_A_8(3,4,kk)- &
               F_A_8(1,1,kk)*F_A_8(4,2,kk)*F_A_8(2,4,kk)*F_A_8(3,3,kk) &
               -F_A_8(2,1,kk)*F_A_8(1,2,kk)*F_A_8(3,3,kk)*F_A_8(4,4,kk)+ &
               F_A_8(2,1,kk)*F_A_8(1,2,kk)*F_A_8(3,4,kk)*F_A_8(4,3,kk)+F_A_8(2,1,kk) &
              *F_A_8(3,2,kk)*F_A_8(1,3,kk)*F_A_8(4,4,kk)-F_A_8(2,1,kk)*F_A_8(3,2,kk) &
              *F_A_8(1,4,kk)*F_A_8(4,3,kk)-F_A_8(2,1,kk)*F_A_8(4,2,kk) &
              *F_A_8(1,3,kk)*F_A_8(3,4,kk)+F_A_8(2,1,kk)*F_A_8(4,2,kk)*F_A_8(1,4,kk)* &
               F_A_8(3,3,kk)+F_A_8(3,1,kk)*F_A_8(1,2,kk)*F_A_8(2,3,kk) &
               *F_A_8(4,4,kk)-F_A_8(3,1,kk)*F_A_8(1,2,kk)*F_A_8(2,4,kk)*F_A_8(4,3,kk) &
              -F_A_8(3,1,kk)*F_A_8(2,2,kk)*F_A_8(1,3,kk)*F_A_8(4,4,kk) &
              +F_A_8(3,1,kk)*F_A_8(2,2,kk)*F_A_8(1,4,kk)*F_A_8(4,3,kk)+F_A_8(3,1,kk)* &
               F_A_8(4,2,kk)*F_A_8(1,3,kk)*F_A_8(2,4,kk)-F_A_8(3,1,kk) &
               *F_A_8(4,2,kk)*F_A_8(1,4,kk)*F_A_8(2,3,kk)-F_A_8(4,1,kk)*F_A_8(1,2,kk) &
              *F_A_8(2,3,kk)*F_A_8(3,4,kk)+F_A_8(4,1,kk)*F_A_8(1,2,kk) &
              *F_A_8(2,4,kk)*F_A_8(3,3,kk)+F_A_8(4,1,kk)*F_A_8(2,2,kk)*F_A_8(1,3,kk) &
              *F_A_8(3,4,kk)-F_A_8(4,1,kk)*F_A_8(2,2,kk)*F_A_8(1,4,kk) &
              *F_A_8(3,3,kk)-F_A_8(4,1,kk)*F_A_8(3,2,kk)*F_A_8(1,3,kk)*F_A_8(2,4,kk) &
               +F_A_8(4,1,kk)*F_A_8(3,2,kk)*F_A_8(1,4,kk)*F_A_8(2,3,kk)

       F_B_8(1,1,kk)= (F_A_8(2,2,kk)*F_A_8(3,3,kk)*F_A_8(4,4,kk)-F_A_8(2,2,kk)* &
         F_A_8(3,4,kk)*F_A_8(4,3,kk)-F_A_8(3,2,kk) &
        *F_A_8(2,3,kk)*F_A_8(4,4,kk)+F_A_8(3,2,kk)*F_A_8(2,4,kk)* &
         F_A_8(4,3,kk)+F_A_8(4,2,kk)*F_A_8(2,3,kk)*F_A_8(3,4,kk)-F_A_8(4,2,kk) &
        *F_A_8(2,4,kk)*F_A_8(3,3,kk))/det_8

       F_B_8(1,2,kk)=-(F_A_8(1,2,kk)*F_A_8(3,3,kk)*F_A_8(4,4,kk)-F_A_8(1,2,kk)* &
         F_A_8(3,4,kk)*F_A_8(4,3,kk)-F_A_8(3,2,kk) &
        *F_A_8(1,3,kk)*F_A_8(4,4,kk)+F_A_8(3,2,kk)*F_A_8(1,4,kk)*F_A_8(4,3,kk)+ &
        F_A_8(4,2,kk)*F_A_8(1,3,kk)*F_A_8(3,4,kk)-F_A_8(4,2,kk) &
        *F_A_8(1,4,kk)*F_A_8(3,3,kk))/det_8

       F_B_8(1,3,kk)=(F_A_8(1,2,kk)*F_A_8(2,3,kk)*F_A_8(4,4,kk)-F_A_8(1,2,kk)* &
        F_A_8(2,4,kk)*F_A_8(4,3,kk)-F_A_8(2,2,kk) &
        *F_A_8(1,3,kk)*F_A_8(4,4,kk)+F_A_8(2,2,kk)*F_A_8(1,4,kk)*F_A_8(4,3,kk)+ &
         F_A_8(4,2,kk)*F_A_8(1,3,kk)*F_A_8(2,4,kk)-F_A_8(4,2,kk) &
        *F_A_8(1,4,kk)*F_A_8(2,3,kk))/det_8

       F_B_8(1,4,kk)=-(F_A_8(1,2,kk)*F_A_8(2,3,kk)*F_A_8(3,4,kk)-F_A_8(1,2,kk)* &
         F_A_8(2,4,kk)*F_A_8(3,3,kk)-F_A_8(2,2,kk) &
        *F_A_8(1,3,kk)*F_A_8(3,4,kk)+F_A_8(2,2,kk)*F_A_8(1,4,kk)*F_A_8(3,3,kk)+ &
         F_A_8(3,2,kk)*F_A_8(1,3,kk)*F_A_8(2,4,kk)-F_A_8(3,2,kk) &
        *F_A_8(1,4,kk)*F_A_8(2,3,kk)) /det_8

       F_B_8(2,1,kk)= -(F_A_8(2,1,kk)*F_A_8(3,3,kk)*F_A_8(4,4,kk)-F_A_8(2,1,kk)* &
         F_A_8(3,4,kk)*F_A_8(4,3,kk)-F_A_8(3,1,kk) &
        *F_A_8(2,3,kk)*F_A_8(4,4,kk)+F_A_8(3,1,kk)*F_A_8(2,4,kk)*F_A_8(4,3,kk)+ &
         F_A_8(4,1,kk)*F_A_8(2,3,kk)*F_A_8(3,4,kk)-F_A_8(4,1,kk) &
        *F_A_8(2,4,kk)*F_A_8(3,3,kk))/det_8

       F_B_8(2,2,kk)=(F_A_8(1,1,kk)*F_A_8(3,3,kk)*F_A_8(4,4,kk)-F_A_8(1,1,kk)* &
         F_A_8(3,4,kk)*F_A_8(4,3,kk)-F_A_8(3,1,kk) &
        *F_A_8(1,3,kk)*F_A_8(4,4,kk)+F_A_8(3,1,kk)*F_A_8(1,4,kk)*F_A_8(4,3,kk)+ &
         F_A_8(4,1,kk)*F_A_8(1,3,kk)*F_A_8(3,4,kk)-F_A_8(4,1,kk) &
        *F_A_8(1,4,kk)*F_A_8(3,3,kk))/det_8

       F_B_8(2,3,kk)=-(F_A_8(1,1,kk)*F_A_8(2,3,kk)*F_A_8(4,4,kk)-F_A_8(1,1,kk)* &
         F_A_8(2,4,kk)*F_A_8(4,3,kk)-F_A_8(2,1,kk) &
        *F_A_8(1,3,kk)*F_A_8(4,4,kk)+F_A_8(2,1,kk)*F_A_8(1,4,kk)*F_A_8(4,3,kk) &
        +F_A_8(4,1,kk)*F_A_8(1,3,kk)*F_A_8(2,4,kk)-F_A_8(4,1,kk) &
        *F_A_8(1,4,kk)*F_A_8(2,3,kk))/det_8

       F_B_8(2,4,kk) = (F_A_8(1,1,kk)*F_A_8(2,3,kk)*F_A_8(3,4,kk)-F_A_8(1,1,kk)* &
         F_A_8(2,4,kk)*F_A_8(3,3,kk)-F_A_8(2,1,kk) &
        *F_A_8(1,3,kk)*F_A_8(3,4,kk)+F_A_8(2,1,kk)*F_A_8(1,4,kk)*F_A_8(3,3,kk) &
        +F_A_8(3,1,kk)*F_A_8(1,3,kk)*F_A_8(2,4,kk)-F_A_8(3,1,kk) &
        *F_A_8(1,4,kk)*F_A_8(2,3,kk))/det_8

       F_B_8(3,1,kk) = (F_A_8(2,1,kk)*F_A_8(3,2,kk)*F_A_8(4,4,kk)-F_A_8(2,1,kk)* &
         F_A_8(3,4,kk)*F_A_8(4,2,kk)-F_A_8(3,1,kk) &
        *F_A_8(2,2,kk)*F_A_8(4,4,kk)+F_A_8(3,1,kk)*F_A_8(2,4,kk)*F_A_8(4,2,kk)+ &
         F_A_8(4,1,kk)*F_A_8(2,2,kk)*F_A_8(3,4,kk)-F_A_8(4,1,kk) &
        *F_A_8(2,4,kk)*F_A_8(3,2,kk))/det_8


       F_B_8(3,2,kk) = -(F_A_8(1,1,kk)*F_A_8(3,2,kk)*F_A_8(4,4,kk)-F_A_8(1,1,kk)* &
         F_A_8(3,4,kk)*F_A_8(4,2,kk)-F_A_8(3,1,kk) &
        *F_A_8(1,2,kk)*F_A_8(4,4,kk)+F_A_8(3,1,kk)*F_A_8(1,4,kk)*F_A_8(4,2,kk)+ &
         F_A_8(4,1,kk)*F_A_8(1,2,kk)*F_A_8(3,4,kk)-F_A_8(4,1,kk) &
        *F_A_8(1,4,kk)*F_A_8(3,2,kk))/det_8

       F_B_8(3,3,kk) = (F_A_8(1,1,kk)*F_A_8(2,2,kk)*F_A_8(4,4,kk)-F_A_8(1,1,kk)* &
         F_A_8(2,4,kk)*F_A_8(4,2,kk)-F_A_8(2,1,kk) &
        *F_A_8(1,2,kk)*F_A_8(4,4,kk)+F_A_8(2,1,kk)*F_A_8(1,4,kk)*F_A_8(4,2,kk)+ &
         F_A_8(4,1,kk)*F_A_8(1,2,kk)*F_A_8(2,4,kk)-F_A_8(4,1,kk) &
        *F_A_8(1,4,kk)*F_A_8(2,2,kk))/det_8

       F_B_8(3,4,kk) = -(F_A_8(1,1,kk)*F_A_8(2,2,kk)*F_A_8(3,4,kk)-F_A_8(1,1,kk)* &
         F_A_8(2,4,kk)*F_A_8(3,2,kk)-F_A_8(2,1,kk) &
        *F_A_8(1,2,kk)*F_A_8(3,4,kk)+F_A_8(2,1,kk)*F_A_8(1,4,kk)*F_A_8(3,2,kk)+ &
         F_A_8(3,1,kk)*F_A_8(1,2,kk)*F_A_8(2,4,kk)-F_A_8(3,1,kk) &
        *F_A_8(1,4,kk)*F_A_8(2,2,kk))/det_8

       F_B_8(4,1,kk) = -(F_A_8(2,1,kk)*F_A_8(3,2,kk)*F_A_8(4,3,kk)-F_A_8(2,1,kk)* &
         F_A_8(3,3,kk)*F_A_8(4,2,kk)-F_A_8(3,1,kk) &
        *F_A_8(2,2,kk)*F_A_8(4,3,kk)+F_A_8(3,1,kk)*F_A_8(2,3,kk)*F_A_8(4,2,kk)+ &
         F_A_8(4,1,kk)*F_A_8(2,2,kk)*F_A_8(3,3,kk)-F_A_8(4,1,kk) &
        *F_A_8(2,3,kk)*F_A_8(3,2,kk))/det_8

       F_B_8(4,2,kk) = (F_A_8(1,1,kk)*F_A_8(3,2,kk)*F_A_8(4,3,kk)-F_A_8(1,1,kk)* &
        F_A_8(3,3,kk)*F_A_8(4,2,kk)-F_A_8(3,1,kk) &
        *F_A_8(1,2,kk)*F_A_8(4,3,kk)+F_A_8(3,1,kk)*F_A_8(1,3,kk)*F_A_8(4,2,kk) &
        +F_A_8(4,1,kk)*F_A_8(1,2,kk)*F_A_8(3,3,kk)-F_A_8(4,1,kk) &
        *F_A_8(1,3,kk)*F_A_8(3,2,kk))/det_8


       F_B_8(4,3,kk) = -(F_A_8(1,1,kk)*F_A_8(2,2,kk)*F_A_8(4,3,kk)-F_A_8(1,1,kk)* &
         F_A_8(2,3,kk)*F_A_8(4,2,kk)-F_A_8(2,1,kk) &
        *F_A_8(1,2,kk)*F_A_8(4,3,kk)+F_A_8(2,1,kk)*F_A_8(1,3,kk)*F_A_8(4,2,kk)+ &
         F_A_8(4,1,kk)*F_A_8(1,2,kk)*F_A_8(2,3,kk)-F_A_8(4,1,kk) &
        *F_A_8(1,3,kk)*F_A_8(2,2,kk))/det_8


       F_B_8(4,4,kk) = (F_A_8(1,1,kk)*F_A_8(2,2,kk)*F_A_8(3,3,kk)-F_A_8(1,1,kk)* &
         F_A_8(2,3,kk)*F_A_8(3,2,kk)-F_A_8(2,1,kk) &
        *F_A_8(1,2,kk)*F_A_8(3,3,kk)+F_A_8(2,1,kk)*F_A_8(1,3,kk)*F_A_8(3,2,kk) &
        +F_A_8(3,1,kk)*F_A_8(1,2,kk)*F_A_8(2,3,kk)-F_A_8(3,1,kk) &
        *F_A_8(1,3,kk)*F_A_8(2,2,kk))/det_8

        enddo
        endif
!
!     __________________________________________________________________
!     
        return
        end
