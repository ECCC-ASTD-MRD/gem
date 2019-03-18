!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/

        SUBROUTINE RPN_COMM_bcst_world(start,end,mode)
!	Luc Corbeil, 2000-02-04
!	routine qui fait un broadcast selon le mode indique
!
!	mode 0 = broadcast standard, de Pe0 vers MPI_COMM_WOLRD
!	 
!	Sinon, un message d'erreur est retourne.
!
!	ATTENTION !!! WARNING
!	Le code assume que la longueur de "end" est la meme qu'un
!	entier, donc les resultats seront errones si "end" est
!	un real*8, complex*16 par exemple 
!
!	This routine assume that the length of "end" is the same
!	than an integer, so it will return wrong results if
!	"end" is a real*8 or a complex*16 
!
!
        use rpn_comm
        implicit none
        integer start(2), end(2), mode
	integer buff(2),ierr
	integer length
!*
!        include 'mpif.h'

	if (mode.eq.0) then ! broadcast habituel
	  length=(loc(end(2))-loc(start(1)))/(loc(buff(2))-loc(buff(1)))
	  call MPI_BCAST(start,length,MPI_INTEGER,0,pe_defcomm,ierr)

!	elseif () then
!	quand on aura d'autres modes...
!
	else
	  write(rpn_u,*) 'RPN_COMM_bcst_world:: Mode non reconnu: ',mode
	  stop
	endif

	return
	end  SUBROUTINE RPN_COMM_bcst_world
