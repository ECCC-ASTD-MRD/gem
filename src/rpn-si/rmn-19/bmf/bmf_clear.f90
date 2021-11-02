!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2005  Environnement Canada
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
        subroutine bmf_clear      
!
! SUBROUTINE BMF_clear, L. Corbeil
!
! ARGUMENTS
! 
! DESCRIPTION
! Routine that ends the BMF mode started by bmf_init
! and clears the memory used by BMF_gobe
        use bmf_mod
        implicit none
        type(bmf_liste), pointer :: champ_courant
        champ_courant => liste
!	write(*,*) 'clear',associated(liste),associated(champ_courant)
        do
           if(.not.associated(champ_courant)) EXIT
!           write(*,*) 'dealloc',champ_courant%bmf_champ%nom
           liste => champ_courant%champ_suivant
!           write(*,*) 'avant' ,associated(liste)
           deallocate(champ_courant%bmf_champ%tableau)
           deallocate(champ_courant)
           champ_courant => liste
!           write(*,*) associated(champ_courant),associated(liste)
        enddo
        bmf_started=.false.
        bmf_liste_started=.false.
        end subroutine bmf_clear
