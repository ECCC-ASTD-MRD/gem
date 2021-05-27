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
subroutine bmf_splitwrall(nom,gni,gnj,nk,time1,time2, &
                         hgrid,vgrid,dtyp,scat,vecteur)
! SUBROUTINE bmf_splitwrall (nom,gni,gnj,nk,time1,time2, &
! hgrid,vgrid,dtyp,scat,vecteur) L. Corbeil
!
!ARGUMENT
! nom             Variable name
! gni,gnj,nk      Array size
! time1           Timestamp 1 (YYYYMMDD)
! time2           Timestamp 2 (HHMMSSCC)
! hgrid, vgrid    Horizontal and vertical grid descriptor
! dtyp            Data type
! scat            Scatter list indicator
! vecteur         Array to be written
!     
! DESCRIPTION
! Routine which write the array "vecteur" in each files opened 
! by bmf_splitinit with no subdivision at all.
!
  use bmf_modsplit
  implicit none
  integer gni,gnj,nk,time1,time2,hgrid,vgrid,dtyp,scat
  character*4 nom
  integer vecteur((dtyp/10/4)*gni,gnj,nk)
  integer px,py
  integer no_pe,ierr

  integer bmf_write2
  external bmf_write2

  no_pe=0
  call bmf_perturb(nom,vecteur,(dtyp/10/4)*gni,gnj,nk)
  do px=0,bmf_npex-1
  do py=0,bmf_npey-1
      no_pe=no_pe+1
     ierr=bmf_write2(no_pe,nom,gni,1,gni,gnj,1,gnj,nk,1, &
              nk,time1,time2,hgrid,vgrid,dtyp,scat,gni*gnj*nk, &
              vecteur)
     if(ierr.ne.0) then
	write(*,*) "BMF_SPLITWRALL ERROR: error opening ",split_files(no_pe)
     endif
  enddo
  enddo
  return
  end
