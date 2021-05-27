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
 Subroutine bmf_splitoptions(string, array)
 use bmf_modsplit
 character*(*) string
 integer array(*)


 if((string(1:4).eq.'halo').or.(string(1:4).eq.'HALO')) then
    bmf_haloileft=array(1)
    bmf_haloiright=array(2) 
    bmf_halojleft=array(3) 
    bmf_halojright=array(4) 
    bmf_ghaloileft=array(5) 
    bmf_ghaloiright=array(6) 
    bmf_ghalojleft=array(7) 
    bmf_ghalojright=array(8) 
    return
 else if ((string(1:4).eq.'grid').or.(string(1:4).eq.'GRID')) then
    bmf_nig=array(1)
    bmf_njg=array(2)
    return
 else
    write(*,*) 'BMF_SPLITOPTIONS ERROR: unrecognized option'
    stop
 endif
end subroutine bmf_splitoptions



