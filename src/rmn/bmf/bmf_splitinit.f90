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
subroutine bmf_splitinit(npex,npey,path,prefix,date,hour,min,sec,num,numlen,unit)
!SUBROUTINE bmf_splitinit(npex,npey,path,prefix,date,&
! & hour,min,sec,num,numlen,unit), L. Corbeil
!
! ARGUMENTS
!
! IN
! npex        Number of processors for x axis
! npey        Number of processors for y axis
! path        Path to the splitted files   
! prefix      Prefix (character*2)
! date        date (integer, YYYYMMDD)
! hour        hours
! min         minutes
! sec         seconds
! num         multiple of unit;: optional
! numlen      number of digits for num, optional
! unit        1 character to represent type of units, optional
!
! Prefix, date, hour, min, sec, num, numlen and unit are passed directly
! to prog_filename: please refer to librmn documentation.
!     
! DESCRIPTION
! Routine which starts bmf_split mode: it permits to write bmf
! fields by dividing them for a later use in a parallel context.
! Files corresponding to the return value of prog_filename
! will be open and used for writing by bmf_splitwr* routines.  
  use bmf_modsplit
  implicit none
  integer npex,npey,date,hour,min,sec,num,numlen
  character*2 unit, prefix
  character* (*) path
  integer ierr,prog_filename,i,j,indice,petot,fnom,longueur
  external prog_filename,fnom,longueur
  petot=npex*npey
  bmf_npex=npex
  bmf_npey=npey
  if(petot.le.0) then
      write(*,*) 'BMF_SPLITINIT: npex*npey.le.0, abort'
      stop
  endif
  if(.not.allocated(split_files)) then
     allocate(split_files(petot))
     allocate(split_unit(petot))
  else
     write(*,*) 'BMF_SPLITINIT: split mode already started: use SPLITEND first'
  endif
  indice=0
  do i=0,npex-1
      do j=0,npey-1
         indice=indice+1
         ierr=prog_filename(split_files(indice),prefix,date,hour,min,sec, &
                     i,j,num,numlen,unit)
      enddo
  enddo
  if(ierr.ne.0) then
      write(*,*) 'BMF_SPLITINIT: Error prog_filename, abort'
      stop
  endif
  do i=1,petot
!      write(*,*) i,split_unit(i),split_files(i)
      split_files(i)=path(1:longueur(path))//'/'//split_files(i)
      split_unit(i)=0
!      write(*,*) i,split_unit(i),split_files(i)
      ierr=FNOM(split_unit(i),split_files(i),'SEQ/UNF',0)
!      write(*,*) i,split_unit(i),split_files(i)
  enddo
  bmf_haloileft=0
  bmf_haloiright=0
  bmf_halojleft=0
  bmf_haloiright=0
  bmf_ghaloileft=0
  bmf_ghaloiright=0
  bmf_ghalojleft=0
  bmf_ghaloiright=0
  bmf_nig=-1
  bmf_njg=-1
  return 
  end subroutine
