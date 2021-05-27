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
subroutine bmf_splitwritexy2(comm_split_func,nom,gni,gnj,nk,kstart,kend,time1,time2, &
                         hgrid,vgrid,dtyp,scat,vecteur)
! SUBROUTINE bmf_splitwritexy(nom,gni,gnj,nk,kstart,kend,time1,time2,&
!  hgrid,vgrid,dtyp,scat,vecteur) L. Corbeil
! 
! ARGUMENT
! comm_split_func split_function to be called (normaly RPN_COMM_split)
! nom             Variable name
! gni,gnj,nk      Size of the array
! kstart, kend    Starting and ending position of the slice
! time1           Timestamp 1 (YYYYMMDD)
! time2           Timestamp 2 (HHMMSSCC)
! hgrid, vgrid    Horizontal and vertical grid descriptor
! dtyp            Data type
! scat            Scatter list indicator
! vecteur         Array to be written
!     
! DESCRIPTION
! Routine which writes the slice between kstart and kend of the array
! "vecteur" in all files opened by bmf_splitinit. The subdivision is
! guided by parameters gni, gnj, those of bmf_splitinit and those of
! bmf_splithalo.
  use bmf_modsplit
  implicit none
  integer gni,gnj,nk,time1,time2,hgrid,vgrid,dtyp,scat
  integer kstart,kend
  character*4 nom
  integer, target :: vecteur(dtyp/10/4,gni,gnj,kend-kstart+1)
  integer, pointer, dimension(:,:,:,:) :: p_vecteur
  integer px,py,ierr,minx,maxx,nxl,nxlmax,halox,nx0
  integer il,ir,jl,jr,gil,gir,gjl,gjr,vgni,vgnj
  integer miny,maxy,nyl,nylmax,haloy,ny0,no_pe
  integer vgni_bmf,vgnj_bmf
  logical fill
  integer size
  integer  comm_split_func, bmf_write2
  external comm_split_func, bmf_write2
  size=dtyp/10/4
  fill=.false.
  halox=0
  haloy=0
  ierr=0
  no_pe=0
  nullify(p_vecteur)

  il=bmf_haloileft
  ir=bmf_haloiright
  jl=bmf_halojleft
  jr=bmf_haloiright
  gil=bmf_ghaloileft
  gir=bmf_ghaloiright
  gjl=bmf_ghalojleft
  gjr=bmf_ghaloiright

  vgni=gni-gir-gil
  vgnj=gnj-gjl-gjr
  if(bmf_nig.eq.-1) bmf_nig=vgni
  if(bmf_njg.eq.-1) bmf_njg=vgnj
  vgni_bmf=bmf_nig-gir-gil
  vgnj_bmf=bmf_njg-gir-gil

  call bmf_perturb(nom,vecteur,(dtyp/10/4)*gni,gnj,kend-kstart+1)
  if((vgni.gt.bmf_nig).or.(vgnj.gt.bmf_njg)) then
    write(*,*) 'BMF_SPLITWRITEXY: error, trying to split bigger array'
  else if((vgni.lt.bmf_nig-2).or.(vgnj.lt.bmf_njg-2)) then
    write(*,*) 'BMF_SPLITWRITEXY: error, trying to split smaller array'
  else if((vgni.ne.bmf_nig).or.(vgnj.ne.bmf_njg)) then
    write(*,*) 'BMF_SPLITWRITEXY: warning, the array size is sligthly'
    write(*,*) '                less than expected: splitting anyway'
  endif

  do px=0,bmf_npex-1
  do py=0,bmf_npey-1
     no_pe=no_pe+1
     ierr= comm_split_func(px,bmf_npex,vgni_bmf,minx,maxx,nxl,nxlmax, &
                      halox,nx0,fill)
     if(ierr.ne.0) then
          write(*,*) 'BMF_SPLITWRITEXY: error comm_split_func, abort'
          stop
     endif
     if(px.eq.0) then
       if(bmf_npex.eq.1) then
         nxl=gni
       else
         nxl=nxl+gil+ir
       endif
     else if(px.eq.bmf_npex-1) then
       nxl=nxl+gir+il
       nxl=nxl-(bmf_nig-vgni)
       nx0=nx0+gil-il
     else
       nxl=nxl+ir+il
       nx0=nx0+gil-il
     endif
     ierr= comm_split_func(py,bmf_npey,vgnj_bmf,miny,maxy,nyl,nylmax, &
                      haloy,ny0,fill)
     if(ierr.ne.0) then
          write(*,*) 'BMF_SPLITWRITEXY: error comm_split_func, abort'
          stop
     endif

     if(py.eq.0) then
        if(bmf_npey.eq.1) then
          nyl=gnj
        else
          nyl=nyl+gjl+jr
        endif
     else if(py.eq.bmf_npey-1) then
         nyl=nyl+gjr+jl
         nyl=nyl-(bmf_njg-vgnj)
         ny0=ny0+gjl-jl
     else
         nyl=nyl+jr+jl
         ny0=ny0+gjl-jl
     endif
!     allocate(p_vecteur(nxl,nyl,nk))
     p_vecteur => vecteur(1:size,nx0:nx0+nxl-1,ny0:ny0+nyl-1,:)
     
     ierr = bmf_write2(no_pe,nom,nxl,1,nxl,nyl,1,nyl,nk,kstart, &
           kend,time1,time2,hgrid,vgrid,dtyp,scat,nxl*nyl*(kend-kstart+1), &
              p_vecteur)
     if(ierr.ne.0) then
        write(*,*) "BMF_SPLITWRITEXY ERROR: error opening ",split_files(no_pe)
     endif
     nullify(p_vecteur)
  enddo
  enddo
  return
  end
