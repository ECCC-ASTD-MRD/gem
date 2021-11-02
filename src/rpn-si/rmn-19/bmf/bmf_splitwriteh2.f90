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
subroutine bmf_splitwriteh2(comm_split_func,nom,minx,maxx,miny,maxy,nk,ni0,&
     nj0,gni,gnj,time1,time2,dtyp,vecteur,holename)
! SUBROUTINE bmf_splitwriteh(nom,gni,gnj,nk,time1,time2, &
! hgrid,vgrid,dtyp,scat,vecteur,holename) L. Corbeil
!
! ARGUMENT
! comm_split_func split_function to be called (normaly RPN_COMM_split)
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
! Routine which write the array "vecteur" in all files opened by
! bmf_splitinit. The array is splitted according to parameters gni,
! gnj, those of bmf_splitinit and bmf_splithalo. It also handles 
! "holes" defined by bmf_splithole

  use bmf_modsplit
  implicit none
  integer, intent(IN) :: gni,gnj,nk,time1,time2,dtyp
  integer, intent(IN) :: ni0, nj0, minx, maxx, miny,maxy
  character*4, intent(IN) :: nom,holename
  integer, target, intent(IN) :: vecteur((dtyp/10/4),minx:maxx,miny:maxy,nk)
  integer, allocatable, dimension(:,:) :: p_vecteur
  integer px,py,ierr,nxl,nxlmax,halox,nx0
  integer nyl,nylmax,haloy,ny0,no_pe,ind
  logical fill,trouve_trou
  integer j,k,size,hgrid,vgrid,scat,minx2,maxx2,miny2,maxy2
  integer hole_i0, hole_isize, hole_in, hole_t1, hole_t2
  integer hole_j0, hole_jsize, hole_jn
  integer dimx1, dimx2,dimx3, dimx4, dimx5,dimx6
  integer hole(4)

  type(bmf_holelist), pointer :: current_hole

  integer  comm_split_func, bmf_write2
  external comm_split_func, bmf_write2
  size=dtyp/10/4
  fill=.false.
  halox=0
  haloy=0
  ierr=0
  no_pe=0
  vgrid=0
  hgrid=0
  scat=1
!
! Find hole position and size
!
  current_hole => holelist
  trouve_trou=.false.
  if(associated(current_hole)) then
     do
        if(.not.associated(current_hole)) exit
        if(current_hole%bmf_trou%bmf_holename==holename.and.&
             (current_hole%bmf_trou%bmf_holet1.eq.time1).and.&
             (current_hole%bmf_trou%bmf_holet2.eq.time2)) then
           trouve_trou=.true.
           hole_t1=current_hole%bmf_trou%bmf_holet1
           hole_t2=current_hole%bmf_trou%bmf_holet2
           hole_i0=current_hole%bmf_trou%bmf_holei0
           hole_isize=current_hole%bmf_trou%bmf_holeisize
           hole_j0=current_hole%bmf_trou%bmf_holej0
           hole_jsize=current_hole%bmf_trou%bmf_holejsize
           hole_in= hole_i0+hole_isize-1
           hole_jn= hole_j0+hole_jsize-1
           hole(1)=hole_i0
           hole(2)=hole_isize
           hole(3)=hole_j0
           hole(4)=hole_jsize
          exit
        endif
        current_hole => current_hole%trou_suivant
     enddo
     if(.not.trouve_trou) then
        write(*,*) "BMF_SPLITWRITEH: hole not found, returning"
        return
     endif

  endif

  call bmf_perturb(nom,vecteur,(dtyp/10/4)*(maxx-minx+1),(maxy-miny+1),nk)

  do px=0,bmf_npex-1
     ierr= comm_split_func(px,bmf_npex,gni,minx2,maxx2,nxl,nxlmax, &
          halox,nx0,fill)
     if(ierr.ne.0) then
        write(*,*) 'BMF_SPLITWRITEH: error comm_split_func, abort'
        stop
     endif
     do py=0,bmf_npey-1
        no_pe=no_pe+1
        ierr= comm_split_func(py,bmf_npey,gnj,miny2,maxy2,nyl,nylmax, &
             haloy,ny0,fill)
        
        if(ierr.ne.0) then
           write(*,*) 'BMF_SPLITWRITEH: error comm_split_func, abort'
           stop
        endif
        write(*,*) nxlmax*nylmax
        if(.not.allocated(p_vecteur)) then
           allocate(p_vecteur(size,nxlmax*nylmax),stat=ierr)
           if(ierr.ne.0) then
              write(*,*) 'BMF_SPLITWRITEH: error allocate, abort'
              stop
           endif
        endif
        dimx1=max(nx0,ni0)
        dimx2=min(nx0+nxl-1,ni0+hole_i0-1)
        dimx3=max(nx0,ni0+hole_i0-1)
        dimx4=min(nx0+nxl-1,ni0+hole_in-1)
        dimx5=max(nx0,ni0+hole_in)
        dimx6=min(nx0+nxl-1,ni0+gni-1)
        write(*,*) no_pe,'dimx',dimx1,dimx2,dimx3,dimx4,dimx5,dimx6
        
        do k=1,nk
           ind=1
           do j=ny0,ny0+nyl-1
              if((j.ge.nj0).and.(j.lt.nj0+hole_j0-1)) then
                 ! copie toute la rangee
                 if(dimx6-dimx1.gt.0) then
                    write(*,*) 'RAN',j,dimx1,dimx6,size,ind
                    p_vecteur(1:size,ind:ind+dimx6-dimx1)= &
                         vecteur(1:size,dimx1:dimx6,j,k)
                    ind=ind+dimx6-dimx1+1
                    write(*,*) ind
                 endif
              else if((j.ge.nj0+hole_j0-1).and. &
                   (j.lt.nj0+hole_jn)) then
                 ! copie les morceaux aux extremites
                 if(dimx2-dimx1.gt.0) then
                    p_vecteur(1:size,ind:ind+dimx2-dimx1-1)= &
                         vecteur(1:size,dimx1:dimx2-1,j,k)
                    ind=ind+dimx2-dimx1
                 endif
                 write(*,*) 'AA',px,py,dimx6,dimx5
                 if(dimx6-dimx5.gt.0) then
                    p_vecteur(1:size,ind:ind+dimx6-dimx5) = &
                         vecteur(1:size,dimx5:dimx6,j,k)
                    ind=ind+dimx6-dimx5+1
                 endif
              else if((j.ge.nj0+hole_jn).and. &
                   (j.lt.nj0+gnj)) then
                 ! copie toute la rangee
                 if(dimx6-dimx1.gt.0) then
                    p_vecteur(1:size,ind:ind+dimx6-dimx1)= &
                         vecteur(1:size,dimx1:dimx6,j,k)
                    ind=ind+dimx6-dimx1+1
                 endif
              endif
           enddo
! On se positionne au dernier point
           ind=ind-1
           write(*,*)  p_vecteur(1,1:22)
           hgrid=nx0
           vgrid=ny0
           ierr = bmf_write2(no_pe,nom,nxl,ni0,ni0+gni-1,nyl,nj0,nj0+gnj-1,nk,k, &
                k,time1,time2,hgrid,vgrid,dtyp,scat,ind, &
                p_vecteur)
           if(ierr.ne.0) then
              write(*,*) "BMF_SPLITWRITEH ERROR: error writing in ",split_files(no_pe)
           endif
        enddo
     enddo
  enddo
  deallocate(p_vecteur)
  return
end subroutine bmf_splitwriteh2
