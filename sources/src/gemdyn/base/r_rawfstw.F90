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

      subroutine r_rawfstw2 (rs,minx,maxx,miny,maxy,mink,maxk,varname,i0,in,j0,jn,k0,kn, &
                                                      stepno,dt,date0,ig1,ig2,out_file)
      implicit none
#include <arch_specific.hf>
!
      integer minx,maxx,miny,maxy,mink,maxk,stepno,dt,date0,i0,in,j0,jn,k0,kn,ig1,ig2
      real rs(minx:maxx,miny:maxy,mink:maxk)
      character(len=*) varname,out_file
!
      integer fnom,fclos
      integer id_unit,err,ip1,ip2,i,j,nx,ny,cnt
      real wk ((in-i0+1)*(jn-j0+1))
!
      id_unit=0
      err = fnom (id_unit, out_file, 'rnd', 0)
      call fstouv (id_unit,'rnd')
!
      nx = in-i0+1
      ny = jn-j0+1
      ip2 =(stepno*dt)/3600
      do ip1=k0,kn
         cnt=0
         do j=j0,jn
         do i=i0,in
            cnt=cnt+1
            wk(cnt) = rs(i,j,ip1)
         end do
         end do
         if ((ig1 > -1) .and. (ig2 > -1) ) then
            call fstecr (wk,wk,-32,id_unit, date0, int(dt), stepno, nx,ny,1, &
                    ip1,ip2,stepno,'P',varname,'rslt2','Z',ig1,ig2,0,0,1,.false.)
         else
            call fstecr (wk,wk,-32,id_unit, date0, int(dt), stepno, nx,ny,1, &
                    ip1,ip2,stepno,'P',varname,'rslt2','X',1,1,1,1,1,.false.)
         end if
      end do

      call fstfrm (id_unit)
      err = fclos (id_unit)
!
      return
      end

      subroutine r_rawfst_pos (xp,yp,ni,nj,ip1,ip2,ig1,ig2,ig3,ig4,out_file)
      implicit none
#include <arch_specific.hf>
!
      integer ni,nj,ip1,ip2,ig1,ig2,ig3,ig4
      real xp(ni),yp(nj)
      character(len=*) out_file
!
      integer fnom,fclos
      integer id_unit,err
!
      id_unit=0
      err = fnom (id_unit, out_file, 'rnd', 0)
      call fstouv (id_unit,'rnd')
!
      call fstecr (xp,xp,-32,id_unit, 0, 0, 0, ni,1,1, &
                   ip1,ip2,0,'X','>>','POS_X','E', &
                   ig1,ig2,ig3,ig4,1,.false.)
      call fstecr (yp,yp,-32,id_unit, 0, 0, 0, 1,nj,1, &
                   ip1,ip2,0,'X','^^','POS_Y','E', &
                   ig1,ig2,ig3,ig4,1,.false.)

      call fstfrm (id_unit)
      err = fclos (id_unit)
!
      return
      end

