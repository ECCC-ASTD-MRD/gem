!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
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
!**     s/r ez_llwfgdw - converts from grid wind components to std
!       meteorological speed and direction
   subroutine ez_llwfgdw(z1,z2,xlon,li,lj,grtyp,ig1,ig2,ig3,ig4)
   implicit none
   integer li,lj
   real z1(li,lj), z2(li,lj), xlon(li,lj)
   character*1 grtyp
   integer ig1,ig2,ig3,ig4
!
!       AUTHOR     Yves Chartier              April 1991
!
!       REVISION
!       Documentation update             Jul 1994
!
!       LANGUAGE   Fortran 77
!
!       OBJECT (gdw2llw)
!       Converts grid wind components to meteorological wind speed and direction,
!       independent of geographical projection.
!       On output, a 270 degree wind is a westerly wind, meaning U is +ve and
!       V is zero.
!
!       FILES
!
!       ARGUMENTS
!       NAMES          I/O  TYPE  A/S        DESCRIPTION
!       z1           I/O   R     A         On entry, contains the U grid component
!       of wind. On exit, contains wind speed.
!       z2           I/O   R     A         On entry, contains the V component
!       of wind. On exit, contains wind dir.
!       xlon           I     R     A         Array of longitudes
!       li             I     I     S         X dimension of the grid
!       lj             I     I     S         Y dimension of the grid
!       grtyp          I     C     S         Grid type
!       ig1-2-3-4      I     I     S         Grid descriptors
!
!       IMPLICIT
!       include
!
!       MODULES
   external cigaxg
!       *
!
!-------------------------------------------------------------
!
!
!
!       * 1.866025=(1+sin60),   6.371e+6=earth radius in meters.
!
!       rdtodg = 180/pie, dgtord = pie/180
!
#include "pi.cdk"

   real xg1, xg2, xg3, xg4, dgrw
   character*1 gtypout
   real xg(20)
   integer nxg
   real delx,dely,alpha
   real uuu,vvv
   integer i,j,k,n,ntot
   real spd0, dir0
   integer npts,ier
   real x1(2*li*lj),y1(2*li*lj),lat(2*li*lj)


!  les #define qui suivent rendent le code plus lisible
#define uu   z1
#define vv   z2
#define spd  z1
#define dir  z2

   if (grtyp .eq. '!') then
      call ez_lamb_llwfgdw(uu,vv,xlon,li,lj,grtyp,ig1,ig2,ig3,ig4,x1,y1,lat)
      return
   endif


   if (grtyp.eq. 'N')then
      call cigaxg(grtyp,xg1,xg2,xg3,xg4,ig1,ig2,ig3,ig4)
      do j=1,lj
         do i=1,li
            spd0=sqrt(uu(i,j)*uu(i,j)+vv(i,j)*vv(i,j))
            if (spd0.eq. 0.0)then
               dir0= 0.0
            else
              if (uu(i,j).eq. 0.0)then
                 if (vv(i,j).ge. 0.0)then
                    dir0= xlon(i,j)+xg4-90.0
                 else
                    dir0= xlon(i,j)+xg4+90.0
                 endif
              else
                 dir0=xlon(i,j)+xg4-rdtodg*atan2(vv(i,j),uu(i,j))
              endif
            endif
            dir0=amod(amod(dir0,360.0)+360.0,360.0)
            spd(i,j)=spd0
            dir(i,j)=dir0
           enddo
      enddo
      return
   endif

   if (grtyp.eq. 'S')then
      call cigaxg(grtyp,xg1,xg2,xg3,xg4,ig1,ig2,ig3,ig4)
      do j=1,lj
         do i=1,li
            spd0=sqrt(uu(i,j)*uu(i,j)+vv(i,j)*vv(i,j))
            if (spd0.eq. 0.0)then
               dir0= 0.0
            else
               if (uu(i,j).eq. 0.0)then
                  if (vv(i,j).ge. 0.0)then
                     dir0= 90.0 - xlon(i,j)+xg4
                  else
                     dir0= 270.0 - xlon(i,j)+xg4
                 endif
               else
                 dir0=180.0-xlon(i,j)+xg4-rdtodg*atan2(vv(i,j),uu(i,j))
               endif
            endif
            dir0=amod(amod(dir0,360.0)+360.0,360.0)
            spd(i,j)=spd0
            dir(i,j)=dir0
         enddo
      enddo
      return
   endif

   if (grtyp.eq.'A'.or.grtyp.eq.'B'.or.grtyp.eq.'G'.or.grtyp.eq.'L')then
      do j=1,lj
         do i=1,li
            spd0=sqrt(uu(i,j)*uu(i,j)+vv(i,j)*vv(i,j))
            if (spd0.eq. 0.0)then
               dir0= 0.0
            else
               if (uu(i,j).eq. 0.0)then
                  if (vv(i,j).ge. 0.0)then
                     dir0= 180.0
                  else
                     dir0= 0.0
                  endif
               else
                  dir0=270.0 - rdtodg*atan2(vv(i,j),uu(i,j))
               endif
            endif
            dir0=amod(amod(dir0,360.0)+360.0,360.0)
            spd(i,j)=spd0
            dir(i,j)=dir0
         enddo
      enddo
      return
   endif
   write(6, 600) grtyp
 600   format('0',' erreur, bad grid type (llwfgdw) - grtyp = ', a11)
   return
   end

