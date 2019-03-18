#ifdef TEST
	program testit
	real grida(9,9),xa(9),ya(9),gridb(5,5),xb(5),yb(5)
        real xa_8(9),ya_8(9),xb_8(5),yb_8(5)
	real*8 cxa(9),cxb(9),cxc(9),cxd(9)
	real*8 cya(9),cyb(9),cyc(9),cyd(9)
	integer indxx(9),indxy(9)
        character*8 hint

        hint = 'LAG2D'
	do j=1,9
	do i=1,9
	  grida(i,j)=0
	  xa(i)=i+2
          xa_8(i)=xa(i)
	  ya(j)=j+2
          ya_8(j)=ya(j)
	  indxx(i)=-1
	  indxy(i)=-1
	enddo
	enddo
	do j=9,1,-1
	  PRINT 101,(GRIDA(I,J),i=1,9)
	ENDDO
	print *,'xa=',xa
	print *,'ya=',ya
	do j=1,5
	do i=1,5
	  gridb(i,j)=i*2+1
	  xb(i)=i*2+1
	  yb(j)=j*2+1
          xb_8(i)=xb(i)
          yb_8(j)=yb(j)
	enddo
	enddo
	do j=5,1,-1
	  PRINT 101,(GRIDB(I,J),i=1,5)
	ENDDO
	print *,'xb=',xb
	print *,'yb=',yb
	call grid_to_grid_coef(xa_8,9,xb_8,5,indxx,cxa,cxb,cxc,cxd,hint)
	print *,'indxx=',indxx
	print 101,cxa
	print 101,cxb
	print 101,cxc
	print 101,cxd
	call grid_to_grid_coef(ya_8,9,yb_8,5,indxy,cya,cyb,cyc,cyd,hint)
	print *,'indxy=',indxy
	print 101,cya
	print 101,cyb
	print 101,cyc
	print 101,cyd
	call grid_to_grid(grida,9,9,gridb,5,5,xa,ya,xb,yb)
	print *,' out of grid_to_grid'
	do j=9,1,-1
	  PRINT 101,(GRIDA(I,J),i=1,9)
	ENDDO
      call grid_to_grid_interp(grida,9,9,gridb,5,5,indxx,indxy,
     %                    cxa,cxb,cxc,cxd,cya,cyb,cyc,cyd,hint)
	print *,' out of grid_to_grid_interp'
	do j=9,1,-1
	  PRINT 101,(GRIDA(I,J),i=1,9)
	ENDDO
101	format(2x,11f6.1)
	stop
	end
*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
***S/R grid_to_grid - SPLINE INTERPOLATION.
*
      SUBROUTINE grid_to_grid(FI,IFI,JFI,F,IF,JF,XI,YI,X,Y)
*
*AUTHOR   - D. ROBERTSON (original routine D1INT2)
*
*REVISION  M. VALIN - 2004 rearranged original code of D1INT2 to make it stand alone
*
*ARGUMENTS
*   OUT   - FI   - INTERPOLATED FIELD.
*   IN    - IFI  - X-DIMENSION OF FIELD FI.
*         - JFI  - Y-DIMENSION OF FIELD FI.
*         - F    - ORIGINAL FIELD.
*         - IF   - X-DIMENSION OF FIELD F.
*         - JF   - Y-DIMENSION OF FILED F.
*         - XI   - X LOCATION OF THE VALUES OF THE INTERPOLATED FIELD FI.
*         - YI   - Y LOCATION OF THE VALUES OF THE INTERPOLATED FIELD FI.
*         - X    - X LOCATION OF THE VALUES OF THE ORIGINAL FIELD F.
*         - Y    - Y LOCATION OF THE VALUES OF THE ORIGINAL FIELD F.
*
*NOTES    - this will only work if grid spacing is CONSTANT
*           SPLINE interpolation is used
*
*---------------------------------------------------------------------------
*
      implicit none
      integer if,jf,ifi,jfi,k,k1,kk
      integer i,j,l,l1,ll
      real we,we1,we2,wn,ww,wm,wz,wz1,wz2,z,ze,zf,zg,zl,zz
      REAL FI     (IFI,JFI)
      REAL F      (IF ,JF )
      REAL FX     (IF ,JF )
      REAL FY     (IF ,JF )
      REAL FXY    (IF ,JF )
      REAL XI     (IFI)
      REAL YI     (JFI)
      REAL X      (IF )
      REAL Y      (JF )
*         - FX   - ARRAY OF SIZE (IF,JF) THAT CONTAINS COMPUTED
*                  PARTIAL DERIVATIVE OF FIELD F WITH RESPECT TO X.
*         - FY   - PARTIAL DERIVATIVE WITH RESPECT TO Y.
*         - FXY  - PARTIAL SECOND DERIVATIVE WITH RESPECT TO X AND Y.
*         - HX   - GRID-LENGTH along X
*         - HY   - GRID-LENGTH along Y
*         - ZA   - WORKING VECTOR OF LENGTH (JFI).
*         - ZB   - WORKING VECTOR OF LENGTH (JFI).
*         - ZC   - WORKING VECTOR OF LENGTH (JFI).
*         - ZD   - WORKING VECTOR OF LENGTH (JFI).
      REAL hx,hx2
      REAL hy,hy2
      REAL ZA     (JFI)
      REAL ZB     (JFI)
      REAL ZC     (JFI)
      REAL ZD     (JFI)
*
*----------------------------------------------------------------------
*
      hx=1.0/(x(2)-x(1))
      hx2=.5*hx
      hy=1.0/(y(2)-y(1))
      hy2=.5*hy
!  compute df/dx
      do j=1,jf
        fx(1,j)=(f(2,j)-f(1,j))*hx
        fx(if,j)=(f(if,j)-f(if-1,j))*hx
        do i=2,if-1
          fx(i,j)=(f(i+1,j)-f(i-1,j))*hx2
        enddo
      enddo
!  compute df/dy
      do j=2,jf-1
        do i=1,if
          fy(i,1)=(f(i,2)-f(i,1))*hy
          fy(i,jf)=(f(i,jf)-f(i,jf-1))*hy
        enddo
        do i=1,if
          fy(i,j)=(f(i,j+1)-f(i,j-1))*hy2
        enddo
      enddo
!  compute df/dxdy
      do j=1,jf
        fxy(1,j)=(fy(2,j)-fy(1,j))*hx
        fxy(if,j)=(fy(if,j)-fy(if-1,j))*hx
        do i=2,if-1
          fxy(i,j)=(fy(i+1,j)-fy(i-1,j))*hx2
        enddo
      enddo
      LL=2
      DO 15 J=1,JFI
      DO 22 L=LL,JF
        L1=L-1
        IF(YI(J).LE.Y(L)) GO TO 23
   22 CONTINUE
*
   23 continue
      LL=L1+1
      WN=YI(J)-Y(L1)
      WE=WN*hy
      WE1=1.-WE
      WE2=WE1*WE1
      WW=2.*WE
      ZA(J)=WE2*WN
      ZB(J)=WE1*WN*WE
      ZC(J)=WE2*(1.+WW)
      ZD(J)=WE*WE*(3.-WW)
   15 CONTINUE
*
      KK=2
      DO 11 I=1,IFI
      DO 12 K=KK,IF
        K1=K-1
        IF(XI(I).LE.X(K)) then
          GO TO 13
        endif
   12 CONTINUE
*
   13 continue
      KK=K1+1
      WM=XI(I)-X(K1)
      WZ=WM*hx
      WZ1=1.-WZ
      WZ2=WZ1*WZ1
      ZZ=2.*WZ
      ZE=WZ2*WM
      ZF=WZ1*WM*WZ
      ZG=WZ2*(1.+ZZ)
      ZL=WZ*WZ*(3.-ZZ)
      LL=2
      DO 115 J=1,JFI
      DO 122 L=LL,JF
        L1=L-1
        IF(YI(J).LE.Y(L)) then
          GO TO 123
        endif
  122 CONTINUE
*
  123 continue
      LL=L1+1
      Z=ZA(J)*(ZE*FXY(K1,L1)-ZF*FXY(K,L1)+ZG*FY(K1,L1)+ZL*FY(K,L1))
      Z=Z-ZB(J)*(ZE*FXY(K1,L)-ZF*FXY(K,L)+ZG*FY(K1,L)+ZL*FY(K,L))
      Z=Z+ZC(J)*(ZE*FX(K1,L1)-ZF*FX(K,L1)+ZG*F(K1,L1)+ZL*F(K,L1))
      Z=Z+ZD(J)*(ZE*FX(K1,L)-ZF*FX(K,L)+ZG*F(K1,L)+ZL*F(K,L))
      FI(I,J)=Z
  115 CONTINUE
*
   11 CONTINUE
*
*--------------------------------------------------------------------
*
      RETURN
      END
#endif
*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
***S/R grid_to_grid_interp 2D interpolator
      subroutine grid_to_grid_interp(FI,IFI,JFI,F,IF,JF,indxx,indxy,
     %           cxa,cxb,cxc,cxd,cya,cyb,cyc,cyd,interp)
      implicit none
      character* (*) interp
      integer if,jf,ifi,jfi
      integer i,j
      REAL FI     (IFI,JFI)
      REAL F      (IF ,JF )
      real*8 cxa(ifi),cxb(ifi),cxc(ifi),cxd(ifi)
      real*8 cya(jfi),cyb(jfi),cyc(jfi),cyd(jfi)
      integer indxx(ifi),indxy(jfi)
*
*ARGUMENTS
*   OUT   - FI   - INTERPOLATED FIELD.
*   IN    - IFI  - X-DIMENSION OF FIELD FI.
*         - JFI  - Y-DIMENSION OF FIELD FI.
*         - F    - ORIGINAL FIELD.
*         - IF   - X-DIMENSION OF FIELD F.
*         - JF   - Y-DIMENSION OF FILED F.
*         - INDXX  PRECOMPUTED TERM FOR INTERPOLATION
*         - INDXY  PRECOMPUTED TERM FOR INTERPOLATION
*         - CXA    PRECOMPUTED TERM FOR INTERPOLATION
*         - CXB    PRECOMPUTED TERM FOR INTERPOLATION
*         - CXC    PRECOMPUTED TERM FOR INTERPOLATION
*         - CXD    PRECOMPUTED TERM FOR INTERPOLATION
*         - CYA    PRECOMPUTED TERM FOR INTERPOLATION
*         - CYB    PRECOMPUTED TERM FOR INTERPOLATION
*         - CYC    PRECOMPUTED TERM FOR INTERPOLATION
*         - CYD    PRECOMPUTED TERM FOR INTERPOLATION
*         - INTERP order of interpolation 'CUB_LAG'/'LINEAR'/'NEAREST'
*
*NOTES
*     grid_to_grid_coef is used to compute interpolation constants that are non
*     data dependent but grid dependent
*     one call to grid_to_grid_coef is needed for each (X and Y) axis
**

      real *8 ta,tb,tc,td
      integer ix,jy

      if (interp.eq.'CUB_LAG') then
        do j=1,jfi
          jy=indxy(j)
          do i=1,ifi
            ix=indxx(i)
            ta=f(ix,jy  )*cxa(i)+f(ix+1,jy  )*cxb(i)+f(ix+2,jy  )
     $                   *cxc(i)+f(ix+3,jy  )*cxd(i)
            tb=f(ix,jy+1)*cxa(i)+f(ix+1,jy+1)*cxb(i)+f(ix+2,jy+1)
     $                   *cxc(i)+f(ix+3,jy+1)*cxd(i)
            tc=f(ix,jy+2)*cxa(i)+f(ix+1,jy+2)*cxb(i)+f(ix+2,jy+2)
     $                   *cxc(i)+f(ix+3,jy+2)*cxd(i)
            td=f(ix,jy+3)*cxa(i)+f(ix+1,jy+3)*cxb(i)+f(ix+2,jy+3)
     $                   *cxc(i)+f(ix+3,jy+3)*cxd(i)
            fi(i,j)=cya(j)*ta+cyb(j)*tb+cyc(j)*tc+cyd(j)*td
          enddo
        enddo
      endif
*
      if (interp.eq.'LINEAR') then
        do j=1,jfi
          jy=indxy(j)
          do i=1,ifi
            ix=indxx(i)
            ta=f(ix,jy  )*cxa(i) + f(ix+1,jy  )*(1.0d0-cxa(i))
            tb=f(ix,jy+1)*cxa(i) + f(ix+1,jy+1)*(1.0d0-cxa(i))
            fi(i,j)=cya(j)*ta + (1.0d0-cya(j))*tb
          enddo
        enddo
      endif
*
      if (interp.eq.'NEAREST') then
        do j=1,jfi
          jy=indxy(j)
          do i=1,ifi
            ix=indxx(i)
            fi(i,j)=f(ix,jy  )
          enddo
        enddo
      endif
*
      return
      end
*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
***S/R grid_to_grid_coef pre-calculate constants for 2D interpolator
      subroutine grid_to_grid_coef(xi,ifi,x,if,index,cxa,cxb,cxc,cxd,
     %                             interp)
      implicit none
      character* (*) interp
      integer if,ifi,index(ifi)
      real*8 xi(ifi),x(if)
      real*8 cxa(ifi),cxb(ifi),cxc(ifi),cxd(ifi)
*ARGUMENTS
*   IN    - XI   - LOCATION OF THE VALUES OF THE INTERPOLATED FIELD
*         - IFI  - DIMENSION OF INTERPOLATED FIELD AND OUTPUTS
*         - X    - LOCATION OF THE VALUES OF THE ORIGINAL FIELD
*         - IF   - DIMENSION OF ORIGINAL FIELD
*         - INTERP order of interpolation 'CUB_LAG'/'LINEAR'/'NEAREST'
*   OUT   - INDEX  PRECOMPUTED TERM FOR INTERPOLATION
*         - CXA    PRECOMPUTED TERM FOR INTERPOLATION
*         - CXB    PRECOMPUTED TERM FOR INTERPOLATION
*         - CXC    PRECOMPUTED TERM FOR INTERPOLATION
*         - CXD    PRECOMPUTED TERM FOR INTERPOLATION
**
      real*8 mid
      parameter (mid = 0.5d0)

      integer ii,i
      real *8 da,db,dc,dd,xa,xb,xc,xd

      real *8 triprd,za,zb,zc,zd,xxi
      triprd(za,zb,zc,zd)=(za-zb)*(za-zc)*(za-zd)
      if (interp.eq.'CUB_LAG') then
        ii=4
        xa=x(1)
        xb=x(2)
        xc=x(3)
        xd=x(4)
        da=1.0/triprd(xa,xb,xc,xd)
        db=1.0/triprd(xb,xa,xc,xd)
        dc=1.0/triprd(xc,xa,xb,xd)
        dd=1.0/triprd(xd,xa,xb,xc)

        do i=1,ifi
          do while( xi(i).gt.xc .and. ii.lt.if)
            ii=ii+1
            xa=xb
            xb=xc
            xc=xd
            xd=x(ii)
            da=1.0/triprd(xa,xb,xc,xd)
            db=1.0/triprd(xb,xa,xc,xd)
            dc=1.0/triprd(xc,xa,xb,xd)
            dd=1.0/triprd(xd,xa,xb,xc)
          end do
          xxi=xi(i)
          index(i)=ii-3
          cxa(i)=da*triprd(xxi,xb,xc,xd)
          cxb(i)=db*triprd(xxi,xa,xc,xd)
          cxc(i)=dc*triprd(xxi,xa,xb,xd)
          cxd(i)=dd*triprd(xxi,xa,xb,xc)
        enddo
      endif
*
      if (interp.eq.'LINEAR') then
        ii=2
        xa=x(1)
        xb=x(2)
        da=xb-xa
        do i=1,ifi
          do while( xi(i).gt.xb .and. ii.lt.if)
            ii=ii+1
            xa=xb
            xb=x(ii)
            da=xb-xa
          end do
          xxi= xi(i)
          index(i)=ii-1
          cxa(i) = (xb-xxi) / da
        enddo
      endif
*
      if (interp.eq.'NEAREST') then
        ii=2
        xa=x(1)
        xb=x(2)
        da=xb-xa
        do i=1,ifi
          do while( xi(i).gt.xb .and. ii.lt.if)
            ii=ii+1
            xa=xb
            xb=x(ii)
            da=xb-xa
          end do
          xxi=xi(i)
          index(i)=ii-1
          db = 1.0d0 - (xb-xxi) / da
          if ( db .gt. mid ) index(i) = ii
        enddo
      endif
*
      return
      end
*
