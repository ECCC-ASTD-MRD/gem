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

!**s/r int_cubvec_lag - YY cubic interpolation for vector fields (local PE)

       subroutine int_cubvec_lag ( FF, Fu_8, Fv_8,Imx1,Imy1, Imx2, Imy2, &
                                 geomgxu_8,geomgy_8, geomgx_8,geomgyv_8, &
                                 Minx,Maxx,Miny,Maxy,Nk,Xi_8,Yi_8,NLEN,  &
                                 s1_8,s2_8)
       implicit none
#include <arch_specific.hf>

       integer Nk,NLEN,Minx,Maxx,Miny,Maxy
       integer Imx1(NLEN),Imy1(NLEN),Imx2(NLEN),Imy2(NLEN)
       real    FF(NLEN*Nk)
       real*8  Fu_8(Minx:Maxx,Miny:Maxy,Nk),Fv_8(Minx:Maxx,Miny:Maxy,Nk), &
               geomgxu_8(Minx:Maxx),geomgy_8(Miny:Maxy), &
               geomgx_8(Minx:Maxx),geomgyv_8(Miny:Maxy), &
               Xi_8(NLEN),Yi_8(NLEN), s1_8(NLEN),s2_8(NLEN)
!
!author
!           Abdessamad Qaddouri - October 2009
!revision
! v4_60 - Qaddouri A.   - initial version
! v4_70 - Desgagne M.   - introduce mono_L
!
       integer k,i,Imx,Jmx,Imy,Jmy
       real*8  WX1,WX2,WX3,WX4
       real*8  WY1,WY2,WY3,WY4
       real*8  X1,XX,X2,X3,X4
       real*8  WWX1,WWX2,WWX3,WWX4
       real*8  WWY1,WWY2,WWY3,WWY4
       real*8  YY,y1,y2,y3,y4,Fi
       real*8  fx1,fx2,fx3,fx4
       real*8  fy1,fy2,fy3,fy4
       real*8  FFu_8, FFv_8
!
!----------------------------------------------------------------------
!
!$omp parallel private (k,i,imx,jmx,imy,jmy, &
!$omp          WX1,WX2,WX3,WX4,      &
!$omp          WY1,WY2,WY3,WY4,      &
!$omp          X1,XX,X2,X3,X4,       &
!$omp          WWX1,WWX2,WWX3,WWX4,  &
!$omp          WWY1,WWY2,WWY3,WWY4,  &
!$omp          YY,y1,y2,y3,y4,Fi,    &
!$omp          fx1,fx2,fx3,fx4,FFu_8,&
!$omp          fy1,fy2,fy3,fy4,FFv_8   )
!$omp do
   Do 200 i=1,NLEN
!For U
       Imx  = iMx1(i)
       Jmx  = iMy1(i)

       X1  = Geomgxu_8(imx)
       X2  = Geomgxu_8(imx+1) - X1
       X3  = Geomgxu_8(imx+2) - X1
       X4  = Geomgxu_8(imx+3) - X1
       XX  = Xi_8(i)   - X1

       WX1  = (XX-X2)*(XX-X3)*(XX-X4)/(-X2*X3*X4)
       WX2  =  XX    *(XX-X3)*(XX-X4)/(X2*(X2-X3)*(X2-X4))
       WX3  =  XX    *(XX-X2)*(XX-X4)/(X3*(X3-X2)*(X3-X4))
       WX4  =  XX    *(XX-X2)*(XX-X3)/(X4*(X4-X2)*(X4-X3))

       Y1  = Geomgy_8(Jmx)
       Y2  = Geomgy_8(Jmx+1) - Y1
       Y3  = Geomgy_8(Jmx+2) - Y1
       Y4  = Geomgy_8(Jmx+3) - Y1
       YY  = Yi_8(i)   - Y1

       WWX1  = (yy-y2)*(yy-y3)*(yy-y4)/(-y2*y3*y4)
       WWX2  =  yy    *(yy-y3)*(yy-y4)/(y2*(y2-y3)*(y2-y4))
       WWX3  =  yy    *(yy-y2)*(yy-y4)/(y3*(y3-y2)*(y3-y4))
       WWX4  =  yy    *(yy-y2)*(yy-y3)/(y4*(y4-y2)*(y4-y3))

!For V
       Imy  = iMx2(i)
       Jmy  = iMy2(i)

       X1  = Geomgx_8(imy)
       X2  = Geomgx_8(imy+1) - X1
       X3  = Geomgx_8(imy+2) - X1
       X4  = Geomgx_8(imy+3) - X1
       XX  = Xi_8(i)   - X1

       WY1  = (XX-X2)*(XX-X3)*(XX-X4)/(-X2*X3*X4)
       WY2  =  XX    *(XX-X3)*(XX-X4)/(X2*(X2-X3)*(X2-X4))
       WY3  =  XX    *(XX-X2)*(XX-X4)/(X3*(X3-X2)*(X3-X4))
       WY4  =  XX    *(XX-X2)*(XX-X3)/(X4*(X4-X2)*(X4-X3))

       Y1  = Geomgyv_8(Jmy)
       Y2  = Geomgyv_8(Jmy+1) - Y1
       Y3  = Geomgyv_8(Jmy+2) - Y1
       Y4  = Geomgyv_8(Jmy+3) - Y1
       YY  = Yi_8(i)   - Y1

       WWY1  = (yy-y2)*(yy-y3)*(yy-y4)/(-y2*y3*y4)
       WWY2  =  yy    *(yy-y3)*(yy-y4)/(y2*(y2-y3)*(y2-y4))
       WWY3  =  yy    *(yy-y2)*(yy-y4)/(y3*(y3-y2)*(y3-y4))
       WWY4  =  yy    *(yy-y2)*(yy-y3)/(y4*(y4-y2)*(y4-y3))

   Do 300 k=1,Nk
!For U
       Fx1 = WX1*Fu_8(imx,jmx,k) + WX2*Fu_8(imx+1,jmx,k) + WX3*Fu_8(imx+2,jmx,k)+ &
             WX4*Fu_8(imx+3,jmx,k)

       Fx2 = WX1*Fu_8(imx,jmx+1,k) + WX2*Fu_8(imx+1,jmx+1,k) + WX3*Fu_8(imx+2,jmx+1,k)+ &
             WX4*Fu_8(imx+3,jmx+1,k)

       Fx3 = WX1*Fu_8(imx,jmx+2,k) + WX2*Fu_8(imx+1,jmx+2,k) + WX3*Fu_8(imx+2,jmx+2,k)+ &
             WX4*Fu_8(imx+3,jmx+2,k)

       Fx4 = WX1*Fu_8(imx,jmx+3,k) + WX2*Fu_8(imx+1,jmx+3,k) + WX3*Fu_8(imx+2,jmx+3,k)+ &
             WX4*Fu_8(imx+3,jmx+3,k)


       FFu_8  = WWX1*fx1+WWX2*fx2+WWX3*fx3+WWX4*fx4


!For V
       Fy1 = WY1*Fv_8(imy,jmy,k) + WY2*Fv_8(imy+1,jmy,k) + WY3*Fv_8(imy+2,jmy,k)+ &
             WY4*Fv_8(imy+3,jmy,k)

       Fy2 = WY1*Fv_8(imy,jmy+1,k) + WY2*Fv_8(imy+1,jmy+1,k) + WY3*Fv_8(imy+2,jmy+1,k)+ &
             WY4*Fv_8(imy+3,jmy+1,k)

       Fy3 = WY1*Fv_8(imy,jmy+2,k) + WY2*Fv_8(imy+1,jmy+2,k) + WY3*Fv_8(imy+2,jmy+2,k)+ &
             WY4*Fv_8(imy+3,jmy+2,k)

       Fy4 = WY1*Fv_8(imy,jmy+3,k) + WY2*Fv_8(imy+1,jmy+3,k) + WY3*Fv_8(imy+2,jmy+3,k)+ &
             WY4*Fv_8(imy+3,jmy+3,k)


       FFv_8  = WWY1*fy1+WWY2*fy2+WWY3*fy3+WWY4*fy4

       FF((i-1)*Nk+k)=real(s1_8(i)*FFu_8 + s2_8(i)*FFv_8)
  300  continue
  200  continue
!$omp enddo
!$omp end parallel
!
!----------------------------------------------------------------------
!
       return
       end
