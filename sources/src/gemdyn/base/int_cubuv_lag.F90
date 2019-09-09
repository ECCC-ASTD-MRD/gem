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

!**s/r int_cubuv_lag - YY cubic interpolation for fields (local PE)

      subroutine int_cubuv_lag ( FF, Fu_8, Fv_8,Imx,Imy, geomgx_8,geomgy_8, &
                                 Minx,Maxx,Miny,Maxy,Nk,Xi_8,Yi_8,NLEN,     &
                                 s1_8,s2_8,s3_8,s4_8)
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Nk,NLEN,Minx,Maxx,Miny,Maxy
      integer, intent(in) :: Imx(NLEN),Imy(NLEN)
      real, intent(out) :: FF(NLEN*Nk*2)
      real*8, intent(in) :: Fu_8(Minx:Maxx,Miny:Maxy,Nk),Fv_8(Minx:Maxx,Miny:Maxy,Nk), &
               geomgx_8(Minx:Maxx),geomgy_8(Miny:Maxy),Xi_8(NLEN),Yi_8(NLEN), &
               s1_8(NLEN),s2_8(NLEN),s3_8(NLEN),s4_8(NLEN)
!
!author
!           Abdessamad Qaddouri - October 2009
!revision
! v4_60 - Qaddouri A.   - initial version
! v4_70 - Desgagne M.   - introduce mono_L
!
      integer :: k,i,Im,Jm
      real*8  :: W1,W2,W3,W4,X1,XX,X2,X3,X4
      real*8  :: WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4,Fi
      real*8  :: fx1,fx2,fx3,fx4
      real*8  :: fy1,fy2,fy3,fy4
      real*8  :: FFu_8, FFv_8
!
!----------------------------------------------------------------------
!
!$omp parallel private (k,i,im,jm, &
!$omp          W1,W2,W3,W4,X1,XX,X2,X3,X4, &
!$omp          WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4,Fi,&
!$omp          fx1,fx2,fx3,fx4,FFu_8,&
!$omp          fy1,fy2,fy3,fy4,FFv_8   )
!$omp do
      do i=1,NLEN
         Im  = iMx(i)
         Jm  = iMy(i)

         X1  = Geomgx_8(im)
         X2  = Geomgx_8(im+1) - X1
         X3  = Geomgx_8(im+2) - X1
         X4  = Geomgx_8(im+3) - X1
         XX  = Xi_8(i)   - X1

         W1  = (XX-X2)*(XX-X3)*(XX-X4)/(-X2*X3*X4)
         W2  =  XX    *(XX-X3)*(XX-X4)/(X2*(X2-X3)*(X2-X4))
         W3  =  XX    *(XX-X2)*(XX-X4)/(X3*(X3-X2)*(X3-X4))
         W4  =  XX    *(XX-X2)*(XX-X3)/(X4*(X4-X2)*(X4-X3))
         Y1  = Geomgy_8(Jm)
         Y2  = Geomgy_8(Jm+1) - Y1
         Y3  = Geomgy_8(Jm+2) - Y1
         Y4  = Geomgy_8(Jm+3) - Y1
         YY  = Yi_8(i)   - Y1

         WW1  = (yy-y2)*(yy-y3)*(yy-y4)/(-y2*y3*y4)
         WW2  =  yy    *(yy-y3)*(yy-y4)/(y2*(y2-y3)*(y2-y4))
         WW3  =  yy    *(yy-y2)*(yy-y4)/(y3*(y3-y2)*(y3-y4))
         WW4  =  yy    *(yy-y2)*(yy-y3)/(y4*(y4-y2)*(y4-y3))

         do k=1,nk
            !For U
            Fx1 = W1*Fu_8(Im,Jm,k) + W2*Fu_8(Im+1,Jm,k) + W3*Fu_8(Im+2,Jm,k)+ &
                  W4*Fu_8(Im+3,jm,k)

            Fx2 = W1*Fu_8(Im,Jm+1,k) + W2*Fu_8(Im+1,jm+1,k) + W3*Fu_8(Im+2,Jm+1,k)+ &
                  W4*Fu_8(Im+3,Jm+1,k)

            Fx3 = W1*Fu_8(Im,Jm+2,k) + W2*Fu_8(Im+1,jm+2,k) + W3*Fu_8(Im+2,Jm+2,k)+ &
                  W4*Fu_8(Im+3,Jm+2,k)

            Fx4 = W1*Fu_8(Im,jm+3,k) + W2*Fu_8(Im+1,jm+3,k) + W3*Fu_8(Im+2,Jm+3,k)+ &
                  W4*Fu_8(Im+3,Jm+3,k)

            FFu_8  = WW1*fx1+WW2*fx2+WW3*fx3+WW4*fx4

            !For V
            Fy1 = W1*Fv_8(Im,Jm,k) + W2*Fv_8(Im+1,Jm,k) + W3*Fv_8(Im+2,Jm,k)+ &
                  W4*Fv_8(Im+3,jm,k)

            Fy2 = W1*Fv_8(Im,Jm+1,k) + W2*Fv_8(Im+1,jm+1,k) + W3*Fv_8(Im+2,Jm+1,k)+ &
                  W4*Fv_8(Im+3,Jm+1,k)

            Fy3 = W1*Fv_8(Im,Jm+2,k) + W2*Fv_8(Im+1,jm+2,k) + W3*Fv_8(Im+2,Jm+2,k)+ &
                  W4*Fv_8(Im+3,Jm+2,k)

            Fy4 = W1*Fv_8(Im,jm+3,k) + W2*Fv_8(Im+1,jm+3,k) + W3*Fv_8(Im+2,Jm+3,k)+ &
                  W4*Fv_8(Im+3,Jm+3,k)


            FFv_8  = WW1*fy1+WW2*fy2+WW3*fy3+WW4*fy4

            FF((i-1)*Nk*2+ k*2-1)=real(s1_8(i)*FFu_8 + s2_8(i)*FFv_8)
            FF((i-1)*Nk*2+ k*2)=real(s3_8(i)*FFu_8 + s4_8(i)*FFv_8)
         end do
      end do
!$omp end do
!$omp end parallel
!
!----------------------------------------------------------------------
!
      return
      end
