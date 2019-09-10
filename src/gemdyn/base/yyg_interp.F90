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

!**s/r yyg_int_cub88 - YY cubic interpolation for scalar

      subroutine yyg_int_cub88 ( F_dest, F_src, Imx,Imy, geomgx,geomgy,&
                                  Minx,Maxx,Miny,Maxy,Nk,Xi,Yi,NLEN )
      use glb_ld
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Nk,NLEN,Minx,Maxx,Miny,Maxy
      integer, dimension(NLEN), intent(in) :: Imx, Imy
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_src
      real(kind=REAL64), dimension(NLEN*Nk)    , intent(out):: F_dest
      real(kind=REAL64), dimension(Minx:Maxx), intent(in) :: geomgx
      real(kind=REAL64), dimension(Miny:Maxy), intent(in) :: geomgy
      real(kind=REAL64), dimension(NLEN)     , intent(in) :: Xi, Yi
!
!author
!           Abdessamad Qaddouri - October 2009
!revision
! v4_60 - Qaddouri A.   - initial version
!
       integer k,i,j,Im, Jm
       real(kind=REAL64)  W1,W2,W3,W4,X1,XX,X2,X3,X4
       real(kind=REAL64)  WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4
       real(kind=REAL64)  fx1,fx2,fx3,fx4
!
!----------------------------------------------------------------------
!
!$omp parallel private (k,i,j,im,jm, &
!$omp          W1,W2,W3,W4,X1,XX,X2,X3,X4, &
!$omp          WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4,&
!$omp          fx1,fx2,fx3,fx4)
!$omp do
   Do i=1,NLEN
       Im  = iMx(i)
       Jm  = iMy(i)

       X1  = Geomgx(im)
       X2  = Geomgx(im+1) - X1
       X3  = Geomgx(im+2) - X1
       X4  = Geomgx(im+3) - X1
       XX  = Xi(i)   - X1

       W1  = (XX-X2)*(XX-X3)*(XX-X4)/(-X2*X3*X4)
       W2  =  XX    *(XX-X3)*(XX-X4)/(X2*(X2-X3)*(X2-X4))
       W3  =  XX    *(XX-X2)*(XX-X4)/(X3*(X3-X2)*(X3-X4))
       W4  =  XX    *(XX-X2)*(XX-X3)/(X4*(X4-X2)*(X4-X3))

       Y1  = Geomgy(Jm)
       Y2  = Geomgy(Jm+1) - Y1
       Y3  = Geomgy(Jm+2) - Y1
       Y4  = Geomgy(Jm+3) - Y1
       YY  = Yi(i)   - Y1

       WW1  = (yy-y2)*(yy-y3)*(yy-y4)/(-y2*y3*y4)
       WW2  =  yy    *(yy-y3)*(yy-y4)/(y2*(y2-y3)*(y2-y4))
       WW3  =  yy    *(yy-y2)*(yy-y4)/(y3*(y3-y2)*(y3-y4))
       WW4  =  yy    *(yy-y2)*(yy-y3)/(y4*(y4-y2)*(y4-y3))

       Do k=1,Nk

          Fx1 = W1*dble(F_src(Im  ,Jm  ,k)) + W2*dble(F_src(Im+1,Jm  ,k)) &
               +W3*dble(F_src(Im+2,Jm  ,k)) + W4*dble(F_src(Im+3,jm  ,k))

          Fx2 = W1*dble(F_src(Im  ,Jm+1,k)) + W2*dble(F_src(Im+1,Jm+1,k)) &
               +W3*dble(F_src(Im+2,Jm+1,k)) + W4*dble(F_src(Im+3,Jm+1,k))

          Fx3 = W1*dble(F_src(Im  ,Jm+2,k)) + W2*dble(F_src(Im+1,Jm+2,k)) &
               +W3*dble(F_src(Im+2,Jm+2,k)) + W4*dble(F_src(Im+3,Jm+2,k))

          Fx4 = W1*dble(F_src(Im  ,Jm+3,k)) + W2*dble(F_src(Im+1,Jm+3,k)) &
               +W3*dble(F_src(Im+2,Jm+3,k)) + W4*dble(F_src(Im+3,Jm+3,k))

          F_dest((i-1)*Nk+k)  = WW1*fx1+WW2*fx2+WW3*fx3+WW4*fx4
       end do
      end do
!$omp enddo
!$omp end parallel
!
!----------------------------------------------------------------------
!
       return
       end subroutine yyg_int_cub88

!**s/r yyg_int_cub48 - YY cubic interpolation for scalar

      subroutine yyg_int_cub48 ( F_dest, F_src, Imx,Imy, geomgx,geomgy,&
                                  Minx,Maxx,Miny,Maxy,Nk,Xi,Yi,NLEN )
      use glb_ld
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Nk,NLEN,Minx,Maxx,Miny,Maxy
      integer, dimension(NLEN), intent(in) :: Imx, Imy
      real  , dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_src
      real(kind=REAL64), dimension(NLEN*Nk)    , intent(out):: F_dest
      real(kind=REAL64), dimension(Minx:Maxx), intent(in) :: geomgx
      real(kind=REAL64), dimension(Miny:Maxy), intent(in) :: geomgy
      real(kind=REAL64), dimension(NLEN)     , intent(in) :: Xi, Yi
!
!author
!           Abdessamad Qaddouri - October 2009
!revision
! v4_60 - Qaddouri A.   - initial version
!
       integer k,i,j,Im, Jm
       real(kind=REAL64)  W1,W2,W3,W4,X1,XX,X2,X3,X4
       real(kind=REAL64)  WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4
       real(kind=REAL64)  fx1,fx2,fx3,fx4
!
!----------------------------------------------------------------------
!
!$omp parallel private (k,i,j,im,jm, &
!$omp          W1,W2,W3,W4,X1,XX,X2,X3,X4, &
!$omp          WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4,&
!$omp          fx1,fx2,fx3,fx4)
!$omp do
   Do i=1,NLEN
       Im  = iMx(i)
       Jm  = iMy(i)

       X1  = Geomgx(im)
       X2  = Geomgx(im+1) - X1
       X3  = Geomgx(im+2) - X1
       X4  = Geomgx(im+3) - X1
       XX  = Xi(i)   - X1

       W1  = (XX-X2)*(XX-X3)*(XX-X4)/(-X2*X3*X4)
       W2  =  XX    *(XX-X3)*(XX-X4)/(X2*(X2-X3)*(X2-X4))
       W3  =  XX    *(XX-X2)*(XX-X4)/(X3*(X3-X2)*(X3-X4))
       W4  =  XX    *(XX-X2)*(XX-X3)/(X4*(X4-X2)*(X4-X3))

       Y1  = Geomgy(Jm)
       Y2  = Geomgy(Jm+1) - Y1
       Y3  = Geomgy(Jm+2) - Y1
       Y4  = Geomgy(Jm+3) - Y1
       YY  = Yi(i)   - Y1

       WW1  = (yy-y2)*(yy-y3)*(yy-y4)/(-y2*y3*y4)
       WW2  =  yy    *(yy-y3)*(yy-y4)/(y2*(y2-y3)*(y2-y4))
       WW3  =  yy    *(yy-y2)*(yy-y4)/(y3*(y3-y2)*(y3-y4))
       WW4  =  yy    *(yy-y2)*(yy-y3)/(y4*(y4-y2)*(y4-y3))

       Do k=1,Nk

          Fx1 = W1*dble(F_src(Im  ,Jm  ,k)) + W2*dble(F_src(Im+1,Jm  ,k)) &
               +W3*dble(F_src(Im+2,Jm  ,k)) + W4*dble(F_src(Im+3,jm  ,k))

          Fx2 = W1*dble(F_src(Im  ,Jm+1,k)) + W2*dble(F_src(Im+1,Jm+1,k)) &
               +W3*dble(F_src(Im+2,Jm+1,k)) + W4*dble(F_src(Im+3,Jm+1,k))

          Fx3 = W1*dble(F_src(Im  ,Jm+2,k)) + W2*dble(F_src(Im+1,Jm+2,k)) &
               +W3*dble(F_src(Im+2,Jm+2,k)) + W4*dble(F_src(Im+3,Jm+2,k))

          Fx4 = W1*dble(F_src(Im  ,Jm+3,k)) + W2*dble(F_src(Im+1,Jm+3,k)) &
               +W3*dble(F_src(Im+2,Jm+3,k)) + W4*dble(F_src(Im+3,Jm+3,k))

          F_dest((i-1)*Nk+k)  = WW1*fx1+WW2*fx2+WW3*fx3+WW4*fx4
       end do
      end do
!$omp enddo
!$omp end parallel
!
!----------------------------------------------------------------------
!
       return
       end subroutine yyg_int_cub48

!**s/r yyg_int_cub - YY cubic interpolation for scalar

      subroutine yyg_int_cub ( F_dest, F_src, Imx,Imy, geomgx,geomgy,   &
                               Minx,Maxx,Miny,Maxy,Nk,Xi,Yi,NLEN,mono_l )
      use glb_ld
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: mono_L
      integer, intent(in) :: Nk,NLEN,Minx,Maxx,Miny,Maxy
      integer, dimension(NLEN), intent(in) :: Imx, Imy
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in) :: F_src
      real, dimension(NLEN*Nk)    , intent(out):: F_dest
      real(kind=REAL64), dimension(Minx:Maxx), intent(in) :: geomgx
      real(kind=REAL64), dimension(Miny:Maxy), intent(in) :: geomgy
      real(kind=REAL64), dimension(NLEN)     , intent(in) :: Xi, Yi

!author
!           Abdessamad Qaddouri - October 2009
!revision
! v4_60 - Qaddouri A.   - initial version

       integer k,i,j,Im, Jm
       real(kind=REAL64)  W1,W2,W3,W4,X1,XX,X2,X3,X4
       real(kind=REAL64)  WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4
       real(kind=REAL64)  fx1,fx2,fx3,fx4,prmax,prmin
!
!----------------------------------------------------------------------
!
!$omp parallel private (k,i,j,im,jm, &
!$omp     W1,W2,W3,W4,X1,XX,X2,X3,X4, &
!$omp     WW1,WW2,WW3,WW4,YY,y1,y2,y3,y4,&
!$omp     fx1,fx2,fx3,fx4,prmax,prmin)
!$omp do
       Do i=1,NLEN
          Im  = iMx(i)
          Jm  = iMy(i)

          X1  = Geomgx(im)
          X2  = Geomgx(im+1) - X1
          X3  = Geomgx(im+2) - X1
          X4  = Geomgx(im+3) - X1
          XX  = Xi(i)   - X1

          W1  = (XX-X2)*(XX-X3)*(XX-X4)/(-X2*X3*X4)
          W2  =  XX    *(XX-X3)*(XX-X4)/(X2*(X2-X3)*(X2-X4))
          W3  =  XX    *(XX-X2)*(XX-X4)/(X3*(X3-X2)*(X3-X4))
          W4  =  XX    *(XX-X2)*(XX-X3)/(X4*(X4-X2)*(X4-X3))

          Y1  = Geomgy(Jm)
          Y2  = Geomgy(Jm+1) - Y1
          Y3  = Geomgy(Jm+2) - Y1
          Y4  = Geomgy(Jm+3) - Y1
          YY  = Yi(i)   - Y1

          WW1  = (yy-y2)*(yy-y3)*(yy-y4)/(-y2*y3*y4)
          WW2  =  yy    *(yy-y3)*(yy-y4)/(y2*(y2-y3)*(y2-y4))
          WW3  =  yy    *(yy-y2)*(yy-y4)/(y3*(y3-y2)*(y3-y4))
          WW4  =  yy    *(yy-y2)*(yy-y3)/(y4*(y4-y2)*(y4-y3))

          Do k=1,Nk

             if (mono_l) then
                prmax=max(F_src(Im,Jm,k),F_src(Im+1,Jm,k),F_src(Im+2,Jm,k),F_src(Im+3,jm,k)   ,&
                F_src(Im,Jm+1,k),F_src(Im+1,jm+1,k),F_src(Im+2,Jm+1,k),F_src(Im+3,Jm+1,k),&
                F_src(Im,Jm+2,k),F_src(Im+1,jm+2,k),F_src(Im+2,Jm+2,k),F_src(Im+3,Jm+2,k),&
                F_src(Im,jm+3,k),F_src(Im+1,jm+3,k),F_src(Im+2,Jm+3,k),F_src(Im+3,Jm+3,k))
                prmin=min(F_src(Im,Jm,k),F_src(Im+1,Jm,k),F_src(Im+2,Jm,k),F_src(Im+3,jm,k)   ,&
                F_src(Im,Jm+1,k),F_src(Im+1,jm+1,k),F_src(Im+2,Jm+1,k),F_src(Im+3,Jm+1,k),&
                F_src(Im,Jm+2,k),F_src(Im+1,jm+2,k),F_src(Im+2,Jm+2,k),F_src(Im+3,Jm+2,k),&
                F_src(Im,jm+3,k),F_src(Im+1,jm+3,k),F_src(Im+2,Jm+3,k),F_src(Im+3,Jm+3,k))
             end if

             Fx1 = W1*dble(F_src(Im  ,Jm  ,k)) + W2*dble(F_src(Im+1,Jm  ,k)) &
                  +W3*dble(F_src(Im+2,Jm  ,k)) + W4*dble(F_src(Im+3,jm  ,k))

             Fx2 = W1*dble(F_src(Im  ,Jm+1,k)) + W2*dble(F_src(Im+1,Jm+1,k)) &
                  +W3*dble(F_src(Im+2,Jm+1,k)) + W4*dble(F_src(Im+3,Jm+1,k))

             Fx3 = W1*dble(F_src(Im  ,Jm+2,k)) + W2*dble(F_src(Im+1,Jm+2,k)) &
                  +W3*dble(F_src(Im+2,Jm+2,k)) + W4*dble(F_src(Im+3,Jm+2,k))

             Fx4 = W1*dble(F_src(Im  ,Jm+3,k)) + W2*dble(F_src(Im+1,Jm+3,k)) &
                  +W3*dble(F_src(Im+2,Jm+3,k)) + W4*dble(F_src(Im+3,Jm+3,k))

             F_dest((i-1)*Nk+k)= WW1*fx1 + WW2*fx2 + WW3*fx3 + WW4*fx4
             if (mono_L) F_dest((i-1)*Nk+k)= max(real(prmin), min(real(prmax),F_dest((i-1)*Nk+k)))

          end do
       end do
!$omp enddo
!$omp end parallel
!
!----------------------------------------------------------------------
!
       return
       end subroutine yyg_int_cub

!**s/r yyg_int_lin - YY linear interpolation for scalar

       subroutine yyg_int_lin ( FF,F, Imx,Imy, Geomgx,Geomgy,&
                                minx,maxx,miny,maxy,Nk,Xi,Yi,NLEN )
       use geomh
      use, intrinsic :: iso_fortran_env
       implicit none
#include <arch_specific.hf>

       integer minx,miny,maxx,maxy,Nk,NLEN
       integer Imx(NLEN),Imy(NLEN)
       real FF(NLEN*Nk),F(minx:maxx,miny:maxy,Nk)
       real(kind=REAL64) geomgx(Minx:Maxx),geomgy(Miny:Maxy)
       real(kind=REAL64) Xi(NLEN),Yi(NLEN)

!author
!           Abdessamad Qaddouri - October 2009

       integer i,j,k,Im, Jm
       real(kind=REAL64) betax,betax1,betay,betay1
!
!     ---------------------------------------------------------------
!
!$omp parallel private (i,j,k,im,jm, &
!$omp     betax,betax1,betay,betay1)
!$omp do

       Do i=1,NLEN

          Im = Imx(i)
          Jm = Imy(i)
          betax= (Xi(i)-Geomgx(Im))*geomh_inv_hx_8
          betax1= (1.0d0-betax)
          betay=(Yi(i)-Geomgy(Jm))*geomh_inv_hy_8
          betay1=1.0d0-betay
          Do k=1,Nk
             FF((i-1)*Nk+k) = betay1*(betax1*F(Im,Jm  ,k)+betax*F(Im+1,Jm  ,k))&
                             +betay *(betax1*F(Im,Jm+1,k)+betax*F(Im+1,Jm+1,k))
          end do
       end do

!$omp enddo
!$omp end parallel
!
!     ---------------------------------------------------------------
!
      return
      end subroutine yyg_int_lin

!**s/r yyg_int_near - YY nearest interpolation for scalar

       subroutine yyg_int_near ( FF,F, Imx,Imy, Geomgx,Geomgy,&
                                minx,maxx,miny,maxy,Nk,Xi,Yi,NLEN )
       use geomh
      use, intrinsic :: iso_fortran_env
       implicit none
#include <arch_specific.hf>

       integer Minx,Maxx,Miny,Maxy,Nk,NLEN
       integer Imx(NLEN),Imy(NLEN)
       real FF(NLEN*Nk),F(Minx:Maxx,Miny:Maxy,Nk)
       real(kind=REAL64) geomgx(Minx:Maxx),geomgy(Miny:Maxy)
       real(kind=REAL64) Xi(NLEN),Yi(NLEN)

!
!author
!           Abdessamad Qaddouri - October 2009

       integer i,j,k,Im,Jm
       real(kind=REAL64)  betax,betax1,betay,betay1
!
!     ---------------------------------------------------------------
!
!$omp parallel private (i,j,k,im,jm,       &
!$omp      betax,betax1,betay,betay1)
!$omp do
       Do i=1,NLEN
          Im = imx(i)
          Jm = imy(i)

          betax= ( Xi(i) - Geomgx(Im))*geomh_inv_hx_8
          betax1= (1.0-betax)
          if (betax <= betax1) then
             betax=0.0d0
             betax1=1.0d0
          else
             betax=1.0d0
             betax1=0.0d0
          end if

          betay=(Yi(i) - Geomgy(Jm))*geomh_inv_hy_8
          betay1=1.0-betay
          if (betay <= betay1) then
             betay=0.0d0
             betay1=1.0d0
          else
             betay=1.0d0
             betay1=0.0d0
          end if
          Do k=1,Nk
             FF((i-1)*Nk+k)= betay1*(betax1*F(Im,Jm,k)+betax*F(Im+1,Jm,k))+ &
             betay*(betax1*F(Im,Jm+1,k)+betax*F(Im+1,Jm+1,k))
          end do
       end do

!$omp enddo
!$omp end parallel
!
!     ---------------------------------------------------------------
!
      return
      end subroutine yyg_int_near
