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
!**s/r int_lin_lag - to do YY linear interpolation for scalar
       subroutine int_lin_lag3(FF,F,Imx,Imy,Geomgx,Geomgy,minx,maxx,miny,maxy,Nk,Xi,Yi,NLEN)


       use geomh
       implicit none
#include <arch_specific.hf>

       integer minx,miny,maxx,maxy,Nk,NLEN
       integer Imx(NLEN),Imy(NLEN)
       real FF(NLEN*Nk)
       real*8 F(minx:maxx,miny:maxy,Nk),geomgx(Minx:Maxx),geomgy(Miny:Maxy)
       real*8 Xi(NLEN),Yi(NLEN)
!
!author
!           Abdessamad Qaddouri - October 2009
!  PLEASE consult Abdessamad or Vivian before modifying this routine.
!revision
!  v4_8    V.Lee correction in interpolation (MPI precision sensitive)
!
       integer i,j,k,Im, Jm
       real*8 FF_8,betax,betax1,betay,betay1


!$omp parallel private (i,j,k,im,jm, &
!$omp     FF_8,betax,betax1,betay,betay1)
!$omp do

   Do 200 i=1,NLEN

       Im = Imx(i)
       Jm = Imy(i)
       betax= (Xi(i)-Geomgx(Im))*geomh_inv_hx_8
       betax1= (1.0d0-betax)
       betay=(Yi(i)-Geomgy(Jm))*geomh_inv_hy_8
       betay1=1.0d0-betay
   Do 300 k=1,Nk
       FF_8 = betay1*(betax1*F(Im,Jm,k)+betax*F(Im+1,Jm,k))+ &
                    betay*(betax1*F(Im,Jm+1,k)+betax*F(Im+1,Jm+1,k))
       FF((i-1)*Nk+k)=real(FF_8)
   300 continue
   200 continue
!$omp enddo
!$omp end parallel
      return
      end

      subroutine int_lin_lag2(FF,F,Imx,Imy,Geomgx,Geomgy,minx,maxx,miny,maxy,Xi,Yi)

       use geomh
       implicit none
#include <arch_specific.hf>
!
!author
!           Abdessamad Qaddouri - October 2009
!
!  v4_8    V.Lee correction in interpolation (MPI precision sensitive)
!

       integer Imx,Imy,minx,miny,maxx,maxy
       integer Im, Jm
       real*8 F(minx:maxx,miny:maxy),geomgx(Minx:Maxx),geomgy(Miny:Maxy)
       real*8 FF,Xi,Yi
       real*8 betax,betax1,betay,betay1

       Im = Imx
       Jm = Imy
           betax= (Xi-Geomgx(Im))*geomh_inv_hx_8
           betax1= (1.0-betax)
           betay=(Yi-Geomgy(Jm))*geomh_inv_hy_8
           betay1=1.0-betay
           FF= betay1*(betax1*F(Im,Jm)+betax*F(Im+1,Jm))+ &
                    betay*(betax1*F(Im,Jm+1)+betax*F(Im+1,Jm+1))
      return
      end

      subroutine int_lin_lag(FF,F,Imx,Imy,Ni,Nj,Nx,Ny,Xi, &
                                            Yi,x,y)

       use geomh
       implicit none
!
!author
!           Abdessamad Qaddouri - October 2009
!

       integer Ni,Nj,Imx(Ni,Nj),Imy(Ni,Nj),Nx,Ny
       integer i,j
       integer Im, Jm
       real*8 FF(Ni,Nj)
       real*8 F(Nx,Ny),x(*),y(*)
       real*8 Xi(Ni,Nj),Yi(Ni,Nj)
       real*8 betax,betax1,betay,betay1

       do j=1, Nj
        do i = 1, Ni
           Im = imx(i,j)
           Jm = imy(i,j)
           betax= ( Xi(i,j)-x(Im))*geomh_inv_hx_8
           betax1= (1.0-betax)
           betay=(Yi(i,j)-y(Jm))*geomh_inv_hy_8
           betay1=1.0-betay
           FF(i,j)= betay1*(betax1*F(Im,Jm)+betax*F(Im+1,Jm))+ &
                    betay*(betax1*F(Im,Jm+1)+betax*F(Im+1,Jm+1))
        end do
      end do
      return
      end
