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
!**s/r int_near_lag - to do YY nearest interpolation for scalar
!
       subroutine int_near_lag3(FF,F,Imx,Imy,Geomgx,Geomgy,Minx,Maxx,Miny,Maxy,Nk,Xi,Yi,NLEN)

       use geomh
       implicit none
#include <arch_specific.hf>

       integer Minx,Maxx,Miny,Maxy,Nk,NLEN
       integer Imx(NLEN),Imy(NLEN)
       real FF(NLEN*Nk)
       real*8 F(Minx:Maxx,Miny:Maxy,Nk),geomgx(Minx:Maxx),geomgy(Miny:Maxy)
       real*8 Xi(NLEN),Yi(NLEN)

!
!author
!           Abdessamad Qaddouri - October 2009
!  PLEASE consult Abdessamad or Vivian before modifying this routine.
!
!
!revision
!  v4_8    V.Lee correction in interpolation (MPI precision sensitive)
!
!
       integer i,j,k,Im,Jm
       real*8  FF_8,betax,betax1,betay,betay1
!

!$omp parallel private (i,j,k,im,jm,       &
!$omp      FF_8,betax,betax1,betay,betay1)
!$omp do
   Do 200 i=1,NLEN
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
   Do 300 k=1,Nk
       FF_8= betay1*(betax1*F(Im,Jm,k)+betax*F(Im+1,Jm,k))+ &
           betay*(betax1*F(Im,Jm+1,k)+betax*F(Im+1,Jm+1,k))
       FF((i-1)*Nk+k)=real(FF_8)
   300 continue
   200 continue
!$omp enddo
!$omp end parallel

      return
      end

      subroutine int_near_lag2(FF,F,Imx,Imy,Geomgx,Geomgy,Minx,Maxx,Miny,Maxy,Xi,Yi)

       use geomh
       implicit none
#include <arch_specific.hf>
!
!author
!           Abdessamad Qaddouri - October 2009
!
       integer Imx,Imy,Minx,Maxx,Miny,Maxy
       integer Im, Jm
       real*8 F(Minx:Maxx,Miny:Maxy),geomgx(Minx:Maxx),geomgy(Miny:Maxy)
       real*8 FF,Xi,Yi
       real*8 betax,betax1,betay,betay1
!
       Im = imx
       Jm = imy
           betax= ( Xi-Geomgx(Im))*geomh_inv_hx_8
           betax1= (1.0-betax)
           if (betax <= betax1) then
             betax=0.0d0
             betax1=1.0d0
           else
             betax=1.0d0
             betax1=0.0d0
           end if
           betay=(Yi-Geomgy(Jm))*geomh_inv_hy_8
           betay1=1.0-betay
           if (betay <= betay1) then
             betay=0.0d0
             betay1=1.0d0
           else
             betay=1.0d0
             betay1=0.0d0
           end if
           FF= betay1*(betax1*F(Im,Jm)+betax*F(Im+1,Jm))+ &
                    betay*(betax1*F(Im,Jm+1)+betax*F(Im+1,Jm+1))
      return
      end
