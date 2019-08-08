!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------


subroutine verint_lnp (F_dest,F_dstlev,F_src,F_srclev,n,nks,nkd)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   !@object cubic vertical interpolation in ln(p)
   !@arguments
   ! F_dest        O           Interpolated field (destination)
   ! F_src         I           source       field (source)
   ! F_srclev      I           source levels
   ! F_dstlev      I           destination levels
   integer n, nkd, nks
   real F_dest(n,nkd), F_dstlev(n,nkd)
   real F_src (n,nks), F_srclev(n,nks)
   !@author - Methot/Laroche - April 97 - v1_01
   !@revision
   ! v2_10 - L. Corbeil          - rewrited for optimization,
   ! v2_10                         removed e_vcubique
   ! v2_30 - L. Corbeil          - renamed vte_ (no more called in gemntr)
   ! v2_31 - Lee V.              - F_srclev,F_dstlev is calculated outside

   include "rmnlib_basics.inc"

   integer i,k,iter,niter,lev,lev_lin
   integer top(n),bot(n),topcub(n),botcub(n),ref(n)
   real(REAL64) :: deltalev,prxd,prda,prdb,prsaf,prsbf,prsad,prsbd
   !     ---------------------------------------------------------------
   !***  General case, we'll care about VT later
   !
   !     First, we find the level we are by squeezing the destination
   !     between increasing bot() and decreasing top(). We need log_2_(nks)
   !     iteration to squeeze it completely (each integer between 1 and
   !     nks can be expressed as a sum of power of 2 such as sum c_i 2^i
   !     where c_i = 0 or 1 and i lower or equal than log_2_(nks)
   !
   !     WARNING
   !  niter calculation is ok for nks.lt. 2097152: should be ok for now...
   !  (Maybe the grid will be that precise in 2010!) (I don't even bother
   !  to add an if statement, for performance purpose...)
   if (real(int(log(real(nks))/log(2.0))) == &
        log(real(nks))/log(2.0)) then
      niter=int(log(real(nks))/log(2.0))
   else
      niter=int(log(real(nks))/log(2.0))+1
   endif

   !     squeeze...

   do k=1,nkd
      do i=1,n
         top(i)=nks
         bot(i)=1
      enddo
      do iter=1,niter
         do i=1,n
            !     divide by two (the old fashioned way...)
            ref(i)=ishft(top(i)+bot(i),-1)
            !     adjust top or bot
            if(F_dstlev(i,k) < F_srclev(i,ref(i))) then
               top(i)=ref(i)
            else
               bot(i)=ref(i)
            endif
         enddo
      enddo
      !     adjusting top and bot to ensure we can perform cubic interpolation
      do i=1,n
         botcub(i)=max(2,bot(i))
         topcub(i)=min(nks-1,top(i))
      enddo
      !     cubic or linear interpolation
      do i=1,n
         lev=botcub(i)
         lev_lin=bot(i)
         deltalev=(F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin))

         !     Interpolation: if not enough points to perform cubic interpolation
         !                    we use linear interpolation or persistency

         if (lev.ne.lev_lin .or. topcub(i).ne.top(i)) then

            !     persistancy of this interval

            if (F_dstlev(i,k) <= F_srclev(i,1)) then
               F_dest(i,k) = F_src(i,1)
            else if (F_dstlev(i,k) >= F_srclev(i,nks)) then
               F_dest(i,k) = F_src(i,nks)
            else
               !     linear interpolation
               prxd = (F_dstlev(i,k)-F_srclev(i,lev_lin))/deltalev
               F_dest(i,k) = (1.0-prxd)*F_src(i,lev_lin) &
                    +prxd*F_src(i,lev_lin+1)
            endif
         else
            !     cubic interpolation
            prxd = (F_dstlev(i,k)-F_srclev(i,lev_lin))/ &
                 deltalev
            prda = ((F_src(i,lev_lin+1)-F_src(i,lev_lin-1))/ &
                 (F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin-1))* &
                 deltalev)
            prdb = ((F_src(i,lev_lin+2)-F_src(i,lev_lin))/ &
                 (F_srclev(i,lev_lin+2)-F_srclev(i,lev_lin))* &
                 deltalev)
            prsaf= (1.0+2.0*prxd)*(1.0-prxd)*(1.0-prxd)
            prsbf= (3.0-2.0*prxd)*prxd*prxd
            prsad= prxd*(1.0-prxd)*(1.0-prxd)
            prsbd= (1.0-prxd)*prxd*prxd
            F_dest(i,k) = F_src(i,lev_lin  )*prsaf &
                 +F_src(i,lev_lin+1)*prsbf+prda*prsad &
                 -prdb*prsbd
         endif
      enddo
   enddo
   !     ---------------------------------------------------------------
   return
end subroutine verint_lnp

