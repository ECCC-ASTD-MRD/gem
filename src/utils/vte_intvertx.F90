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

!/@*
subroutine vte_intvertx3(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,F_type)
   implicit none
!!!#include <arch_specific.hf>
   !@ojective cubic vertical interpolation from eta/sigma to hybrid
   !@arguments
   ! F_dch                   Interpolated field
   ! F_sch                   source       field
   ! F_srclev                source levels
   ! F_dstlev                destination levels
   ! F_var_S                 name of the variable
   ! F_type                  order of interpolation ('linear' or 'cubic')
   integer,intent(in) :: n, nkd, nks
   real :: F_dch(n,nkd), F_sch(n,nks), F_srclev(n,nks),F_dstlev(n,nkd)
   character(len=*) :: F_var_S,F_type
   !---------------------------------------------------------------
   if (any(F_type(1:1) == (/'l','L'/))) then
      call vte_intvertx4(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,nkd)
   else
      call vte_intvertx4(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,0)
   endif
   !---------------------------------------------------------------
   return
end subroutine vte_intvertx3


!/@*
subroutine vte_intvertx4(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,F_nlinbot)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use tdpack
   implicit none
!!!#include <arch_specific.hf>
   !@ojective cubic vertical interpolation from eta/sigma to hybrid
   !@arguments
   ! F_dch                   Interpolated field
   ! F_sch                   source       field
   ! F_srclev                source levels
   ! F_dstlev                destination levels
   ! F_var_S                 name of the variable
   ! F_nlinbot               number of levels at bottom linearly interpolated
   integer,intent(in) :: n, nkd, nks,F_nlinbot
   real :: F_dch(n,nkd), F_sch(n,nks), F_srclev(n,nks),F_dstlev(n,nkd)
   character(len=*) :: F_var_S
   !@author - Methot/Laroche - April 97 - v1_01
   !@revision
   ! v2_10 - L. Corbeil          - rewrited for optimization,
   ! v2_10                         removed e_vcubique
   ! v2_30 - L. Corbeil          - renamed vte_ (no more called in gemntr)
   ! v2_31 - Lee V.              - F_srclev,F_dstlev is calculated outside
   ! v3_32 - Plante A.           - correction of F_dstlev(n,nks) to F_dstlev(n,nkd)
   !@description
   !     General case, we'll care about VT later
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
   !*@/
   include "rmnlib_basics.inc"

   integer :: i,k,iter,niter,lev,lev_lin,k_surf,k_ciel
   integer :: top(n),bot(n),topcub(n),botcub(n),ref(n)
   real(REAL64) :: deltalev,prxd,prda,prdb,prsaf,prsbf,prsad,prsbd
   logical :: ascending_L
   !---------------------------------------------------------------
   ascending_L = (F_srclev(1,1) < F_srclev(1,nks))
   k_surf = 1
   k_ciel = nks
   if (.not.ascending_L) then
      k_surf = nks
      k_ciel = 1
   endif
   if (real(int(log(real(nks))/log(2.0))).eq. &
        log(real(nks))/log(2.0)) then
      niter=int(log(real(nks))/log(2.0))
   else
      niter=int(log(real(nks))/log(2.0))+1
   endif

   !- squeeze...

   do k=1,nkd
      top(1:n)=nks
      bot(1:n)=1
      if (ascending_L) then
         do iter=1,niter
            do i=1,n
               !     divide by two (the old fashioned way...)
               ref(i)=ishft(top(i)+bot(i),-1)
               !     adjust top or bot
               if(F_dstlev(i,k).lt.F_srclev(i,ref(i))) then
                  top(i)=ref(i)
               else
                  bot(i)=ref(i)
               endif
            enddo
         enddo
      else
         do iter=1,niter
            do i=1,n
               !     divide by two (the old fashioned way...)
               ref(i)=ishft(top(i)+bot(i),-1)
               !     adjust top or bot
               if(F_dstlev(i,k).gt.F_srclev(i,ref(i))) then
                  top(i)=ref(i)
               else
                  bot(i)=ref(i)
               endif
            enddo
         enddo
      endif
      !- adjusting top and bot to ensure we can perform cubic interpolation
      do i=1,n
         botcub(i)=max(2,bot(i))
         topcub(i)=min(nks-1,top(i))
      enddo
      !- cubic or linear interpolation
      do i=1,n
         lev=botcub(i)
         lev_lin=bot(i)
         deltalev=(F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin))

         !- Interpolation: either if not enough points to perform cubic interp
         !                 or if linear interpolation is requested use linear
         !                 interpolation and constant extrapolation
         if((lev.ne.lev_lin).or.(topcub(i).ne.top(i)).or.k<=F_nlinbot) then
            !- persistancy of this interval
            if(F_dstlev(i,k).le.F_srclev(i,k_surf)) then
               F_dch(i,k) = F_sch(i,k_surf)
            else if(F_dstlev(i,k).ge.F_srclev(i,k_ciel)) then
               F_dch(i,k) = F_sch(i,k_ciel)
            else
               !- linear interpolation
               prxd=(F_dstlev(i,k)-F_srclev(i,lev_lin))/deltalev
               F_dch(i,k) = (1.0-prxd)*F_sch(i,lev_lin) &
                    +prxd*F_sch(i,lev_lin+1)
            endif
         else
            !- cubic interpolation
            prxd = (F_dstlev(i,k)-F_srclev(i,lev_lin))/ &
                 deltalev
            prda = ((F_sch(i,lev_lin+1)-F_sch(i,lev_lin-1))/ &
                 (F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin-1))* &
                 deltalev)
            prdb = ((F_sch(i,lev_lin+2)-F_sch(i,lev_lin))/ &
                 (F_srclev(i,lev_lin+2)-F_srclev(i,lev_lin))* &
                 deltalev)
            prsaf= (1.0+2.0*prxd)*(1.0-prxd)*(1.0-prxd)
            prsbf= (3.0-2.0*prxd)*prxd*prxd
            prsad= prxd*(1.0-prxd)*(1.0-prxd)
            prsbd= (1.0-prxd)*prxd*prxd
            F_dch(i,k) = F_sch(i,lev_lin  )*prsaf &
                 +F_sch(i,lev_lin+1)*prsbf+prda*prsad &
                 -prdb*prsbd
         endif
      enddo
   enddo

   !- special case for VT
   if(F_var_S == 'VT') then
      do k=1,nkd
         do i=1,n
            if(F_srclev(i,k_ciel).lt.F_dstlev(i,k)) then
               F_dch(i,k) = F_sch(i,k_ciel) * exp (  &
                    rgasd_8*stlo_8*(F_dstlev(i,k)-F_srclev(i,k_ciel)) )
            endif
         enddo
      enddo
   endif
   !     ---------------------------------------------------------------
   return
end subroutine vte_intvertx4


!/@*
subroutine vte_intvertx_isodst(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,F_nlinbot)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use tdpack
   implicit none
!!!#include <arch_specific.hf>
   !@ojective cubic vertical interpolation from eta/sigma to hybrid
   !@arguments
   ! F_dch                   Interpolated field
   ! F_sch                   source       field
   ! F_srclev                source levels
   ! F_dstlev                destination levels
   ! F_var_S                 name of the variable
   ! F_nlinbot               number of levels at bottom linearly interpolated
   integer,intent(in) :: n, nkd, nks,F_nlinbot
   real :: F_dch(n,nkd), F_sch(n,nks), F_srclev(n,nks),F_dstlev(nkd)
   character(len=*) :: F_var_S
   !@author - Methot/Laroche - April 97 - v1_01
   !@revision
   ! v2_10 - L. Corbeil          - rewrited for optimization,
   ! v2_10                         removed e_vcubique
   ! v2_30 - L. Corbeil          - renamed vte_ (no more called in gemntr)
   ! v2_31 - Lee V.              - F_srclev,F_dstlev is calculated outside
   ! v3_32 - Plante A.           - correction of F_dstlev(n,nks) to F_dstlev(n,nkd)
   !@description
   !     General case, we'll care about VT later
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
   !*@/
   include "rmnlib_basics.inc"

   integer :: i,k,iter,niter,lev,lev_lin,k_surf,k_ciel
   integer :: top(n),bot(n),topcub(n),botcub(n),ref(n)
   real(REAL64) :: deltalev,prxd,prda,prdb,prsaf,prsbf,prsad,prsbd
   logical :: ascending_L
   !---------------------------------------------------------------
   ascending_L = (F_srclev(1,1) < F_srclev(1,nks))
   k_surf = 1
   k_ciel = nks
   if (.not.ascending_L) then
      k_surf = nks
      k_ciel = 1
   endif
   if (real(int(log(real(nks))/log(2.0))).eq. &
        log(real(nks))/log(2.0)) then
      niter=int(log(real(nks))/log(2.0))
   else
      niter=int(log(real(nks))/log(2.0))+1
   endif

   !- squeeze...

   do k=1,nkd
      top(1:n)=nks
      bot(1:n)=1
      if (ascending_L) then
         do iter=1,niter
            do i=1,n
               !     divide by two (the old fashioned way...)
               ref(i)=ishft(top(i)+bot(i),-1)
               !     adjust top or bot
               if(F_dstlev(k).lt.F_srclev(i,ref(i))) then
                  top(i)=ref(i)
               else
                  bot(i)=ref(i)
               endif
            enddo
         enddo
      else
         do iter=1,niter
            do i=1,n
               !     divide by two (the old fashioned way...)
               ref(i)=ishft(top(i)+bot(i),-1)
               !     adjust top or bot
               if(F_dstlev(k).gt.F_srclev(i,ref(i))) then
                  top(i)=ref(i)
               else
                  bot(i)=ref(i)
               endif
            enddo
         enddo
      endif
      !- adjusting top and bot to ensure we can perform cubic interpolation
      do i=1,n
         botcub(i)=max(2,bot(i))
         topcub(i)=min(nks-1,top(i))
      enddo
      !- cubic or linear interpolation
      do i=1,n
         lev=botcub(i)
         lev_lin=bot(i)
         deltalev=(F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin))

         !- Interpolation: either if not enough points to perform cubic interp
         !                 or if linear interpolation is requested use linear
         !                 interpolation and constant extrapolation
         if((lev.ne.lev_lin).or.(topcub(i).ne.top(i)).or.k<=F_nlinbot) then
            !- persistancy of this interval
            if(F_dstlev(k).le.F_srclev(i,k_surf)) then
               F_dch(i,k) = F_sch(i,k_surf)
            else if(F_dstlev(k).ge.F_srclev(i,k_ciel)) then
               F_dch(i,k) = F_sch(i,k_ciel)
            else
               !- linear interpolation
               prxd=(F_dstlev(k)-F_srclev(i,lev_lin))/deltalev
               F_dch(i,k) = (1.0-prxd)*F_sch(i,lev_lin) &
                    +prxd*F_sch(i,lev_lin+1)
            endif
         else
            !- cubic interpolation
            prxd = (F_dstlev(k)-F_srclev(i,lev_lin))/ &
                 deltalev
            prda = ((F_sch(i,lev_lin+1)-F_sch(i,lev_lin-1))/ &
                 (F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin-1))* &
                 deltalev)
            prdb = ((F_sch(i,lev_lin+2)-F_sch(i,lev_lin))/ &
                 (F_srclev(i,lev_lin+2)-F_srclev(i,lev_lin))* &
                 deltalev)
            prsaf= (1.0+2.0*prxd)*(1.0-prxd)*(1.0-prxd)
            prsbf= (3.0-2.0*prxd)*prxd*prxd
            prsad= prxd*(1.0-prxd)*(1.0-prxd)
            prsbd= (1.0-prxd)*prxd*prxd
            F_dch(i,k) = F_sch(i,lev_lin  )*prsaf &
                 +F_sch(i,lev_lin+1)*prsbf+prda*prsad &
                 -prdb*prsbd
         endif
      enddo
   enddo

   !- special case for VT
   if(F_var_S == 'VT') then
      do k=1,nkd
         do i=1,n
            if(F_srclev(i,k_ciel).lt.F_dstlev(k)) then
               F_dch(i,k) = F_sch(i,k_ciel) * exp (  &
                    rgasd_8*stlo_8*(F_dstlev(k)-F_srclev(i,k_ciel)) )
            endif
         enddo
      enddo
   endif
   !     ---------------------------------------------------------------
   return
end subroutine vte_intvertx_isodst



