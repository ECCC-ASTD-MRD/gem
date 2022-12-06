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

module vertical_interpolation
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use clib_itf_mod, only: clib_toupper
  use tdpack, only: GRAV_8, RGASD_8, STLO_8
  implicit none
!!!#include <arch_specific.hf>
  private
  public :: vertint2

  include "rmnlib_basics.inc"

contains

      subroutine vertint2 ( F_dch,F_dstlev,nkd, F_sch,F_srclev,nks , &
                            Minx,Maxx,Miny,Maxy,F_i0,F_in,F_j0,F_jn, &
                            varname, inttype, Schumann_list , levtype)
      implicit none

!!!#include <arch_specific.hf>

      character(len=*), optional, intent(IN) :: &
                        varname, inttype, Schumann_list, levtype
      integer         ,           intent(IN) :: &
                        Minx,Maxx,Miny,Maxy,F_i0,F_in,F_j0,F_jn,nkd, nks
      real            ,           intent(IN) :: &
                        F_sch   (Minx:Maxx,Miny:Maxy,nks), &
                        F_srclev(Minx:Maxx,Miny:Maxy,nks), &
                        F_dstlev(Minx:Maxx,Miny:Maxy,nkd)
      real            ,           intent(OUT) :: &
                        F_dch   (Minx:Maxx,Miny:Maxy,nkd)

   !@ojective vertical interpolation logp to logp

   !@arguments
   ! F_dch                   Interpolated field
   ! F_sch                   source       field
   ! F_srclev                source levels (logp)
   ! F_dstlev                destination levels (logp)
   ! F_inttype_S             order of interpolation ('linear' or 'cubic')
   ! F_Schumann_list_S       Schumann extrapolation

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
   !  niter calculation is ok for nks < 2097152: should be ok for now...
   !  (Maybe the grid will be that precise in 2010!) (I don't even bother
   !  to add an if statement, for performance purpose...)

      integer, parameter :: MAXLIST=20
      character(len=1) :: levtype_S
      character(len=12)  Schumann_list_S(MAXLIST)
      character(len=521) tmp_list,name
      logical :: ascending_L, flag
      integer :: i,j,k,iter,niter,lev,lev_lin,k_surf, &
                 k_ciel,nlinbot,n_schum,istat,deb,fin,coma
      integer, dimension (Minx:Maxx,Miny:Maxy) :: top, bot, &
                                                  topcub, botcub, ref
      real(REAL64) :: deltalev,prxd,prda,prdb,prsaf,prsbf, &
                 prsad,prsbd
!
!-------------------------------------------------------------------
!
      if (nks<2) then
         do k=1,nkd
            F_dch(:,:,k) = F_sch(:,:,nks)
         end do
         return
      endif
      
      levtype_S='P'
      if(present(levtype))then
         select case(trim(levtype))
         case('P')
            levtype_S='P'
         case('H')
            levtype_S='H'
         end select
      end if
      nlinbot=0
      if (present(inttype)) then
         if (any(trim(inttype) == (/'linear','LINEAR'/))) nlinbot=nkd
      end if

      n_schum= 3
      Schumann_list_S(1) = 'TT'
      Schumann_list_S(2) = 'VT'
      Schumann_list_S(3) = 'PW_TT:P'
      Schumann_list_S(4:MAXLIST) = '#@!$%&'

      if (present(Schumann_list)) then
         tmp_list= Schumann_list
         istat = clib_toupper(tmp_list)
         deb=1 ; fin=len(trim(tmp_list))
 55      tmp_list= tmp_list(deb:fin)
         if ( ( fin > 0 ) .and. ( deb <= fin ) )then
            coma= index (tmp_list,',') - 1
            if (coma < 0) coma= fin
            name= tmp_list(1:coma)
            deb= coma+2
            coma= index (name,':')
            if ((coma > 0).and.(name(1:coma-1)=='REMOVE')) then
               name= name(coma+1:)
               flag= .false. ; n_schum= n_schum - 1
               do i=1,n_schum
                  if (Schumann_list_S(i) == name) flag= .true.
                  if (flag) Schumann_list_S(i)= Schumann_list_S(i+1)
               end do
               Schumann_list_S(n_schum+1:)= '#@!$%&'
            else
               if ( n_schum + 1 <= MAXLIST ) then
                  n_schum= n_schum + 1
                  Schumann_list_S(n_schum)= name
               end if
            end if
            goto 55
         end if
      end if

      ascending_L = (F_srclev(F_i0,F_j0,1) < F_srclev(F_i0,F_j0,nks))
      k_surf = 1
      k_ciel = nks
      if (.not.ascending_L) then
         k_surf = nks
         k_ciel = 1
      end if
      if (real(int(log(real(nks))/log(2.0))) == &
      log(real(nks))/log(2.0)) then
      niter=int(log(real(nks))/log(2.0))
      else
         niter=int(log(real(nks))/log(2.0))+1
      end if

!$omp parallel private(i,j,k,iter,lev,lev_lin,top,bot,topcub,botcub,ref, &
!$omp                  deltalev,prxd,prda,prdb,prsaf,prsbf,prsad,prsbd)  &
!$omp    shared(ascending_L,nlinbot,k_surf,k_ciel,niter)

!$omp do
      do k=1,nkd
         top=nks
         bot=1
         if (ascending_L) then
            do iter=1,niter
               do j= F_j0, F_jn
               do i= F_i0, F_in
!     divide by two (the old fashioned way...)
                  ref(i,j)=ishft(top(i,j)+bot(i,j),-1)
!     adjust top or bot
                  if(F_dstlev(i,j,k) < F_srclev(i,j,ref(i,j))) then
                     top(i,j)=ref(i,j)
                  else
                     bot(i,j)=ref(i,j)
                  end if
               end do
               end do
            end do
         else
            do iter=1,niter
               do j= F_j0, F_jn
               do i= F_i0, F_in
!     divide by two (the old fashioned way...)
                  ref(i,j)=ishft(top(i,j)+bot(i,j),-1)
!     adjust top or bot
                  if(F_dstlev(i,j,k) > F_srclev(i,j,ref(i,j))) then
                     top(i,j)=ref(i,j)
                  else
                     bot(i,j)=ref(i,j)
                  end if
               end do
               end do
            end do
         end if
!- adjusting top and bot to ensure we can perform cubic interpolation
         do j= F_j0, F_jn
         do i= F_i0, F_in
            botcub(i,j)=max(    2,bot(i,j))
            topcub(i,j)=min(nks-1,top(i,j))
         end do
         end do
!- cubic or linear interpolation
         do j= F_j0, F_jn
         do i= F_i0, F_in
            lev=botcub(i,j)
            lev_lin=bot(i,j)
            deltalev=(F_srclev(i,j,lev_lin+1)-F_srclev(i,j,lev_lin))

!- Interpolation: either if not enough points to perform cubic interp
!                 or if linear interpolation is requested use linear
!                 interpolation and constant extrapolation
            if((lev /= lev_lin).or.(topcub(i,j) /= top(i,j)).or.k<=nlinbot) then
!- persistancy of this interval
               if(F_dstlev(i,j,k) <= F_srclev(i,j,k_surf)) then
                  F_dch(i,j,k) = F_sch(i,j,k_surf)
               else if(F_dstlev(i,j,k) >= F_srclev(i,j,k_ciel)) then
                  F_dch(i,j,k) = F_sch(i,j,k_ciel)
               else
!- linear interpolation
                  prxd=(F_dstlev(i,j,k)-F_srclev(i,j,lev_lin))/deltalev
                  F_dch(i,j,k) = (1.0-prxd)*F_sch(i,j,lev_lin  ) &
                  +prxd *F_sch(i,j,lev_lin+1)
               end if
            else
!- cubic interpolation
               prxd = (F_dstlev(i,j,k)-F_srclev(i,j,lev_lin)) / deltalev
               prda = ((F_sch   (i,j,lev_lin+1)-F_sch   (i,j,lev_lin-1))/ &
                       (F_srclev(i,j,lev_lin+1)-F_srclev(i,j,lev_lin-1))* &
                        deltalev)
               prdb = ((F_sch    (i,j,lev_lin+2)-F_sch   (i,j,lev_lin))/ &
                       (F_srclev (i,j,lev_lin+2)-F_srclev(i,j,lev_lin))* &
                        deltalev)
               prsaf= (1.0+2.0*prxd)*(1.0-prxd)*(1.0-prxd)
               prsbf= (3.0-2.0*prxd)*prxd*prxd
               prsad= prxd*(1.0-prxd)*(1.0-prxd)
               prsbd= (1.0-prxd)*prxd*prxd
               F_dch(i,j,k) =  F_sch(i,j,lev_lin  )*prsaf            &
                             + F_sch(i,j,lev_lin+1)*prsbf+prda*prsad &
                             - prdb*prsbd
            end if
         end do
         end do
      end do
!$omp end do

      if (present(varname)) then
      if (any(Schumann_list_S(1:n_schum)==trim(varname))) then
      if (levtype_S == 'P')then
!$omp do
         do k=1,nkd
            do j= F_j0, F_jn
            do i= F_i0, F_in
               if(F_srclev(i,j,k_ciel) < F_dstlev(i,j,k)) then
                  F_dch(i,j,k) = F_sch(i,j,k_ciel) * exp (  &
                  rgasd_8*stlo_8*(F_dstlev(i,j,k)-F_srclev(i,j,k_ciel)) )
               end if
            end do
            end do
         end do
!$omp end do
      elseif(levtype_S == 'H')then
!$omp do
         do k=1,nkd
            do j= F_j0, F_jn
            do i= F_i0, F_in
               if(F_srclev(i,j,k_surf) > F_dstlev(i,j,k)) then
                  F_dch(i,j,k) = F_sch(i,j,k_surf) + grav_8*stlo_8*(F_srclev(i,j,k_surf)-F_dstlev(i,j,k))
               end if
            end do
            end do
         end do
!$omp end do
      end if
      end if
      end if

!$omp end parallel

!
!-------------------------------------------------------------------
!
      return
      end subroutine vertint2

end module vertical_interpolation
