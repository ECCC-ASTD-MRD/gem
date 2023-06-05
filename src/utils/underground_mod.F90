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

!/@*
module underground_mod
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   private
   public :: vt_on_lieb_levels, mslp, gz_vt_on_pres
   !@Author 2018, Andre Plante
   !*@/

   include "rmnlib_basics.inc"

contains

   !/@*
   subroutine vt_on_lieb_levels(ttx, vt, gz, fis0, wlao,&
        lieb_levels, lieb_conv, lieb_maxite, &
        Minx,Maxx,Miny,Maxy,ni,nj, nkund, nk, &
        liebxch_iter, G_halox, G_haloy, G_periodx, G_periody )
      !@Object 
      !@Arguments
      integer, intent(in) :: lieb_maxite, &
           Minx,Maxx,Miny,Maxy,ni,nj, nkund, nk
      real, intent(inout) :: ttx (Minx:Maxx,Miny:Maxy,nkund)
      real, intent(in) :: &
           vt  (Minx:Maxx,Miny:Maxy,nk   ), &
           gz  (Minx:Maxx,Miny:Maxy,nk   ), &
           fis0(Minx:Maxx,Miny:Maxy), &
           wlao(Minx:Maxx,Miny:Maxy),&
           lieb_levels(nkund), lieb_conv
      integer, intent(in), optional :: liebxch_iter, G_halox, G_haloy
      logical, intent(in), optional :: G_periodx, G_periody
      !*@/

      real :: mask(Minx:Maxx,Miny:Maxy,nkund)
      integer l_liebxch_iter, l_halox, l_haloy
      logical :: l_periodx, l_periody
      !----------------------------------------------------------------------

      l_liebxch_iter=1; l_halox=0; l_haloy=0
      l_periodx=.false.; l_periody=.false.
      if(present(liebxch_iter)) l_liebxch_iter = liebxch_iter
      if(present(G_halox)) l_halox = G_halox
      if(present(G_haloy)) l_haloy = G_haloy
      if(present(G_periodx)) l_periodx = G_periodx
      if(present(G_periody)) l_periody = G_periody

      call vt_trial(ttx, mask, vt, gz, fis0, wlao,&
           lieb_levels, Minx,Maxx,Miny,Maxy,ni,nj, nkund, nk )

      call fill_halo (ttx,Minx,Maxx,Miny,Maxy,ni,nj,nkund)

      call liebman_comm (ttx,mask,lieb_conv,lieb_maxite,&
           Minx,Maxx,Miny,Maxy,ni,nj,nkund, &
           l_liebxch_iter, l_halox, l_haloy ,&
           l_periodx, l_periody )
      return
   end subroutine vt_on_lieb_levels


   !/@*
   subroutine vt_trial(ttx, mask, vt, gz, fis0, wlao,&
        lieb_levels, Minx,Maxx,Miny,Maxy,ni,nj,nkund,nk )

      ! Obtain vt on liebman levels
      use tdpack, only: grav_8, stlo_8, pi_8
      implicit none
!!!#include <arch_specific.hf>
      !@Object 
      !@Arguments
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,ni,nj, nkund, nk
      real, intent(inout) :: ttx (Minx:Maxx,Miny:Maxy,nkund), &
           mask(Minx:Maxx,Miny:Maxy,nkund  )
      real, intent(in) :: &
           vt  (Minx:Maxx,Miny:Maxy,nk   ), &
           gz  (Minx:Maxx,Miny:Maxy,nk   ), &
           fis0(Minx:Maxx,Miny:Maxy), &
           wlao(Minx:Maxx,Miny:Maxy),&
           lieb_levels(nkund)
      !@Author Michel Desgagne - Fall 2011 (out_liebman.F90)
      !@Revision
      !     Andre Plante    - June 2018, generalize for cases with 3D vt,gz or just
      !                                  surface vt and gz (nk = 1)
      !*@/
      integer :: i,j,k,kk,kgrnd
      real :: grad, htx(nkund)
      external :: fil_halo
      !----------------------------------------------------------------------
!$omp parallel private (grad, kgrnd) shared (htx)
!$omp do

      do k=1,nkund

         !        Store fictitious height level in htx
         htx(k) = lieb_levels(k) * grav_8

         do j=1,nj
            do i=1,ni

               ! Determine if fictitious level is above or below ground
               ttx(i,j,k) = fis0(i,j) - htx(k)

               if ( ttx(i,j,k) <= 0 .and. nk > 1 )then

                  ! Fictitious level is above ground and vt is 3D.
                  ! Temperature is obtained by linear INTerpolation.
                  ! Identify above ground grid point.
                  mask(i,j,k) = 0.0

                  do kk= nk, 1, -1
                     kgrnd = kk
                     ttx(i,j,k) = gz (i,j,kk) - htx(k)
                     if ( ttx(i,j,k) > 0. ) goto 10
                  enddo
10                continue
                  if( kgrnd == nk )then
                     grad = 0.
                  else
                     grad = - (vt(i,j,kgrnd) - vt(i,j,kgrnd+1) ) / &
                          (gz(i,j,kgrnd) - gz(i,j,kgrnd+1) )
                  endif
                  ttx(i,j,k) = vt (i,j,kgrnd) + grad * ttx(i,j,k)

               else

                  ! Fictitious level is under ground or
                  ! at/above ground with only surface data (nk = 1).
                  ! Temperature is obtained by linear EXTrapolation

                  ! Identify under and above/at ground grid point
                  if( ttx(i,j,k) > 0. )then
                     mask(i,j,k) = 1.0
                  else
                     mask(i,j,k) = 0.0
                  endif

                  if ( abs( wlao(i,j)*180./pi_8 ) >= 49. ) then
                     ttx(i,j,k) = vt(i,j,Nk) +  .0005 * ttx(i,j,k)
                  else
                     ttx(i,j,k) = vt(i,j,Nk) + stlo_8 * ttx(i,j,k)
                  endif

               endif

            end do
         end do

      end do

!$omp enddo
!$omp end parallel
      !----------------------------------------------------------------------
      return
   end subroutine vt_trial


   !/@*
   subroutine liebman_comm (F_field,  F_mask, F_conv, F_maxite, &
        Minx, Maxx, Miny, Maxy, ni, nj, Nk, &
        liebxch_iter, G_halox, G_haloy, G_periodx, G_periody)
      implicit none
!!!#include <arch_specific.hf>
      !@Object 
      !@Arguments
      integer, intent(in) :: F_maxite, Minx, Maxx, Miny, Maxy, NK,&
           ni, nj
      real, intent(in) :: F_conv
      real, intent(inout) :: F_field(Minx:Maxx,Miny:Maxy,NK)
      real, intent(in) ::  F_mask   (Minx:Maxx,Miny:Maxy,NK)
      integer, intent(in), optional :: liebxch_iter, G_halox, G_haloy
      logical, intent(in), optional :: G_periodx, G_periody
      !@author Michel Desgagne - after version v1_03 of liebman.ftn
      !*@/
      integer i,j,k,ite,i0,in,ic,j0,jn,jc,count,ier
      integer l_liebxch_iter, l_halox, l_haloy
      real    prfact,prmax(Nk),prmaxall(Nk),prmod
      logical :: comm_L, l_periodx, l_periody
      !----------------------------------------------------------------------
      l_liebxch_iter=1; l_halox=0; l_haloy=0
      l_periodx=.false.; l_periody=.false.

      if(present(liebxch_iter)) l_liebxch_iter = liebxch_iter
      if(present(G_halox)) l_halox = G_halox
      if(present(G_haloy)) l_haloy = G_haloy
      if(present(G_periodx)) l_periodx = G_periodx
      if(present(G_periody)) l_periody = G_periody

      comm_L = l_halox > 0 .or. l_haloy > 0
      if(comm_L) call rpn_comm_xch_halo(F_field, Minx,Maxx,Miny,Maxy,ni,nj,Nk, &
           l_halox,l_haloy,l_periodx,l_periody,ni,0)

      prfact   = 1.75 * 0.25
      prmaxall = 1000.

!$omp parallel shared (i0,in,ic,j0,jn,jc,prmax,prmaxall,prfact, &
!$omp                  count)          &
!$omp          private (i,j,k,ite,prmod)

!$omp single
      i0 = 1 ; in = ni ; ic = 1
      j0 = 1 ; jn = nj ; jc = 1
!$omp end single

      do ite=1,F_maxite
!$omp do
         do k=1,Nk
            prmax(k) = 0.0
         end do
!$omp enddo

!$omp do
         do k=1,Nk
            if ( prmaxall(k) > F_conv ) then
               do j=j0,jn,jc
                  do i=i0,in,ic
                     prmod = prfact * F_mask(i,j,k)               * &
                          (F_field(i-1,j,k) + F_field(i+1,j,k) + &
                          F_field(i,j-1,k) + F_field(i,j+1,k) - &
                          4.*F_field(i,j,k))
                     prmax(k) = max ( prmax(k), abs(prmod) )
                     F_field(i,j,k) = F_field(i,j,k) + prmod
                  enddo
               enddo
            endif
         enddo
!$omp enddo

!$omp single
         if ((i0 == 1).and.(j0 == 1)) then
            i0=ni ; in=1 ; ic=-1
         endif
         if ((i0 == ni).and.(j0 == 1)) then
            j0=nj ; jn=1 ; jc=-1
         endif
         if ((i0 == ni).and.(j0 == nj)) then
            i0 = 1 ; in = ni ; ic =  1
            j0 = nj ; jn = 1 ; jc = -1
         endif
         if ((i0 == 1).and.(j0 == nj)) then
            i0 = 1 ; in = ni ; ic = 1
            j0 = 1 ; jn = nj ; jc = 1
         endif
         if(comm_L)then
            call rpn_comm_ALLREDUCE (prmax,prmaxall,Nk,"MPI_REAL","MPI_MAX", &
                 "grid",ier)
         else
            prmaxall = prmax
         endif
         count = 0
         do k=1,Nk
            if ( prmaxall(k) < F_conv ) count = count + 1
         end do
!$omp end single

         if ( count == Nk ) exit

!$omp single
         if (mod(ite,liebxch_iter) == 0 .and. comm_L ) &
              call rpn_comm_xch_halo( F_field, Minx,Maxx,Miny,Maxy,ni,nj,Nk, &
              l_halox,l_haloy,l_periodx,l_periody,ni,0)
!$omp end single

      end do
!$omp end parallel
      !----------------------------------------------------------------------
      return
   end subroutine liebman_comm


   !/@*
   subroutine mslp( F_pnm,   F_vts,   F_fis, F_lnps, F_la, &
        F_vtund, F_lieb_levels, F_und, &
        Minx,Maxx,Miny,Maxy,ni,nj, Nk)
      use tdpack, only: grav_8, rgasd_8, stlo_8, pi_8
      implicit none
!!!#include <arch_specific.hf>
      !@Arguments
      !  Name        I/O                 Description
      !----------------------------------------------------------------
      ! F_pnm        O    - MSL pressure
      ! F_vts        I    - surface virtual temperature
      ! F_fis        I    - surface geopotential height
      ! F_lnps       I    - surface log of hydrostatic pressure
      ! F_la         I    - geographical latitude (radian)
      ! F_vtund      I    - virtual temperatures for underground extrapolation
      ! F_und        I    - number of virtual temperature levels for underground
      !                     extrapolation
      !                   = 0 if no underground temperature is used and the
      !                       the traditional scheme will be used
      integer Minx,Maxx,Miny,Maxy,ni,nj, Nk
      real F_pnm(Minx:Maxx,Miny:Maxy), F_vts(Minx:Maxx,Miny:Maxy,Nk)
      real F_fis(Minx:Maxx,Miny:Maxy,Nk)
      real F_lnps(Minx:Maxx,Miny:Maxy,Nk), F_la (Minx:Maxx,Miny:Maxy)
      integer F_und
      real    F_vtund(Minx:Maxx,Miny:Maxy,F_und),F_lieb_levels(F_und)
      !@author andre methot - alain patoine - after pnm1
      !      !@revision
      ! v2_00 - Lee V.            - initial MPI version (from pnm2 v1_03)
      ! v3_00 - Desgagne & Lee    - Lam configuration
      !@object
      !******************************************************************************
      !                                                                             *
      ! The hypsometric equation is used:                                           *
      !                                                                             *
      !                        /    \                                               *
      !                        | p  |                                               *
      !                   _    |  t |                                               *
      ! fi  - fi   = - R  T ln |----|                                          (1)  *
      !   t     b       d      | p  |                                               *
      !                        |  b |                                               *
      !                        \    /                                               *
      !                                                                             *
      ! Here the subscript t and b stand respectively for top and bottom of the     *
      ! considered layer.                                                           *
      !                                               dT                            *
      ! We consider a constant temperature lapse rate --- = - L                     *
      !                                           _   dfi                           *
      ! (e.g. L = STLO) and use the definition of T:                                *
      !                                                                             *
      !         /                \                                                  *
      !         |   fi  - fi     |                                                  *
      ! _       |     t     b    |                                                  *
      ! T = - L |----------------| ,                                           (2)  *
      !         |    / T   /   \ |                                                  *
      !         | ln |  t / T  | |                                                  *
      !         |    \   /   b / |                                                  *
      !         \                /                                                  *
      !                                                                             *
      ! into expression (1) and get an expression for p :                           *
      !                                                b                            *
      !                                                                             *
      !             /    / T   /   \ \                                              *
      !             | ln |  b / T  | |                                              *
      !             |    \   /   t / |                                              *
      ! p  = p  exp | -------------- |                                              *
      !  b    t     |      R  L      |                                              *
      !             |       d        |                                              *
      !             \                /                                              *
      !                                                                             *
      ! In the case where L -> 0, we have to revert to expression (1) in which      *
      !        _                                                                    *
      ! we use T = T .                                                              *
      !             t                                                               *
      !                                                                             *
      ! At points where we want to use underground temperatures for calculation,    *
      ! we recursively compute the pressure at the bottom of each layer.            *
      !                                                                             *
      ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      !                                                                             *
      ! The temperature lapse rate in each virtual layer is either computed using   *
      ! the provided temperatures or assumed to be the Shumman-Newel lapse rate     *
      ! for the layer near ground except, when the temperatures exceeds given       *
      ! critical values:                                                            *
      !                                                                             *
      ! T  = 301.75 - lat / 4                                                       *
      !  c                                                                          *
      !                                                                             *
      ! If T  is lower than T , the algorithm ensures that the bottom temperature   *
      !     t                c                                                      *
      ! is not greater than T .                                                     *
      !                      c                                                      *
      ! Else, if T  exceeds T , then the bottom temperature is set to:              *
      !           t          c                                                      *
      !                                                                             *
      !                              2                                              *
      ! T  = T  - 0.005 * ( T  - T  ) .                                             *
      !  b    c              t    c                                                 *
      !                                                                             *
      !******************************************************************************
      !
      !@notes
      !   All fields in arguments are assumed to be workable on the same grid
      !   (fni x fnj). This grid could be the staggered or the non staggered.
      !
      !   It is important that the data stored in F_vtund and fiund be ordered
      !   in the proper manner:
      !   F_vtund(i,j,1) and F_fiund(1) --> highest level
      !   F_vtund(i,j,2) and F_fiund(2) --> second highest level
      !   ......................................and so on
      !*@/
      integer i, j, k, pnund,   pn1
      real    prl, prvtc, prsmall
      real    prlptop, prvttop, prfitop
      real    prlpbot, prvtbot, prfibot
      real, dimension(F_und) :: fiund
      !
      !
      do k=1, F_und
         fiund(k) = F_lieb_levels(k) * grav_8
      end do

      prsmall = .001
      !
      do j= 1, nj
         do i= 1, ni
            !
            !        calculation of critical temperature
            !                       --------------------
            !
            prvtc = 301.75 - abs( (F_la(i,j) * 180.) / ( 4. * pi_8) )
            !
            !
            do pnund=1,F_und+1
               if ( pnund > F_und ) go to 30
               if ( F_fis(i,j,nk) > fiund(pnund) ) go to 30
            enddo
            !
30          continue
            !
            prlptop = F_lnps(i,j,nk)
            prvttop = F_vts(i,j,nk)
            prfitop = F_fis(i,j,nk)
            !
            do pn1=pnund,F_und
               !
               if ( prvttop <= prvtc ) then
                  prvtbot = min( F_vtund(i,j,pn1),  prvtc )
               else
                  prvtbot = prvtc - 0.005 * ( prvttop - prvtc ) **2
               endif
               !
               prfibot  = fiund (pn1)
               !
               if ( abs(prvtbot-prvttop) <= prsmall ) then
                  prlpbot = prlptop + (prfitop-prfibot)/(rgasd_8*prvttop)
               else
                  prl     = - ( prvttop - prvtbot ) / ( prfitop - prfibot )
                  prlpbot = prlptop + (log(prvtbot/prvttop)) / (rgasd_8*prl)
               endif
               !
               prlptop = prlpbot
               prvttop = F_vtund(i,j,pn1)
               prfitop = prfibot
               !
            end do
            !
            if ( prvttop <= prvtc ) then
               prvtbot = min( 1.0d0*prvttop + stlo_8 * 1.0d0*prfitop, 1.0d0*prvtc)
            else
               prvtbot = prvtc - 0.005 * ( prvttop - prvtc ) **2
            endif
            !
            !        calculation of MSL pressure
            !                       ------------
            !
            if ((abs(prvtbot-prvttop) <= prsmall) .or. (prfitop <= 0.0)) then
               F_pnm(i,j) = exp (prlptop+prfitop/(rgasd_8*prvttop))
            else
               prl = - ( prvttop - prvtbot ) / ( prfitop )
               F_pnm(i,j)=exp (prlptop+(log(prvtbot/prvttop))/(rgasd_8*prl))
            endif
            !
         end do
      end do
      !
      return
   end subroutine mslp


   !/@*
   subroutine gz_vt_on_pres(F_gzout, F_vtout,  F_pres, Nkout, &
        F_gzin,  F_vtin,   F_wlnph, F_la, &
        F_vtund, F_zund,   F_nundr, &
        F_cubzt_L, F_linbot, &
        Minx,Maxx,Miny,Maxy,ni,nj,F_Nk)
      use tdpack, only: grav_8, rgasd_8, stlo_8, pi_8
      implicit none
!!!#include <arch_specific.hf>
      !@Object interpolation of geopotential and virtual temperature on
      !        given pressure levels
      !@Arguments
      logical F_cubzt_L
      integer F_nundr, F_linbot
      integer Minx,Maxx,Miny,Maxy,ni,nj,F_Nk,Nkout
      real    F_pres(Nkout)
      real    F_gzout(Minx:Maxx,Miny:Maxy,Nkout),F_vtout(Minx:Maxx,Miny:Maxy,Nkout)
      real    F_gzin (Minx:Maxx,Miny:Maxy,F_Nk), F_vtin (Minx:Maxx,Miny:Maxy,F_Nk)
      real    F_wlnph (Minx:Maxx,Miny:Maxy,F_Nk), F_la   (Minx:Maxx,Miny:Maxy)
      real    F_vtund(Minx:Maxx,Miny:Maxy,F_nundr), F_zund(F_nundr)
      !*@/
      integer i, j, k, kk, pnk, pnkm, pnindex(ni), pnund,   pn1
      real    prlprso
      real    prd, pre, prr
      real    prfm0, prfm1, prfm2, prfm3, prfl2
      real    prl, prsmall
      real    prlptop, prvttop, prfitop
      real    prlpbot, prvtbot, prfibot
      real    logpres(Nkout)
      real, dimension(F_nundr) :: gzund
      real(REAL64) :: invprd
      !-------------------------------------------------------------------
      prsmall= .001
      logpres(1:Nkout) = log(F_pres(1:Nkout))

      do k=1,F_nundr
         gzund(k) = F_zund(k) * grav_8
      end do

!$omp parallel private(i,k,kk,pnk,pnkm,pnindex,prlprso, &
!$omp    prd,pre,prr,prfm0,prfm1,prfm2,prfm3,prfl2,prl, &
!$omp    pnund,pn1,prlptop,prvttop,prfitop, &
!$omp    prlpbot,prvtbot,prfibot,invprd) &
!$omp          shared (logpres)

!$omp do
      do j  = 1, nj
         do kk = 1, Nkout
            pnindex = 0
            prlprso = logpres(kk)

            do k= 1, F_nk
               do i= 1, ni
                  if ( prlprso > F_wlnph(i,j,k) ) pnindex(i) = k
               enddo
            enddo

            do 300 i= 1, ni
               !**************************************************************
               ! If:    output pressure   <   hydrostatic pressure on the
               !                              first level of the model
               ! Then:  upward extrapolation
               !**************************************************************

               if ( pnindex(i) == 0 ) then
                  prd = prlprso - F_wlnph(i,j,1)
                  F_vtout(i,j,kk) = F_vtin(i,j,1) + prd &
                       * (F_vtin(i,j,1)-F_vtin(i,j,2)) &
                       / (F_wlnph(i,j,1)-F_wlnph(i,j,2))

                  F_gzout(i,j,kk) = F_gzin(i,j,1) - prd * rgasd_8 &
                       * (F_vtin(i,j,1) + F_vtout(i,j,kk)) * 0.5

                  !***********************************************************
                  ! If:    output pressure   >   hydrostatic pressure on the
                  !                              last level of the model
                  ! Then:  downward extrapolation
                  ! The hypsometric equation is used:
                  !                         /    \
                  !                         | p  |
                  !                    _    |  t |
                  !  fi  - fi   = - R  T ln |----|   (1)
                  !    t     b       d      | p  |
                  !                         |  b |
                  !                         \    /
                  !  Here the subscript t and b stand respectively for top
                  !  and bottom of theconsidered layer.
                  !                                                dT
                  !  We consider a constant temperature lapse rate --- = - L
                  !                                            _   dfi
                  !  (e.g. L = STLO) and use the definition of T:
                  !          /                \
                  !          |   fi  - fi     |
                  !  _       |     t     b    |
                  !  T = - L |----------------| ,   (2)
                  !          |    / T   /   \ |
                  !          | ln |  t / T  | |
                  !          |    \   /   b / |
                  !          \                /
                  !  into expression (1) and get an expression for T :
                  !                                                 b
                  !               /                     \
                  !               |         / p   /   \ |
                  !  T  = T   exp | R  L ln |  b / p  | |   (3)
                  !   b    t      |  d      \   /   t / |
                  !               \                     /
                  !  Then, we use the definition of L, to get an expression for fi :
                  !                                                               b
                  !              /         \
                  !              | T  - T  |
                  !              \  t    b /
                  !  fi  = fi  + -----------    (4)
                  !    b     t        L
                  !
                  ! In the case where L -> 0, we have to revert to expression (1) in which
                  !        _
                  ! we use T = T .
                  !             t
                  !
                  ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                  !
                  ! At points where we want to use underground temperatures for extrapolation,
                  ! we first determine the layer bottom pressure using (3) rearranged:
                  !
                  !             /    / T   /   \ \
                  !             | ln |  b / T  | |
                  !             |    \   /   t / |
                  ! p  = p  exp | -------------- |
                  !  b    t     |      R  L      |
                  !             |       d        |
                  !             \                /
                  !
                  ! In the case where L -> 0, we have to revert to the expression (1) in which
                  !        _
                  ! we use T = T .
                  !             t
                  !
                  ! Then, if the layer bottom pressure is larger than the destination pressure,
                  ! we proceed with calculation (3) and (4). Otherwise, we update the variables
                  ! at top and bottom for next layer calculation and iterate.
                  !
                  !***********************************************************
               else if ( pnindex(i) == F_nk ) then

                  do pnund=1,F_nundr
                     if ( F_gzin(i,j,F_nk) > gzund(pnund) ) go to 30
                  enddo
30                prlptop = F_wlnph(i,j,F_nk)
                  prvttop = F_vtin (i,j,F_nk)
                  prfitop = F_gzin (i,j,F_nk)

                  do pn1=pnund,F_nundr

                     prvtbot = F_vtund (i,j,pn1)
                     prfibot = gzund (pn1)

                     if ( abs(prvtbot-prvttop) <= prsmall ) then
                        prlpbot = prlptop + (prfitop-prfibot)/(rgasd_8*prvttop)
                        if ( prlpbot >= prlprso ) then
                           F_vtout(i,j,kk) = prvttop
                           F_gzout(i,j,kk) = prfitop + &
                                rgasd_8*prvttop*(prlptop-prlpbot)
                           go to 300
                        endif
                     else
                        prl     = - ( prvttop - prvtbot ) / ( prfitop - prfibot )
                        prlpbot = prlptop + (log(prvtbot/prvttop)) / (rgasd_8*prl)
                        if ( prlpbot >= prlprso ) then
                           F_vtout(i,j,kk) = prvttop * &
                                exp ( rgasd_8 * prl * (prlprso-prlptop))
                           F_gzout(i,j,kk) = prfitop + &
                                (prvttop-F_vtout(i,j,kk)) / prl
                           go to 300
                        endif
                     endif

                     prlptop = prlpbot
                     prvttop = prvtbot
                     prfitop = prfibot
                  end do

                  prl = stlo_8
                  if ( abs (F_la(i,j)*180./pi_8) >= 49.0 ) prl = .0005
                  F_vtout(i,j,kk) = prvttop * &
                       exp ( rgasd_8 * prl * (prlprso-prlptop))
                  F_gzout(i,j,kk) = prfitop + (prvttop-F_vtout(i,j,kk)) / prl

                  !***********************************************************
                  !  Else, interpolate between appropriate levels
                  !***********************************************************
               else
                  !        ***************************************************
                  !        * NOTE ABOUT "F_linbot"
                  !        * this parameter is used to force a linear interpolation in a
                  !        * certain number of layers (equal to F_linbot) close to the bottom
                  !        * of the model even if F_cubzt_L is .true.
                  !        * it has no effect if F_cubzt_L is .false.
                  !        ***************************************************
                  pnkm = pnindex(i)
                  pnk  = pnindex(i) + 1
                  prd  = F_wlnph(i,j,pnk) - F_wlnph(i,j,pnkm)
                  invprd = 1.0/prd
                  pre = prlprso - 0.5 * ( F_wlnph(i,j,pnk) + F_wlnph(i,j,pnkm) )

                  if ( F_cubzt_L .and. ( pnk < F_nk+1-F_linbot ) ) then
                     prr = 0.125 * prd * prd - 0.5 * pre * pre
                     prfm0 = 0.5 * ( F_gzin(i,j,pnk) + F_gzin(i,j,pnkm) )
                     prfm1 = ( F_gzin(i,j,pnk) - F_gzin(i,j,pnkm) ) * invprd
                     prfm2 = - rgasd_8 &
                          * ( F_vtin(i,j,pnk) - F_vtin(i,j,pnkm) ) * invprd
                     prfm3 = - rgasd_8 * ( F_vtin(i,j,pnk) + F_vtin(i,j,pnkm) )
                     prfm3 = ( prfm3 - prfm1 - prfm1 ) * invprd * invprd
                     prfl2 = prfm2 + 2.0 * pre * prfm3
                     F_gzout(i,j,kk) = prfm0 + pre * prfm1 - prr * prfl2
                     F_vtout(i,j,kk) = prfm1 + pre * prfl2 - 2.0 * prr * prfm3
                     F_vtout(i,j,kk) = - F_vtout(i,j,kk) / rgasd_8
                  else
                     prfm0 = 0.5 * ( F_gzin(i,j,pnk) + F_gzin(i,j,pnkm) )
                     prfm1 = ( F_gzin(i,j,pnk) - F_gzin(i,j,pnkm) ) * invprd
                     F_gzout(i,j,kk)= prfm0 + pre * prfm1
                     prfm0 = 0.5 * ( F_vtin(i,j,pnk) + F_vtin(i,j,pnkm) )
                     prfm1 = ( F_vtin(i,j,pnk) - F_vtin(i,j,pnkm) ) * invprd
                     F_vtout(i,j,kk)= prfm0 + pre * prfm1
                  endif
               endif
300         end do
         end do
      end do
!$omp enddo
!$omp end parallel
      !-------------------------------------------------------------------
      return
   end subroutine gz_vt_on_pres

end module underground_mod

