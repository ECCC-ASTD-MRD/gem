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

!**s/r ac_posi - Find the subgrid Grdc_gid,Grdc_gjd,Grdc_gif,Grdc_gjf
!                of the current model grid configuration to
!                extract for the next cascade grid

      subroutine ac_posi (xlon,ylat,dimgx,dimgy,prout)
      use HORgrid_options
      use gem_options
      use grdc_options
      use glb_ld
      use lun
      use tr3d
      use rstr
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: prout
      integer, intent(in) :: dimgx,dimgy
      real, intent(in) :: xlon(dimgx), ylat(dimgy)

      logical :: flag_hu,printout_L
      integer :: i,k,cnt,iref,jref,pil
      real :: latr,lonr
      real(kind=REAL64) :: x0, xl, y0, yl
      real(kind=REAL64), dimension(:), allocatable :: ac_xp, ac_yp
      real(kind=REAL64), dimension(dimgx) :: xpx
      real(kind=REAL64), dimension(dimgy) :: ypx
!
!---------------------------------------------------------------------
!

      Grdc_gid= 0 ; Grdc_gjd= 0
      Grdc_gif= 0 ; Grdc_gjf= 0

      printout_L= prout .and. (.not.Rstri_rstn_L)
      if (printout_L) write(output_unit, 1000)

      if ( (Grdc_ndt < 0) .or. (Grd_yinyang_S == 'YAN') .or. &
           (Grdc_ni == 0) .or. (Grdc_nj == 0) .or.           &
           ((Grdc_dx < 0.).and.(Grdc_dy < 0.)) ) then
         Grdc_ndt = -1
         goto 9999
      end if

      if (Grdc_fullgrid_L) then

         Grdc_gid= 1 ; Grdc_gif= dimgx
         Grdc_gjd= 1 ; Grdc_gjf= dimgy

      else

         if (Grdc_dy < 0) Grdc_dy = Grdc_dx
         if (Grdc_dx < 0) Grdc_dx = Grdc_dy

         xpx = dble(xlon)
         ypx = dble(ylat)

         pil = Grdc_maxcfl + Grd_bsc_base + Grd_bsc_ext1

         lonr = Grdc_lonr
         latr = Grdc_latr

         if ((Grdc_iref==-1) .and. (Grdc_jref==-1)) then
            iref = Grdc_ni / 2 + pil
            jref = Grdc_nj / 2 + pil

            if (mod(Grdc_ni,2)==0) then
               lonr = Grdc_lonr - dble(Grdc_dx)/2.d0
            else
               iref = iref + 1
            end if

            if (mod(Grdc_nj,2)==0) then
               latr = Grdc_latr - dble(Grdc_dy)/2.d0
            else
               jref = Grdc_nj / 2 + pil + 1
            end if
         else
            iref = Grdc_iref + pil
            jref = Grdc_jref + pil
         end if

         Grdc_ni = Grdc_ni + 2*pil
         Grdc_nj = Grdc_nj + 2*pil
         allocate(ac_xp(Grdc_ni),ac_yp(Grdc_nj))

         x0   = lonr - dble(iref-1+G_halox)    * dble(Grdc_dx)
         y0   = latr - dble(jref-1+G_haloy)    * dble(Grdc_dy)
         xl   = x0   + dble(Grdc_ni+2*G_halox) * dble(Grdc_dx)
         yl   = y0   + dble(Grdc_nj+2*G_haloy) * dble(Grdc_dy)

         if (x0 < 0.0) x0= x0 + 360.0d0
         if (xl < 0.0) xl= xl + 360.0d0

         call set_gemHgrid4 ( ac_xp, ac_yp, Grdc_ni, Grdc_nj, &
                              Grdc_dx, Grdc_dy, x0,xl,y0,yl, .false. )

         if (ac_xp(Grdc_ni) <= xpx(1    )) Grdc_ndt= -1
         if (ac_xp(1      ) >= xpx(dimgx)) Grdc_ndt= -1
         if (ac_yp(Grdc_nj) <= ypx(1    )) Grdc_ndt= -1
         if (ac_yp(1      ) >= ypx(dimgy)) Grdc_ndt= -1

         if (Grdc_ndt > 0) then

            do i=1,dimgx
               if (xpx(i) <= ac_xp(1)      ) Grdc_gid=i
               if (xpx(i) <= ac_xp(Grdc_ni)) Grdc_gif=i
            end do

            do i=1,dimgy
               if (ypx(i) <= ac_yp(1)      ) Grdc_gjd=i
               if (ypx(i) <= ac_yp(Grdc_nj)) Grdc_gjf=i
            end do

            Grdc_gid = Grdc_gid - 2
            Grdc_gjd = Grdc_gjd - 2
            Grdc_gif = Grdc_gif + 3
            Grdc_gjf = Grdc_gjf + 3

            if (printout_L) then
               if (Grdc_gid < 1)     write(output_unit,1015) 'WEST'
               if (Grdc_gjd < 1)     write(output_unit,1015) 'SOUTH'
               if (Grdc_gif > dimgx) write(output_unit,1015) 'EAST'
               if (Grdc_gjf > dimgy) write(output_unit,1015) 'NORTH'
            end if

            Grdc_gid = max(Grdc_gid, 1    )
            Grdc_gjd = max(Grdc_gjd, 1    )
            Grdc_gif = min(Grdc_gif, dimgx)
            Grdc_gjf = min(Grdc_gjf, dimgy)

         end if

         if ( ((Grdc_gif-Grdc_gid+1) < 2) .or. &
              ((Grdc_gjf-Grdc_gjd+1) < 2) ) then
              Grdc_ndt = -1
         end if

      end if

 9999 if (Grdc_ndt > 0) then

         if (Grdc_trnm_S(1) == '@#$%') then
            do i=1,Tr3d_ntr
               Grdc_trnm_S(i) = Tr3d_name_S(i)
            end do
            Grdc_ntr = Tr3d_ntr
         else
            cnt = 0 ;  flag_hu = .false.
            do k=1,max_trnm
               if (Grdc_trnm_S(k) == '@#$%') exit
               flag_hu= (trim(Grdc_trnm_S(k)) == 'HU')
               do i=1,Tr3d_ntr
                  if (trim(Grdc_trnm_S(k)) == trim(Tr3d_name_S(i))) then
                     cnt=cnt+1
                     Grdc_trnm_S(cnt) = Tr3d_name_S(i)
                     exit
                  end if
               end do
            end do

            if (.not.flag_hu) then
               cnt=cnt+1
               Grdc_trnm_S(cnt) = 'HU'
            end if
            Grdc_ntr = cnt
         end if

         if ((prout).and.(.not.Rstri_rstn_L)) then
            write (output_unit,1005) Grdc_gid,Grdc_gif,Grdc_gjd,Grdc_gjf,Grdc_ndt,Grdc_start,Grdc_end
            write (output_unit,1006) Grdc_gif-Grdc_gid+1,xpx(Grdc_gid),xpx(Grdc_gif),&
                           Grdc_gjf-Grdc_gjd+1,ypx(Grdc_gjd),ypx(Grdc_gjf)
            write (output_unit,1001)
            write (output_unit,'(5(x,a))') Grdc_trnm_S(1:Grdc_ntr)
         end if

      else

         if (prout) write (output_unit,1004)
         Grdc_gid= 0 ; Grdc_gjd= 0
         Grdc_gif= 0 ; Grdc_gjf= 0

      end if

      if (printout_L) write (output_unit,1200)

      if (allocated(ac_xp)) deallocate(ac_xp)
      if (allocated(ac_yp)) deallocate(ac_yp)

 1000 format (/,22('#'),' SELF CASCADE CONFIGURATION ',22('#'))
 1001 format ( ' Cascade grid: Tracers to be written for cascade run are: ')
 1004 format (' NO SELF CASCADE DATA will be produced')
 1005 format (' SELF CASCADE OUTPUT WILL BE PRODUCED ################'/&
              ' Cascade grid: Grdc_gid,Grdc_gif=',2I5, ';    Grdc_gjd,Grdc_gjf=',2I5/&
              ' Cascade grid: Grdc_ndt,Grdc_start,Grdc_end=',3i5)
 1006 format ( ' Cascade grid output dimensions and coverage: ',&
                 /1X,' NI=',I5,' FROM x0=',F11.5,' TO xl=',F11.5,' DEGREES' &
                 /1X,' NJ=',I5,' FROM y0=',F11.5,' TO yl=',F11.5,' DEGREES')
 1015 format (' Output grid too large: will be clipped at ',a,' boundary of current grid')
 1200 format (72('#'))
!
!--------------------------------------------------------------------
      return
      end

