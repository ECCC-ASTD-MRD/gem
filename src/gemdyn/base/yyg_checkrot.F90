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

!**s/r yyg_checkrot

      integer function yyg_checkrot ()
      use HORgrid_options
      use lun
      implicit none
#include <arch_specific.hf>

!author
!     V. Lee/A. Qaddouri - April 2011
!
!revision
! v4_40 - Qaddouri/Lee     - initial version
!
!-------------------------------------------------------------------
!
      yyg_checkrot = -1

      if (Grd_xlat1 < 0.0 .and. Grd_xlat2 < 0.0 .and. &
                                 Grd_xlon2 > Grd_xlon1) then
          print *,'ERROR: Grd_xlat1,Grdxlat2 < 0.0 and Grd_xlon2 > Grd_xlon1'
          if (Lun_out > 0) then
             write(Lun_out,1001)  Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2
             write(Lun_out,8000)
          end if
          return
      else if (Grd_xlat1 >= 0.0 .and. Grd_xlat2 >= 0.0 .and. &
                                 Grd_xlon1 > Grd_xlon2) then
          print *,'ERROR: Grd_xlat1,Grd_xlat2 >= 0.0 and Grd_xlon1 > Grd_xlon2'
          if (Lun_out > 0) then
             write(Lun_out,1001)  Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2
             write(Lun_out,8000)
          end if
          return
      end if

      yyg_checkrot = 1
!
 1001 format(/,' WRONG YIN GRID CONFIGURATION --- ABORT ---'/, &
               ' Grd_xlat1,Grd_xlon1,Grd_xlat2,Grd_xlon2:'/4f10.3/)
 8000 format (/,'========= ERROR IN S/R yyg_checkrot ============='/)
!
!-------------------------------------------------------------------
!
      return
      end
