!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

module chm_phyvar_mod
!
   use phybus, only: me_moins, p0_plus, sigm, sigt, o3lplus,      &
                     tplus, huplus, qcplus, uplus, vplus, wplus,  &
                     dlat, dlon, dxdy, flusolis, glsea, cang, h,  &
                     qdiag, tdiag, tsurf, udiag, vdiag, ilmo, ue, &
                     rainrate, gzmom, gztherm, kt, ftot, f12,     &
                     rnflx, snoflx, fevp, qrkfc, kfcrf, kfcsf,    &
                     kfmrf, kfmsf, mg, alvis, lwc, cldrad, o3s,   &
                     shear2, rig, ncplus, o3ce

   save

   integer(kind=4) psn
   integer(kind=4) snodp
   integer(kind=4) urban
   integer(kind=4) vegf
   integer(kind=4) wsoil

end module chm_phyvar_mod
