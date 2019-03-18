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

module cstv
   implicit none
   public
   save

!
! Contains some constant parameters during the integration.
!
!________________________________________________________________________________________________________
!                                                                                                        |
! SOME CONSTANT PARAMETERS DURING THE INTEGRATION                                                        |
!_______________________________________________________________________________________________________ |
!                    |                                                                                   |
!     NAME           |                         DESCRIPTION                                               |
!--------------------|---------------------------------------------------------------------------------- |
!                    |                                                                                   |
! Cstv_dtA_8         | timestep used by the model:traditional advection =>Cstv_dtA_8 = Cstv_dt_8 * 0.5d0 |
!                    |                             consistent advection =>Cstv_dtA_8 = Cstv_tau_8        |
! Cstv_dtD_8         | Cstv_dtD_8 = Cstv_dt_8 - Cstv_dtA_8                                               |
! Cstv_Tau_8         | effective timestep: tau=dt*bA                                                     |
! Cstv_Tau_m_8       | effective timestep for mom. eqns: tau_m=dt*bA_m                                   |
! Cstv_Tau_nh_8      | effective timestep for nonhydro. eqns: tau_nh=dt*bA_nh                            |
! Cstv_invT_8        | inverse of tau                                                                    |
! Cstv_invT_m_8      | inverse of tau_m                                                                  |
! Cstv_invT_nh_8     | inverse of tau_nh                                                                 |
! Cstv_Beta_8        | ratio Beta=(1-bA)/bA                                                              |
! Cstv_Beta_m_8      | ratio Beta=(1-bA_m)/bA_m                                                          |
! Cstv_Beta_nh_8     | ratio Beta=(1-bA_nh)/bA_nh                                                        |
! Cstv_Sstar_8       | S*                                                                                |
! Cstv_pref_8        | a reference pressure : 100000 PASCALS                                             |
! Cstv_ptop_8        | pressure at the top                                                               |
! Cstv_Ztop_8        | ln(ptop)                                                                          |
! Cstv_Zsrf_8        | ln(pSref)                                                                         |
! Cstv_hco0_8        | Helmholtz constant                                                                |
! Cstv_hco1_8        | Helmholtz constant                                                                |
! Cstv_hco2_8        | Helmholtz constant                                                                |
! Cstv_bar0_8        | value for barotropic case, ZERO otherwise                                         |
! Cstv_bar1_8        | ZERO  for barotropic case, one  otherwise                                         |
!--------------------|------------------------------------------------------------------------------------
!
      real*8 :: Cstv_dtA_8   ,Cstv_dtD_8    ,Cstv_dtzA_8   ,&
                Cstv_Tau_8   ,Cstv_invT_8   ,Cstv_Beta_8   ,&
                Cstv_tau_m_8 ,Cstv_invT_m_8 ,Cstv_Beta_m_8 ,&
                Cstv_tau_nh_8,Cstv_invT_nh_8,Cstv_Beta_nh_8,&
                Cstv_Sstar_8 ,Cstv_pref_8   ,Cstv_dtzD_8   ,&
                Cstv_ptop_8  ,Cstv_Ztop_8   ,Cstv_Zsrf_8   ,&
                Cstv_hco0_8  ,Cstv_hco1_8   ,Cstv_hco2_8   ,&
                Cstv_bar0_8  ,Cstv_bar1_8   ,Cstv_dt_8

end module cstv
