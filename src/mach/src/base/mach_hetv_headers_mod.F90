!begin trap head
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

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_hetv_headers_mod.ftn90
! Creation       : S. Menard, H. Landry, Juillet 2008
! Description    : Modules defining explicit interfaces for mach_hetv* subroutines
!
! Extra info     :
!
!============================================================================

module mach_hetv_headers_mod
   interface
!end trap head

subroutine mach_hetv_activity(ne, t, m_h, m_hso4, m_nh4, m_no3, m_so4, &
                              g_h_hso4, g_h2_so4, g_h_no3, g_nh4_no3,  &
                              g_nh42_so4, g_nh4_hso4, g_nh43_hso42)
   integer(kind=4),           intent (in) :: ne
   real(kind=8),              intent (in) :: t           (ne)
   real(kind=8),              intent (in) :: m_h         (ne)
   real(kind=8),              intent (in) :: m_hso4      (ne)
   real(kind=8),              intent (in) :: m_nh4       (ne)
   real(kind=8),              intent (in) :: m_no3       (ne)
   real(kind=8),              intent (in) :: m_so4       (ne)
   real(kind=8),              intent(out) :: g_h_hso4    (ne)
   real(kind=8),              intent(out) :: g_h2_so4    (ne)
   real(kind=8),    optional, intent(out) :: g_h_no3     (ne)
   real(kind=8),    optional, intent(out) :: g_nh4_no3   (ne)
   real(kind=8),    optional, intent(out) :: g_nh42_so4  (ne)
   real(kind=8),    optional, intent(out) :: g_nh4_hso4  (ne)
   real(kind=8),    optional, intent(out) :: g_nh43_hso42(ne)
end subroutine mach_hetv_activity

subroutine mach_hetv_case1(npts, nr, ne1, so4_i, no3_i, nh4_i, hso4_i,    &
                           hno3_i, h_i, nh3_i, ambis_i, lwn_i, t_i, rh_i, &
                           ts_i, ta_i, tn_i, k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne1
   integer(kind=4), intent   (in) :: ncas    (npts)
   real(kind=8),    intent   (in) :: k0      (nr)
   real(kind=8),    intent   (in) :: p1      (nr)
   real(kind=8),    intent   (in) :: p2      (nr)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: ambis_i (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
   real(kind=8),    intent   (in) :: t_i     (npts)
   real(kind=8),    intent   (in) :: rh_i    (npts)
end subroutine mach_hetv_case1

subroutine mach_hetv_case10(npts, nr, ne10, so4_i, no3_i, nh4_i, hso4_i, &
                            hno3_i, h_i, nh3_i, amsul_i, amnit_i,        &
                            mdrh_amnit_amsul_i, drh_amnit_i, lwn_i, t_i, &
                            rh_i, rho_i, ts_i, ta_i, tn_i, k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne10
   integer(kind=4), intent   (in) :: ncas              (npts)
   real(kind=8),    intent   (in) :: k0                (nr)
   real(kind=8),    intent   (in) :: p1                (nr)
   real(kind=8),    intent   (in) :: p2                (nr)
   real(kind=8),    intent(inout) :: so4_i             (npts)
   real(kind=8),    intent(inout) :: no3_i             (npts)
   real(kind=8),    intent(inout) :: nh4_i             (npts)
   real(kind=8),    intent(inout) :: hso4_i            (npts)
   real(kind=8),    intent(inout) :: hno3_i            (npts)
   real(kind=8),    intent(inout) :: h_i               (npts)
   real(kind=8),    intent(inout) :: nh3_i             (npts)
   real(kind=8),    intent(inout) :: amsul_i           (npts)
   real(kind=8),    intent(inout) :: amnit_i           (npts)
   real(kind=8),    intent(inout) :: lwn_i             (npts)
   real(kind=8),    intent   (in) :: t_i               (npts)
   real(kind=8),    intent   (in) :: rh_i              (npts)
   real(kind=8),    intent   (in) :: rho_i             (npts)
   real(kind=8),    intent   (in) :: ta_i              (npts)
   real(kind=8),    intent   (in) :: ts_i              (npts)
   real(kind=8),    intent   (in) :: tn_i              (npts)
   real(kind=8),    intent   (in) :: mdrh_amnit_amsul_i(npts)
   real(kind=8),    intent   (in) :: drh_amnit_i       (npts)
end subroutine mach_hetv_case10

subroutine mach_hetv_case11(npts, nr, ne11, nc, so4_i, no3_i, nh4_i, hso4_i, &
                            hno3_i, h_i, nh3_i, amsul_i, lwn_i, t_i, rh_i,   &
                            rho_i, ts_i, ta_i, tn_i, k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne11
   integer(kind=4), intent   (in) :: nc
   integer(kind=4), intent   (in) :: ncas   (npts)
   real(kind=8),    intent   (in) :: k0     (nr)
   real(kind=8),    intent   (in) :: p1     (nr)
   real(kind=8),    intent   (in) :: p2     (nr)
   real(kind=8),    intent(inout) :: so4_i  (npts)
   real(kind=8),    intent(inout) :: no3_i  (npts)
   real(kind=8),    intent(inout) :: nh4_i  (npts)
   real(kind=8),    intent(inout) :: hso4_i (npts)
   real(kind=8),    intent(inout) :: hno3_i (npts)
   real(kind=8),    intent(inout) :: h_i    (npts)
   real(kind=8),    intent(inout) :: nh3_i  (npts)
   real(kind=8),    intent(inout) :: amsul_i(npts)
   real(kind=8),    intent(inout) :: lwn_i  (npts)
   real(kind=8),    intent   (in) :: t_i    (npts)
   real(kind=8),    intent   (in) :: rh_i   (npts)
   real(kind=8),    intent   (in) :: rho_i  (npts)
   real(kind=8),    intent   (in) :: ta_i   (npts)
   real(kind=8),    intent   (in) :: ts_i   (npts)
   real(kind=8),    intent   (in) :: tn_i   (npts)
end subroutine mach_hetv_case11

subroutine mach_hetv_case12(npts, nr, ne12, so4_i, no3_i, nh4_i, hso4_i, hno3_i,   &
                            h_i, nh3_i, lwn_i, t_i, rh_i, rho_i, ta_i, ts_i, tn_i, &
                            k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne12
   integer(kind=4), intent   (in) :: ncas(npts)
   real(kind=8),    intent   (in) :: k0(nr)
   real(kind=8),    intent   (in) :: p1(nr)
   real(kind=8),    intent   (in) :: p2(nr)
   real(kind=8),    intent(inout) :: so4_i      (npts)
   real(kind=8),    intent(inout) :: no3_i      (npts)
   real(kind=8),    intent(inout) :: nh4_i      (npts)
   real(kind=8),    intent(inout) :: hso4_i     (npts)
   real(kind=8),    intent(inout) :: hno3_i     (npts)
   real(kind=8),    intent(inout) :: h_i        (npts)
   real(kind=8),    intent(inout) :: nh3_i      (npts)
   real(kind=8),    intent(inout) :: lwn_i      (npts)
   real(kind=8),    intent   (in) :: t_i        (npts)
   real(kind=8),    intent   (in) :: rh_i       (npts)
   real(kind=8),    intent   (in) :: rho_i      (npts)
   real(kind=8),    intent   (in) :: ta_i       (npts)
   real(kind=8),    intent   (in) :: ts_i       (npts)
   real(kind=8),    intent   (in) :: tn_i       (npts)
end subroutine mach_hetv_case12

subroutine mach_hetv_case2(npts, nr, ne2, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, lwn_i, t_i, rh_i, ts_i, ta_i, tn_i,     &
                           k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne2
   integer(kind=4), intent   (in) :: ncas (npts)
   real(kind=8),    intent   (in) :: k0   (nr)
   real(kind=8),    intent   (in) :: p1   (nr)
   real(kind=8),    intent   (in) :: p2   (nr)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
   real(kind=8),    intent   (in) :: t_i     (npts)
   real(kind=8),    intent   (in) :: rh_i    (npts)
end subroutine mach_hetv_case2

subroutine mach_hetv_case3(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i, h_i, &
                           nh3_i, ambis_i, leto_i, lwn_i, ts_i, ta_i, tn_i, ncas)
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ncas    (npts)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: ambis_i (npts)
   real(kind=8),    intent(inout) :: leto_i  (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
end subroutine mach_hetv_case3

subroutine mach_hetv_case4(npts, nr, ne4, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, ambis_i, leto_i, mdrh_leto_ambis_i,     &
                           drh_ambis_i, lwn_i, t_i, rh_i, ts_i, ta_i, tn_i,    &
                           k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne4
   integer(kind=4), intent   (in) :: ncas             (npts)
   real(kind=8),    intent   (in) :: k0               (nr)
   real(kind=8),    intent   (in) :: p1               (nr)
   real(kind=8),    intent   (in) :: p2               (nr)
   real(kind=8),    intent(inout) :: so4_i            (npts)
   real(kind=8),    intent(inout) :: no3_i            (npts)
   real(kind=8),    intent(inout) :: nh4_i            (npts)
   real(kind=8),    intent(inout) :: hso4_i           (npts)
   real(kind=8),    intent(inout) :: hno3_i           (npts)
   real(kind=8),    intent(inout) :: h_i              (npts)
   real(kind=8),    intent(inout) :: nh3_i            (npts)
   real(kind=8),    intent(inout) :: ambis_i          (npts)
   real(kind=8),    intent(inout) :: leto_i           (npts)
   real(kind=8),    intent(inout) :: lwn_i            (npts)
   real(kind=8),    intent   (in) :: t_i              (npts)
   real(kind=8),    intent   (in) :: rh_i             (npts)
   real(kind=8),    intent   (in) :: ta_i             (npts)
   real(kind=8),    intent   (in) :: ts_i             (npts)
   real(kind=8),    intent   (in) :: tn_i             (npts)
   real(kind=8),    intent   (in) :: mdrh_leto_ambis_i(npts)
   real(kind=8),    intent   (in) :: drh_ambis_i      (npts)
end subroutine mach_hetv_case4

subroutine mach_hetv_case5(npts, nr, ne5, nc, so4_i, no3_i, nh4_i, hso4_i, &
                           hno3_i, h_i, nh3_i, leto_i, lwn_i, t_i, rh_i,   &
                           ts_i, ta_i, tn_i, k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne5
   integer(kind=4), intent   (in) :: nc
   integer(kind=4), intent   (in) :: ncas    (npts)
   real(kind=8),    intent   (in) :: k0      (nr)
   real(kind=8),    intent   (in) :: p1      (nr)
   real(kind=8),    intent   (in) :: p2      (nr)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: leto_i  (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
   real(kind=8),    intent   (in) :: t_i     (npts)
   real(kind=8),    intent   (in) :: rh_i    (npts)
end subroutine mach_hetv_case5

subroutine mach_hetv_case6(npts, so4_i, no3_i, nh4_i, hso4_i,   &
                           hno3_i, h_i, nh3_i, amsul_i, leto_i, &
                           lwn_i, ts_i, ta_i, tn_i, ncas)
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ncas    (npts)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: amsul_i (npts)
   real(kind=8),    intent(inout) :: leto_i  (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
end subroutine mach_hetv_case6

subroutine mach_hetv_case7(npts, nr, ne7, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, amsul_i, leto_i, mdrh_leto_amsul_i,     &
                           drh_leto_i, lwn_i, t_i, rh_i, ts_i, ta_i, tn_i,     &
                           k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne7
   integer(kind=4), intent   (in) :: ncas             (npts)
   real(kind=8),    intent   (in) :: k0               (nr)
   real(kind=8),    intent   (in) :: p1               (nr)
   real(kind=8),    intent   (in) :: p2               (nr)
   real(kind=8),    intent(inout) :: so4_i            (npts)
   real(kind=8),    intent(inout) :: no3_i            (npts)
   real(kind=8),    intent(inout) :: nh4_i            (npts)
   real(kind=8),    intent(inout) :: hso4_i           (npts)
   real(kind=8),    intent(inout) :: hno3_i           (npts)
   real(kind=8),    intent(inout) :: h_i              (npts)
   real(kind=8),    intent(inout) :: nh3_i            (npts)
   real(kind=8),    intent(inout) :: amsul_i          (npts)
   real(kind=8),    intent(inout) :: leto_i           (npts)
   real(kind=8),    intent(inout) :: lwn_i            (npts)
   real(kind=8),    intent   (in) :: t_i              (npts)
   real(kind=8),    intent   (in) :: rh_i             (npts)
   real(kind=8),    intent   (in) :: ta_i             (npts)
   real(kind=8),    intent   (in) :: ts_i             (npts)
   real(kind=8),    intent   (in) :: tn_i             (npts)
   real(kind=8),    intent   (in) :: mdrh_leto_amsul_i(npts)
   real(kind=8),    intent   (in) :: drh_leto_i       (npts)
end subroutine mach_hetv_case7

subroutine mach_hetv_case8(npts, nr, ne8, nc, so4_i, no3_i, nh4_i, hso4_i, &
                           hno3_i, h_i, nh3_i, amsul_i, lwn_i, t_i, rh_i,  &
                           ts_i, ta_i, tn_i, k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne8
   integer(kind=4), intent   (in) :: nc
   integer(kind=4), intent   (in) :: ncas    (npts)
   real(kind=8),    intent   (in) :: k0      (nr)
   real(kind=8),    intent   (in) :: p1      (nr)
   real(kind=8),    intent   (in) :: p2      (nr)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: amsul_i (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
   real(kind=8),    intent   (in) :: t_i     (npts)
   real(kind=8),    intent   (in) :: rh_i    (npts)
end subroutine mach_hetv_case8

subroutine mach_hetv_case9(npts, ne9, hno3, nh3, amsul, amnit, &
                           kamnit, ta, ts, tn)
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne9
   real(kind=8),    intent(inout) :: hno3  (npts)
   real(kind=8),    intent(inout) :: nh3   (npts)
   real(kind=8),    intent(inout) :: amsul (npts)
   real(kind=8),    intent(inout) :: amnit (npts)
   real(kind=8),    intent   (in) :: ta    (npts)
   real(kind=8),    intent   (in) :: ts    (npts)
   real(kind=8),    intent   (in) :: tn    (npts)
   real(kind=8),    intent   (in) :: kamnit(npts)
end subroutine mach_hetv_case9

subroutine mach_hetv_corrhno3(npts, nr, ne, nc, so4_i, no3_i, nh4_i, hso4_i, hno3_i,  &
                              h_i, lwn_i, t_i, rh_i, rho_i, k0, p1, p2, ncas)
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: nc
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne
   integer(kind=4), intent   (in) :: ncas(npts)
   real(kind=8),    intent   (in) :: k0(nr)
   real(kind=8),    intent   (in) :: p1(nr)
   real(kind=8),    intent   (in) :: p2(nr)
   real(kind=8),    intent(inout) :: so4_i (npts)
   real(kind=8),    intent(inout) :: no3_i (npts)
   real(kind=8),    intent(inout) :: nh4_i (npts)
   real(kind=8),    intent(inout) :: hso4_i(npts)
   real(kind=8),    intent(inout) :: hno3_i(npts)
   real(kind=8),    intent(inout) :: h_i   (npts)
   real(kind=8),    intent(inout) :: lwn_i (npts)
   real(kind=8),    intent   (in) :: t_i   (npts)
   real(kind=8),    intent   (in) :: rh_i  (npts)
   real(kind=8),    intent   (in) :: rho_i (npts)
end subroutine mach_hetv_corrhno3

subroutine mach_hetv_main_2cases(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                                 h_i, nh3_i, lwn_i, t_i, rh_i, rho_i,       &
                                 case_number)
   integer(kind=4), intent   (in) :: npts
   real(kind=8),    intent(inout) :: so4_i      (npts)
   real(kind=8),    intent(inout) :: no3_i      (npts)
   real(kind=8),    intent(inout) :: nh4_i      (npts)
   real(kind=8),    intent(inout) :: hso4_i     (npts)
   real(kind=8),    intent(inout) :: hno3_i     (npts)
   real(kind=8),    intent(inout) :: h_i        (npts)
   real(kind=8),    intent(inout) :: nh3_i      (npts)
   real(kind=8),    intent(inout) :: lwn_i      (npts)
   real(kind=8),    intent   (in) :: t_i        (npts)
   real(kind=8),    intent   (in) :: rh_i       (npts)
   real(kind=8),    intent   (in) :: rho_i      (npts)
   real(kind=4),    intent  (out) :: case_number(npts)
end subroutine mach_hetv_main_2cases

subroutine mach_hetv_main_12cases(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i,  &
                                  h_i, nh3_i, amsul_i, ambis_i, amnit_i,      & 
                                  leto_i, lwn_i, t_i, rh_i, rho_i, case_number)
   integer(kind=4), intent   (in) :: npts
   real(kind=8),    intent(inout) :: so4_i      (npts)
   real(kind=8),    intent(inout) :: no3_i      (npts)
   real(kind=8),    intent(inout) :: nh4_i      (npts)
   real(kind=8),    intent(inout) :: hso4_i     (npts)
   real(kind=8),    intent(inout) :: hno3_i     (npts)
   real(kind=8),    intent(inout) :: h_i        (npts)
   real(kind=8),    intent(inout) :: nh3_i      (npts)
   real(kind=8),    intent  (out) :: amsul_i    (npts)
   real(kind=8),    intent  (out) :: ambis_i    (npts)
   real(kind=8),    intent  (out) :: amnit_i    (npts)
   real(kind=8),    intent  (out) :: leto_i     (npts)
   real(kind=8),    intent  (out) :: lwn_i      (npts)
   real(kind=8),    intent   (in) :: t_i        (npts)
   real(kind=8),    intent   (in) :: rh_i       (npts)
   real(kind=8),    intent   (in) :: rho_i      (npts)
   real(kind=4),    intent  (out) :: case_number(npts)
end subroutine mach_hetv_main_12cases

subroutine mach_hetv_hetchem(gascon, aerocon, tempk, zpres, aeronum, &
                             rhrow, rhorow, ibulk, jlat, kount, pni, pnk)
   use mach_cam_utils_mod,    only: isize, icom, maxns
   integer(kind=4), intent   (in) :: ibulk
   integer(kind=4), intent   (in) :: jlat
   integer(kind=4), intent   (in) :: kount
   integer(kind=4), intent   (in) :: pni, pnk
   real(kind=4),    intent(inout) :: gascon (pni, pnk, maxns)
   real(kind=4),    intent(inout) :: aerocon(pni, pnk, icom, isize)
   real(kind=4),    intent   (in) :: tempk  (pni, pnk)
   real(kind=4),    intent   (in) :: zpres  (pni, pnk)
   real(kind=4),    intent   (in) :: aeronum(pni, pnk, isize)
   real(kind=4),    intent   (in) :: rhrow  (pni, pnk)
   real(kind=4),    intent   (in) :: rhorow (pni, pnk)
end subroutine mach_hetv_hetchem

subroutine mach_hetv_hetvcall(npts, ghetf, hetf, t_i, rh_i, rho_i, jlat, kount)
   use mach_hetv_mod,         only: maxhet, maxghet
   use mach_cam_utils_mod,    only: isize
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: jlat
   integer(kind=4), intent   (in) :: kount
   real(kind=8),    intent   (in) :: t_i   (npts)
   real(kind=8),    intent   (in) :: rh_i  (npts)
   real(kind=8),    intent   (in) :: rho_i (npts)
   real(kind=8),    intent(inout) :: ghetf (npts, maxghet)
   real(kind=8),    intent(inout) :: hetf  (npts, maxhet, isize)
end subroutine mach_hetv_hetvcall

subroutine mach_hetv_poly3v (a1, a2, a3, root, islv, ne)
   integer(kind=4), intent (in) :: ne
   real(kind=8),    intent (in) :: a1  (ne)
   real(kind=8),    intent (in) :: a2  (ne)
   real(kind=8),    intent (in) :: a3  (ne)
   real(kind=8),    intent(out) :: root(ne)
   integer(kind=4), intent(out) :: islv(ne)
end subroutine mach_hetv_poly3v

subroutine mach_hetv_rebin(het, dhet_chem, tempk, zpres, binnum, kount, jlat, &
                           pni, pnk)
   use mach_hetv_mod,      only: maxhet
   use mach_cam_utils_mod, only: isize
    integer(kind=4), intent   (in) :: kount
    integer(kind=4), intent   (in) :: jlat
    integer(kind=4), intent   (in) :: pni, pnk
    real(kind=8),    intent(inout) :: het      (pni, pnk, maxhet, isize)
    real(kind=8),    intent(inout) :: dhet_chem(pni, pnk, maxhet)
    real(kind=4),    intent   (in) :: tempk    (pni, pnk)
    real(kind=4),    intent   (in) :: zpres    (pni, pnk)
    real(kind=4),    intent   (in) :: binnum   (pni, pnk, isize)
end subroutine mach_hetv_rebin

subroutine mach_hetv_water(ne, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)
   
   integer(kind=4), intent (in) :: ne
   real(kind=8),    intent (in) :: so4 (ne)
   real(kind=8),    intent (in) :: h   (ne)
   real(kind=8),    intent (in) :: no3 (ne)
   real(kind=8),    intent (in) :: nh4 (ne)
   real(kind=8),    intent (in) :: hso4(ne)
   real(kind=8),    intent (in) :: awu (ne)
   real(kind=8),    intent (in) :: law (ne)
   real(kind=8),    intent (in) :: lwo (ne)
   real(kind=8),    intent(out) :: lwn (ne)
end subroutine mach_hetv_water

end interface
end module mach_hetv_headers_mod
