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
! Fichier/File   : mach_hetv_activity.ftn90
! Creation       : P. Makar, V. Bouchet, S. Gravel, B. Pabla, S. Menard
! Description    : Calculates activity coefficients for gridpoints in the
!                  range 1 to ne.  The system includes the following
!                  includes the following aqueous phase reactions:
!
!                  HNO3(g,eq)    <=> H+ + NO3-
!                  HSO4-         <=> H+ + SO4=
!                  NH4HSO4(s)    <=> NH4+ + HSO4-
!                  (NH4)2SO4(s)  <=> 2 NH4+ + SO4=
!                  (NH4)3HSO4(s) <=> 3 NH4+ + HSO4- + SO4=
!                  NH4NO3(aq)    <=> NH4+ + NO3-
!                  NH3(aq) + H+  <=> NH4+
!                  H2O(aq)       <=> H+ + OH-
!
!                  The ions are H+, NO3-, SO4=, HSO4-, OH-, and NH4+
!
! Extra info     : Reference:  Kim, Y.P., J.H. Seinfeld and P. Saxeena, "Atmospheric
!                  Gas-Aerosol Equilibrium I.  Thermodynamic Model",  Aerosol Sci. and
!                  Technology, 19: 157-181, 1993.   Pilinis, C. and J.H. Seinfeld,
!                  "Continued development of a general equilibrium model for inorganic
!                  multicomponent atmospheric aerosols", Atm. Env. 21, 2453-2466, 1987.
!                  Equation numbers in brackets in the code refer to equation numbers
!                  in the 1987 paper.
!
!                  Temperature correction is from Meissner, H.P. and Peppas, N.A.,
!                  "Activity coefficients - aqueous solutions of polybasic acids and their
!                  salts", AIChE Journal, Vol 19, No. 4, pp 806-809, equations 1 and 2.
!                  Temperature correction implemented only for the KM method;
!                  variable 'method' is set to = "KM".
!
!                  Cap on activity coefficients imposed to prevent underflow/overflow
!                  errors in calculating net equilibrium
!
!                  Pitzer method binary activity coefficient values:
!                  HNO3, H2SO4, H-HSO4, NH4NO3, (NH4)2SO4, NH4HSO4 values from
!                  Kim, Y.P., J.H. Seinfeld, P. Saxena, Atmospheric Gas-Aerosol
!                  Equilibrium I.  Thermodynamic Model., Aerosol Sci. Tech.
!                  vol 19:157-181, (1993).
!                  Note that the values for gamma_h_OH and gamma_nh4_oh are
!                  both assumed to be 1.0, due to lack of data to the contrary.
!                  This should not be a bad approximation, since these are expected
!                  to be minor components of the net system.
!
!                  Kusik and Meissner method binary activity coefficients
!                  values from Kim, Y.P., J.H. Seinfeld, P. Saxena, Atmospheric Gas-Aerosol
!                  Equilibrium I.  Thermodynamic Model., Aerosol Sci. Tech.
!                  vol 19:157-181, (1993).
!                  Note that the values for gamma_h_OH and gamma_nh4_oh are
!                  both assumed to be 1.0, due to lack of data to the contrary.
!                  This should not be a bad approximation, since these are expectedv
!                  to be minor components of the net system.


!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_activity(ne, t, m_h, m_hso4, m_nh4, m_no3, m_so4, &
                              g_h_hso4, g_h2_so4, g_h_no3, g_nh4_no3,  &
                              g_nh42_so4, g_nh4_hso4, g_nh43_hso42)
!!if_off
   use mach_hetv_mod, only: method
   implicit none
!!if_on
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
!!if_off
!
! local variables
!
   real(kind=8), parameter     :: ag = 0.511d0  !(gk/mole)**0.5
   real(kind=8), parameter     :: gmin = 1.0d-05, gmax = 1.0d+05 !cap on activity coefficients
   integer                     :: i
   real(kind=8), parameter     :: ionmin = 1.0d-30  ! minimum allowable ionic concentrati
   real(kind=8)                :: fac, tc
   real(kind=8)                :: ttt, fgamma
   real(kind=8), dimension(ne) :: g0_h_no3, g0_h2_so4, g0_h_hso4, g0_nh4_no3, g0_nh42_so4, g0_nh4_hso4
   real(kind=8), dimension(ne) :: f_nh4, f_h, f_no3, f_hso4, f_so4
!
!  pitzer method binary activity coefficient parameters:
   real(kind=8), parameter :: beta0_h_no3 = 0.1119d0, beta0_h2_so4 = -0.09330d0, beta0_h_hso4 = 0.25713d0
   real(kind=8), parameter :: beta0_nh4_no3 = -0.0154d0, beta0_nh42_so4 = 0.04763d0, beta0_nh4_hso4 =  0.04494d0
   real(kind=8), parameter :: beta1_h_no3 = 0.3206d0, beta1_h2_so4 = 0.32381d0, beta1_h_hso4 = 0.35308d0
   real(kind=8), parameter :: beta1_nh4_no3 = 0.1120d0, beta1_nh42_so4 = 0.44459d0, beta1_nh4_hso4 = 0.23594d0
   real(kind=8)            :: b_h_no3, b_h2_so4, b_h_hso4
   real(kind=8)            :: b_nh4_no3, b_nh42_so4, b_nh4_hso4
   real(kind=8), parameter :: c_h_no3 = 0.0010d0 * 1.5d0, c_h2_so4 = 0.021162d0 * 1.5d0, c_h_hso4 = -0.002830d0 * 1.5d0
   real(kind=8), parameter :: c_nh4_no3 = -0.00003d0 * 1.5d0, c_nh42_so4 = -0.001311d0 * 1.5d0, c_nh4_hso4 = -0.002920d0 * 1.5d0
   real(kind=8)            :: m_hno3, m_h2so4, m_hhso4, m_nh4no3, m_nh42so4
   real(kind=8)            :: m_nh4hso4
!  kusik and meissner binary coefficient parameters:
   real(kind=8), parameter :: q_nh42so4 = -0.25d0, q_nh4no3 = -1.15d0
   real(kind=8), parameter :: q_h2so4 = -0.1d0, q_h_hso4 = 8.0d0, q_hno3 = 2.6d0
   real(kind=8)            :: cfac, tfac
   real(kind=8), parameter, dimension(63) :: pten = &
                              (/1.0d-31, 1.0d-30, 1.0d-29, 1.0d-28, 1.0d-27,  &
                                1.0d-26, 1.0d-25, 1.0d-24, 1.0d-23, 1.0d-22,  &
                                1.0d-21, 1.0d-20, 1.0d-19, 1.0d-18, 1.0d-17,  &
                                1.0d-16, 1.0d-15, 1.0d-14, 1.0d-13, 1.0d-12,  &
                                1.0d-11, 1.0d-10, 1.0d-09, 1.0d-08, 1.0d-07,  &
                                1.0d-06, 1.0d-05, 1.0d-04, 1.0d-03, 1.0d-02,  &
                                1.0d-01, 1.0d+00, 1.0d+01, 1.0d+02, 1.0d+03,  &
                                1.0d+04, 1.0d+05, 1.0d+06, 1.0d+07, 1.0d+08,  &
                                1.0d+09, 1.0d+10, 1.0d+11, 1.0d+12, 1.0d+13,  &
                                1.0d+14, 1.0d+15, 1.0d+16, 1.0d+17, 1.0d+18,  &
                                1.0d+19, 1.0d+20, 1.0d+21, 1.0d+22, 1.0d+23,  &
                                1.0d+24, 1.0d+25, 1.0d+26, 1.0d+27, 1.0d+28,  &
                                1.0d+29, 1.0d+30, 1.0d+31/)
                                 !
   real(kind=8), parameter, dimension(14) :: pcoef1 = & ! Taylor series power coefficients
                        (/1.147074559772972471d-11, 1.60590438368216146d-10, &
                           2.08767569878680990d-09, 2.50521083854417188d-08, &
                           2.75573192239858907d-07, 2.75573192239858907d-06, &
                           2.48015873015873016d-05, 1.98412698412698413d-04, &
                           1.38888888888888889d-03, 8.33333333333333333d-03, &
                           4.16666666666666667d-02, 1.66666666666666667d-01, &
                           0.50000000000000000d0,   1.00000000000000000d0 /)
   real(kind=8), parameter, dimension(14) :: pcoef2 = & ! Taylor series power coefficients
                        (/1.35086294762237060d-06,  8.21341253543939690d-06, &
                          4.63715166425722500d-05,  2.41666725544247160d-04, &
                          1.15449977899843590d-03,  5.01392883377544360d-03, &
                          1.95976946264785410d-02,  6.80893650744371080d-02, &
                          2.06995848696868180d-01,  5.39382929195581620d-01, &
                          1.17125514891226730d0,    2.03467859229347690d0,   &
                          2.65094905523919970d0,    2.30258509299404590d0/)
   real(kind=8)                :: fion, lnion, p_hno3, p_h2so4_b, p_nh4no3, p_h_hso4, p_nh42so4,  &
                                  x_h_no3, x_h2_so4, x_h_hso4, x_nh4_no3, x_nh42_so4, x_nh4_hso4, &
                                  y_no3_h, y_so4_h2, y_hso4_h, y_oh_h, y_no3_nh4, y_so4_nh42,     &
                                  y_hso4_nh4, y_oh_nh4
   real(kind=8)                :: g_h_no3_p, g_h_hso4_p, g_nh42_so4_p, g_nh4_hso4_p,   &
                                  g_nh4_no3_p, g_h2_so4_p, g_h_no3_r, g_h_hso4_r,      &
                                  g_nh42_so4_r, g_nh4_hso4_r, g_nh4_no3_r, g_h2_so4_r, &
                                  g_h_no3_v, g_h_hso4_v, g_nh42_so4_v, g_nh4_hso4_v,   &
                                  g_nh4_no3_v, g_h2_so4_v
   real(kind=8), dimension(ne) :: ion, sqion, m_oh
!
!----------------------------------------------------------------------
!Name      |Description                           |T  |Dimensions |In/|
!          |                                      |Y  |           |Out|
!          |                                      |P  |           |Loc|
!          |                                      |E  |           |cal|
!----------------------------------------------------------------------
! ag       | Debye-Huckel constant for activity   | R |           | L |
!          |coefficient; 0.511 (kg/mole)**.5      |   |           |   |
! b_h_hso4 | Pitzer method "Bgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation. 1=H, 2=HSO4. See         |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! b_h_no3  | Pitzer method "Bgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation.  1=H, 2=NO3.  See        |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! b_h2_so4 | Pitzer method "Bgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation.  1=H2, 2=SO4.  See       |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! b_nh4_hso| Pitzer method "Bgamma12" parameter   | R |           | L |
!4         |for binary activity coefficient       |   |           |   |
!          |calculation.  1=NH4, 2=HSO4.  See     |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! b_nh4_no3| Pitzer method "Bgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation.  1=NH4, 2=NO3.  See      |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! b_nh42_so| Pitzer method "Bgamma12" parameter   | R |           | L |
!4         |for binary activity coefficient       |   |           |   |
!          |calculation.  1=(NH4)2, 2=SO4.  See   |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! beta0_h_h| Pitzer method beta(0) H_HSO4         | R |           | L |
!so4       |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta0_h_n| Pitzer method beta(0) HNO3 parameter | R |           | L |
!o3        |for binary activity coefficient       |   |           |   |
!          |calculation.  See Kim et al (1993)    |   |           |   |
!          |Table 8.                              |   |           |   |
! beta0_h2_| Pitzer method beta(0) H2SO4          | R |           | L |
!so4       |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta0_nh4| Pitzer method beta(0) NH4HSO4        | R |           | L |
!_hso4     |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta0_nh4| Pitzer method beta(0) NH4NO3         | R |           | L |
!_no3      |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta0_nh4| Pitzer method beta(0) (NH4)2SO4      | R |           | L |
!2_so4     |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta1_h_h| Pitzer method beta(1) H_HSO4         | R |           | L |
!so4       |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta1_h_n| Pitzer method beta(1) HNO3 parameter | R |           | L |
!o3        |for binary activity coefficient       |   |           |   |
!          |calculation.  See Kim et al (1993)    |   |           |   |
!          |Table 8.                              |   |           |   |
! beta1_h2_| Pitzer method beta(1) H2SO4          | R |           | L |
!so4       |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta1_nh4| Pitzer method beta(1) NH4HSO4        | R |           | L |
!_hso4     |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta1_nh4| Pitzer method beta(1) NH4NO3         | R |           | L |
!_no3      |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! beta1_nh4| Pitzer method beta(1) (NH4)2SO4      | R |           | L |
!2_so4     |parameter for binary activity         |   |           |   |
!          |coefficient calculation.  See Kim et  |   |           |   |
!          |al (1993) Table 8.                    |   |           |   |
! c_h_hso4 | Pitzer method "Cgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation.  1=H, 2=HSO4.  See       |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! c_h_no3  | Pitzer method "Cgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation.  1=H, 2=NO3.  See        |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! c_h2_so4 | Pitzer method "Cgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation.  1=H2, 2=SO4.  See       |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! c_nh4_hso| Pitzer method "Cgamma12" parameter   | R |           | L |
!4         |for binary activity coefficient       |   |           |   |
!          |calculation.  1=NH4, 2=HSO4.  See     |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! c_nh4_no3| Pitzer method "Cgamma12" parameter   | R |           | L |
!          |for binary activity coefficient       |   |           |   |
!          |calculation.  1=NH4, 2=NO3.  See      |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! c_nh42_so| Pitzer method "Cgamma12" parameter   | R |           | L |
!4         |for binary activity coefficient       |   |           |   |
!          |calculation.  1=(NH4)2, 2=SO4.  See   |   |           |   |
!          |formulae following Table 8, Kim et    |   |           |   |
!          |al.(1993).                            |   |           |   |
! f_h      | Bromley multicomponent activity      | R | ne        | L |
!          |coefficient F1 (or F2) parameter (eqn |   |           |   |
!          |12, Kim et al., 1993). 1= H           |   |           |   |
! f_hso4   | Bromley multicomponent activity      | R | ne        | L |
!          |coefficient F1 (or F2) parameter (eqn |   |           |   |
!          |12, Kim et al., 1993). 2=HSO4         |   |           |   |
! f_nh4    | Bromley multicomponent activity      | R | ne        | L |
!          |coefficient F1 (or F2) parameter (eqn |   |           |   |
!          |12, Kim et al., 1993). 1= NH4         |   |           |   |
! f_no3    | Bromley multicomponent activity      | R | ne        | L |
!          |coefficient F1 (or F2) parameter (eqn |   |           |   |
!          |12, Kim et al., 1993). 2=NO3          |   |           |   |
! f_oh     | Bromley multicomponent activity      | R | ne        | L |
!          |coefficient F1 (or F2) parameter (eqn |   |           |   |
!          |12, Kim et al., 1993). 2=OH           |   |           |   |
! f_so4    | Bromley multicomponent activity      | R | ne        | L |
!          |coefficient F1 (or F2) parameter (eqn |   |           |   |
!          |12, Kim et al., 1993). 2=SO4          |   |           |   |
! fac      | Pitzer method binary activity        | R |           | L |
!          |coefficient Bgamma12 factor; used in  |   |           |   |
!          |Bgamma12 calculation:  Table 8, Kim   |   |           |   |
!          |et al. (1993).                        |   |           |   |
! fgamma   | Pitzer binary activity coefficient   | R | ne        | L |
!          |f(gamma) parameter, Table 8, Kim et   |   |           |   |
!          |al (1993).                            |   |           |   |
! g_h_hso4 | Activity coefficient of H-HSO4.      | R | ne        | O |
! g_h_no3  | Activity coefficient of HNO3.        | R | ne        | O |
! g_h_oh   | Activity coefficient of water.       | R | ne        | O |
! g_h2_so4 | Activity coefficient of H2SO4.       | R | ne        | O |
! g_nh4_hso| Activity coefficient of NH4HSO4.     | R | ne        | O |
!4         |                                      |   |           |   |
! g_nh4_no3| Activity coefficient of NH4NO3.      | R | ne        | O |
! g_nh4_oh | Activity coefficient of NH4OH.       | R |           | O |
! g_nh42_so| Activity coefficient of (NH4)2SO4.   | R | ne        | O |
!4         |                                      |   |           |   |
! g_nh43_hs| Activity coefficient of              | R | ne        | O |
!o4        |(NH4)3H(SO4)2.                        |   |           |   |
! g0_h_hso4| Binary activity coefficient for      | R | ne        | L |
!          |H-HSO4.                               |   |           |   |
! g0_h_no3 | Binary activity coefficient for      | R | ne        | L |
!          |HNO3.                                 |   |           |   |
! g0_h2_so4| Binary activity coefficient for      | R | ne        | L |
!          |H2SO4.                                |   |           |   |
! g0_nh4_hs| Binary activity coefficient for      | R | ne        | L |
!o4        |NH4HSO4.                              |   |           |   |
! g0_nh4_no| Binary activity coefficient for      | R | ne        | L |
!3         |NH4NO3.                               |   |           |   |
! g0_nh42_s| Binary activity coefficient for      | R | ne        | L |
!o4        |(NH4)2SO4.                            |   |           |   |
! i        | Loop index.                          | I |           | L |
! ion      | Ionic strength of solution; 0.5*sum  | R | ne        | L |
!          |of {molalities*square of number of    |   |           |   |
!          |unit charges) for each ion.           |   |           |   |
! m_h      | Molality (moles/kg H2O) of H+.       | R | ne        | I |
! m_h2so4  | Molality of H2SO4 required to match  | R |           | L |
!          |the ionic strength of the current     |   |           |   |
!          |mixture.                              |   |           |   |
! m_hhso4  | Molality of H_HSO4 required to match | R |           | L |
!          |the ionic strength of the current     |   |           |   |
!          |mixture.                              |   |           |   |
! m_hno3   | Molality of HNO3 required to match   | R |           | L |
!          |the ionic strength of the current     |   |           |   |
!          |mixture.                              |   |           |   |
! m_hso4   | Molality (moles/kg H2O) of HSO4-.    | R | ne        | I |
! m_nh4    | Molality (moles/kg H2O) of NH4+.     | R | ne        | I |
! m_nh42so4| Molality of (NH4)2SO4 required to    | R |           | L |
!          |match the ionic strength of the       |   |           |   |
!          |current mixture.                      |   |           |   |
! m_nh4hso4| Molality of NH4HSO4 required to      | R |           | L |
!          |match the ionic strength of the       |   |           |   |
!          |current mixture.                      |   |           |   |
! m_nh4no3 | Molality of NH4NO3 required to match | R |           | L |
!          |the ionic strength of the current     |   |           |   |
!          |mixture.                              |   |           |   |
! m_no3    | Molality (moles/kg H2O) of NO3-.     | R | ne        | I |
! m_oh     | Molality of (moles/kg H2O) OH-.      | R | ne        | I |
! m_so4    | Molality (moles/kg H2O) of SO4=.     | R | ne        | I |
! method   | Method to use for binary activity    | C | C*20      | I |
!          |coefficient calculation.  Current     |   |           |   |
!          |choices are "Pitzer" and "KM" (for    |   |           |   |
!          |Kusik and Meissner).                  |   |           |   |
! ne       | Number of gridpoints for which       | I |          | I |
!          |activity coefficients are desired.    |   |           |   |
! q_h_hso4 | Kusik and Meissner binary activity   | R |           | L |
!          |coefficient q parameter for H_HSO4.   |   |           |   |
!          |Table 10, Kim et al.(1993).           |   |           |   |
! q_h2so4  | Kusik and Mei
!
!  The conditions for which this case is called are as follows:
!       TA/TS < 1, rh < drh_ambis
!  The reactions representing this case are as follows(bracketed numbers
!  indicate that the solution is done in stages, as numbered):
!  (1)
!     HSO4 <=> H + SO4 , kHSO4
!     NH4HSO4 = TA - NH4
!     TS = SO4 + HSO4 + NH4HSO4
!     H + NH4 = 2 SO4 + HSO4
!  (2)
!     NH4HSO4 <=> NH4 + HSO4, kNH4HSO4
!  The system of equations to be solved in this case are as follows:
!   Equations
!  (1)
!     kHSO4 (HSO4)- (H)(SO4) = 0
!     TA - NH4 - NH4HSO4 = 0,
!     TS - SO4 - HSO4 - NH4HSO4 = 0,
!     H + NH4 - 2 SO4 - HSO4 = 0.
!  (2)
!     kNH4HSO4 - (NH4)(HSO4) = 0
!  The solution to the system of equations:
!
!  (1)
!     Let
!     b = kHSO4 - TA + TS,
!     c = kHSO4 (TA - NH4 - TS)
!     SO4 = 0.5*(-b + sqrt(b**2 -4c)),
!     NH4HSO4 = TA - NH4,
!     H = SO4 + TS - TA,
!     HSO4 = TS - SO4 - TA + NH4
!  (2)
!     kNH4HSO4 - (NH4)(HSO4) = 0
!  Variable, Range
!     NH4, 0 to TA
!  -----------------------------------------------------------------
!ssner binary activity   | R |           | L |
!          |coefficient q parameter for H2SO4.    |   |           |   |
!          |Table 10, Kim et al.(1993).           |   |           |   |
! q_hno3   | Kusik and Meissner binary activity   | R |           | L |
!          |coefficient q parameter for HNO3.     |   |           |   |
!          |Table 10, Kim et al.(1993).           |   |           |   |
! q_nh42so4| Kusik and Meissner binary activity   | R |           | L |
!          |coefficient q parameter for           |   |           |   |
!          |(NH4)2SO4.  Table 10, Kim et          |   |           |   |
!          |al.(1993).                            |   |           |   |
! q_nh4no3 | Kusik and Meissner binary activity   | R |           | L |
!          |coefficient q parameter for NH4NO3.   |   |           |   |
!          |Table 10, Kim et al.(1993).           |   |           |   |
! sqion    | Square root of ionic strength (see   | R | ne         | L |
!          |"ion").                               |   |           |   |
! x_h_hso4 | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H, 2=HSO4        |   |           |   |
! x_h_no3  | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H, 2=NO3         |   |           |   |
! x_h_oh   | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H, 2=OH          |   |           |   |
! x_h2_so4 | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H2, 2=SO4        |   |           |   |
! x_nh4_hso| Bromley multicomponent activity      | R | ne         | L |
!4         |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=NH4, 2=HSO4      |   |           |   |
! x_nh4_no3| Bromley multicomponent activity      | R | ne         | L |
!          |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=NH4, 2=NO3       |   |           |   |
! x_nh4_oh | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=NH4, 2=OH        |   |           |   |
! x_nh42_so| Bromley multicomponent activity      | R | ne         | L |
!4         |coefficient X12 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=(NH4)2, 2=SO4    |   |           |   |
! y_hso4_h | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H, 2=HSO4        |   |           |   |
! y_hso4_nh| Bromley multicomponent activity      | R | ne         | L |
!4         |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=NH4, 2=HSO4      |   |           |   |
! y_no3_h  | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H, 2=NO3         |   |           |   |
! y_no3_nh4| Bromley multicomponent activity      | R | ne         | L |
!          |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=NH4, 2=NO3       |   |           |   |
! y_oh_h   | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H, 2=OH          |   |           |   |
! y_oh_nh4 | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=NH4, 2=OH        |   |           |   |
! y_so4_h2 | Bromley multicomponent activity      | R | ne         | L |
!          |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=H2, 2=SO4        |   |           |   |
! y_so4_nh4| Bromley multicomponent activity      | R | ne         | L |
!2         |coefficient Y21 parameter (eqn 12,    |   |           |   |
!          |Kim et al., 1993). 1=(NH4)2, 2=SO4    |   |           |   |
!----------------------------------------------------------------------
!
!  First, calculate the ionic strength, and fgamma:
!  Note:  the value of m_oh input here is actually oh*gamma(oh)*gamma(h);
!  the actual value of [OH] can't be determined from the system of
!  equations used for SO4-NO3-NH4.  It will be assumed to be zero
!  here, but, if other OH containing salts are added to the system
!  at some later date, these could be used to determine the OH
!  concentration.
!IBM* UNROLL(8)
   do i = 1, ne
      m_oh(i) = ionmin
!
!  Note:  minimum values of H, NO3, SO4, HSO4 and NH4
!  all assumed to be ionmin here.  The routine is sometimes
!  called for dry aerosols; this adjustment is necessary to
!  prevent the ionic concentration from going to infinity
!  in these cases.

      ion(i) = 0.5d0 * ( max(m_h(i), ionmin) + max(m_no3(i), ionmin) + max(m_so4(i) * 4.0d0, ionmin * 4.0d0) +  &
                       max(m_hso4(i), ionmin) + max(m_oh(i), ionmin) + max(m_nh4(i), ionmin))
      sqion(i) = sqrt(ion(i))
   end do

   if (method == 'Pitzer') then

!  Pitzer's method binary activity coefficients:
!  Calculate B values for all possible ion combinations:

      do i = 1, ne
!  (eq 17)
         fgamma = -0.392d0 * (sqion(i) / (1.0d0 + 1.2d0 * sqion(i)) +        &
                     1.666666666666667D0 * log(1.0d0 + 1.2d0 * sqion(i)))
         fac = (1.0d0 - exp(-2.0d0 * sqion(i))*(1.0d0 + 2.0d0 * ( sqion(i)-ion(i)))) &
                  / (4.0d0 * ion(i))
!  (eq 18)
         b_h_no3 = 2.0d0 * (beta0_h_no3 +  beta1_h_no3 * fac)
         b_h2_so4 = 2.0d0 * (beta0_h2_so4 +  beta1_h2_so4 * fac)
         b_h_hso4 = 2.0d0 * (beta0_h_hso4 +  beta1_h_hso4 * fac)
         b_nh4_no3 = 2.0d0 * (beta0_nh4_no3 +  beta1_nh4_no3 * fac)
         b_nh42_so4 = 2.0d0 * (beta0_nh42_so4 +  beta1_nh42_so4 * fac)
         b_nh4_hso4 = 2.0d0 * (beta0_nh4_hso4 +  beta1_nh4_hso4 * fac)
!  Molalities of electrolytes for given ionic strength:
!  egs.  2 H+ + SO4=:  molality of H2SO4 same as SO4
!  I = 0.5* { [H+] + [SO4=]*4 }; [H+] = 2 [SO4=]
!  so I = 0.5*6*[SO4=], molality of H2SO4 = 1/3 I.
         m_hno3 = ion(i)
         m_h2so4 = 0.33333333333333D0 * ion(i)
         m_hhso4 = ion(i)
         m_nh4no3 = ion(i)
         m_nh42so4 = 0.33333333333333D0 * ion(i)
         m_nh4hso4 = ion(i)
!  (eq 16)  Pitzer's method binary activity coefficients
         g0_h_no3(i) = exp(fgamma + m_hno3 * (b_h_no3 + m_hno3 *  c_h_no3))
         g0_h2_so4(i) = exp(2.0d0 * fgamma + m_h2so4 * (1.33333333333333D0 * b_h2_so4     &
                     + m_h2so4 * (2.0d0 ** 2.5d0) * 0.33333333333333D0 * c_h2_so4))
         g0_h_hso4(i) = exp(fgamma + m_hhso4 * (b_h_hso4                        &
                     + m_hhso4 * c_h_hso4))
         g0_nh4_no3(i) = exp(fgamma + m_nh4no3 * (b_nh4_no3                     &
                     + m_nh4no3 *c_nh4_no3))
         g0_nh42_so4(i) = exp(2.0d0 * fgamma+m_nh42so4*(1.333333333333333D0*    &
         b_nh42_so4 +m_nh42so4*(2.0d0**2.5d0)*0.33333333333333D0*c_nh42_so4))
         g0_nh4_hso4(i) = exp(fgamma+m_nh4hso4*(b_nh4_hso4 &
                        + m_nh4hso4 * c_nh4_hso4))
      end do
   end if

!  Kusik and Meissner Method:

   if (method == 'KM') then

!  Generate base 10 log part of Gamma0 functions:

!IBM* UNROLL(8)
      do i = 1, ne
         fion = 1.0d0 + 0.1d0 * ion(i)
         lnion = log(fion)

!  Determine residual power for series approximation:
!  Note that in some cases this will (later) require
!  inverting the power portion of the formula.

         p_hno3 = lnion * 0.4d0
         p_hno3 = (1.0d0 + 0.581d0 * (fion * fion * fion /                   &
                  ((((((((((((((pcoef1(1)  * p_hno3 + pcoef1(2))  * p_hno3   &
                             + pcoef1(3))  * p_hno3 + pcoef1(4))  * p_hno3   &
                             + pcoef1(5))  * p_hno3 + pcoef1(6))  * p_hno3   &
                             + pcoef1(7))  * p_hno3 + pcoef1(8))  * p_hno3   &
                             + pcoef1(9))  * p_hno3 + pcoef1(10)) * p_hno3   &
                             + pcoef1(11)) * p_hno3 + pcoef1(12)) * p_hno3   &
                             + pcoef1(13)) * p_hno3 + pcoef1(14)) * p_hno3   &
                             + 1.0d0) - 1.0d0))
         p_hno3    = min(max(1.0d-32, p_hno3),    1.0d+32)

!  Taylor series for q_h2so4 = 0.7:
!
!         p_h2so4_a(i) = (1.D0 + 0.7045D0*( fion/(
!     2    (((((((((((((pcoef1(14)*p_h2so4_a(i)  !  1/14!
!     3     + pcoef1(13))*p_h2so4_a(i)            !  1/13!
!     4     + pcoef1(12))*p_h2so4_a(i)            !  1/12!
!     5     + pcoef1(11))*p_h2so4_a(i)            !  1/11!
!     6     + pcoef1(10))*p_h2so4_a(i)            !  1/10!
!     7     + pcoef1(9))*p_h2so4_a(i)            !  1/9!
!     8     + pcoef1(8))*p_h2so4_a(i)            !  1/8!
!     9     + pcoef1(7))*p_h2so4_a(i)            !  1/7!
!     1     + pcoef1(6))*p_h2so4_a(i)            !  1/6!
!     2     + pcoef1(5))*p_h2so4_a(i)            !  1/5!
!     3     + pcoef1(4))*p_h2so4_a(i)            !  1/4!
!     4     + pcoef1(3))*p_h2so4_a(i)            !  1/3!
!     5     + pcoef1(2))*p_h2so4_a(i)            !  1/2!
!     6     + pcoef1(1))*p_h2so4_a(i)            !  1/1!
!     7     + 1.0d0   ) - 1.0d0 ))
!
!  Taylor series for q_h2so4 = -0.1:
!
         p_h2so4_b = lnion * 0.1d0
         p_h2so4_b = (1.0D0 + 0.7565D0 * (1.0D0 /                            &
              ((((((((((((((pcoef1(1)  * p_h2so4_b + pcoef1(2))  * p_h2so4_b &
                         + pcoef1(3))  * p_h2so4_b + pcoef1(4))  * p_h2so4_b &
                         + pcoef1(5))  * p_h2so4_b + pcoef1(6))  * p_h2so4_b &
                         + pcoef1(7))  * p_h2so4_b + pcoef1(8))  * p_h2so4_b &
                         + pcoef1(9))  * p_h2so4_b + pcoef1(10)) * p_h2so4_b &
                         + pcoef1(11)) * p_h2so4_b + pcoef1(12)) * p_h2so4_b &
                         + pcoef1(13)) * p_h2so4_b + pcoef1(14)) * p_h2so4_b &
                         + 1.0d0) - 1.0d0))
         p_h2so4_b = min(max(1.0d-32, p_h2so4_b), 1.0d+32)
!
         p_h_hso4 = (1.0D0 + 0.23D0 * (fion * fion * fion * fion * fion *    &
                    fion * fion * fion - 1.0D0))
         p_h_hso4  = min(max(1.0d-32, p_h_hso4),  1.0d+32)
!
         p_nh4no3 = lnion * 0.15d0
         p_nh4no3 = (1.0d0 + 0.82475d0 * (1.0d0 / fion /                     &
             ((((((((((((((pcoef1(1)  * p_nh4no3  + pcoef1(2))  * p_nh4no3   &
                        + pcoef1(3))  * p_nh4no3  + pcoef1(4))  * p_nh4no3   &
                        + pcoef1(5))  * p_nh4no3  + pcoef1(6))  * p_nh4no3   &
                        + pcoef1(7))  * p_nh4no3  + pcoef1(8))  * p_nh4no3   &
                        + pcoef1(9))  * p_nh4no3  + pcoef1(10)) * p_nh4no3   &
                        + pcoef1(11)) * p_nh4no3  + pcoef1(12)) * p_nh4no3   &
                        + pcoef1(13)) * p_nh4no3  + pcoef1(14)) * p_nh4no3   &
                        + 1.0d0) - 1.0d0))
         p_nh4no3  = min(max(1.0d-32, p_nh4no3),  1.0d+32)
!
         p_nh42so4 = lnion * 0.25d0
         p_nh42so4 = (1.0d0 + 0.76625d0 * (1.0d0 /                           &
              ((((((((((((((pcoef1(1)  * p_nh42so4 + pcoef1(2))  * p_nh42so4 &
                         + pcoef1(3))  * p_nh42so4 + pcoef1(4))  * p_nh42so4 &
                         + pcoef1(5))  * p_nh42so4 + pcoef1(6))  * p_nh42so4 &
                         + pcoef1(7))  * p_nh42so4 + pcoef1(8))  * p_nh42so4 &
                         + pcoef1(9))  * p_nh42so4 + pcoef1(10)) * p_nh42so4 &
                         + pcoef1(11)) * p_nh42so4 + pcoef1(12)) * p_nh42so4 &
                         + pcoef1(13)) * p_nh42so4 + pcoef1(14)) * p_nh42so4 &
                         + 1.0d0) - 1.0d0))
         p_nh42so4 = min(max(1.0d-32, p_nh42so4), 1.0d+32)

!  Calculate log10(Gamma*), pg 167, Kim, Seinfeld and Saxena, 1993:
         fion = exp(-0.023d0 * ion(i) * ion(i) * ion(i))
         g0_h_no3(i) = (-0.5107d0 * sqion(i) / (1.0d0 +          &
                       (1.0d0 + 0.055d0 * q_hno3 * fion) * sqion(i)))
         g0_h2_so4(i) = (-0.5107d0 * sqion(i) / (1.0d0 +         &
                       (1.0d0 + 0.055d0 * q_h2so4 * fion) * sqion(i)))
         g0_h_hso4(i) = (-0.5107d0 * sqion(i) / (1.0d0 +         &
                       (1.0d0 + 0.055d0 * q_h_hso4 * fion) * sqion(i)))
         g0_nh4_no3(i) = (-0.5107d0 * sqion(i) / (1.0d0 +        &
                       (1.0d0 + 0.055d0 * q_nh4no3 * fion) * sqion(i)))
         g0_nh42_so4(i) = (-0.5107d0 * sqion(i) / (1.0d0 +       &
                       (1.0d0 + 0.055d0 * q_nh42so4 * fion) * sqion(i)))
!
         g0_h_no3(i)    = g0_h_no3(i)    + log10(p_hno3)
         g0_h2_so4(i)   = g0_h2_so4(i)   + log10(p_h2so4_b)
         g0_h_hso4(i)   = g0_h_hso4(i)   + log10(p_h_hso4)
         g0_nh4_no3(i)  = g0_nh4_no3(i)  + log10(p_nh4no3)
         g0_nh42_so4(i) = g0_nh42_so4(i) + log10(p_nh42so4)

!  Correct Gamma0 at 25C to the value at the gridpoint temperature:
!  base 10 log of correction factor:

         tc =   t(i) - 273.15d0
         cfac = ((5.0D-03 * tc - 0.125D0) *  &
         (3.9D-02 * ion(i) ** 0.92 - 0.41D0 * sqion(i) / (1.0D0 + sqion(i))))
         tfac = 1.125D0 - 5D-03 * tc
         g0_h_no3(i)    = cfac + tfac * g0_h_no3(i)
         g0_h2_so4(i)   = cfac + tfac * g0_h2_so4(i)
         g0_h_hso4(i)   = cfac + tfac * g0_h_hso4(i)
         g0_nh4_no3(i)  = cfac + tfac * g0_nh4_no3(i)
         g0_nh42_so4(i) = cfac + tfac * g0_nh42_so4(i)

!  calculate binary activity coefficients, k&m:
         g0_h2_so4(i) = g0_h2_so4(i) + g0_h2_so4(i)
         g0_nh42_so4(i) = g0_nh42_so4(i) + g0_nh42_so4(i)
         g0_nh4_hso4(i) = 0.5d0 * (g0_nh42_so4(i) + g0_h2_so4(i))
      end do
   end if

!  (eq 14, 15) X and Y factors
   do i = 1, ne
      ttt = 1.0d0 / ion(i)
      x_h_no3 = max(m_h(i), ionmin) * ttt
      y_no3_h = max(m_no3(i), ionmin) * ttt
      x_h2_so4 = 2.25d0 * max(m_h(i), ionmin) * ttt
      y_so4_h2 = 2.25d0 * max(m_so4(i), ionmin) * ttt
      x_h_hso4 = max(m_h(i), ionmin) * ttt
      y_hso4_h = max(m_hso4(i), ionmin)*ttt
!      x_h_oh = max(m_h(i), ionmin) * ttt
      y_oh_h = max(m_oh(i), ionmin) * ttt
      x_nh4_no3 = max(m_nh4(i), ionmin) * ttt
      y_no3_nh4 = max(m_no3(i), ionmin) * ttt
      x_nh42_so4 = 2.25d0 * max(m_nh4(i), ionmin) *ttt
      y_so4_nh42 = 2.25d0 * max(m_so4(i), ionmin) *ttt
      x_nh4_hso4 = max(m_nh4(i), ionmin) * ttt
      y_hso4_nh4 = max(m_hso4(i), ionmin) * ttt
!      x_nh4_oh = max(m_nh4(i), ionmin) * ttt
      y_oh_nh4 = max(m_oh(i), ionmin) * ttt
!
!  Cations (F1's of Pilinis and Seinfeld):
!  Note:  g0_nh4_oh and g_h_oh are both 1.0, so
!  their log10 is zero, dropping one of the terms
!  of the summation.
      ttt = 1.0D0 / (1.0D0 + sqion(i))
      f_nh4(i) = y_no3_nh4 * (g0_nh4_no3(i)) + y_so4_nh42 * (g0_nh42_so4(i))    &
               + y_hso4_nh4 * (g0_nh4_hso4(i)) + ag * sqion(i) * ttt *       &
               (y_no3_nh4 + 2.0D0 * y_so4_nh42 + y_hso4_nh4 + y_oh_nh4)

      f_h(i) = y_no3_h * (g0_h_no3(i)) + y_so4_h2 * (g0_h2_so4(i))              &
               + y_hso4_h * (g0_h_hso4(i)) + ag * sqion(i) * ttt *           &
               (y_no3_h + 2.0D0 * y_so4_h2 + y_hso4_h + y_oh_h)

!  Anions (F2's of Pilinis and Seinfeld):
      f_no3(i) = x_h_no3 * (g0_h_no3(i)) + x_nh4_no3 * (g0_nh4_no3(i))          &
               + ag * sqion(i) * ttt * (x_h_no3 + x_nh4_no3)

      f_hso4(i) = x_h_hso4 * (g0_h_hso4(i)) + x_nh4_hso4 * (g0_nh4_hso4(i))     &
               + ag * sqion(i) * ttt * (x_h_hso4 + x_nh4_hso4)

  !    f_oh(i) = ag * sqion(i) * ttt * (x_h_oh + x_nh4_oh)

      f_so4(i) = x_h2_so4 * (g0_h2_so4(i)) + x_nh42_so4 * (g0_nh42_so4(i))      &
               + ag * sqion(i) * ttt * 2.0D0 *( x_h2_so4 + x_nh42_so4)
   end do
!
!  Activity Coefficients:
   do i = 1, ne
      ttt = -ag * sqion(i) / (1.0d0 + sqion(i))
!  H+, HSO4-:
      g_h_hso4(i) = ttt + 0.5d0 * (f_h(i) + f_hso4(i))
!  2 H+, SO4=:
      g_h2_so4(i) = 2.0d0 * (ttt + 0.333333333333333d0 * (f_h(i) + f_so4(i) * 0.5d0))
!
!  Limit the range of activity coefficients to be between
!  -30 and + 30.  This is required (a) as a lower/upper
!  number limit exceedence prevention step and (b) in order
!  to make use of the following approximation for 10**p,
!  where the value of p is subdivided into an integer and
!  residual part.  ISORROPIA allowed a lookup table search
!  as a low accuracy way of getting the residual power;
!  here, a higher accuracy series is used.
!
!  H+, HSO4-:
      g_h_hso4(i) = min(max(g_h_hso4(i), -30.0d0), 30.0d0)
!  2 H+, SO4=:
      g_h2_so4(i) = min(max(g_h2_so4(i), -30.0d0), 30.0d0)
!
!  Note that the algorithm employed here makes use of
!  sign functions in order to ensure that the smallest
!  possible value is used for the operand of the series
!  approximation to the base ten power.  Thus 10**6.78
!  becomes 10**7 x 10**(-0.22), 10**6.32 becomes
!  10**6 x 10**0.32, 10**(-6.78) becomes 10.**(-7) x 10.**(0.22)
!  and 10.**(-6.32) becomes 10.**(-6) x 10.**(-0.32), with
!  the second of the two powers in each case being approximated
!  by a 14th order series for 10**p, -1<p<1.  By keeping
!  -0.5 < p < 0.5, the accuracy of the series improves from
!  six figures to 10 figures.
!
!  Separate the power into adjusted modulus and remainder.
!  H+, HSO4-:
      g_h_hso4_p = dble(int(g_h_hso4(i)))
      g_h_hso4_r = g_h_hso4(i) - g_h_hso4_p
      g_h_hso4_v = sign(1.0d0, g_h_hso4_r) * 0.5d0 *  &
      (1.0d0 + sign(1.0d0, (g_h_hso4_r * g_h_hso4_r - 0.25d0)))
      g_h_hso4_p = g_h_hso4_p + g_h_hso4_v
      g_h_hso4_r = g_h_hso4_r - g_h_hso4_v
!
!  2 H+, SO4=:
      g_h2_so4_p = dble(int(g_h2_so4(i)))
      g_h2_so4_r = g_h2_so4(i) - g_h2_so4_p
      g_h2_so4_v = sign(1.0d0, g_h2_so4_r) * 0.5d0 *  &
      (1.0d0 + sign(1.0d0, (g_h2_so4_r * g_h2_so4_r - 0.25d0)))
      g_h2_so4_p = g_h2_so4_p + g_h2_so4_v
      g_h2_so4_r = g_h2_so4_r - g_h2_so4_v
!
!  Use power series + indirect addressing to determine value of
!  power.  Coefficients are (ln(10)^14)/14!, (ln(10)^13)/13!, etc.,
!  in descending order.
!  H+, HSO4-:
      g_h_hso4(i) = pten(int(g_h_hso4_p) + 32) *                              &
             ((((((((((((((pcoef2(1)  * g_h_hso4_r + pcoef2(2))  * g_h_hso4_r &
                        + pcoef2(3))  * g_h_hso4_r + pcoef2(4))  * g_h_hso4_r &
                        + pcoef2(5))  * g_h_hso4_r + pcoef2(6))  * g_h_hso4_r &
                        + pcoef2(7))  * g_h_hso4_r + pcoef2(8))  * g_h_hso4_r &
                        + pcoef2(9))  * g_h_hso4_r + pcoef2(10)) * g_h_hso4_r &
                        + pcoef2(11)) * g_h_hso4_r + pcoef2(12)) * g_h_hso4_r &
                        + pcoef2(13)) * g_h_hso4_r + pcoef2(14)) * g_h_hso4_r &
                        + 1.0d0)

!  2 H+, SO4=:
      g_h2_so4(i) = pten(int(g_h2_so4_p) + 32) *                              &
            ((((((((((((((pcoef2(1)  * g_h2_so4_r + pcoef2(2))  * g_h2_so4_r  &
                       + pcoef2(3))  * g_h2_so4_r + pcoef2(4))  * g_h2_so4_r  &
                       + pcoef2(5))  * g_h2_so4_r + pcoef2(6))  * g_h2_so4_r  &
                       + pcoef2(7))  * g_h2_so4_r + pcoef2(8))  * g_h2_so4_r  &
                       + pcoef2(9))  * g_h2_so4_r + pcoef2(10)) * g_h2_so4_r  &
                       + pcoef2(11)) * g_h2_so4_r + pcoef2(12)) * g_h2_so4_r  &
                       + pcoef2(13)) * g_h2_so4_r + pcoef2(14)) * g_h2_so4_r  &
                       + 1.0d0)

!  Impose maximum and minimum limits on activity coefficients.
!  This becomes necessary for cases in which the activities
!  when taken to powers in the equilibrium expressions, might
!  otherwise result in overflow or underflows.  Greatest
!  power used is ^5, and maxima and minima of 1E5, 1d-5 are
!  employed, limiting the effects to 1E25, 1d-25 (for letovicite).

      g_h_hso4(i) = min(gmax, max(gmin, g_h_hso4(i)))
      g_h2_so4(i) = min(gmax, max(gmin, g_h2_so4(i)))
   end do
!
   if (present(g_h_no3)) then
      do i = 1, ne
         ttt = -ag * sqion(i) / (1.0d0 + sqion(i))
!  H+, NO3-:
         g_h_no3(i) = ttt + 0.5d0 * (f_h(i) + f_no3(i))
!  NH4+, NO3-:
         g_nh4_no3(i) = ttt + 0.5d0 * (f_nh4(i) + f_no3(i))
!
!  H+, NO3-:
         g_h_no3(i) = min(max(g_h_no3(i), -30.0d0), 30.0d0)
!  NH4+, NO3-:
         g_nh4_no3(i) = min(max(g_nh4_no3(i), -30.0d0), 30.0d0)
!
!  H+, NO3-:
         g_h_no3_p = dble(int(g_h_no3(i)))
         g_h_no3_r = g_h_no3(i) - g_h_no3_p
         g_h_no3_v = sign(1.0d0, g_h_no3_r) * 0.5d0 *  &
                     (1.0d0 + sign(1.0d0, (g_h_no3_r * g_h_no3_r - 0.25d0)))
         g_h_no3_p = g_h_no3_p + g_h_no3_v
         g_h_no3_r = g_h_no3_r - g_h_no3_v
         g_h_no3(i) = pten(int(g_h_no3_p) + 32) *                             &
             ((((((((((((((pcoef2(1)  * g_h_no3_r + pcoef2(2))  * g_h_no3_r   &
                        + pcoef2(3))  * g_h_no3_r + pcoef2(4))  * g_h_no3_r   &
                        + pcoef2(5))  * g_h_no3_r + pcoef2(6))  * g_h_no3_r   &
                        + pcoef2(7))  * g_h_no3_r + pcoef2(8))  * g_h_no3_r   &
                        + pcoef2(9))  * g_h_no3_r + pcoef2(10)) * g_h_no3_r   &
                        + pcoef2(11)) * g_h_no3_r + pcoef2(12)) * g_h_no3_r   &
                        + pcoef2(13)) * g_h_no3_r + pcoef2(14)) * g_h_no3_r   &
                        + 1.0d0)
         g_h_no3(i) = min(gmax, max(gmin, g_h_no3(i)))
!
!  NH4+, NO3-:                                              APPEND MODE
         g_nh4_no3_p = dble(int(g_nh4_no3(i)))
         g_nh4_no3_r = g_nh4_no3(i) - g_nh4_no3_p
         g_nh4_no3_v = sign(1.0d0, g_nh4_no3_r) * 0.5d0 *  &
                 (1.0d0 + sign(1.0d0, (g_nh4_no3_r * g_nh4_no3_r - 0.25d0)))
         g_nh4_no3_p = g_nh4_no3_p + g_nh4_no3_v
         g_nh4_no3_r = g_nh4_no3_r - g_nh4_no3_v
         g_nh4_no3(i) = pten(int(g_nh4_no3_p) + 32) *                         &
          ((((((((((((((pcoef2(1)   * g_nh4_no3_r + pcoef2(2))  * g_nh4_no3_r &
                      + pcoef2(3))  * g_nh4_no3_r + pcoef2(4))  * g_nh4_no3_r &
                      + pcoef2(5))  * g_nh4_no3_r + pcoef2(6))  * g_nh4_no3_r &
                      + pcoef2(7))  * g_nh4_no3_r + pcoef2(8))  * g_nh4_no3_r &
                      + pcoef2(9))  * g_nh4_no3_r + pcoef2(10)) * g_nh4_no3_r &
                      + pcoef2(11)) * g_nh4_no3_r + pcoef2(12)) * g_nh4_no3_r &
                      + pcoef2(13)) * g_nh4_no3_r + pcoef2(14)) * g_nh4_no3_r &
                      + 1.0d0)
         g_nh4_no3(i) = min(gmax, max(gmin, g_nh4_no3(i)))

      end do
   end if
!
!  NH4+, HSO4-:
   if (present(g_nh4_hso4) .or. present(g_nh43_hso42)) then
      do i = 1, ne
         ttt = -ag * sqion(i) / (1.0d0 + sqion(i))
         g_nh4_hso4(i) = ttt + 0.5d0 * (f_nh4(i) + f_hso4(i))
         g_nh4_hso4(i) = min(max(g_nh4_hso4(i), -30.0d0), 30.0d0)
!
         g_nh4_hso4_p = dble(int(g_nh4_hso4(i)))
         g_nh4_hso4_r = g_nh4_hso4(i) - g_nh4_hso4_p
         g_nh4_hso4_v = sign(1.0d0, g_nh4_hso4_r) * 0.5d0 *  &
             (1.0d0 + sign(1.0d0, (g_nh4_hso4_r * g_nh4_hso4_r - 0.25d0)))
         g_nh4_hso4_p = g_nh4_hso4_p + g_nh4_hso4_v
         g_nh4_hso4_r = g_nh4_hso4_r - g_nh4_hso4_v
         g_nh4_hso4(i) = pten(int(g_nh4_hso4_p) + 32) *                       &
         ((((((((((((((pcoef2(1)  * g_nh4_hso4_r + pcoef2(2))  * g_nh4_hso4_r &
                    + pcoef2(3))  * g_nh4_hso4_r + pcoef2(4))  * g_nh4_hso4_r &
                    + pcoef2(5))  * g_nh4_hso4_r + pcoef2(6))  * g_nh4_hso4_r &
                    + pcoef2(7))  * g_nh4_hso4_r + pcoef2(8))  * g_nh4_hso4_r &
                    + pcoef2(9))  * g_nh4_hso4_r + pcoef2(10)) * g_nh4_hso4_r &
                    + pcoef2(11)) * g_nh4_hso4_r + pcoef2(12)) * g_nh4_hso4_r &
                    + pcoef2(13)) * g_nh4_hso4_r + pcoef2(14)) * g_nh4_hso4_r &
                    + 1.0d0)
         g_nh4_hso4(i) = min(gmax, max(gmin, g_nh4_hso4(i)))
      end do
   end if
!
!  (NH4+)2, SO4=:
   if (present(g_nh42_so4) .or. present(g_nh43_hso42)) then
      do i = 1, ne
         ttt = -ag * sqion(i) / (1.0d0 + sqion(i))
         g_nh42_so4(i) = 2.0d0 * (ttt + 0.3333333333333333d0 * (f_nh4(i) + f_so4(i) * 0.5d0))
         g_nh42_so4(i) = min(max(g_nh42_so4(i), -30.0d0), 30.0d0)
!
         g_nh42_so4_p = dble(int(g_nh42_so4(i)))
         g_nh42_so4_r = g_nh42_so4(i) - g_nh42_so4_p
         g_nh42_so4_v = sign(1.0d0, g_nh42_so4_r) * 0.5d0 *  &
             (1.0d0 + sign(1.0d0, (g_nh42_so4_r * g_nh42_so4_r - 0.25d0)))
         g_nh42_so4_p = g_nh42_so4_p + g_nh42_so4_v
         g_nh42_so4_r = g_nh42_so4_r - g_nh42_so4_v
         g_nh42_so4(i) = pten(int(g_nh42_so4_p) + 32) *                       &
        ((((((((((((((pcoef2(1)  * g_nh42_so4_r + pcoef2(2))  * g_nh42_so4_r  &
                   + pcoef2(3))  * g_nh42_so4_r + pcoef2(4))  * g_nh42_so4_r  &
                   + pcoef2(5))  * g_nh42_so4_r + pcoef2(6))  * g_nh42_so4_r  &
                   + pcoef2(7))  * g_nh42_so4_r + pcoef2(8))  * g_nh42_so4_r  &
                   + pcoef2(9))  * g_nh42_so4_r + pcoef2(10)) * g_nh42_so4_r  &
                   + pcoef2(11)) * g_nh42_so4_r + pcoef2(12)) * g_nh42_so4_r  &
                   + pcoef2(13)) * g_nh42_so4_r + pcoef2(14)) * g_nh42_so4_r  &
                   + 1.0d0)
         g_nh42_so4(i) = min(gmax, max(gmin, g_nh42_so4(i)))
      end do
   end if
!

!  Final adjustments for constant and letovicite values:
!   do i = 1, ne
!!  H+, OH-
!      g_h_oh(i) = 1.0D0
!!  NH4+, OH-:
!      g_nh4_oh(i) = 1.0D0
!   end do
   if (present(g_nh43_hso42)) then
!  (NH4)3HSO4:
!      g_nh43_hso42 = min(gmax, max(gmin, g_nh43_hso42))
      g_nh43_hso42 = sqrt(g_nh4_hso4 * g_nh42_so4) ! Array operations
   end if

   return
end subroutine mach_hetv_activity
