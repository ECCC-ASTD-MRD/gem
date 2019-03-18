!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms 
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer 
!version 3 or (at your option) any later version that should be found at: 
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html 
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software; 
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec), 
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
!*** S/P MO_GWSPACTRUM
MODULE mo_gwspectrum

  REAL*8  :: kstar=5.E-5
  INTEGER :: naz=8
  REAL*8  :: m_min=5.23598E-4
  REAL*8  :: slope=1.
  REAL*8  :: f1=1.5
  REAL*8  :: f2=0.3
  REAL*8  :: f3=1.0
  REAL*8  :: f5=2.0
  REAL*8  :: f6=1.0
  REAL*8  :: alt_cutoff=105.E3
  INTEGER :: icutoff=0
  REAL*8  :: smco=2.0
  INTEGER :: nsmax=5

END MODULE mo_gwspectrum
