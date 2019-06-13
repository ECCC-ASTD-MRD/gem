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
!**S/P getbpm
      function getbpm()
      implicit none
!!!#include <arch_specific.hf>
      integer getbpm
!Author   R. Benoit (Sept 93)
!
!Arguments
!         None
!
!Object   Find number of bits per integer word
!
!*
!
      integer i, m
!
      i=1
      m=1
 1    i=2*i+1
      m=m+1
      if (i.gt.0) go to 1
!
!      print *,' getbpm.  sign reversal at m=', m
      getbpm=m
      return
      end

